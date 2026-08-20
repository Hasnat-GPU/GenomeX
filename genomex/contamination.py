"""Is this assembly one organism, or more than one?

Four independent signals, deliberately chosen so that no single one can carry a
verdict on its own:

  1. Tetranucleotide frequency (TNF).  Genome-wide 4-mer composition is
     species-characteristic and roughly constant along a chromosome, so a contig
     whose TNF vector sits far from the assembly's core composition came from
     somewhere else.  (Teeling et al. 2004; the signal every binner uses.)
  2. GC content.  Cruder than TNF and correlated with it, but interpretable and
     hard to argue with when it is extreme.
  3. Duplicated single-copy markers.  If a marker that must occur once per genome
     occurs on two unrelated contigs, either the assembly is redundant or it
     contains two genomes.  This is how CheckM detects contamination.
  4. Bimodality.  A deterministic 2-means split in TNF-PCA space: if the
     assembly really is a mixture, the split is stable, both bins are large, and
     they differ in GC.  If it is one organism, the split is arbitrary.

Signals are reported separately and then combined, so a reader can disagree with
the verdict while still trusting the evidence.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .composition import (
    N_CANON,
    NullCurve,
    base_tiles,
    chi2_divergence,
    counts_from_codes,
    frequencies,
    fwer_threshold,
    gc_percent_from_freq,
    kmer_codes,
    null_curve,
)
from .fasta import Assembly, Contig

# 4-mers: 256 raw, 136 canonical after merging each with its reverse complement.
K = 4
_COMPLEMENT = {0: 3, 1: 2, 2: 1, 3: 0}  # A<->T, C<->G with A=0,C=1,G=2,T=3

_LUT = np.full(256, 255, dtype=np.uint8)
for _base, _code in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    _LUT[_base] = _code


def _revcomp_index(idx: int) -> int:
    """Reverse-complement of a 4-mer expressed as a base-4 integer."""
    out = 0
    for i in range(K):
        base = (idx >> (2 * i)) & 3
        out = (out << 2) | _COMPLEMENT[base]
    return out


_CANON_GROUP = {}
_canon_order: list[int] = []
for _i in range(4**K):
    _rc = _revcomp_index(_i)
    _key = min(_i, _rc)
    if _key not in _CANON_GROUP:
        _CANON_GROUP[_key] = len(_canon_order)
        _canon_order.append(_key)
N_CANON = len(_canon_order)  # 136

#: Marker count of `bacteria_odb10`, the set the duplication cutoffs below were
#: originally written against.
_REFERENCE_MARKER_SET = 124

#: Verdict cutoffs on the *rate* of cross-contig marker duplication, not on its
#: count. A mixture holding a foreign genome at mass fraction f displaces roughly
#: f x N of the N single-copy markers onto foreign contigs, so the rate estimates
#: f and does not depend on N; the raw count scales with the marker set and means
#: nothing on its own.
#:
#: Counting is what a lineage-specific set breaks. Swapping bacteria_odb10 (124
#: markers) for burkholderiales_odb10 (688) multiplies every count by 5.5 without
#: adding one base of contamination, and finished reference genomes carrying an
#: ordinary 1% of cross-replicon paralogs jumped straight to `likely`.
#:
#: The two values are the previous cutoffs of 5 and 2 divided by 124. That is a
#: change of units, not a retuning: on bacteria_odb10 the rule is unchanged.
DUP_RATE_LIKELY = 5 / _REFERENCE_MARKER_SET      # 4.03%
DUP_RATE_POSSIBLE = 2 / _REFERENCE_MARKER_SET    # 1.61%
_RAW_TO_CANON = np.array(
    [_CANON_GROUP[min(i, _revcomp_index(i))] for i in range(4**K)], dtype=np.int32
)


def tetranucleotide_freq(seq: str) -> np.ndarray:
    """Canonical 4-mer frequency vector (length 136) for one sequence."""
    codes = _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]
    if codes.size < K:
        return np.zeros(N_CANON, dtype=np.float64)
    valid = codes != 255
    windows = np.lib.stride_tricks.sliding_window_view(codes, K)
    window_ok = np.lib.stride_tricks.sliding_window_view(valid, K).all(axis=1)
    if not window_ok.any():
        return np.zeros(N_CANON, dtype=np.float64)
    w = windows[window_ok].astype(np.int32)
    raw = (w[:, 0] << 6) | (w[:, 1] << 4) | (w[:, 2] << 2) | w[:, 3]
    counts = np.bincount(_RAW_TO_CANON[raw], minlength=N_CANON).astype(np.float64)
    total = counts.sum()
    return counts / total if total else counts


def _robust_z(values: np.ndarray) -> np.ndarray:
    """Median/MAD z-score. Robust to the very outliers we are hunting for."""
    med = float(np.median(values))
    mad = float(np.median(np.abs(values - med)))
    scale = 1.4826 * mad
    if scale < 1e-12:
        std = float(values.std())
        scale = std if std > 1e-12 else 1.0
    return (values - med) / scale


def _two_means(points: np.ndarray, iterations: int = 50) -> np.ndarray:
    """Deterministic 2-means: seed at the extremes of the first component.

    No RNG, so a run is reproducible byte-for-byte -- which matters more here
    than the marginal quality of a random restart.
    """
    if len(points) < 2:
        return np.zeros(len(points), dtype=int)
    axis = points[:, 0]
    centroids = np.array([points[axis.argmin()], points[axis.argmax()]], dtype=float)
    labels = np.zeros(len(points), dtype=int)
    for _ in range(iterations):
        d = np.linalg.norm(points[:, None, :] - centroids[None, :, :], axis=2)
        new_labels = d.argmin(axis=1)
        if np.array_equal(new_labels, labels):
            break
        labels = new_labels
        for c in (0, 1):
            if (labels == c).any():
                centroids[c] = points[labels == c].mean(axis=0)
    return labels


@dataclass
class ContigVerdict:
    name: str
    length: int
    gc_percent: float
    tnf_distance: float
    tnf_z: float
    gc_z: float
    duplicated_markers: int          # duplicated markers whose other copies sit elsewhere
    bin_label: int
    flags: list[str] = field(default_factory=list)
    suspicion: float = 0.0
    # core                  composition matches the assembly and no marker conflict
    # contaminant_candidate compositionally foreign; the real contamination signal
    # replicon_candidate    large and distinct but no displaced markers: plasmid/chromid
    # atypical_host_region  distinct composition, but holds the assembly's only
    #                       copies of core markers -- this organism's own
    #                       chromosome; prophage or island, not contamination
    # marker_conflict       shares duplicated core markers but composition is typical
    call: str = "core"

    @property
    def suspect(self) -> bool:
        return bool(self.flags)


@dataclass
class ContaminationResult:
    verdict: str                      # clean | possible | likely
    reasons: list[str]
    suspect_bp: int
    suspect_fraction: float
    n_suspect_contigs: int
    contigs: list[ContigVerdict]
    bins: dict
    params: dict

    def summary(self) -> dict:
        return {
            "verdict": self.verdict,
            "reasons": self.reasons,
            "suspect_contigs": self.n_suspect_contigs,
            "suspect_bp": self.suspect_bp,
            "suspect_fraction_percent": round(100 * self.suspect_fraction, 2),
            "replicon_candidates": [
                {"contig": c.name, "length": c.length, "gc_percent": c.gc_percent}
                for c in self.contigs if c.call == "replicon_candidate"
            ],
            "atypical_host_regions": [
                {"contig": c.name, "length": c.length, "gc_percent": c.gc_percent}
                for c in self.contigs if c.call == "atypical_host_region"
            ],
            "marker_conflict_contigs": [
                c.name for c in self.contigs if c.call == "marker_conflict"
            ],
            "bins": self.bins,
            "params": self.params,
        }

    def suspect_contig_names(self) -> set[str]:
        """Contigs believed to be foreign.

        Excludes plasmid/chromid candidates and atypical host regions alike:
        the genes on both are this organism's real biology, and the comparative
        step writes off everything named here as an artefact."""
        return {c.name for c in self.contigs if c.call == "contaminant_candidate"}

    def flagged_contig_names(self) -> set[str]:
        """Every contig with any composition flag, replicon candidates included."""
        return {c.name for c in self.contigs if c.suspect}


def detect_contamination(
    asm: Assembly,
    *,
    duplicated_marker_contigs: dict[str, list[str]] | None = None,
    contig_marker_counts: dict[str, int] | None = None,
    total_markers: int | None = None,
    min_contig_length: int = 3000,
    tnf_z_threshold: float = 4.0,
    gc_z_threshold: float = 4.0,
    marker_dup_threshold: int = 2,
    min_contigs_for_zscores: int = 8,
    tnf_distance_threshold: float = 0.02,
    gc_absolute_threshold: float = 5.0,
    replicon_min_length: int = 50_000,
    alpha: float = 0.05,
) -> ContaminationResult:
    """Score every contig, then combine per-contig evidence into one verdict.

    Two corrections that real assemblies force:

    * Below `min_contigs_for_zscores` sequences there is no population to compute
      a robust z-score against -- a finished 3-replicon genome would otherwise
      produce z=74 from a MAD of nearly zero.  Absolute thresholds are used there
      instead, and the verdict stays deliberately weaker.
    * A large, compositionally distinct contig carrying no displaced single-copy
      markers is at least as likely to be a plasmid or a chromid as a
      contaminant.  Composition alone cannot separate those, so it is called a
      `replicon_candidate` and excluded from the contamination fraction.
    * A compositionally distinct group holding the assembly's *only* copies of
      part of the single-copy core is this organism's own chromosome: a second
      organism arrives with its own copies, which duplicate the host's.  It is
      called an `atypical_host_region` -- a prophage or island candidate -- and
      likewise excluded.  Calling it contamination instead turned every outlier
      the family-wise threshold admits into a genome-level verdict, which cost
      12.5% false positives on shreddings of finished genomes; see
      `docs/benchmark-fragmentation.md`.
    """
    duplicated_marker_contigs = duplicated_marker_contigs or {}

    # Only *cross-contig* duplication is contamination evidence at contig level:
    # a marker duplicated twice on one contig is a paralog, not a second organism.
    dup_per_contig: dict[str, int] = {}
    for _busco, contigs in duplicated_marker_contigs.items():
        distinct = set(contigs)
        if len(distinct) < 2:
            continue
        for name in distinct:
            dup_per_contig[name] = dup_per_contig.get(name, 0) + 1

    scored: list[Contig] = [c for c in asm.contigs if c.length >= min_contig_length]
    total_bp = asm.total_bp

    if len(scored) < 3:
        return ContaminationResult(
            verdict="undetermined",
            reasons=[
                f"only {len(scored)} contigs of at least {min_contig_length} bp -- "
                "composition-based detection needs more sequence to establish a baseline"
            ],
            suspect_bp=0,
            suspect_fraction=0.0,
            n_suspect_contigs=0,
            contigs=[],
            bins={},
            params={"min_contig_length": min_contig_length},
        )

    # Length-conditioned scoring. Every contig is scored against a null built
    # from windows of its own length tiled inside this assembly, so a 3 kb contig
    # is compared with native 3 kb sequence rather than with a 4 Mb chromosome.
    codes = [kmer_codes(c.seq) for c in scored]
    tile_sets = [base_tiles(cd) for cd in codes]
    totals = [counts_from_codes(cd) for cd in codes]
    curve = null_curve(tile_sets, totals)

    tnf = np.vstack([frequencies(t) for t in totals])
    lengths = np.array([c.length for c in scored], dtype=float)
    gc = np.array([100.0 * c.gc for c in scored], dtype=float)

    # Core composition = length-weighted centroid. Long contigs define "self".
    weights = lengths / lengths.sum()
    centroid = (tnf * weights[:, None]).sum(axis=0)

    # Correlation distance to the core: scale-free and the standard TNF measure.
    tc = tnf - tnf.mean(axis=1, keepdims=True)
    cc = centroid - centroid.mean()
    denom = np.linalg.norm(tc, axis=1) * np.linalg.norm(cc)
    corr = np.divide((tc * cc).sum(axis=1), denom, out=np.zeros(len(scored)), where=denom > 0)
    tnf_distance = 1.0 - corr

    # Divergence of each contig from the assembly excluding itself, standardised
    # by what native sequence of that length does. `_robust_z` is kept for the
    # degenerate case where a genome cannot supply enough windows for a curve.
    def _score(exclude: set[int]) -> tuple[np.ndarray, np.ndarray, np.ndarray, NullCurve]:
        """Score every contig against the assembly with `exclude` held out.

        Held-out contigs are still scored -- they are simply not allowed to define
        what native looks like. Both the reference composition and the null curve
        exclude them, because a contaminant that contributes to its own reference
        hides, and one that contributes to the null widens it for everyone.
        """
        keep = [i for i in range(len(scored)) if i not in exclude]
        ref_pool = np.sum([totals[i] for i in keep], axis=0) if keep else np.sum(totals, axis=0)
        curve_ = null_curve(
            [tile_sets[i] for i in keep], [totals[i] for i in keep]
        ) if keep else curve0
        chi2_ = np.zeros(len(scored))
        tnf_z_ = np.zeros(len(scored))
        gc_z_ = np.zeros(len(scored))
        for i, c in enumerate(scored):
            # Leave-one-out only for contigs still inside the reference pool.
            ref_counts = ref_pool - totals[i] if i in keep else ref_pool
            ref = frequencies(ref_counts) if ref_counts.sum() > 0 else centroid
            chi2_[i] = chi2_divergence(tnf[i], ref)
            mu_i, sd_i = curve_.at(c.length)
            tnf_z_[i] = (math.log(chi2_[i]) - mu_i) / sd_i if chi2_[i] > 0 else 0.0
            gc_z_[i] = (gc[i] - gc_percent_from_freq(ref)) / curve_.gc_sigma_at(c.length)
        return chi2_, tnf_z_, gc_z_, curve_

    curve0 = curve
    chi2, tnf_z, gc_z, curve = _score(set())

    # One refit, deterministic, no convergence loop. The first pass is computed
    # against a reference that includes whatever foreign sequence is present; at
    # appreciable contamination that drags the reference toward the contaminant
    # and makes the host look divergent. Re-scoring against the survivors fixes it.
    first_pass_z = fwer_threshold(len(scored), alpha=alpha)
    outliers = {i for i in range(len(scored)) if tnf_z[i] > first_pass_z}
    refit_applied = bool(outliers) and len(outliers) < len(scored) - 1
    if refit_applied:
        chi2, tnf_z, gc_z, curve = _score(outliers)

    # PCA on centred TNF for the bimodality check and for plotting.
    centred = tnf - tnf.mean(axis=0)
    try:
        _u, _s, vt = np.linalg.svd(centred, full_matrices=False)
        pcs = centred @ vt[:2].T
    except np.linalg.LinAlgError:  # pragma: no cover - degenerate input
        pcs = np.zeros((len(scored), 2))
    labels = _two_means(pcs)

    bin_stats = {}
    for c in (0, 1):
        mask = labels == c
        if mask.any():
            bin_stats[f"bin{c}"] = {
                "contigs": int(mask.sum()),
                "bp": int(lengths[mask].sum()),
                "bp_fraction_percent": round(100 * float(lengths[mask].sum()) / total_bp, 2),
                "mean_gc_percent": round(float(np.average(gc[mask], weights=lengths[mask])), 2),
            }
    gc_gap = 0.0
    if "bin0" in bin_stats and "bin1" in bin_stats:
        gc_gap = abs(bin_stats["bin0"]["mean_gc_percent"] - bin_stats["bin1"]["mean_gc_percent"])
    minor_fraction = (
        min(b["bp_fraction_percent"] for b in bin_stats.values()) if len(bin_stats) == 2 else 0.0
    )
    bins = {
        **bin_stats,
        "gc_gap_percent": round(gc_gap, 2),
        "minor_bin_fraction_percent": round(minor_fraction, 2),
        "bimodal": bool(
            gc_gap >= 3.0
            and minor_fraction >= 5.0
            and len(scored) >= min_contigs_for_zscores
        ),
    }

    # One cutoff, stated as an error rate rather than a constant. A fixed z of 4.0
    # meant a different family-wise error rate on a 3-contig genome than on a
    # 1243-contig one, so the same assembly shredded finer was held to a laxer
    # standard per contig and a stricter one overall.
    z_star = fwer_threshold(len(scored), alpha=alpha)
    use_zscores = not curve.theoretical
    core_gc = float(np.average(gc, weights=lengths))

    # Group compositionally flagged contigs into candidate replicons before
    # judging them. A plasmid is one biological object whether the assembler
    # emitted it as a single 785 kb contig or as thirty 27 kb ones, so the
    # replicon test has to be applied to the group's mass, not to each fragment's
    # length. Judging per contig meant that shredding a megaplasmid turned one
    # `replicon_candidate` into thirty `contaminant_candidate`s -- fragmentation
    # changing the biology.
    #
    # Two flagged contigs join the same group when they diverge from each other
    # no more than native sequence of their lengths does, measured against the
    # same null curve used for everything else.
    flagged_idx = [
        i for i in range(len(scored))
        if tnf_z[i] > z_star or abs(gc_z[i]) > z_star
    ]
    parent: dict[int, int] = {i: i for i in flagged_idx}

    def _find(a: int) -> int:
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for ai, i in enumerate(flagged_idx):
        for j in flagged_idx[ai + 1:]:
            shorter = min(scored[i].length, scored[j].length)
            mu_ij, sd_ij = curve.at(shorter)
            d_ij = chi2_divergence(tnf[i], tnf[j])
            if d_ij <= 0:
                z_ij = -math.inf
            else:
                z_ij = (math.log(d_ij) - mu_ij) / sd_ij
            if z_ij <= z_star:
                ri, rj = _find(i), _find(j)
                if ri != rj:
                    parent[ri] = rj

    groups: dict[int, list[int]] = {}
    for i in flagged_idx:
        groups.setdefault(_find(i), []).append(i)
    group_bp = {root: sum(scored[i].length for i in members) for root, members in groups.items()}
    group_dups = {
        root: sum(dup_per_contig.get(scored[i].name, 0) for i in members)
        for root, members in groups.items()
    }
    group_of = {i: _find(i) for i in flagged_idx}

    # Expected core markers per base, used to ask whether a flagged group looks
    # like chromosomal sequence or like an extrachromosomal element.
    # The rate is estimated from markers actually located on contigs, never from
    # the lineage size: if a scan found nothing, there is no rate, and "carries no
    # markers" then says nothing about whether a group is extrachromosomal.
    marker_rate: float | None = None
    if contig_marker_counts:
        observed_total = sum(contig_marker_counts.values())
        if observed_total and total_bp:
            marker_rate = observed_total / total_bp
    group_carries_markers: dict[int, bool] = {}
    group_markers: dict[int, tuple[int, float]] = {}
    for root, members in groups.items():
        observed = sum(contig_marker_counts.get(scored[i].name, 0) for i in members) \
            if contig_marker_counts else 0
        if marker_rate is None:
            group_carries_markers[root] = False
            group_markers[root] = (observed, 0.0)
            continue
        expected = marker_rate * sum(scored[i].length for i in members)
        group_markers[root] = (observed, expected)
        group_carries_markers[root] = observed >= max(1.0, 0.5 * expected)

    verdicts: list[ContigVerdict] = []
    for i, c in enumerate(scored):
        flags: list[str] = []
        if tnf_z[i] > z_star:
            flags.append(f"composition_outlier(z={tnf_z[i]:.1f})")
        if abs(gc_z[i]) > z_star:
            flags.append(f"gc_outlier(z={gc_z[i]:.1f})")
        composition_flagged = bool(flags)
        ndup = dup_per_contig.get(c.name, 0)
        if ndup >= marker_dup_threshold:
            flags.append(f"displaced_markers(n={ndup})")
        root = group_of.get(i)
        group_mass = group_bp.get(root, 0) if root is not None else 0
        group_dup = group_dups.get(root, 0) if root is not None else 0
        if not composition_flagged and ndup >= marker_dup_threshold:
            # Marker duplication without any compositional support. Real evidence
            # at genome level, but not grounds for calling this specific contig
            # foreign -- in a multipartite genome the chromosome and the chromid
            # legitimately share paralogs of a few core markers.
            call = "marker_conflict"
        elif not composition_flagged:
            call = "core"
        elif group_dup >= marker_dup_threshold:
            call = "contaminant_candidate"
        elif marker_rate is None:
            # No marker scan was supplied, so "carries no core markers" is not a
            # fact about this group, merely an absence of evidence. Report the
            # anomaly rather than explaining it away as a plasmid.
            call = "contaminant_candidate"
        elif group_carries_markers.get(root, False):
            # Carries the single-copy core at roughly the genome-wide rate, and
            # not one of those markers is duplicated anywhere else in the
            # assembly -- so these are this organism's *only* copies of them. A
            # group holding the sole copy of part of the core is this organism's
            # own chromosome, however foreign its composition looks. A second
            # organism brings its own copies of a universal single-copy set and
            # they duplicate the host's; that is the branch above, and it is the
            # branch that has ground truth behind it.
            #
            # Calling this contamination instead cost specificity outright. The
            # FWER threshold admits a compositional outlier on some share of
            # perfectly clean assemblies by construction, and this arm turned
            # every one of them into a genome-level verdict. Measured over 320
            # shreddings of eight finished genomes: 40 false positives, 12.5%,
            # rising from 1.6% at 10 pieces to 26.6% at 1000. Zero afterwards,
            # with the outliers still flagged and still reported.
            # It also discarded the biology: a prophage or a genomic island is
            # precisely a chromosomal region with atypical composition, and the
            # comparative step drops genes sitting on suspect contigs.
            n_obs, n_exp = group_markers.get(root, (0, 0.0))
            call = "atypical_host_region"
            flags.append(
                f"sole_copy_core_markers(n={n_obs}, expected={n_exp:.1f}, displaced=0)"
            )
        elif group_mass >= replicon_min_length:
            # A compositionally coherent group of this mass, carrying no displaced
            # core markers, is as likely a plasmid or chromid as a foreign
            # organism. Mass is measured over the group, so the call survives
            # fragmentation.
            call = "replicon_candidate"
            n_frag = len(groups.get(root, []))
            flags.append(
                "distinct_replicon_or_contaminant"
                + (f"(group={n_frag} contigs, {group_mass // 1000} kb)" if n_frag > 1 else "")
            )
        else:
            call = "contaminant_candidate"

        suspicion = (
            max(0.0, float(tnf_z[i])) / z_star
            + abs(float(gc_z[i])) / z_star
            + ndup / max(1, marker_dup_threshold)
        )
        verdicts.append(
            ContigVerdict(
                name=c.name,
                length=c.length,
                gc_percent=round(100 * c.gc, 2),
                tnf_distance=round(float(tnf_distance[i]), 5),
                tnf_z=round(float(tnf_z[i]), 2),
                gc_z=round(float(gc_z[i]), 2),
                duplicated_markers=ndup,
                bin_label=int(labels[i]),
                flags=flags,
                suspicion=round(suspicion, 3),
                call=call,
            )
        )

    suspect = [v for v in verdicts if v.call == "contaminant_candidate"]
    replicons = [v for v in verdicts if v.call == "replicon_candidate"]
    host_regions = [v for v in verdicts if v.call == "atypical_host_region"]
    suspect_bp = sum(v.length for v in suspect)
    suspect_fraction = suspect_bp / total_bp if total_bp else 0.0

    n_dup_markers = len(duplicated_marker_contigs)
    cross_contig_dups = sum(
        1 for contigs in duplicated_marker_contigs.values() if len(set(contigs)) > 1
    )

    reasons: list[str] = []
    if suspect:
        reasons.append(
            f"{len(suspect)} contigs ({round(100 * suspect_fraction, 2)}% of assembly) "
            "deviate from core composition"
        )
    if replicons:
        reasons.append(
            f"{len(replicons)} large contig(s) ({round(100 * sum(v.length for v in replicons) / total_bp, 2)}% "
            "of assembly) have distinct composition but carry no displaced core markers -- "
            "consistent with a plasmid or a second chromosome, not resolvable by composition alone"
        )
    if host_regions:
        reasons.append(
            f"{len(host_regions)} contig(s) "
            f"({round(100 * sum(v.length for v in host_regions) / total_bp, 2)}% of assembly) "
            "have distinct composition but hold this assembly's only copies of core "
            "single-copy markers -- this organism's own chromosome, so a prophage or "
            "acquired island rather than contamination"
        )
    if cross_contig_dups:
        share = (
            f" ({round(100 * cross_contig_dups / total_markers, 2)}% of the "
            f"{total_markers}-marker set)" if total_markers else ""
        )
        reasons.append(
            f"{cross_contig_dups} single-copy markers duplicated across "
            f"different contigs{share}"
        )
    if bins["bimodal"]:
        reasons.append(
            f"composition is bimodal: minor bin holds {bins['minor_bin_fraction_percent']}% "
            f"of bases at a {bins['gc_gap_percent']} point GC offset"
        )
    if not use_zscores:
        reasons.append(
            f"only {len(scored)} contigs scored -- too few for outlier statistics, so "
            "absolute composition thresholds were used and this verdict is weak evidence"
        )

    # The bimodality field no longer votes. A 2-means split always returns two
    # bins, so on any fragmented assembly some split cleared the old thresholds;
    # it fired on 57 of 72 published genomes CheckM2 calls clean. It is retained
    # as reported evidence, not as a trigger.
    # Duplication is judged as a rate when the marker-set size is known, and as a
    # raw count when it is not. A caller that cannot say how many markers were
    # searched has not supplied enough information to form a rate, and inventing
    # a denominator would be worse than falling back to the older, set-specific
    # behaviour -- which the fallback reproduces exactly on bacteria_odb10.
    if total_markers:
        dup_rate = cross_contig_dups / total_markers
        # `>=`, so that on a 124-marker set the rule is bit-for-bit the old one:
        # 5/124 >= 5/124 must hold, or the change of units quietly moves the line.
        dup_likely = dup_rate >= DUP_RATE_LIKELY
        dup_possible = dup_rate >= DUP_RATE_POSSIBLE
    else:
        dup_rate = None
        dup_likely = cross_contig_dups >= 5
        dup_possible = cross_contig_dups >= 2

    if suspect_fraction > 0.05 or dup_likely:
        verdict = "likely"
    elif suspect_fraction > 0.01 or dup_possible:
        verdict = "possible"
    else:
        verdict = "clean"
        reasons = reasons or ["no compositional outliers and no cross-contig marker duplication"]

    return ContaminationResult(
        verdict=verdict,
        reasons=reasons,
        suspect_bp=suspect_bp,
        suspect_fraction=suspect_fraction,
        n_suspect_contigs=len(suspect),
        contigs=verdicts,
        bins=bins,
        params={
            "min_contig_length": min_contig_length,
            "zscore_statistics_used": use_zscores,
            "alpha": alpha,
            "fwer_z_threshold": round(z_star, 3),
            "null_curve": curve.summary(),
            "flagged_groups": len(groups),
            "marker_rate_per_mb": round(1e6 * marker_rate, 2) if marker_rate else None,
            "refit_applied": refit_applied,
            "refit_excluded_contigs": len(outliers) if refit_applied else 0,
            "tnf_z_threshold": tnf_z_threshold,
            "gc_z_threshold": gc_z_threshold,
            "tnf_distance_threshold": tnf_distance_threshold,
            "gc_absolute_threshold": gc_absolute_threshold,
            "marker_dup_threshold": marker_dup_threshold,
            "replicon_candidates": [v.name for v in replicons],
            "atypical_host_regions": [v.name for v in host_regions],
            "contigs_scored": len(scored),
            "contigs_skipped_short": len(asm.contigs) - len(scored),
            "duplicated_markers_total": n_dup_markers,
            "duplicated_markers_cross_contig": cross_contig_dups,
            "marker_set_size": total_markers,
            # Share of the single-copy core sitting on more than one contig.
            # Rises with the foreign fraction and is comparable across marker
            # sets, which is what the verdict needs. It is emphatically not a
            # contamination percentage: on `bench/mixture_ladder.py` a 5%
            # mixture produced 4.8-11.6% here, and clean finished genomes sit
            # near 1% from ordinary cross-replicon paralogy alone. Report the
            # verdict, not this number.
            "cross_contig_duplication_percent": (
                round(100 * dup_rate, 2) if dup_rate is not None else None
            ),
            "duplication_rate_thresholds_percent": [
                round(100 * DUP_RATE_POSSIBLE, 2), round(100 * DUP_RATE_LIKELY, 2)
            ],
        },
    )


def write_contig_table(result: ContaminationResult, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = [
        "contig", "length", "gc_percent", "tnf_distance", "tnf_z", "gc_z",
        "displaced_markers", "bin", "call", "suspicion", "flags",
    ]
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(header) + "\n")
        for c in sorted(result.contigs, key=lambda x: -x.suspicion):
            fh.write(
                "\t".join(
                    [
                        c.name, str(c.length), f"{c.gc_percent}", f"{c.tnf_distance}",
                        f"{c.tnf_z}", f"{c.gc_z}", str(c.duplicated_markers),
                        str(c.bin_label), c.call, f"{c.suspicion}", ";".join(c.flags) or "-",
                    ]
                )
                + "\n"
            )
    return path
