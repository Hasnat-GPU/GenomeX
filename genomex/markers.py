"""Single-copy ortholog recovery -- completeness and, crucially, duplication.

This is a BUSCO-compatible scan, not BUSCO itself: it runs HMMER directly
against the OrthoDB v10 HMM profiles that BUSCO ships, and applies BUSCO's own
per-marker score and length cutoffs. Numbers land close to BUSCO's but are not
guaranteed identical -- BUSCO additionally re-predicts genes per candidate
region (and uses metaeuk/miniprot for eukaryotes). Reports must say so.

Two outputs matter downstream:
  completeness  -- how much of the expected single-copy core is present
  duplication   -- markers found more than once, and *which contigs* carry them.
                   Widespread duplication spread across unrelated contigs is the
                   canonical signature of a mixed (contaminated) assembly.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .runtime import Runtime


@dataclass
class Lineage:
    path: Path
    name: str
    n_markers: int
    scores: dict[str, float]
    lengths: dict[str, tuple[float, float]]  # busco_id -> (sigma, mean_length)

    @classmethod
    def load(cls, path: str | Path) -> "Lineage":
        path = Path(path)
        hmm_dir = path / "hmms"
        if not hmm_dir.is_dir():
            raise FileNotFoundError(
                f"{path} does not look like a BUSCO lineage directory (no hmms/). "
                "Download one from https://busco-data.ezlab.org/v5/data/lineages/"
            )
        scores: dict[str, float] = {}
        for line in (path / "scores_cutoff").read_text().splitlines():
            if line.strip():
                bid, score = line.split()[:2]
                scores[bid] = float(score)
        lengths: dict[str, tuple[float, float]] = {}
        for line in (path / "lengths_cutoff").read_text().splitlines():
            if line.strip():
                parts = line.split()
                # columns: busco_id, 0, sigma, mean_length
                bid, sigma, mean_len = parts[0], float(parts[2]), float(parts[3])
                lengths[bid] = (sigma if sigma > 0 else 1.0, mean_len)
        cfg: dict[str, str] = {}
        cfg_path = path / "dataset.cfg"
        if cfg_path.exists():
            for line in cfg_path.read_text().splitlines():
                if "=" in line:
                    k, v = line.split("=", 1)
                    cfg[k.strip()] = v.strip()
        return cls(
            path=path,
            name=cfg.get("name", path.name),
            n_markers=int(cfg.get("number_of_BUSCOs", len(scores))),
            scores=scores,
            lengths=lengths,
        )

    def combined_hmm(self, workdir: str | Path) -> Path:
        """Concatenate the per-marker HMMs into one searchable profile database."""
        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)
        combined = workdir / f"{self.name}.all.hmm"
        if combined.exists() and combined.stat().st_size > 0:
            return combined
        with open(combined, "wb") as out:
            for hmm in sorted((self.path / "hmms").glob("*.hmm")):
                out.write(hmm.read_bytes())
        return combined


@dataclass
class MarkerHit:
    busco_id: str
    protein_id: str
    contig: str
    score: float
    matched_aa: int
    status: str  # complete | fragmented


@dataclass
class MarkerResult:
    lineage: str
    n_markers: int
    single: list[str] = field(default_factory=list)
    duplicated: list[str] = field(default_factory=list)
    fragmented: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    hits: list[MarkerHit] = field(default_factory=list)

    @property
    def complete(self) -> int:
        return len(self.single) + len(self.duplicated)

    @property
    def completeness(self) -> float:
        return round(100.0 * self.complete / self.n_markers, 2) if self.n_markers else 0.0

    @property
    def duplication_percent(self) -> float:
        return round(100.0 * len(self.duplicated) / self.n_markers, 2) if self.n_markers else 0.0

    def contig_marker_counts(self) -> dict[str, int]:
        counts: dict[str, int] = defaultdict(int)
        for h in self.hits:
            if h.status == "complete":
                counts[h.contig] += 1
        return dict(counts)

    def duplicated_marker_contigs(self) -> dict[str, list[str]]:
        """busco_id -> contigs carrying each complete copy, for duplicated markers."""
        dup = set(self.duplicated)
        per: dict[str, list[str]] = defaultdict(list)
        for h in self.hits:
            if h.status == "complete" and h.busco_id in dup:
                per[h.busco_id].append(h.contig)
        return dict(per)

    def summary(self) -> dict:
        n = self.n_markers or 1
        return {
            "lineage": self.lineage,
            "markers_total": self.n_markers,
            "complete": self.complete,
            "single_copy": len(self.single),
            "duplicated": len(self.duplicated),
            "fragmented": len(self.fragmented),
            "missing": len(self.missing),
            "completeness_percent": self.completeness,
            "duplication_percent": self.duplication_percent,
            "busco_style_string": (
                "C:{c}%[S:{s}%,D:{d}%],F:{f}%,M:{m}%,n:{n}".format(
                    c=self.completeness,
                    s=round(100.0 * len(self.single) / n, 2),
                    d=self.duplication_percent,
                    f=round(100.0 * len(self.fragmented) / n, 2),
                    m=round(100.0 * len(self.missing) / n, 2),
                    n=self.n_markers,
                )
            ),
            "method": "HMMER hmmsearch vs BUSCO odb10 profiles (BUSCO-compatible, not BUSCO)",
        }


#: BUSCO keeps only matches scoring within this fraction of a marker's best hit.
#: Without it, every weak hit above a permissive score cutoff counts as another
#: copy and the duplication rate -- which the contamination module consumes --
#: is badly inflated. (BUSCO 5 hmmer.py, _remove_low_scoring_matches.)
BITSCORE_RETENTION = 0.85


def sum_hmm_len(intervals: list[tuple[int, int]]) -> int:
    """Profile positions covered by a set of (hmm_from, hmm_to) domain coordinates.

    Measured on the *HMM* axis, not the sequence envelope: how much of the
    profile the protein covers is what the length cutoffs are calibrated against,
    and an envelope can extend well past the modelled region.

    The merge follows BUSCO's own implementation, which absorbs a region only
    when its start falls inside an existing one -- deliberately reproduced,
    including its quirks, so the two tools count the same way.
    """
    used: list[list[int]] = []
    for region in sorted(intervals, key=lambda x: x[0]):
        for u in used:
            if u[0] <= region[0] <= u[1]:
                if region[1] > u[1]:
                    u[1] = region[1]
                break
        else:
            used.append(list(region))
    return sum(hi - lo + 1 for lo, hi in used)


def _union_length(intervals: list[tuple[int, int]]) -> int:
    """Deprecated alias for :func:`sum_hmm_len`."""
    return sum_hmm_len(intervals)


def _resolve_cross_marker_claims(scored: list[MarkerHit]) -> list[MarkerHit]:
    """A protein belongs to one marker: the one it scores highest against.

    Two of BUSCO's rules, which act across markers rather than within one.
    (BUSCO 5 hmmer.py, `_remove_lower_ranked_duplicates` and
    `_remove_remaining_duplicate_matches`, both run by `_remove_duplicates`
    before the 85%-retention step -- hence the order here.)

    * **Rank first.** A protein good enough to be a complete copy of something is
      not also a fragment of something else, so its fragmented claims go.
    * **Then best score.** Among the claims left at its rank, the protein is kept
      only under the marker it scores highest against.

    Large paralogous families are what this is for. In `fungi_odb10`, three
    AAA-ATPases each clear the score cutoff of the Pex1 marker *and* of the
    AAA-ATPase marker, scoring higher against the latter. Without this rule they
    count as extra copies of Pex1 and the marker is reported duplicated when it
    is single -- measured on S. cerevisiae, the only marker of 758 on which
    GenomeX and BUSCO 5.8.3 disagreed.

    It is nearly inert on bacteria, which is why the 496/496 bacterial agreement
    held without it: re-deriving all 72 genomes changes 1 marker of 9,002 under
    `bacteria_odb10` and 5 of 49,963 under `burkholderiales_odb10`. Bacterial
    core markers are largely ribosomal and rarely share a protein between
    families the way eukaryotic ATPase families do.

    Ties are broken on `busco_id` so the result does not depend on dict order.
    """
    by_protein: dict[str, list[MarkerHit]] = defaultdict(list)
    for h in scored:
        by_protein[h.protein_id].append(h)

    keep: list[MarkerHit] = []
    for _protein, claims in by_protein.items():
        complete = [h for h in claims if h.status == "complete"]
        ranked = complete or claims
        winner = max(ranked, key=lambda h: (h.score, h.busco_id))
        keep.append(winner)
    return keep


def parse_domtbl(
    path: str | Path,
    lineage: Lineage,
    gene_contig: dict[str, str],
    *,
    bitscore_retention: float = BITSCORE_RETENTION,
) -> MarkerResult:
    """Apply BUSCO's score, length and bitscore-retention rules to a domtblout table.

    Three steps, in BUSCO's order:

    1. drop hits below the marker's score cutoff;
    2. classify each surviving hit complete or fragmented from its HMM-profile
       coverage against the marker's length cutoff (BUSCO's odb10 rule is
       ``zeta = (length - size) / sigma``, fragmented when ``zeta > 2``);
    3. resolve claims across markers, so each protein belongs to exactly one
       (see :func:`_resolve_cross_marker_claims`);
    4. per marker, drop hits scoring below ``bitscore_retention`` of that
       marker's best hit -- only then count copies.

    Step 4 is what separates a genuine second copy from a distant paralog that
    happened to clear a permissive cutoff. Without it a marker whose best hit
    scores 297 also counts hits at 17 and 19 as copies. Step 3 catches the other
    direction: a protein that is really *another* marker's, counted here as a
    second copy because it cleared this marker's cutoff too.
    """
    # (busco_id, protein_id) -> (full_sequence_score, [(hmm_from, hmm_to), ...])
    acc: dict[tuple[str, str], tuple[float, list[tuple[int, int]]]] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 21:
                continue
            protein_id, busco_id = f[0], f[3]
            full_score = float(f[7])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            key = (busco_id, protein_id)
            if key not in acc:
                acc[key] = (full_score, [])
            acc[key][1].append((hmm_from, hmm_to))

    scored: list[MarkerHit] = []
    for (busco_id, protein_id), (full_score, coords) in acc.items():
        cutoff = lineage.scores.get(busco_id)
        if cutoff is None or full_score < cutoff:
            continue
        sigma, mean_len = lineage.lengths.get(busco_id, (1.0, 0.0))
        matched = sum_hmm_len(coords)
        zeta = (mean_len - matched) / sigma if sigma else 0.0
        status = "fragmented" if zeta > 2 else "complete"
        scored.append(
            MarkerHit(
                busco_id=busco_id,
                protein_id=protein_id,
                contig=gene_contig.get(protein_id, protein_id.rsplit("_", 1)[0]),
                score=full_score,
                matched_aa=matched,
                status=status,
            )
        )

    # One protein may only be claimed by one marker. BUSCO does this *before*
    # the retention step, so the order here is its order.
    scored = _resolve_cross_marker_claims(scored)

    # Retain only matches within `bitscore_retention` of each marker's best hit.
    best: dict[str, float] = {}
    for h in scored:
        best[h.busco_id] = max(best.get(h.busco_id, 0.0), h.score)
    hits = [h for h in scored if h.score >= bitscore_retention * best[h.busco_id]]

    by_marker: dict[str, list[MarkerHit]] = defaultdict(list)
    for h in hits:
        by_marker[h.busco_id].append(h)

    res = MarkerResult(lineage=lineage.name, n_markers=lineage.n_markers, hits=hits)
    for busco_id in lineage.scores:
        marker_hits = by_marker.get(busco_id, [])
        n_complete = sum(1 for h in marker_hits if h.status == "complete")
        if n_complete == 1:
            res.single.append(busco_id)
        elif n_complete > 1:
            res.duplicated.append(busco_id)
        elif any(h.status == "fragmented" for h in marker_hits):
            res.fragmented.append(busco_id)
        else:
            res.missing.append(busco_id)
    return res


MARKER_TABLE_COLUMNS = ["busco_id", "status", "protein_id", "contig", "score", "matched_aa"]


def write_marker_table(result: MarkerResult, path: str | Path) -> Path:
    """Per-marker classification, one row per marker.

    Aggregate percentages hide disagreement: two tools can both report 100%
    complete while differing on which markers are duplicated. This table is what
    makes a marker-by-marker comparison against BUSCO possible.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    status_of = {}
    for bid in result.single:
        status_of[bid] = "Complete"
    for bid in result.duplicated:
        status_of[bid] = "Duplicated"
    for bid in result.fragmented:
        status_of[bid] = "Fragmented"
    for bid in result.missing:
        status_of[bid] = "Missing"

    hits_by_marker: dict[str, list[MarkerHit]] = defaultdict(list)
    for h in result.hits:
        hits_by_marker[h.busco_id].append(h)

    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(MARKER_TABLE_COLUMNS) + "\n")
        for bid in sorted(status_of):
            status = status_of[bid]
            hits = sorted(hits_by_marker.get(bid, []), key=lambda h: -h.score)
            if not hits:
                fh.write(f"{bid}\t{status}\t-\t-\t-\t-\n")
                continue
            for h in hits:
                fh.write(
                    f"{bid}\t{status}\t{h.protein_id}\t{h.contig}\t{h.score:.1f}\t{h.matched_aa}\n"
                )
    return path


def scan_markers(
    proteins_faa: str | Path,
    lineage: Lineage,
    outdir: str | Path,
    rt: Runtime,
    *,
    gene_contig: dict[str, str] | None = None,
    force: bool = False,
) -> MarkerResult:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    domtbl = outdir / "markers.domtbl"
    combined = lineage.combined_hmm(outdir.parent / "_hmmdb")
    if force or not domtbl.exists() or domtbl.stat().st_size == 0:
        rt.run(
            [
                "hmmsearch",
                "--domtblout", str(domtbl),
                "--cpu", str(rt.threads),
                "-o", str(outdir / "markers.hmmer.txt"),
                str(combined),
                str(proteins_faa),
            ]
        )
    return parse_domtbl(domtbl, lineage, gene_contig or {})
