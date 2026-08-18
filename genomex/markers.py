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


def _union_length(intervals: list[tuple[int, int]]) -> int:
    """Total aa covered by a set of (from, to) domain envelopes, overlaps merged."""
    if not intervals:
        return 0
    merged: list[list[int]] = []
    for lo, hi in sorted(intervals):
        if merged and lo <= merged[-1][1] + 1:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return sum(hi - lo + 1 for lo, hi in merged)


def parse_domtbl(path: str | Path, lineage: Lineage, gene_contig: dict[str, str]) -> MarkerResult:
    """Apply BUSCO's score and length cutoffs to a HMMER --domtblout table."""
    # (busco_id, protein_id) -> (full_sequence_score, [(env_from, env_to), ...])
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
            env_from, env_to = int(f[19]), int(f[20])
            key = (busco_id, protein_id)
            if key not in acc:
                acc[key] = (full_score, [])
            acc[key][1].append((env_from, env_to))

    hits: list[MarkerHit] = []
    for (busco_id, protein_id), (full_score, envs) in acc.items():
        cutoff = lineage.scores.get(busco_id)
        if cutoff is None or full_score < cutoff:
            continue
        sigma, mean_len = lineage.lengths.get(busco_id, (1.0, 0.0))
        matched = _union_length(envs)
        status = "complete" if matched > (mean_len - 2 * sigma) else "fragmented"
        hits.append(
            MarkerHit(
                busco_id=busco_id,
                protein_id=protein_id,
                contig=gene_contig.get(protein_id, protein_id.rsplit("_", 1)[0]),
                score=full_score,
                matched_aa=matched,
                status=status,
            )
        )

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
