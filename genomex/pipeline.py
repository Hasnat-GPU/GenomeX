"""The unified pipeline: assemblies in, one evidence-linked report out.

Per genome:   stats -> gene calls -> single-copy markers -> contamination
Per pair:     ANI -> orthogroups -> strain-unique genes with causes attached

The link between the two halves is deliberate.  Contamination analysis produces
a set of suspect contigs; the comparative analysis consumes it, so any gene-content
difference that traces back to a suspect contig is labelled as such instead of
being reported as biology.
"""

from __future__ import annotations

import itertools
import json
import time
from dataclasses import dataclass, field
from pathlib import Path

from .compare import (
    PairComparison,
    Pangenome,
    cluster_proteomes,
    compute_ani,
    explain_unique_genes,
    island_enrichment,
    write_unique_gene_table,
)
from .contamination import ContaminationResult, detect_contamination, write_contig_table
from .fasta import Assembly, assembly_stats
from .genes import GeneCalls, predict_genes
from .markers import Lineage, MarkerResult, scan_markers, write_marker_table
from .runtime import Runtime

DEFAULT_LINEAGE = Path.home() / "genomex-work" / "db" / "bacteria_odb10"


@dataclass
class GenomeResult:
    name: str
    path: Path
    assembly: Assembly
    stats: dict
    calls: GeneCalls
    markers: MarkerResult
    contamination: ContaminationResult
    outdir: Path
    seconds: float = 0.0
    sha256: str = ""

    def summary(self) -> dict:
        return {
            "genome": self.name,
            "input": str(self.path),
            "sha256": self.sha256,
            "seconds": round(self.seconds, 1),
            "assembly": self.stats,
            "genes": self.calls.summary(self.stats.get("total_bp", 0)),
            "markers": self.markers.summary(),
            "contamination": self.contamination.summary(),
            "quality_call": self.quality_call(),
        }

    def quality_call(self) -> dict:
        """A single, arguable, explicitly-justified verdict per genome."""
        m = self.markers
        c = self.contamination
        reasons = []
        if m.completeness >= 95:
            completeness_grade = "high"
        elif m.completeness >= 90:
            completeness_grade = "acceptable"
            reasons.append(f"completeness {m.completeness}% is below the 95% comfort line")
        else:
            completeness_grade = "low"
            reasons.append(f"completeness {m.completeness}% -- assembly is missing core genes")

        if c.verdict in ("likely",):
            reasons.append("composition indicates more than one organism")
        elif c.verdict == "possible":
            reasons.append("composition shows minor foreign signal")

        usable = completeness_grade != "low" and c.verdict != "likely"
        return {
            "usable_for_comparative_analysis": usable,
            "completeness_grade": completeness_grade,
            "contamination_verdict": c.verdict,
            "reasons": reasons or ["complete single-copy core, single-organism composition"],
        }


def analyze_genome(
    path: str | Path,
    outdir: str | Path,
    rt: Runtime,
    lineage: Lineage,
    *,
    genetic_code: int = 11,
    min_contig_length: int = 3000,
    keep_sequences: bool = False,
) -> GenomeResult:
    t0 = time.time()
    path = Path(path)
    asm = Assembly.load(path)
    stats = assembly_stats(asm)
    outdir = Path(outdir) / asm.name
    outdir.mkdir(parents=True, exist_ok=True)

    calls = predict_genes(
        path, outdir, rt, genetic_code=genetic_code, total_bp=stats.get("total_bp", 0)
    )
    gene_contig = {g.protein_id: g.contig for g in calls.genes}
    markers = scan_markers(calls.proteins_faa, lineage, outdir, rt, gene_contig=gene_contig)
    write_marker_table(markers, outdir / "markers.tsv")

    contamination = detect_contamination(
        asm,
        duplicated_marker_contigs=markers.duplicated_marker_contigs(),
        contig_marker_counts=markers.contig_marker_counts(),
        total_markers=markers.n_markers,
        min_contig_length=min_contig_length,
    )
    write_contig_table(contamination, outdir / "contigs.tsv")

    if not keep_sequences:
        for c in asm.contigs:
            c.seq = ""  # free ~10 MB per genome; stats and TNF are already computed

    return GenomeResult(
        name=asm.name,
        path=path,
        assembly=asm,
        stats=stats,
        calls=calls,
        markers=markers,
        contamination=contamination,
        outdir=outdir,
        seconds=time.time() - t0,
        sha256=asm.sha256(),
    )


def compare_pair(
    a: GenomeResult,
    b: GenomeResult,
    outdir: str | Path,
    rt: Runtime,
    *,
    pangenome: Pangenome | None = None,
    min_seq_id: float = 0.5,
    coverage: float = 0.8,
) -> PairComparison:
    outdir = Path(outdir) / f"{a.name}__vs__{b.name}"
    outdir.mkdir(parents=True, exist_ok=True)

    ani = compute_ani(a.path, b.path, outdir, rt)

    if pangenome is None:
        clusters = cluster_proteomes(
            {a.name: a.calls.proteins_faa, b.name: b.calls.proteins_faa},
            outdir,
            rt,
            min_seq_id=min_seq_id,
            coverage=coverage,
        )
        pangenome = Pangenome([a.name, b.name], clusters, min_seq_id, coverage)

    # Pairwise membership, evaluated only on these two genomes. When a joint
    # pangenome spans more genomes, a cluster shared by A and C but absent from B
    # is still "unique to A" *for this pair* -- comparing against the whole set
    # would silently answer a different question.
    presence = pangenome.presence()
    shared_reps = [rep for rep, gs in presence.items() if {a.name, b.name} <= gs]
    only_a = [rep for rep, gs in presence.items() if a.name in gs and b.name not in gs]
    only_b = [rep for rep, gs in presence.items() if b.name in gs and a.name not in gs]

    unique_a = explain_unique_genes(
        pangenome.genes_of(a.name, only_a),
        a.calls.genes,
        genes_fna=a.calls.genes_fna,
        suspect_contigs=a.contamination.suspect_contig_names(),
    )
    unique_b = explain_unique_genes(
        pangenome.genes_of(b.name, only_b),
        b.calls.genes,
        genes_fna=b.calls.genes_fna,
        suspect_contigs=b.contamination.suspect_contig_names(),
    )
    write_unique_gene_table(unique_a, outdir / f"unique_to_{a.name}.tsv")
    write_unique_gene_table(unique_b, outdir / f"unique_to_{b.name}.tsv")

    return PairComparison(
        genome_a=a.name,
        genome_b=b.name,
        ani=ani,
        shared_orthogroups=len(shared_reps),
        unique_a=unique_a,
        unique_b=unique_b,
        n_genes_a=a.calls.n_genes,
        n_genes_b=b.calls.n_genes,
        island_stats_a=island_enrichment(
            pangenome.genes_of(a.name, only_a), a.calls.genes
        ),
        island_stats_b=island_enrichment(
            pangenome.genes_of(b.name, only_b), b.calls.genes
        ),
    )


@dataclass
class PipelineResult:
    genomes: list[GenomeResult]
    pairs: list[PairComparison] = field(default_factory=list)
    pangenome: Pangenome | None = None
    runtime: Runtime | None = None
    outdir: Path = Path(".")
    seconds: float = 0.0

    def to_dict(self) -> dict:
        return {
            "genomex_version": __import__("genomex").__version__,
            "seconds": round(self.seconds, 1),
            "genomes": [g.summary() for g in self.genomes],
            "pangenome": self.pangenome.summary() if self.pangenome else None,
            "pairs": [p.summary() for p in self.pairs],
            "provenance": self.runtime.provenance() if self.runtime else {},
        }

    def write_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path


def run_pipeline(
    inputs: list[str | Path],
    outdir: str | Path,
    *,
    lineage_path: str | Path = DEFAULT_LINEAGE,
    pairs: list[tuple[str, str]] | None = None,
    all_pairs: bool = False,
    pangenome: bool = True,
    genetic_code: int = 11,
    min_contig_length: int = 3000,
    min_seq_id: float = 0.5,
    coverage: float = 0.8,
    threads: int | None = None,
    log=print,
) -> PipelineResult:
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rt = Runtime()
    if threads:
        rt.threads = threads
    rt.require("prodigal", "hmmsearch")

    lineage = Lineage.load(lineage_path)
    log(f"lineage: {lineage.name} ({lineage.n_markers} markers)")

    results: list[GenomeResult] = []
    for path in inputs:
        log(f"[{len(results) + 1}/{len(inputs)}] {Path(path).name}")
        r = analyze_genome(
            path, outdir / "genomes", rt, lineage,
            genetic_code=genetic_code, min_contig_length=min_contig_length,
        )
        log(
            f"    {r.stats['n_contigs']} contigs, {r.stats['total_bp'] / 1e6:.2f} Mb, "
            f"{r.calls.n_genes} genes | completeness {r.markers.completeness}% "
            f"dup {r.markers.duplication_percent}% | contamination: {r.contamination.verdict} "
            f"({r.seconds:.0f}s)"
        )
        results.append(r)

    by_name = {r.name: r for r in results}
    pg: Pangenome | None = None
    if pangenome and len(results) >= 2:
        rt.require("mmseqs")
        log(f"clustering proteomes of {len(results)} genomes (MMseqs2)")
        clusters = cluster_proteomes(
            {r.name: r.calls.proteins_faa for r in results},
            outdir / "pangenome", rt, min_seq_id=min_seq_id, coverage=coverage,
        )
        pg = Pangenome([r.name for r in results], clusters, min_seq_id, coverage)
        s = pg.summary()
        log(
            f"    {s['orthogroups_total']} orthogroups | core {s['core']} "
            f"accessory {s['accessory']} strain-unique {s['strain_unique']}"
        )

    todo: list[tuple[str, str]] = []
    if all_pairs:
        todo = [(a.name, b.name) for a, b in itertools.combinations(results, 2)]
    elif pairs:
        todo = list(pairs)
    elif len(results) == 2:
        todo = [(results[0].name, results[1].name)]

    comparisons: list[PairComparison] = []
    if todo:
        rt.require("fastANI", "mmseqs")
    for a_name, b_name in todo:
        if a_name not in by_name or b_name not in by_name:
            log(f"    skipping pair {a_name} vs {b_name}: genome not in this run")
            continue
        log(f"comparing {a_name} vs {b_name}")
        cmp_ = compare_pair(
            by_name[a_name], by_name[b_name], outdir / "pairs", rt,
            pangenome=pg, min_seq_id=min_seq_id, coverage=coverage,
        )
        ani_txt = "n/a (>~20% divergent)" if cmp_.ani.ani is None else f"{cmp_.ani.ani}%"
        log(
            f"    ANI {ani_txt} | shared {cmp_.shared_orthogroups} | "
            f"unique {len(cmp_.unique_a)}/{len(cmp_.unique_b)}"
        )
        comparisons.append(cmp_)

    return PipelineResult(
        genomes=results, pairs=comparisons, pangenome=pg, runtime=rt,
        outdir=outdir, seconds=time.time() - t0,
    )
