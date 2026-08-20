"""The unified pipeline: assemblies in, one evidence-linked report out.

Per genome:   stats -> gene calls -> single-copy markers -> contamination
Per pair:     ANI -> orthogroups -> strain-unique genes with causes attached
Per proteome: single-copy markers only -- see `analyze_proteome`, which exists
              to make the missing two-thirds explicit rather than empty.

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
from .markers import (
    DEFAULT_LINEAGE_NAME,
    Lineage,
    MarkerResult,
    available_lineages,
    lineage_db,
    resolve_lineage,
    scan_markers,
    write_marker_table,
)
from .proteome import NOT_MEASURED, Proteome, assert_nucleotide
from .runtime import Runtime

DEFAULT_LINEAGE = lineage_db() / DEFAULT_LINEAGE_NAME


def lineage_notice(lineage: Lineage) -> list[str]:
    """The scan header, plus what the default set costs when it is the one in use.

    The default is universal by construction, which is why it is the default and
    also why it is the weakest instrument in the box: contamination recall is a
    function of how many single-copy markers a foreign genome can displace. This
    is said at run time rather than only in the docs because a user who never
    learns `--lineage` exists still gets a verdict, and it looks the same as a
    good one.

    It states the cost; it does not act on it. Choosing a set means knowing the
    organism, and a set chosen wrongly reports fake incompleteness.
    """
    lines = [f"lineage: {lineage.name} ({lineage.n_markers} markers)"]
    if DEFAULT_LINEAGE_NAME not in (lineage.name, lineage.path.name):
        return lines

    others = [i for i in available_lineages() if i.path != lineage.path]
    lines.append(
        "    universal set -- detection scales with it: on constructed 2% mixtures a"
    )
    lines.append(
        "    lineage-specific set finds 4 of 4 donor pairs and this one finds 0 of 4."
    )
    lines.append(
        f"    {len(others)} other set(s) installed; run `genomex lineages`."
        if others else
        "    no other sets installed; run `genomex lineages` for where to get them."
    )
    return lines


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

        # "undetermined" is not "clean". A finished single-replicon genome has
        # too few scorable contigs for composition statistics and abstains here;
        # before this branch existed it fell through to `usable: True` with the
        # reason "single-organism composition" -- a positive claim about a
        # measurement that never ran. None of the 72 benchmark genomes reach it,
        # which is exactly why it survived.
        assessed = c.assessed
        if c.verdict == "likely":
            reasons.append("composition indicates more than one organism")
        elif c.verdict == "possible":
            reasons.append("composition shows minor foreign signal")
        elif not assessed:
            why = c.reasons[0] if c.reasons else "no reason recorded"
            reasons.append(f"composition was not assessed: {why}")

        if completeness_grade == "low":
            usable = False
        elif not assessed:
            usable = None  # not False: nothing was found wrong, nothing was checked
        else:
            usable = c.verdict != "likely"

        if not reasons:
            reasons.append("complete single-copy core, single-organism composition")
        return {
            "usable_for_comparative_analysis": usable,
            "completeness_grade": completeness_grade,
            "contamination_verdict": c.verdict,
            "reasons": reasons,
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
    # Before anything measures this file. A proteome reaching assembly_stats is
    # not an error downstream, it is a plausible wrong answer -- see NotAnAssembly.
    assert_nucleotide((c.seq for c in asm.contigs), path)
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


@dataclass
class ProteomeResult:
    """What a protein FASTA can support: completeness and duplication.

    Deliberately not a `GenomeResult` with empty fields. The two carry different
    guarantees, and a type that can be either invites code that forgets which one
    it is holding.
    """

    name: str
    path: Path
    proteome: Proteome
    stats: dict
    markers: MarkerResult
    outdir: Path
    seconds: float = 0.0
    sha256: str = ""

    def summary(self) -> dict:
        return {
            "name": self.name,
            "input": str(self.path),
            "input_type": "proteome",
            "sha256": self.sha256,
            "seconds": round(self.seconds, 1),
            "proteins": self.stats,
            "markers": self.markers.summary(),
            "not_measured": NOT_MEASURED,
            "quality_call": self.quality_call(),
        }

    def quality_call(self) -> dict:
        """Completeness only, and scoped to what it actually measures.

        On an assembly, low completeness means the assembly is missing core
        genes. On a supplied proteome it means the *proteome* is -- which may be
        the assembly, or may be the gene caller that produced the file. GenomeX
        did not call these genes and cannot tell the two apart, so it says both.
        """
        m = self.markers
        reasons: list[str] = []
        if m.completeness >= 95:
            grade = "high"
        elif m.completeness >= 90:
            grade = "acceptable"
            reasons.append(f"completeness {m.completeness}% is below the 95% comfort line")
        else:
            grade = "low"
            reasons.append(
                f"completeness {m.completeness}% -- this proteome is missing core genes. "
                "GenomeX did not call these genes, so it cannot say whether the assembly "
                "lacks them or the gene caller missed them"
            )
        if m.duplication_percent >= 5:
            reasons.append(
                f"{m.duplication_percent}% of markers appear more than once. On an "
                "assembly this is the contamination signal; on a proteome it is equally "
                "consistent with redundant gene models, and there are no contigs to tell "
                "them apart"
            )
        return {
            # Not False: False would read as "we checked and it is unsuitable".
            "usable_for_comparative_analysis": None,
            "completeness_grade": grade,
            "contamination_verdict": "not_measured",
            "reasons": reasons or ["complete single-copy core for this lineage"],
            "scope": "completeness and duplication only -- see not_measured",
        }


def analyze_proteome(
    path: str | Path,
    outdir: str | Path,
    rt: Runtime,
    lineage: Lineage,
) -> ProteomeResult:
    t0 = time.time()
    path = Path(path)
    prot = Proteome.load(path)
    outdir = Path(outdir) / prot.name
    outdir.mkdir(parents=True, exist_ok=True)

    # No contig map, deliberately: `parse_domtbl` then writes the unknown-contig
    # sentinel and sets `contigs_known=False`, which is what makes every
    # contig-shaped accessor downstream refuse instead of guessing.
    markers = scan_markers(path, lineage, outdir, rt)
    write_marker_table(markers, outdir / "markers.tsv")

    return ProteomeResult(
        name=prot.name,
        path=path,
        proteome=prot,
        stats=prot.stats(),
        markers=markers,
        outdir=outdir,
        seconds=time.time() - t0,
        sha256=prot.sha256(),
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
    #: Proteome inputs are kept in their own list rather than mixed into
    #: `genomes`. Every existing consumer of `genomes` -- the report, the
    #: comparative step, the benchmark harnesses -- assumes an assembly behind
    #: each entry. A proteome entry with `assembly` stubbed to zeroes would
    #: satisfy that assumption syntactically and lie to it semantically.
    proteomes: list["ProteomeResult"] = field(default_factory=list)

    def to_dict(self) -> dict:
        return {
            "genomex_version": __import__("genomex").__version__,
            "seconds": round(self.seconds, 1),
            "genomes": [g.summary() for g in self.genomes],
            "proteomes": [p.summary() for p in self.proteomes],
            "pangenome": self.pangenome.summary() if self.pangenome else None,
            "pairs": [p.summary() for p in self.pairs],
            "provenance": self.runtime.provenance() if self.runtime else {},
        }

    def write_json(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")
        return path


def run_proteome_pipeline(
    inputs: list[str | Path],
    outdir: str | Path,
    *,
    lineage_path: str | Path = DEFAULT_LINEAGE,
    threads: int | None = None,
    log=print,
) -> PipelineResult:
    """Marker scan over supplied proteomes. No Prodigal, no contigs, no pairs.

    Requires only hmmsearch. Demanding Prodigal here would refuse to run on a
    machine that has no gene caller and does not need one.
    """
    t0 = time.time()
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    rt = Runtime()
    if threads:
        rt.threads = threads
    rt.require("hmmsearch")

    lineage = Lineage.load(resolve_lineage(lineage_path))
    for line in lineage_notice(lineage):
        log(line)

    results: list[ProteomeResult] = []
    for path in inputs:
        log(f"[{len(results) + 1}/{len(inputs)}] {Path(path).name}")
        r = analyze_proteome(path, outdir / "proteomes", rt, lineage)
        log(
            f"    {r.stats['n_proteins']} proteins, {r.stats['mean_protein_aa']} aa mean "
            f"| completeness {r.markers.completeness}% dup {r.markers.duplication_percent}% "
            f"({r.seconds:.0f}s)"
        )
        results.append(r)

    log(
        "not measured from a proteome: "
        + ", ".join(sorted(NOT_MEASURED))
        + " -- see the report for why"
    )
    return PipelineResult(
        genomes=[], proteomes=results, runtime=rt, outdir=outdir,
        seconds=time.time() - t0,
    )


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

    lineage = Lineage.load(resolve_lineage(lineage_path))
    for line in lineage_notice(lineage):
        log(line)

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
