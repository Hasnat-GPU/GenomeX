"""Why do two isolates from the same environment carry different genes?

The question decomposes into three that can actually be answered from assemblies:

  1. *Are they even the same species?*  fastANI (whole-genome Average Nucleotide
     Identity); the >=95% ANI line is the operational species boundary.
  2. *Which genes differ?*  MMseqs2 clustering of the predicted proteomes into
     orthogroups, then core / accessory / strain-unique partitioning.
  3. *Why?*  For every strain-unique gene, three competing explanations are
     scored against evidence already in hand:
        contamination  -- the gene sits on a contig this pipeline flagged as
                          compositionally foreign, so the "difference" may be an
                          assembly artefact, not biology;
        acquisition    -- the gene sits in a run of consecutive unique genes
                          (a genomic island: prophage, ICE, plasmid remnant),
                          the signature of horizontal transfer;
        divergence     -- an isolated unique gene, consistent with drift, gene
                          loss in the other isolate, or a failed cluster join.

The third step is the point.  A gene-count difference on its own is not a
finding; a gene-count difference partitioned by cause is.
"""

from __future__ import annotations

import random
import re
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

from .fasta import iter_fasta, write_fasta
from .genes import Gene
from .runtime import Runtime

SPECIES_ANI_THRESHOLD = 95.0
ID_SEP = "|"


# --------------------------------------------------------------------------
# Whole-genome relatedness
# --------------------------------------------------------------------------

@dataclass
class AniResult:
    query: str
    reference: str
    ani: float | None
    fragments_mapped: int | None
    fragments_total: int | None
    method: str = "fastANI"

    @property
    def same_species(self) -> bool | None:
        if self.ani is None:
            return None
        return self.ani >= SPECIES_ANI_THRESHOLD

    def summary(self) -> dict:
        return {
            "query": self.query,
            "reference": self.reference,
            "ani_percent": self.ani,
            "fragments_mapped": self.fragments_mapped,
            "fragments_total": self.fragments_total,
            "alignment_fraction_percent": (
                round(100.0 * self.fragments_mapped / self.fragments_total, 2)
                if self.fragments_mapped and self.fragments_total
                else None
            ),
            "same_species_by_ani": self.same_species,
            "species_threshold_percent": SPECIES_ANI_THRESHOLD,
            "method": self.method,
            "note": (
                "fastANI reports no value below roughly 80% ANI; a null here means "
                "the genomes are more distant than the method resolves, not identical."
            ),
        }


def compute_ani(query_fna: str | Path, ref_fna: str | Path, outdir: str | Path, rt: Runtime) -> AniResult:
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    out = outdir / "fastani.txt"
    rt.run(
        ["fastANI", "-q", str(query_fna), "-r", str(ref_fna), "-o", str(out), "-t", str(rt.threads)],
        check=False,
    )
    q, r = Path(query_fna).name, Path(ref_fna).name
    if out.exists() and out.stat().st_size:
        fields = out.read_text().strip().splitlines()[0].split("\t")
        if len(fields) >= 5:
            return AniResult(q, r, float(fields[2]), int(fields[3]), int(fields[4]))
    return AniResult(q, r, None, None, None)


def mash_distances(fna_paths: list[str | Path], outdir: str | Path, rt: Runtime) -> list[dict]:
    """All-vs-all Mash distances -- cheap triage before paying for fastANI."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    sketch = outdir / "all"
    rt.run(["mash", "sketch", "-o", str(sketch), "-p", str(rt.threads)] + [str(p) for p in fna_paths])
    dist_out = outdir / "mash_dist.tsv"
    rt.run(["mash", "dist", "-p", str(rt.threads), f"{sketch}.msh", f"{sketch}.msh"],
           stdout_path=dist_out)
    rows = []
    for line in dist_out.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 5:
            rows.append(
                {
                    "a": Path(f[0]).name,
                    "b": Path(f[1]).name,
                    "mash_distance": float(f[2]),
                    "p_value": float(f[3]),
                    "shared_hashes": f[4],
                    "ani_estimate_percent": round(100 * (1 - float(f[2])), 2),
                }
            )
    return rows


# --------------------------------------------------------------------------
# Orthogroups
# --------------------------------------------------------------------------

def tag_proteins(faa: str | Path, genome: str, out_faa: str | Path) -> Path:
    """Prefix every protein id with its genome so ids stay unique after merging."""
    records = []
    for header, seq in iter_fasta(faa):
        pid = header.split()[0]
        records.append((f"{genome}{ID_SEP}{pid}", seq.rstrip("*")))
    return write_fasta(out_faa, records)


def cluster_proteomes(
    faa_by_genome: dict[str, str | Path],
    outdir: str | Path,
    rt: Runtime,
    *,
    min_seq_id: float = 0.5,
    coverage: float = 0.8,
    force: bool = False,
) -> dict[str, list[str]]:
    """MMseqs2 easy-cluster over the merged proteomes.

    Returns representative -> [tagged member ids].  Defaults (50% identity, 80%
    bidirectional coverage) are the usual bacterial-pangenome operating point:
    tight enough not to merge paralogous families, loose enough to join true
    orthologs across strains.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    merged = outdir / "all_proteins.faa"
    prefix = outdir / "clusters"
    cluster_tsv = Path(f"{prefix}_cluster.tsv")

    if force or not cluster_tsv.exists() or cluster_tsv.stat().st_size == 0:
        records: list[tuple[str, str]] = []
        for genome, faa in faa_by_genome.items():
            for header, seq in iter_fasta(faa):
                pid = header.split()[0]
                records.append((f"{genome}{ID_SEP}{pid}", seq.rstrip("*")))
        write_fasta(merged, records)
        rt.run(
            [
                "mmseqs", "easy-cluster", str(merged), str(prefix), str(outdir / "tmp"),
                "--min-seq-id", str(min_seq_id),
                "-c", str(coverage),
                "--cov-mode", "0",
                "--threads", str(rt.threads),
                "-v", "1",
            ]
        )

    clusters: dict[str, list[str]] = defaultdict(list)
    for line in cluster_tsv.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip():
            continue
        rep, member = line.split("\t")[:2]
        clusters[rep].append(member)
    return dict(clusters)


@dataclass
class Pangenome:
    genomes: list[str]
    clusters: dict[str, list[str]]
    min_seq_id: float = 0.5
    coverage: float = 0.8

    def presence(self) -> dict[str, set[str]]:
        """cluster representative -> set of genomes present in that cluster."""
        return {
            rep: {m.split(ID_SEP, 1)[0] for m in members}
            for rep, members in self.clusters.items()
        }

    def partition(self) -> dict[str, list[str]]:
        n = len(self.genomes)
        core, accessory, unique = [], [], []
        for rep, genomes in self.presence().items():
            if len(genomes) == n:
                core.append(rep)
            elif len(genomes) == 1:
                unique.append(rep)
            else:
                accessory.append(rep)
        return {"core": core, "accessory": accessory, "unique": unique}

    def summary(self) -> dict:
        part = self.partition()
        return {
            "genomes": self.genomes,
            "n_genomes": len(self.genomes),
            "orthogroups_total": len(self.clusters),
            "core": len(part["core"]),
            "accessory": len(part["accessory"]),
            "strain_unique": len(part["unique"]),
            "core_percent": (
                round(100.0 * len(part["core"]) / len(self.clusters), 2) if self.clusters else 0.0
            ),
            "clustering": {
                "tool": "MMseqs2 easy-cluster",
                "min_seq_id": self.min_seq_id,
                "coverage": self.coverage,
                "cov_mode": 0,
            },
        }

    def genes_of(self, genome: str, reps: list[str]) -> list[str]:
        """Protein ids belonging to `genome` within the given clusters."""
        out = []
        for rep in reps:
            for m in self.clusters.get(rep, []):
                g, pid = m.split(ID_SEP, 1)
                if g == genome:
                    out.append(pid)
        return out


# --------------------------------------------------------------------------
# Pairwise difference, with causes attached
# --------------------------------------------------------------------------

@dataclass
class UniqueGene:
    protein_id: str
    contig: str
    start: int
    end: int
    length_aa: int
    gc_percent: float | None
    gc_deviation: float | None
    island_id: int | None
    island_size: int
    on_suspect_contig: bool
    explanation: str

    def as_row(self) -> list[str]:
        return [
            self.protein_id, self.contig, str(self.start), str(self.end),
            str(self.length_aa),
            "" if self.gc_percent is None else f"{self.gc_percent:.2f}",
            "" if self.gc_deviation is None else f"{self.gc_deviation:+.2f}",
            "" if self.island_id is None else str(self.island_id),
            str(self.island_size), "yes" if self.on_suspect_contig else "no",
            self.explanation,
        ]


UNIQUE_GENE_COLUMNS = [
    "protein_id", "contig", "start", "end", "length_aa", "gc_percent",
    "gc_deviation_vs_genome", "island_id", "island_size", "on_suspect_contig",
    "explanation",
]

_NUM = re.compile(r"_(\d+)$")


def _gene_index(protein_id: str) -> int:
    m = _NUM.search(protein_id)
    return int(m.group(1)) if m else -1


def _gene_gc(seq: str) -> float | None:
    s = seq.upper()
    gc = s.count("G") + s.count("C")
    at = s.count("A") + s.count("T")
    return 100.0 * gc / (gc + at) if (gc + at) else None


def _find_islands(
    ordered_ids: list[str], island_min_size: int, island_max_gap: int
) -> list[list[str]]:
    """Runs of unique genes that are consecutive in gene numbering on one contig."""
    islands: list[list[str]] = []
    run: list[str] = []
    for pid in ordered_ids:
        if not run:
            run = [pid]
        elif _gene_index(pid) - _gene_index(run[-1]) <= island_max_gap + 1:
            run.append(pid)
        else:
            if len(run) >= island_min_size:
                islands.append(run)
            run = [pid]
    if len(run) >= island_min_size:
        islands.append(run)
    return islands


def island_enrichment(
    unique_protein_ids: list[str],
    genes: list[Gene],
    *,
    island_min_size: int = 3,
    island_max_gap: int = 1,
    trials: int = 25,
    seed: int = 0,
) -> dict:
    """How many island genes would appear if the unique genes were placed at random?

    When a large share of genes is strain-unique -- routine between genera -- runs
    of three consecutive unique genes arise by chance alone, and calling them
    acquisitions is unfounded.  The observed count is compared against a
    permutation null that keeps each contig's gene count and unique-gene count
    fixed and only reshuffles which genes are unique.  Deterministic by seed.
    """
    by_contig_all: dict[str, list[str]] = defaultdict(list)
    for g in genes:
        by_contig_all[g.contig].append(g.protein_id)
    unique_set = set(unique_protein_ids)

    observed = 0
    expected_total = 0.0
    rng = random.Random(seed)
    for contig, pids in by_contig_all.items():
        ordered = sorted(pids, key=_gene_index)
        uniq_here = [p for p in ordered if p in unique_set]
        if not uniq_here:
            continue
        observed += sum(len(i) for i in _find_islands(uniq_here, island_min_size, island_max_gap))
        k = len(uniq_here)
        for _ in range(trials):
            sample = sorted(rng.sample(ordered, k), key=_gene_index)
            expected_total += sum(
                len(i) for i in _find_islands(sample, island_min_size, island_max_gap)
            )
    expected = expected_total / trials if trials else 0.0
    return {
        "genes_in_islands_observed": observed,
        "genes_in_islands_expected_by_chance": round(expected, 1),
        "enrichment": round(observed / expected, 2) if expected > 0 else None,
        "informative": bool(expected > 0 and observed / expected >= 1.5),
        "null_model": f"{trials} permutations, per-contig gene and unique counts held fixed",
    }


def explain_unique_genes(
    unique_protein_ids: list[str],
    genes: list[Gene],
    *,
    genes_fna: str | Path | None = None,
    suspect_contigs: set[str] | None = None,
    island_min_size: int = 3,
    island_max_gap: int = 1,
) -> list[UniqueGene]:
    """Attach a cause to every strain-unique gene.

    An island is a run of unique genes that are consecutive in the Prodigal gene
    numbering on one contig (allowing `island_max_gap` shared genes in between).
    Runs of at least `island_min_size` are called acquisitions.
    """
    suspect_contigs = suspect_contigs or set()
    by_id = {g.protein_id: g for g in genes}

    gene_gc: dict[str, float] = {}
    genome_gc: float | None = None
    if genes_fna and Path(genes_fna).exists():
        total_gc = total_len = 0
        for header, seq in iter_fasta(genes_fna):
            pid = header.split()[0]
            gc = _gene_gc(seq)
            if gc is not None:
                gene_gc[pid] = gc
                total_gc += seq.upper().count("G") + seq.upper().count("C")
                total_len += sum(seq.upper().count(b) for b in "ACGT")
        genome_gc = 100.0 * total_gc / total_len if total_len else None

    # Group unique genes by contig and find consecutive runs.
    per_contig: dict[str, list[str]] = defaultdict(list)
    for pid in unique_protein_ids:
        g = by_id.get(pid)
        per_contig[g.contig if g else pid.rsplit("_", 1)[0]].append(pid)

    island_of: dict[str, int] = {}
    island_size: dict[int, int] = {}
    island_counter = 0
    for contig, pids in per_contig.items():
        ordered = sorted(pids, key=_gene_index)
        for run in _find_islands(ordered, island_min_size, island_max_gap):
            island_counter += 1
            for p in run:
                island_of[p] = island_counter
            island_size[island_counter] = len(run)

    out: list[UniqueGene] = []
    for pid in unique_protein_ids:
        g = by_id.get(pid)
        contig = g.contig if g else pid.rsplit("_", 1)[0]
        gc = gene_gc.get(pid)
        dev = (gc - genome_gc) if (gc is not None and genome_gc is not None) else None
        iid = island_of.get(pid)
        on_suspect = contig in suspect_contigs
        if on_suspect:
            explanation = "contamination-suspect contig"
        elif iid is not None:
            explanation = "clustered in genomic island (acquisition/HGT candidate)"
        elif dev is not None and abs(dev) >= 5.0:
            explanation = "isolated gene with atypical GC (acquisition candidate)"
        else:
            explanation = "isolated gene, typical composition (divergence or gene loss)"
        out.append(
            UniqueGene(
                protein_id=pid,
                contig=contig,
                start=g.start if g else 0,
                end=g.end if g else 0,
                length_aa=g.length_aa if g else 0,
                gc_percent=gc,
                gc_deviation=dev,
                island_id=iid,
                island_size=island_size.get(iid, 0) if iid else 0,
                on_suspect_contig=on_suspect,
                explanation=explanation,
            )
        )
    return sorted(out, key=lambda u: (u.contig, _gene_index(u.protein_id)))


@dataclass
class PairComparison:
    genome_a: str
    genome_b: str
    ani: AniResult
    shared_orthogroups: int
    unique_a: list[UniqueGene] = field(default_factory=list)
    unique_b: list[UniqueGene] = field(default_factory=list)
    n_genes_a: int = 0
    n_genes_b: int = 0
    island_stats_a: dict = field(default_factory=dict)
    island_stats_b: dict = field(default_factory=dict)

    @staticmethod
    def _breakdown(genes: list[UniqueGene], island_stats: dict | None = None) -> dict:
        counts: dict[str, int] = defaultdict(int)
        for g in genes:
            counts[g.explanation] += 1
        islands = {g.island_id for g in genes if g.island_id}
        return {
            "total": len(genes),
            "by_explanation": dict(counts),
            "islands": len(islands),
            "genes_in_islands": sum(1 for g in genes if g.island_id),
            "on_suspect_contigs": sum(1 for g in genes if g.on_suspect_contig),
            "island_null_test": island_stats or {},
        }

    def summary(self) -> dict:
        return {
            "genome_a": self.genome_a,
            "genome_b": self.genome_b,
            "genes_a": self.n_genes_a,
            "genes_b": self.n_genes_b,
            "ani": self.ani.summary(),
            "shared_orthogroups": self.shared_orthogroups,
            "unique_to_a": self._breakdown(self.unique_a, self.island_stats_a),
            "unique_to_b": self._breakdown(self.unique_b, self.island_stats_b),
        }


def write_unique_gene_table(genes: list[UniqueGene], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("\t".join(UNIQUE_GENE_COLUMNS) + "\n")
        for g in genes:
            fh.write("\t".join(g.as_row()) + "\n")
    return path
