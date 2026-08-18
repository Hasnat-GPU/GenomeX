"""Gene prediction with Prodigal, and parsing of its output.

Prodigal writes protein FASTA headers of the form:

    >contig_3 # 1200 # 2043 # 1 # ID=3_2;partial=00;start_type=ATG;rbs_motif=...

The sequence id is `<contig>_<gene index on that contig>`, which is what lets
downstream modules attribute a marker hit or an orthogroup member back to a
specific contig -- the link that makes contamination visible at gene level.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from .fasta import iter_fasta
from .runtime import Runtime

# Genetic code 11 = bacteria/archaea/plant plastid. 4 = Mycoplasma/Spiroplasma.
DEFAULT_GENETIC_CODE = 11


@dataclass
class Gene:
    protein_id: str      # e.g. NZ_CP012345.1_17
    contig: str
    start: int
    end: int
    strand: int
    partial: str         # '00' complete, '10'/'01'/'11' truncated at an edge
    length_aa: int

    @property
    def nt_length(self) -> int:
        return abs(self.end - self.start) + 1


_HDR = re.compile(
    r"^(?P<pid>\S+)\s+#\s+(?P<start>\d+)\s+#\s+(?P<end>\d+)\s+#\s+(?P<strand>-?1)\s+#\s+(?P<attrs>.*)$"
)


def parse_prodigal_proteins(faa: str | Path) -> list[Gene]:
    genes: list[Gene] = []
    for header, seq in iter_fasta(faa):
        m = _HDR.match(header)
        if not m:
            # Not Prodigal-formatted; fall back to id-only so the pipeline degrades
            # rather than crashing on a user-supplied protein file.
            pid = header.split()[0]
            genes.append(
                Gene(pid, pid.rsplit("_", 1)[0], 0, 0, 0, "00", len(seq.rstrip("*")))
            )
            continue
        pid = m.group("pid")
        attrs = dict(
            kv.split("=", 1) for kv in m.group("attrs").split(";") if "=" in kv
        )
        genes.append(
            Gene(
                protein_id=pid,
                contig=pid.rsplit("_", 1)[0],
                start=int(m.group("start")),
                end=int(m.group("end")),
                strand=int(m.group("strand")),
                partial=attrs.get("partial", "00"),
                length_aa=len(seq.rstrip("*")),
            )
        )
    return genes


@dataclass
class GeneCalls:
    proteins_faa: Path
    genes_fna: Path
    gff: Path
    genes: list[Gene]
    genetic_code: int
    procedure: str

    @property
    def n_genes(self) -> int:
        return len(self.genes)

    def coding_density(self, total_bp: int) -> float:
        if not total_bp:
            return 0.0
        return round(sum(g.nt_length for g in self.genes) / total_bp, 4)

    def summary(self, total_bp: int) -> dict:
        lengths = [g.length_aa for g in self.genes]
        partial = sum(1 for g in self.genes if g.partial != "00")
        return {
            "n_genes": self.n_genes,
            "coding_density": self.coding_density(total_bp),
            "mean_protein_aa": round(sum(lengths) / len(lengths), 1) if lengths else 0.0,
            "partial_genes": partial,
            "genes_per_mb": round(1e6 * self.n_genes / total_bp, 1) if total_bp else 0.0,
            "genetic_code": self.genetic_code,
            "prodigal_procedure": self.procedure,
        }


def predict_genes(
    assembly_path: str | Path,
    outdir: str | Path,
    rt: Runtime,
    *,
    genetic_code: int = DEFAULT_GENETIC_CODE,
    total_bp: int | None = None,
    force: bool = False,
) -> GeneCalls:
    """Run Prodigal. Uses `-p meta` for small/fragmented input, `single` otherwise.

    Prodigal's own guidance is that its self-training mode needs ~100 kb of
    sequence; below that the metagenomic pretrained models are more accurate.
    """
    assembly_path = Path(assembly_path)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    faa = outdir / "proteins.faa"
    fna = outdir / "genes.fna"
    gff = outdir / "genes.gff"

    procedure = "single" if (total_bp or 10**9) >= 100_000 else "meta"

    if force or not faa.exists() or faa.stat().st_size == 0:
        rt.run(
            [
                "prodigal",
                "-i", str(assembly_path),
                "-a", str(faa),
                "-d", str(fna),
                "-f", "gff",
                "-o", str(gff),
                "-p", procedure,
                "-g", str(genetic_code),
                "-q",
            ]
        )

    return GeneCalls(
        proteins_faa=faa,
        genes_fna=fna,
        gff=gff,
        genes=parse_prodigal_proteins(faa),
        genetic_code=genetic_code,
        procedure=procedure,
    )
