"""Proteome input: a protein FASTA in place of an assembly.

`scan_markers` has always taken a proteome -- the assembly path just hands it
Prodigal's output. This module supplies the other half of that route: loading a
user's `.faa`, and stating plainly which of GenomeX's answers it cannot support.

That second half is the whole reason this is a module rather than three lines in
the CLI. A proteome carries no contigs, no coordinates and no nucleotides, so
contamination, ANI and genomic islands are not merely absent -- they must not be
rendered as zeroes. `NOT_MEASURED` below is the machine-readable form of that
refusal, and it distinguishes two things a reader would otherwise conflate:

  impossible    the input does not contain the information, at all
  unimplemented the information is derivable, but this path does not do it

Naming the second honestly matters. Clustering proteomes into orthogroups needs
no nucleotides and MMseqs2 would do it happily; what stops it is that the
comparative step attributes every unique gene to a contig and a coordinate.
Calling that "impossible" would be a claim about biology when it is a fact about
this code.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from pathlib import Path

from .fasta import iter_fasta

#: Residues with no nucleotide meaning under any IUPAC reading. This is the
#: discriminator, and it is deliberately not "how much of the file is ACGTU".
#:
#: Counting ACGTU fails in the direction that matters: a scaffolded assembly
#: padded with N drops below any threshold a proteome clears. A 20%-N scaffold
#: scores 0.80 on plain ACGTU -- under a 0.90 bar it would be rejected as
#: protein. Measured on this alphabet instead: S. cerevisiae proteome 0.3535,
#: its genome 0.0000, a bacterial assembly 0.0000, the 20%-N scaffold 0.0000.
PROTEIN_ONLY = frozenset("EFILPQXZJO")

#: Above this share of protein-only residues the file is protein. The gap
#: between the two populations is 0.35 wide, so the exact value is not
#: load-bearing; it is set low because a false "this is protein" costs an error
#: message and a false "this is nucleotide" costs a fabricated GC content.
PROTEIN_ONLY_FLOOR = 0.05

#: The full nucleotide alphabet including ambiguity codes, used only in the
#: other direction: to say a file is affirmatively DNA rather than merely low in
#: protein-only residues. A 15-residue peptide can contain none of `PROTEIN_ONLY`
#: by chance, and must not be refused for it.
IUPAC_NUCLEOTIDE = frozenset("ACGTURYSWKMBDHVN")
NUCLEOTIDE_CEILING = 0.95

#: Below this many residues the nucleotide test is not applied at all.
#:
#: Sixteen of the twenty-six letters are legal nucleotide codes -- M, K, V and S
#: are ambiguity codes as well as methionine, lysine, valine and serine -- so a
#: short peptide can be 100% nucleotide-legal by chance. E, F, I, L, P and Q are
#: about 30% of residues in real proteins, so the chance a peptide this long
#: contains none of them is ~0.7**200, and any misdirected genome file clears
#: 200 bases trivially. A nucleotide file shorter than this is accepted and
#: scored as protein; it will report 0% completeness, which is the one case here
#: where the input is at fault and the output looks like an error.
MIN_RESIDUES_FOR_ALPHABET_TEST = 200


def protein_only_fraction(seqs) -> tuple[float, int]:
    """Share of residues that cannot be nucleotide, over an iterable of strings."""
    total = hits = 0
    for s in seqs:
        s = s.rstrip("*").upper()
        total += len(s)
        hits += sum(1 for ch in s if ch in PROTEIN_ONLY)
    return (hits / total if total else 0.0), total


#: What an assembly would have answered and a proteome cannot, with the reason
#: attached. Rendered in the report and emitted in the JSON, so a machine reading
#: the output finds an explicit refusal where it would otherwise find nothing.
NOT_MEASURED: dict[str, dict[str, str]] = {
    "assembly_stats": {
        "status": "impossible",
        "reason": "genome size, contig count, N50 and GC are properties of nucleotide "
                  "contigs; a proteome contains none",
    },
    "coding_density": {
        "status": "impossible",
        "reason": "the ratio of coding bases to genome length needs a genome length",
    },
    "contamination": {
        "status": "impossible",
        "reason": "composition analysis compares each contig's tetranucleotide signature "
                  "against the assembly's own; with no contigs there is nothing to place "
                  "a foreign signal on, and no way to say how much of the assembly it is",
    },
    "ani": {
        "status": "impossible",
        "reason": "fastANI aligns nucleotide sequence",
    },
    "genomic_islands": {
        "status": "impossible",
        "reason": "islands are runs of consecutive genes along a contig; a proteome has "
                  "no gene order",
    },
    "orthogroups": {
        "status": "unimplemented",
        "reason": "clustering proteomes needs no nucleotides, but the comparative step "
                  "attributes every unique gene to a contig and a coordinate, which this "
                  "input does not carry",
    },
}


class NotAProteome(ValueError):
    """Raised when a file given as a proteome is not one.

    Worth its own exception because the failure is otherwise silent: HMMER will
    happily scan a six-frame-looking nucleotide file, find almost nothing, and
    report 0.3% completeness -- a number that looks like a finding about the
    genome rather than a mistake about the file.
    """


class NotAnAssembly(ValueError):
    """Raised when a file given as an assembly is a proteome.

    The reciprocal guard, and the more important one. A proteome handed to the
    assembly path is not rejected by anything downstream -- it is *measured*.
    `Assembly.load` makes one contig per protein, and `assembly_stats` then
    reports GC as (Gly + Cys) / (Ala + Cys + Gly + Thr), because those four
    amino acids share their one-letter codes with the four bases.

    On the S. cerevisiae S288C proteome that yields **GC 35.36%** against the
    real genome's 38.15%: not an error value, not out of range, and not
    something a reader can catch. Size comes out as 2.93 Mb (a residue count),
    contigs as 6021, N50 as 0.6 kb. Only the last two look odd, and they look
    like a bad assembly rather than the wrong file.
    """


def assert_nucleotide(seqs, path) -> None:
    """Refuse protein input on the assembly path. See `NotAnAssembly`.

    Fires only on clear evidence -- an ambiguous or tiny file is let through, so
    that an unusual but genuine assembly is never blocked by this check.
    """
    frac, total = protein_only_fraction(seqs)
    if total and frac >= PROTEIN_ONLY_FLOOR:
        raise NotAnAssembly(
            f"{path} looks like protein sequence ({100.0 * frac:.1f}% of residues are "
            "amino-acid-only codes). Assembly analysis on a proteome does not fail -- it "
            "fabricates: GC content comes out of the glycine and cysteine counts, and "
            "lands in the plausible range. Use the `proteins` subcommand."
        )


@dataclass
class Protein:
    name: str
    length_aa: int
    description: str = ""


@dataclass
class Proteome:
    path: Path
    name: str
    proteins: list[Protein] = field(default_factory=list)

    @classmethod
    def load(cls, path: str | Path, name: str | None = None) -> "Proteome":
        path = Path(path)
        proteins: list[Protein] = []
        residues = 0
        protein_only = 0
        nucleotide_like = 0
        for header, seq in iter_fasta(path):
            parts = header.split(None, 1)
            # Trailing '*' is a stop codon marker, not a residue. Prodigal writes
            # one; NCBI does not. Counting it would inflate every length by one
            # on some inputs and not others.
            clean = seq.rstrip("*").upper()
            proteins.append(
                Protein(
                    name=parts[0],
                    length_aa=len(clean),
                    description=parts[1] if len(parts) > 1 else "",
                )
            )
            residues += len(clean)
            protein_only += sum(1 for ch in clean if ch in PROTEIN_ONLY)
            nucleotide_like += sum(1 for ch in clean if ch in IUPAC_NUCLEOTIDE)

        if not proteins:
            raise NotAProteome(f"{path} contains no FASTA records")
        # Both conditions, not either, and only above a sample size where they
        # discriminate -- see MIN_RESIDUES_FOR_ALPHABET_TEST.
        if residues >= MIN_RESIDUES_FOR_ALPHABET_TEST and (
            nucleotide_like / residues >= NUCLEOTIDE_CEILING
            and protein_only / residues < PROTEIN_ONLY_FLOOR
        ):
            raise NotAProteome(
                f"{path} looks like nucleotide sequence "
                f"({100.0 * nucleotide_like / residues:.1f}% of residues are nucleotide "
                f"codes, {100.0 * protein_only / residues:.1f}% amino-acid-only). "
                "Proteome input expects translated protein; use the `qc` or `run` "
                "subcommand for an assembly."
            )

        stem = path.name
        for suffix in (".gz", ".faa", ".fa", ".fasta", ".fas", ".pep"):
            if stem.endswith(suffix):
                stem = stem[: -len(suffix)]
        return cls(path=path, name=name or stem, proteins=proteins)

    def __len__(self) -> int:
        return len(self.proteins)

    @property
    def protein_ids(self) -> list[str]:
        return [p.name for p in self.proteins]

    def sha256(self) -> str:
        h = hashlib.sha256()
        with open(self.path, "rb") as fh:
            for block in iter(lambda: fh.read(1 << 20), b""):
                h.update(block)
        return h.hexdigest()

    def stats(self) -> dict:
        """Everything measurable about the input itself, and nothing more.

        No genes-per-Mb, no coding density: both divide by a genome length this
        input does not have, and both would return 0.0 rather than fail.
        """
        lengths = sorted(p.length_aa for p in self.proteins)
        n = len(lengths)
        return {
            "n_proteins": n,
            "total_aa": sum(lengths),
            "mean_protein_aa": round(sum(lengths) / n, 1),
            "median_protein_aa": lengths[n // 2],
            "shortest_protein_aa": lengths[0],
            "longest_protein_aa": lengths[-1],
            "source": "user-supplied protein FASTA (not called by GenomeX)",
        }
