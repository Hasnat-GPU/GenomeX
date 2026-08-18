"""Synthetic genomes with known ground truth.

Contamination detection is easy to make *look* like it works -- flag a few
contigs, print a verdict, nobody checks. These generators exist so the tests can
check: a genome built from one composition must come back clean, and a genome
built from two must come back with the planted contigs named.

Composition differences are produced by a first-order Markov chain with a
species-specific transition matrix, so the planted contigs differ in k-mer
composition the way a real foreign organism does -- not merely in GC.
"""

from __future__ import annotations

import random
from dataclasses import dataclass

BASES = "ACGT"


def _transition_matrix(gc: float, seed: int) -> dict[str, list[float]]:
    """Per-species base transition probabilities with the requested GC target."""
    rng = random.Random(seed)
    at, gcp = (1.0 - gc) / 2, gc / 2
    base_probs = {"A": at, "C": gcp, "G": gcp, "T": at}
    matrix: dict[str, list[float]] = {}
    for prev in BASES:
        # Species-specific tilt: same GC, different di/tetra-nucleotide structure.
        weights = []
        for b in BASES:
            tilt = 1.0 + 0.6 * (rng.random() - 0.5)
            weights.append(max(1e-6, base_probs[b] * tilt))
        total = sum(weights)
        matrix[prev] = [w / total for w in weights]
    return matrix


def make_sequence(
    length: int, gc: float, species_seed: int, stream_seed: int = 0, *, start: str = "A"
) -> str:
    """One organism == one transition matrix.

    `species_seed` picks the composition and must be shared by every contig of the
    same organism -- otherwise each contig is its own species in k-mer terms and a
    contamination test built on such a fixture proves nothing.  `stream_seed` only
    varies the sampling, the way two contigs of one chromosome differ.
    """
    matrix = _transition_matrix(gc, species_seed)
    rng = random.Random(species_seed * 7919 + stream_seed * 104729 + length)
    out = [start]
    for _ in range(length - 1):
        out.append(rng.choices(BASES, weights=matrix[out[-1]], k=1)[0])
    return "".join(out)


@dataclass
class SyntheticGenome:
    name: str
    records: list[tuple[str, str]]
    foreign_contigs: set[str]

    @property
    def total_bp(self) -> int:
        return sum(len(s) for _, s in self.records)


def clean_genome(
    name: str = "host",
    *,
    n_contigs: int = 12,
    contig_len: int = 40_000,
    gc: float = 0.62,
    seed: int = 1,
) -> SyntheticGenome:
    """One organism: one composition, no planted contigs."""
    records = [
        (f"{name}_contig{i + 1}", make_sequence(contig_len, gc, seed, i))
        for i in range(n_contigs)
    ]
    return SyntheticGenome(name=name, records=records, foreign_contigs=set())


def contaminated_genome(
    name: str = "mixed",
    *,
    n_host: int = 12,
    n_foreign: int = 3,
    contig_len: int = 40_000,
    host_gc: float = 0.62,
    foreign_gc: float = 0.35,
    seed: int = 1,
) -> SyntheticGenome:
    """Host plus a planted foreign organism. `foreign_contigs` is the answer key."""
    records = [
        (f"{name}_host{i + 1}", make_sequence(contig_len, host_gc, seed, i))
        for i in range(n_host)
    ]
    foreign = [
        (f"{name}_foreign{i + 1}", make_sequence(contig_len, foreign_gc, 5000 + seed, i))
        for i in range(n_foreign)
    ]
    records.extend(foreign)
    return SyntheticGenome(
        name=name, records=records, foreign_contigs={h for h, _ in foreign}
    )


def gene_pair_genomes(
    *,
    shared: int = 40,
    unique_a: int = 6,
    unique_b: int = 9,
    gene_len_aa: int = 220,
    seed: int = 11,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]], dict]:
    """Two protein sets with a known shared core and known private genes.

    Returns (proteins_a, proteins_b, truth) where truth records the expected
    core/unique counts that compare.classify_orthogroups must reproduce.
    """
    rng = random.Random(seed)
    aa = "ACDEFGHIKLMNPQRSTVWY"

    def protein() -> str:
        return "M" + "".join(rng.choice(aa) for _ in range(gene_len_aa - 1))

    core = [protein() for _ in range(shared)]
    a_only = [protein() for _ in range(unique_a)]
    b_only = [protein() for _ in range(unique_b)]

    prot_a = [(f"A_contig1_{i + 1}", s) for i, s in enumerate(core + a_only)]
    prot_b = [(f"B_contig1_{i + 1}", s) for i, s in enumerate(core + b_only)]
    truth = {"shared": shared, "unique_a": unique_a, "unique_b": unique_b}
    return prot_a, prot_b, truth
