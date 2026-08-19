"""Fragmentation invariance: the ground truth that costs nothing.

Shred one finished genome into 10, 50, 100, 300 and 1000 contigs. Same DNA, same
organism, one composition -- only the contig count changes. A contamination
detector must return the same verdict all the way down the ladder.

This is the defect stated in units the repository owns. The bimodality rule fires
on 0/9 finished genomes and 23/24 assemblies with 100-300 contigs, on genomes
CheckM2 calls clean; if firing rate is a function of contig count, the detector is
measuring assembly fragmentation and calling it contamination.

Crucially this measurement is independent of CheckM2. Any threshold tuned after
seeing the 72-genome comparison is fitted to it; this ladder can be used to
develop against without that hazard, because the right answer is known a priori
from how the input was constructed.

Usage:
    python bench/fragmentation_ladder.py --genomes A.fna B.fna --out ladder.tsv
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomex.contamination import detect_contamination  # noqa: E402
from genomex.fasta import Assembly, Contig  # noqa: E402
from genomex.genes import parse_prodigal_proteins  # noqa: E402

RUNGS = (10, 50, 100, 300, 1000)


def load_marker_positions(run_dir: Path, genome: str) -> list[tuple[str, int]]:
    """(contig, midpoint) for every core marker located in a completed run.

    Marker evidence is half of what the detector reasons with -- without it no
    plasmid call can be justified -- so a ladder run without markers measures a
    crippled detector. Scanning once on the intact assembly and mapping the
    coordinates onto each shredding gives every rung exactly the marker set the
    pipeline would compute, at no extra cost and with no ground truth leaked:
    shredding only cuts, so a marker's position is unchanged.
    """
    gdir = run_dir / "genomes" / genome
    faa, tsv = gdir / "proteins.faa", gdir / "markers.tsv"
    if not (faa.exists() and tsv.exists()):
        return []
    by_id = {g.protein_id: g for g in parse_prodigal_proteins(faa)}
    out: list[tuple[str, int]] = []
    for line in tsv.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("	")
        if len(f) < 4 or f[1] == "Missing" or f[2] == "-":
            continue
        gene = by_id.get(f[2])
        if gene and gene.start:
            out.append((gene.contig, (gene.start + gene.end) // 2))
    return out


def marker_counts_for(pieces: list[tuple[Contig, str, int, int]],
                      markers: list[tuple[str, int]]) -> dict[str, int]:
    """Assign each marker to the fragment whose parent span contains it."""
    counts: dict[str, int] = {}
    for contig, parent, lo, hi in pieces:
        n = sum(1 for pc, pos in markers if pc == parent and lo <= pos < hi)
        if n:
            counts[contig.name] = n
    return counts


def shred_with_spans(asm: Assembly, n_pieces: int, seed: int = 0):
    """Shred, returning both the assembly and each piece's parent span."""
    rng = random.Random(seed * 1000 + n_pieces)
    total = asm.total_bp
    pieces: list[Contig] = []
    spans: list[tuple[Contig, str, int, int]] = []
    for c in asm.contigs:
        share = max(1, round(n_pieces * c.length / total))
        if share == 1 or c.length < 2000:
            piece = Contig(name=f"{c.name}_p1", seq=c.seq)
            pieces.append(piece)
            spans.append((piece, c.name, 0, c.length))
            continue
        cuts = sorted(rng.sample(range(1000, c.length - 1000), min(share - 1, c.length // 1000 - 2)))
        bounds = [0, *cuts, c.length]
        for i, (lo, hi) in enumerate(zip(bounds, bounds[1:]), start=1):
            piece = Contig(name=f"{c.name}_p{i}", seq=c.seq[lo:hi])
            pieces.append(piece)
            spans.append((piece, c.name, lo, hi))
    return Assembly(path=asm.path, name=f"{asm.name}@{n_pieces}", contigs=pieces), spans


def shred(asm: Assembly, n_pieces: int, seed: int = 0) -> Assembly:
    """Cut an assembly into ~n_pieces contigs at reproducible random breakpoints.

    Breakpoints are drawn per replicon in proportion to its length, so a genome
    with a chromosome and a plasmid keeps both represented at every rung.
    """
    rng = random.Random(seed * 1000 + n_pieces)
    total = asm.total_bp
    pieces: list[Contig] = []
    for c in asm.contigs:
        share = max(1, round(n_pieces * c.length / total))
        if share == 1 or c.length < 2000:
            pieces.append(Contig(name=f"{c.name}_p1", seq=c.seq))
            continue
        cuts = sorted(rng.sample(range(1000, c.length - 1000), min(share - 1, c.length // 1000 - 2)))
        bounds = [0, *cuts, c.length]
        for i, (lo, hi) in enumerate(zip(bounds, bounds[1:]), start=1):
            pieces.append(Contig(name=f"{c.name}_p{i}", seq=c.seq[lo:hi]))
    return Assembly(path=asm.path, name=f"{asm.name}@{n_pieces}", contigs=pieces)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomes", nargs="+", required=True, type=Path)
    ap.add_argument("--out", type=Path, default=Path("bench/fragmentation_ladder.tsv"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument(
        "--run-dir", type=Path, default=None,
        help="a completed GenomeX --outdir; supplies the marker evidence the "
             "detector needs to tell a plasmid from a contaminant",
    )
    args = ap.parse_args()

    rows = [
        "genome\trequested_pieces\tscored_contigs\tverdict\tsuspect_contigs\t"
        "suspect_fraction_pct\tbimodal\tmean_contig_kb"
    ]
    flips = 0
    per_genome: dict[str, list[str]] = {}

    for path in args.genomes:
        asm = Assembly.load(path)
        verdicts = []
        print(f"\n{asm.name}: {len(asm.contigs)} contigs, {asm.total_bp / 1e6:.2f} Mb")
        markers = load_marker_positions(args.run_dir, asm.name) if args.run_dir else []
        if args.run_dir and not markers:
            print(f"  warning: no marker scan found for {asm.name}; "
                  f"the plasmid/contaminant distinction is disabled")
        for n in RUNGS:
            frag, spans = shred_with_spans(asm, n, seed=args.seed)
            counts = marker_counts_for(spans, markers) if markers else None
            res = detect_contamination(frag, contig_marker_counts=counts)
            scored = res.params.get("contigs_scored", 0)
            mean_kb = frag.total_bp / max(1, len(frag.contigs)) / 1000
            verdicts.append(res.verdict)
            rows.append(
                f"{asm.name}\t{n}\t{scored}\t{res.verdict}\t{res.n_suspect_contigs}\t"
                f"{round(100 * res.suspect_fraction, 2)}\t{res.bins.get('bimodal')}\t{mean_kb:.1f}"
            )
            print(
                f"  {n:>5} pieces -> {scored:>4} scored | {res.verdict:<12} "
                f"suspect {res.n_suspect_contigs:>3} ({100 * res.suspect_fraction:5.2f}%) "
                f"bimodal={res.bins.get('bimodal')}"
            )
        per_genome[asm.name] = verdicts
        if len(set(verdicts)) > 1:
            flips += 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"\n{'=' * 62}")
    print(f"genomes with a verdict that changes down the ladder: {flips}/{len(per_genome)}")
    for name, vs in per_genome.items():
        mark = "FLIPS" if len(set(vs)) > 1 else "stable"
        print(f"  {mark:<7} {name:<22} {' -> '.join(vs)}")
    print(f"\nwritten: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
