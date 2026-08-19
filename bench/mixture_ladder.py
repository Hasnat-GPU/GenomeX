"""Known-ratio mixtures: ground truth for contamination *recall*.

`fragmentation_ladder.py` answers "does the verdict survive shredding?" -- a
specificity instrument. It has nothing to say about recall, because nothing
foreign is ever added. The CheckM2 comparison does measure recall, but over the
72-genome collection only four genomes exceed 5% contamination, so recall there
is an average of four coin flips, and CheckM2 is a held-out test set that should
be spent once rather than developed against.

This harness supplies the missing instrument. Take a finished genome, shred it,
splice in fragments of a *different* genome until they account for a known
fraction of the total, and ask the detector to name them. The truth is exact:
every contig is labelled by construction. Nothing is inferred, nothing is
borrowed from another tool, and the answer cannot drift when a threshold moves.

Two design choices matter for honesty:

* **Donor fragments are the same mean length as host fragments.** Real
  contaminants are often shorter -- lower coverage assembles worse -- and length
  would then be a giveaway the length-conditioned null could exploit. Matching
  the lengths removes that shortcut and measures composition alone.
* **Donor fragments carry their share of the single-copy core.** They are drawn
  from a real genome, so a 5% mixture displaces roughly 5% of the markers onto
  foreign sequence. That is precisely what should distinguish a contaminant from
  a plasmid, and withholding it would test a detector nobody runs.

Donor distance is the second axis. A same-genus donor at nearly identical GC is
the case that matters -- cross-contamination in a sequencing run is usually
between related isolates on the same plate -- and it is the case GC screening
cannot touch.

Usage:
    python bench/mixture_ladder.py \\
        --data-dir ~/genomex-work/data/Bacteria \\
        --run-dir  ~/genomex-work/runs/sweep \\
        --pair GCF_000020045.1:GCF_001449005.1 \\
        --out bench/mixture_ladder.tsv
"""

from __future__ import annotations

import argparse
import random
import sys
from dataclasses import dataclass
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomex.contamination import detect_contamination  # noqa: E402
from genomex.fasta import Assembly, Contig  # noqa: E402
from genomex.genes import parse_prodigal_proteins  # noqa: E402

from fragmentation_ladder import shred_with_spans  # noqa: E402

#: Donor mass as a percentage of host mass. 0 is the control: the same host, the
#: same shredding, nothing added. Any flag there is a false positive with no
#: possible innocent reading.
RUNGS = (0.0, 1.0, 2.0, 5.0, 10.0, 20.0)

#: Host fragment count. 300 pieces of an 8-9 Mb genome is ~29 kb each -- an
#: ordinary draft contig, well above the 3 kb scoring floor, and fine-grained
#: enough that a 1% mixture is three fragments rather than one. Coarser
#: fragments cannot express the low rungs: at 55 kb a 1% target rounds to 0.7%
#: or 1.3% and the ladder's x-axis stops meaning what it says.
HOST_PIECES = 300

DONOR_PREFIX = "DONOR__"


@dataclass(frozen=True)
class MarkerHit:
    busco_id: str
    contig: str
    midpoint: int


def load_markers(run_dir: Path, genome: str) -> list[MarkerHit]:
    """Every located core marker, with the identity needed to detect duplication.

    Positions come from one scan of the intact assembly. Shredding only cuts, so
    a marker's coordinate is unchanged by it; mapping the coordinates forward
    gives each mixture exactly the marker table the pipeline would compute, with
    no ground truth leaked and no HMM search per rung.
    """
    gdir = run_dir / "genomes" / genome
    faa, tsv = gdir / "proteins.faa", gdir / "markers.tsv"
    if not (faa.exists() and tsv.exists()):
        return []
    by_id = {g.protein_id: g for g in parse_prodigal_proteins(faa)}
    out: list[MarkerHit] = []
    for line in tsv.read_text(encoding="utf-8").splitlines()[1:]:
        f = line.split("\t")
        if len(f) < 4 or f[1] == "Missing" or f[2] == "-":
            continue
        gene = by_id.get(f[2])
        if gene and gene.start:
            out.append(MarkerHit(f[0], gene.contig, (gene.start + gene.end) // 2))
    return out


def assign_markers(
    pieces: list[tuple[Contig, str, int, int]],
    markers: list[MarkerHit],
    rename: dict[str, str] | None = None,
) -> tuple[dict[str, int], dict[str, list[str]]]:
    """Map markers onto fragments: (counts per contig, contigs per busco id)."""
    counts: dict[str, int] = {}
    where: dict[str, list[str]] = {}
    for contig, parent, lo, hi in pieces:
        name = (rename or {}).get(contig.name, contig.name)
        for m in markers:
            if m.contig == parent and lo <= m.midpoint < hi:
                counts[name] = counts.get(name, 0) + 1
                where.setdefault(m.busco_id, []).append(name)
    return counts, where


@dataclass
class Mixture:
    assembly: Assembly
    donor_contigs: set[str]
    donor_bp: int
    marker_counts: dict[str, int]
    duplicated: dict[str, list[str]]
    achieved_pct: float


def build_mixture(
    host: Assembly,
    donor: Assembly,
    host_markers: list[MarkerHit],
    donor_markers: list[MarkerHit],
    pct: float,
    *,
    seed: int = 0,
    host_pieces: int = HOST_PIECES,
) -> Mixture:
    """Shred host and donor to matched fragment lengths, splice to `pct` by mass."""
    host_frag, host_spans = shred_with_spans(host, host_pieces, seed=seed)
    counts, where = assign_markers(host_spans, host_markers)
    pieces = list(host_frag.contigs)
    donor_names: set[str] = set()
    donor_bp = 0

    # `pct` is the donor's share of the *mixture*, which is what a contamination
    # percentage means, so the donor mass solves d/(h+d) = pct/100. Taking pct of
    # the host instead silently delivers 16.7% when the row says 20%, and every
    # recall figure plotted against that axis is then mislabelled.
    target = int(round(host.total_bp * pct / (100.0 - pct))) if pct < 100 else 0
    if target > 0:
        mean_piece = max(1, host.total_bp // max(1, len(host_frag.contigs)))
        n_donor = max(2, round(donor.total_bp / mean_piece))
        donor_frag, donor_spans = shred_with_spans(donor, n_donor, seed=seed + 1)
        rename = {c.name: f"{DONOR_PREFIX}{c.name}" for c, _, _, _ in donor_spans}

        # Only fragments of a host-like size are eligible. `shred_with_spans`
        # leaves a small replicon whole when its proportional share rounds to
        # one, and a whole 300 kb plasmid dropped into a 1% mixture would both
        # overshoot the target and hand the detector a length cue that no real
        # contaminant provides.
        usable = [s for s in donor_spans if 3000 <= s[0].length <= 2 * mean_piece]
        order = list(range(len(usable)))
        random.Random(seed * 7919 + int(pct * 100)).shuffle(order)
        chosen = []
        for i in order:
            if donor_bp >= target:
                break
            contig, parent, lo, hi = usable[i]
            if donor_bp + contig.length > target * 1.15:
                continue  # would overshoot; a later, smaller fragment may fit
            chosen.append(usable[i])
            donor_bp += contig.length
        # Selection order is randomised, but the pieces are re-sorted before use
        # so the assembly is a deterministic function of (host, donor, pct).
        chosen.sort(key=lambda s: s[0].name)
        for contig, parent, lo, hi in chosen:
            pieces.append(Contig(name=rename[contig.name], seq=contig.seq))
            donor_names.add(rename[contig.name])
        d_counts, d_where = assign_markers(chosen, donor_markers, rename)
        counts.update(d_counts)
        for busco, names in d_where.items():
            where.setdefault(busco, []).extend(names)

    mixed = Assembly(
        path=host.path, name=f"{host.name}+{donor.name}@{pct:g}pct", contigs=pieces
    )
    return Mixture(
        assembly=mixed,
        donor_contigs=donor_names,
        donor_bp=donor_bp,
        marker_counts=counts,
        duplicated=where,
        achieved_pct=100.0 * donor_bp / max(1, mixed.total_bp),
    )


def evaluate(mix: Mixture, total_markers: int) -> dict:
    res = detect_contamination(
        mix.assembly,
        contig_marker_counts=mix.marker_counts,
        duplicated_marker_contigs=mix.duplicated,
        total_markers=total_markers,
    )
    by_len = {c.name: c.length for c in mix.assembly.contigs}
    suspect = res.suspect_contig_names()
    replicon = {c.name for c in res.contigs if c.call == "replicon_candidate"}
    truth = mix.donor_contigs

    hit_bp = sum(by_len[n] for n in suspect & truth)
    suspect_bp = sum(by_len[n] for n in suspect)
    n_dup = sum(1 for v in mix.duplicated.values() if len(set(v)) >= 2)
    return {
        "verdict": res.verdict,
        "donor_pct": round(mix.achieved_pct, 2),
        "donor_contigs": len(truth),
        "recall_ctg": round(len(suspect & truth) / len(truth), 3) if truth else None,
        "recall_bp": round(hit_bp / mix.donor_bp, 3) if mix.donor_bp else None,
        "precision_bp": round(hit_bp / suspect_bp, 3) if suspect_bp else None,
        "host_false_pos": len(suspect - truth),
        "donor_as_replicon": len(replicon & truth),
        "duplicated_markers": n_dup,
        "duplication_pct": round(100 * n_dup / max(1, total_markers), 2),
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-dir", type=Path, required=True)
    ap.add_argument("--run-dir", type=Path, required=True)
    ap.add_argument(
        "--pair", action="append", required=True, metavar="HOST:DONOR",
        help="accessions, e.g. GCF_000020045.1:GCF_001449005.1; repeatable",
    )
    ap.add_argument("--out", type=Path, default=Path("bench/mixture_ladder.tsv"))
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--total-markers", type=int, default=124,
                    help="size of the marker set the run-dir was scanned with; "
                         "bacteria_odb10 is 124, burkholderiales_odb10 is 688")
    ap.add_argument("--host-pieces", type=int, default=HOST_PIECES)
    args = ap.parse_args()

    cols = [
        "host", "donor", "requested_pct", "donor_pct", "verdict", "donor_contigs",
        "recall_ctg", "recall_bp", "precision_bp", "host_false_pos",
        "donor_as_replicon", "duplicated_markers", "duplication_pct",
    ]
    rows = ["\t".join(cols)]
    caught: list[tuple[float, bool]] = []

    for spec in args.pair:
        h_acc, _, d_acc = spec.partition(":")
        host = Assembly.load(args.data_dir / f"{h_acc}.fna")
        donor = Assembly.load(args.data_dir / f"{d_acc}.fna")
        hm = load_markers(args.run_dir, h_acc)
        dm = load_markers(args.run_dir, d_acc)
        if not hm or not dm:
            print(f"skip {spec}: no marker scan in {args.run_dir}")
            continue
        print(f"\n{h_acc} ({host.total_bp/1e6:.2f} Mb) <- {d_acc} "
              f"({donor.total_bp/1e6:.2f} Mb), markers {len(hm)}/{len(dm)}")
        print(f"  {'donor%':>7} {'verdict':<12} {'recall_bp':>9} {'prec_bp':>8} "
              f"{'fp':>4} {'as_repl':>8} {'dup':>4}")
        for pct in RUNGS:
            mix = build_mixture(host, donor, hm, dm, pct, seed=args.seed,
                                host_pieces=args.host_pieces)
            m = evaluate(mix, args.total_markers)
            rows.append("\t".join(str(x) for x in [
                h_acc, d_acc, pct, m["donor_pct"], m["verdict"], m["donor_contigs"],
                m["recall_ctg"], m["recall_bp"], m["precision_bp"],
                m["host_false_pos"], m["donor_as_replicon"],
                m["duplicated_markers"], m["duplication_pct"],
            ]))
            print(f"  {m['donor_pct']:>7.2f} {m['verdict']:<12} "
                  f"{str(m['recall_bp']):>9} {str(m['precision_bp']):>8} "
                  f"{m['host_false_pos']:>4} {m['donor_as_replicon']:>8} "
                  f"{m['duplicated_markers']:>4}")
            if pct > 0:
                caught.append((pct, m["verdict"] != "clean"))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(rows) + "\n", encoding="utf-8")

    print(f"\n{'=' * 62}")
    for level in sorted({p for p, _ in caught}):
        hits = [ok for p, ok in caught if p == level]
        print(f"  {level:>5.1f}% donor: verdict non-clean on {sum(hits)}/{len(hits)} pairs")
    print(f"\nwritten: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
