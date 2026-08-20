"""Per-marker comparison of GenomeX's marker scan against BUSCO, on a proteome.

`compare_to_busco.py` compares whole GenomeX runs, which means it can only be
pointed at lineages the pipeline can run end to end -- today, bacteria. This
script isolates the marker scan itself: give it a protein FASTA and a lineage,
and it scores that file both ways and diffs marker by marker.

That isolation is the point. It answers "does the classification generalise to
this lineage?" without needing a gene caller for the lineage's domain, which is
what made it possible to validate `fungi_odb10` while GenomeX still has no
eukaryotic gene calling at all.

It is also how the cross-marker rule in `genomex.markers` was found. Scoring the
S. cerevisiae S288C proteome against fungi_odb10 gave 757 of 758 markers
identical to BUSCO 5.8.3; the one disagreement was a real missing rule, inert on
the bacterial collection and decisive here (`docs/decisions.md` #15).

Usage:
    micromamba run -n busco busco -i proteome.faa -l /path/to/fungi_odb10 \\
        -o scer_prot --out_path ~/genomex-work/busco_fungi -m proteins --offline

    python bench/compare_proteome_to_busco.py \\
        --proteins ~/genomex-work/data/Fungi/GCF_000146045.2.faa \\
        --lineage  ~/genomex-work/db/fungi_odb10 \\
        --busco    ~/genomex-work/busco_fungi/scer_prot \\
        --workdir  ~/genomex-work/fungal_probe
"""

from __future__ import annotations

import argparse
import collections
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from genomex.markers import Lineage, scan_markers, write_marker_table  # noqa: E402
from genomex.runtime import Runtime  # noqa: E402

TAB = "\t"


def load_busco_full_table(run_dir: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """BUSCO's per-marker calls: full_table.tsv, columns Busco id / Status / Sequence."""
    matches = list(run_dir.glob("run_*/full_table.tsv"))
    if not matches:
        raise SystemExit(f"no run_*/full_table.tsv under {run_dir}")
    status: dict[str, str] = {}
    seqs: dict[str, list[str]] = collections.defaultdict(list)
    for line in matches[0].read_text(encoding="utf-8").splitlines():
        if line.startswith("#") or not line.strip():
            continue
        f = line.split(TAB)
        status[f[0]] = f[1]
        if len(f) > 2 and f[2] != "-":
            seqs[f[0]].append(f[2])
    return status, seqs


def load_genomex_table(path: Path) -> tuple[dict[str, str], dict[str, list[str]]]:
    """GenomeX writes one row per hit, every row carrying the marker's status.
    The vocabulary is BUSCO's, so the two tables compare directly."""
    status: dict[str, str] = {}
    seqs: dict[str, list[str]] = collections.defaultdict(list)
    for line in path.read_text(encoding="utf-8").splitlines()[1:]:
        if not line.strip():
            continue
        f = line.split(TAB)
        status[f[0]] = f[1]
        if f[2] != "-":
            seqs[f[0]].append(f[2])
    return status, seqs


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--proteins", required=True, type=Path, help="protein FASTA")
    ap.add_argument("--lineage", required=True, type=Path, help="BUSCO odb10 lineage dir")
    ap.add_argument("--busco", required=True, type=Path,
                    help="a completed BUSCO -m proteins output dir for the same file")
    ap.add_argument("--workdir", type=Path, default=Path("bench/proteome_scan"))
    ap.add_argument("--threads", type=int, default=6)
    args = ap.parse_args(argv)

    lin = Lineage.load(args.lineage)
    out = args.workdir / args.proteins.stem
    out.mkdir(parents=True, exist_ok=True)
    res = scan_markers(args.proteins, lin, out, Runtime(threads=args.threads))
    write_marker_table(res, out / "markers.tsv")

    b_status, b_seqs = load_busco_full_table(args.busco)
    g_status, g_seqs = load_genomex_table(out / "markers.tsv")

    print(f"lineage {lin.name}: {lin.n_markers} markers")
    print(f"  GenomeX {res.summary()['busco_style_string']}")
    print(f"  BUSCO   {dict(sorted(collections.Counter(b_status.values()).items()))}")
    print(f"  GenomeX {dict(sorted(collections.Counter(g_status.values()).items()))}")

    disagree = sorted(b for b in b_status if b_status[b] != g_status.get(b))
    # Aggregate percentages are the weaker check: two tools can agree on the
    # headline and still differ on *which* markers are duplicated, which is the
    # number the contamination module consumes.
    named = sorted(
        b for b in b_status
        if b_status[b] in ("Complete", "Duplicated")
        and sorted(b_seqs.get(b, [])) != sorted(g_seqs.get(b, []))
    )
    for b in disagree:
        print(f"  DIFFERS {b}: BUSCO={b_status[b]} {sorted(b_seqs.get(b, []))} "
              f"| GenomeX={g_status.get(b)} {sorted(g_seqs.get(b, []))}")
    for b in named:
        if b not in disagree:
            print(f"  SAME STATUS, DIFFERENT PROTEIN {b}: "
                  f"BUSCO {sorted(b_seqs.get(b, []))} | GenomeX {sorted(g_seqs.get(b, []))}")

    agree = len(b_status) - len(disagree)
    print(f"\n{agree}/{len(b_status)} markers identical in status, "
          f"{len(b_status) - len(named)}/{len(b_status)} also naming the same protein(s)")
    return 0 if not disagree and not named else 1


if __name__ == "__main__":
    raise SystemExit(main())
