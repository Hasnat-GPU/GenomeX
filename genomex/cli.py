"""Command line entry point.

    python -m genomex run A.fna B.fna --outdir runs/demo
    python -m genomex run *.fna --outdir runs/all --all-pairs
    python -m genomex qc A.fna --outdir runs/qc
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .pipeline import DEFAULT_LINEAGE, run_pipeline
from .report import write_html_report
from .runtime import ToolMissing


def _add_common(p: argparse.ArgumentParser) -> None:
    p.add_argument("inputs", nargs="+", help="assembly FASTA files (.fna/.fa/.fasta, optionally .gz)")
    p.add_argument("--outdir", "-o", default="runs/latest", help="output directory")
    p.add_argument(
        "--lineage",
        default=str(DEFAULT_LINEAGE),
        help="BUSCO odb10 lineage directory (default: %(default)s)",
    )
    p.add_argument("--genetic-code", type=int, default=11, help="Prodigal translation table (11 bacteria)")
    p.add_argument("--min-contig-length", type=int, default=3000,
                   help="shortest contig entering composition analysis")
    p.add_argument("--threads", type=int, default=None)
    p.add_argument("--quiet", action="store_true")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="genomex",
        description="Unified genome QC: completeness, contamination, and comparative gene content.",
    )
    sub = p.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="full pipeline: QC every genome, then compare")
    _add_common(run)
    run.add_argument("--pair", action="append", nargs=2, metavar=("A", "B"),
                     help="compare this specific pair (repeatable); names are input file stems")
    run.add_argument("--all-pairs", action="store_true", help="compare every pair")
    run.add_argument("--no-pangenome", action="store_true",
                     help="skip the joint clustering step (pairs still cluster on demand)")
    run.add_argument("--min-seq-id", type=float, default=0.5, help="MMseqs2 clustering identity")
    run.add_argument("--coverage", type=float, default=0.8, help="MMseqs2 clustering coverage")

    qc = sub.add_parser("qc", help="per-genome QC only: no clustering, no ANI")
    _add_common(qc)

    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    log = (lambda *a: None) if args.quiet else print

    inputs = [Path(x) for x in args.inputs]
    missing = [str(x) for x in inputs if not x.exists()]
    if missing:
        print(f"error: input not found: {', '.join(missing)}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    try:
        result = run_pipeline(
            inputs,
            outdir,
            lineage_path=args.lineage,
            pairs=[tuple(p) for p in (getattr(args, "pair", None) or [])] or None,
            all_pairs=getattr(args, "all_pairs", False),
            pangenome=(args.command == "run" and not getattr(args, "no_pangenome", False)),
            genetic_code=args.genetic_code,
            min_contig_length=args.min_contig_length,
            min_seq_id=getattr(args, "min_seq_id", 0.5),
            coverage=getattr(args, "coverage", 0.8),
            threads=args.threads,
            log=log,
        )
    except ToolMissing as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 3

    json_path = result.write_json(outdir / "genomex.json")
    html_path = write_html_report(result, outdir / "report.html")
    if result.runtime:
        result.runtime.write_provenance(outdir / "provenance.json")
    log(f"\nreport: {html_path}")
    log(f"json:   {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
