"""Command line entry point.

    python -m genomex run A.fna B.fna --outdir runs/demo
    python -m genomex run *.fna --outdir runs/all --all-pairs
    python -m genomex qc A.fna --outdir runs/qc
    python -m genomex proteins P.faa --lineage fungi_odb10 --outdir runs/prot
    python -m genomex lineages
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .markers import (
    DEFAULT_LINEAGE_NAME,
    LINEAGE_DB_ENV,
    available_lineages,
    lineage_db,
)
from .pipeline import run_pipeline, run_proteome_pipeline
from .proteome import NotAnAssembly, NotAProteome
from .report import write_html_report
from .runtime import ToolMissing


def _add_base(p: argparse.ArgumentParser, inputs_help: str) -> None:
    p.add_argument("inputs", nargs="+", help=inputs_help)
    p.add_argument("--outdir", "-o", default="runs/latest", help="output directory")
    p.add_argument(
        "--lineage",
        default=DEFAULT_LINEAGE_NAME,
        metavar="NAME|DIR",
        help=(
            "BUSCO odb10 marker set, by name or by directory (default: %(default)s). "
            "A name is looked up under $GENOMEX_DB; run `genomex lineages` to see "
            "what is installed and what the choice costs"
        ),
    )
    p.add_argument("--threads", type=int, default=None)
    p.add_argument("--quiet", action="store_true")


def _add_common(p: argparse.ArgumentParser) -> None:
    """Flags for the assembly path. Kept off `proteins`, where a translation
    table and a contig-length floor would both be accepted and ignored."""
    _add_base(p, "assembly FASTA files (.fna/.fa/.fasta, optionally .gz)")
    p.add_argument("--genetic-code", type=int, default=11, help="Prodigal translation table (11 bacteria)")
    p.add_argument("--min-contig-length", type=int, default=3000,
                   help="shortest contig entering composition analysis")


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

    prot = sub.add_parser(
        "proteins",
        help="completeness and duplication from a protein FASTA (no contamination, no comparison)",
        description=(
            "Score a supplied proteome against a BUSCO odb10 lineage. This is the "
            "whole of what a proteome supports: with no contigs there is no "
            "contamination analysis, no ANI and no gene order, and the report says "
            "so rather than showing them empty. Use it for lineages GenomeX cannot "
            "call genes for -- eukaryotes -- or when you already have a proteome you "
            "trust."
        ),
    )
    _add_base(prot, "protein FASTA files (.faa/.fa/.fasta, optionally .gz)")

    sub.add_parser(
        "lineages",
        help="list the installed BUSCO marker sets and what choosing one costs",
        description=(
            "The marker set decides what a completeness number means and how much "
            "contamination is detectable at all, and it was previously reachable "
            "only by knowing both that --lineage existed and where the directory "
            "lived. This lists what is installed. It does not recommend one: "
            "picking a set means knowing the organism, and the wrong set reports "
            "incompleteness that is not there."
        ),
    )

    return p


def _list_lineages() -> int:
    db = lineage_db()
    installed = available_lineages(db)
    if not installed:
        print(f"no marker sets found under {db}\n")
        print("Download one from https://busco-data.ezlab.org/v5/data/lineages/")
        print(f"and unpack it there, or point {LINEAGE_DB_ENV} at where yours live.")
        return 1

    width = max(len(i.name) for i in installed)
    print(f"marker sets under {db}\n")
    for i in installed:
        tag = "  <- default" if i.name == DEFAULT_LINEAGE_NAME else ""
        print(f"  {i.name:<{width}}  {i.n_markers:>4} markers{tag}")
    print(f"\nUse any of them by name:  --lineage {installed[-1].name}\n")
    print("What the choice costs, measured (docs/benchmark-mixture.md):")
    print("  A universal set is smaller, so a foreign genome displaces fewer of its")
    print("  markers. On constructed mixtures at 2% donor, a lineage-specific set")
    print("  detects 4 of 4 pairs and bacteria_odb10 detects 0 of 4. Against CheckM2")
    print("  at >=5%, recall is 0.00 at 124 markers and 0.25 at 688.")
    print("\n  The other direction is worse: a set from the wrong clade reports core")
    print("  genes missing that the organism never had. GenomeX will not infer a")
    print("  lineage from a genome, and neither should a script wrapping it.")
    return 0


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.command == "lineages":
        return _list_lineages()

    log = (lambda *a: None) if args.quiet else print

    inputs = [Path(x) for x in args.inputs]
    missing = [str(x) for x in inputs if not x.exists()]
    if missing:
        print(f"error: input not found: {', '.join(missing)}", file=sys.stderr)
        return 2

    outdir = Path(args.outdir)
    if args.command == "proteins":
        try:
            result = run_proteome_pipeline(
                inputs, outdir, lineage_path=args.lineage,
                threads=args.threads, log=log,
            )
        except NotAProteome as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 4
        except FileNotFoundError as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        except ToolMissing as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 3
        return _emit(result, outdir, log)

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
    except NotAnAssembly as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 4
    except FileNotFoundError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    except ToolMissing as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 3

    return _emit(result, outdir, log)


def _emit(result, outdir: Path, log) -> int:
    json_path = result.write_json(outdir / "genomex.json")
    html_path = write_html_report(result, outdir / "report.html")
    if result.runtime:
        result.runtime.write_provenance(outdir / "provenance.json")
    log(f"\nreport: {html_path}")
    log(f"json:   {json_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
