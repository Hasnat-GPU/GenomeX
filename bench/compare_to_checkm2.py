"""Compare GenomeX against CheckM2 over a whole collection.

Two different comparisons live here, and conflating them would be dishonest.

**Deterministic quantities** -- genome size, GC, N50, contig count, CDS count,
coding density -- are computed from the same FASTA by both tools. These must
agree to within rounding. Any real gap is a bug in GenomeX, and this is the
strongest independent check on `fasta.py` and `genes.py` that exists here.

**Judgement calls** -- completeness and contamination -- are produced by
different methods on different bases. CheckM2 uses a machine-learning model over
a reference proteome database; GenomeX uses BUSCO odb10 single-copy markers plus
composition. They should broadly agree, but exact equality is not expected and
would be suspicious. What matters is whether GenomeX flags the genomes CheckM2
considers contaminated, and stays quiet on the ones it does not.

Usage:
    python bench/compare_to_checkm2.py --genomex ~/genomex-work/runs/sweep \\
        --checkm2 ~/genomex-work/checkm2_out/quality_report.tsv \\
        --out docs/benchmark-contamination.md
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

#: MIMAG-style thresholds. High-quality drafts are <5% contamination;
#: medium-quality allows <10%.
CONTAMINATION_THRESHOLDS = (5.0, 10.0)

#: Deterministic quantities and the relative tolerance each must agree within.
STAT_TOLERANCE = {
    "genome_size": 0.0,       # identical bases in, identical total out
    "total_contigs": 0.0,
    "gc_percent": 0.5,        # CheckM2 rounds; allow half a point
    "contig_n50": 0.0,
    "total_cds": 0.02,        # both call genes with Prodigal, settings may differ
    "coding_density": 0.05,
}


def load_checkm2(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row.get("Name") or row.get("name")
            if not name:
                continue
            out[name] = row
    return out


def load_genomex(run_dir: Path) -> dict[str, dict]:
    data = json.loads((run_dir / "genomex.json").read_text(encoding="utf-8"))
    return {g["genome"]: g for g in data["genomes"]}


def _f(value, default=None):
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def stat_pairs(gx: dict, cm: dict) -> dict[str, tuple[float, float]]:
    """The quantities both tools compute from the FASTA alone."""
    a, genes = gx["assembly"], gx["genes"]
    return {
        "genome_size": (a["total_bp"], _f(cm.get("Genome_Size"))),
        "total_contigs": (a["n_contigs"], _f(cm.get("Total_Contigs"))),
        "gc_percent": (a["gc_percent"], _f(cm.get("GC_Content", 0)) * 100
                       if _f(cm.get("GC_Content"), 1) <= 1 else _f(cm.get("GC_Content"))),
        "contig_n50": (a["n50"], _f(cm.get("Contig_N50"))),
        "total_cds": (genes["n_genes"], _f(cm.get("Total_Coding_Sequences"))),
        "coding_density": (genes["coding_density"], _f(cm.get("Coding_Density"))),
    }


def relative_gap(x: float, y: float) -> float:
    if x is None or y is None:
        return math.nan
    if x == y:
        return 0.0
    denom = max(abs(x), abs(y), 1e-9)
    return abs(x - y) / denom


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomex", required=True, type=Path, help="a GenomeX --outdir")
    ap.add_argument("--checkm2", required=True, type=Path, help="CheckM2 quality_report.tsv")
    ap.add_argument("--out", type=Path, default=Path("docs/benchmark-contamination.md"))
    args = ap.parse_args(argv)

    gx_all = load_genomex(args.genomex)
    cm_all = load_checkm2(args.checkm2)
    shared = sorted(set(gx_all) & set(cm_all))
    if not shared:
        print(f"no genomes in common.\n  genomex: {len(gx_all)}\n  checkm2: {len(cm_all)}")
        return 1

    lines: list[str] = ["# GenomeX vs CheckM2", ""]
    lines.append(
        "CheckM2 estimates completeness and contamination with a machine-learning "
        "model over a reference proteome database. GenomeX uses BUSCO odb10 "
        "single-copy markers plus tetranucleotide composition. The two are not "
        "expected to produce identical numbers, and identical numbers would be "
        "suspicious. This page separates what must agree from what need not."
    )
    lines.append("")
    # Which marker set was used is not a detail. bacteria_odb10 has 124 markers,
    # so one duplicated marker is 0.81% of the core and a 5% contaminant may
    # displace none at all; burkholderiales_odb10 has 688 and quantises at 0.15%.
    # A contamination result quoted without its lineage cannot be compared to
    # another one.
    lineages = sorted({
        g["markers"].get("lineage") for g in gx_all.values() if g["markers"].get("lineage")
    })
    sizes = sorted({
        g["markers"].get("markers_total") for g in gx_all.values()
        if g["markers"].get("markers_total")
    })
    scanned = (
        f", scanned against `{', '.join(lineages)}` "
        f"({', '.join(str(s) for s in sizes)} markers)" if lineages else ""
    )
    lines.append(
        f"{len(shared)} genomes compared{scanned}. "
        "Generated by `bench/compare_to_checkm2.py`."
    )
    lines.append("")

    # ---- deterministic quantities ----------------------------------------
    lines.append("## Quantities that must agree")
    lines.append("")
    lines.append(
        "Both tools read the same FASTA, so these are not estimates. A gap here "
        "is a bug, not a difference of method."
    )
    lines.append("")
    lines.append("| quantity | genomes | max relative gap | tolerance | verdict |")
    lines.append("|---|---|---|---|---|")
    stat_failures: dict[str, list[tuple[str, float, float]]] = defaultdict(list)
    for stat, tol in STAT_TOLERANCE.items():
        gaps = []
        for g in shared:
            x, y = stat_pairs(gx_all[g], cm_all[g])[stat]
            gap = relative_gap(x, y)
            if not math.isnan(gap):
                gaps.append(gap)
                if gap > tol:
                    stat_failures[stat].append((g, x, y))
        worst = max(gaps) if gaps else math.nan
        ok = not stat_failures[stat]
        lines.append(
            f"| {stat} | {len(gaps)} | {worst:.4f} | {tol} | "
            f"{'agree' if ok else f'**{len(stat_failures[stat])} exceed tolerance**'} |"
        )
    lines.append("")
    if any(stat_failures.values()):
        lines.append("<details><summary>Genomes exceeding tolerance</summary>")
        lines.append("")
        lines.append("| quantity | genome | GenomeX | CheckM2 |")
        lines.append("|---|---|---|---|")
        for stat, rows in stat_failures.items():
            for g, x, y in rows[:40]:
                lines.append(f"| {stat} | {g} | {x} | {y} |")
        lines.append("")
        lines.append("</details>")
        lines.append("")

    # ---- completeness ----------------------------------------------------
    diffs = []
    for g in shared:
        a = gx_all[g]["markers"]["completeness_percent"]
        b = _f(cm_all[g].get("Completeness"))
        if b is not None:
            diffs.append((g, a, b, a - b))
    lines.append("## Completeness")
    lines.append("")
    if diffs:
        deltas = [d for _, _, _, d in diffs]
        mean_abs = sum(abs(d) for d in deltas) / len(deltas)
        lines.append(
            f"Mean absolute difference **{mean_abs:.2f} points** over {len(deltas)} genomes "
            f"(GenomeX minus CheckM2: min {min(deltas):+.2f}, max {max(deltas):+.2f})."
        )
        lines.append("")
        lines.append("| difference | genomes |")
        lines.append("|---|---|")
        buckets = Counter()
        for d in deltas:
            if abs(d) <= 1:
                buckets["within 1 point"] += 1
            elif abs(d) <= 5:
                buckets["1-5 points"] += 1
            elif abs(d) <= 10:
                buckets["5-10 points"] += 1
            else:
                buckets["over 10 points"] += 1
        for k in ("within 1 point", "1-5 points", "5-10 points", "over 10 points"):
            if buckets[k]:
                lines.append(f"| {k} | {buckets[k]} |")
        lines.append("")
        worst = sorted(diffs, key=lambda r: -abs(r[3]))[:10]
        lines.append("Largest differences:")
        lines.append("")
        lines.append("| genome | GenomeX | CheckM2 | difference |")
        lines.append("|---|---|---|---|")
        for g, a, b, d in worst:
            lines.append(f"| {g} | {a}% | {b}% | {d:+.2f} |")
        lines.append("")

    # ---- contamination ---------------------------------------------------
    lines.append("## Contamination")
    lines.append("")
    lines.append(
        "GenomeX returns a category with evidence, CheckM2 a percentage, so the "
        "honest comparison is whether the category tracks the percentage."
    )
    lines.append("")
    lines.append("| GenomeX verdict | genomes | CheckM2 contamination: median | min | max |")
    lines.append("|---|---|---|---|---|")
    by_verdict: dict[str, list[float]] = defaultdict(list)
    for g in shared:
        v = gx_all[g]["contamination"]["verdict"]
        c = _f(cm_all[g].get("Contamination"))
        if c is not None:
            by_verdict[v].append(c)
    for v in ("clean", "possible", "likely", "undetermined"):
        vals = sorted(by_verdict.get(v, []))
        if not vals:
            continue
        med = vals[len(vals) // 2]
        lines.append(f"| {v} | {len(vals)} | {med:.2f}% | {vals[0]:.2f}% | {vals[-1]:.2f}% |")
    lines.append("")

    for thr in CONTAMINATION_THRESHOLDS:
        tp = fp = fn = tn = 0
        missed: list[tuple[str, float]] = []
        for g in shared:
            c = _f(cm_all[g].get("Contamination"))
            if c is None:
                continue
            flagged = gx_all[g]["contamination"]["verdict"] in ("possible", "likely")
            contaminated = c >= thr
            if contaminated and flagged:
                tp += 1
            elif contaminated and not flagged:
                fn += 1
                missed.append((g, c))
            elif not contaminated and flagged:
                fp += 1
            else:
                tn += 1
        lines.append(f"### Against a CheckM2 threshold of {thr}%")
        lines.append("")
        lines.append("| | CheckM2 contaminated | CheckM2 clean |")
        lines.append("|---|---|---|")
        lines.append(f"| **GenomeX flagged** | {tp} | {fp} |")
        lines.append(f"| **GenomeX quiet** | {fn} | {tn} |")
        lines.append("")
        recall = tp / (tp + fn) if (tp + fn) else float("nan")
        prec = tp / (tp + fp) if (tp + fp) else float("nan")
        lines.append(
            f"Recall {recall:.2f} ({tp}/{tp + fn} caught), precision {prec:.2f} "
            f"({tp}/{tp + fp} flags justified at this threshold)."
        )
        if missed:
            lines.append("")
            lines.append("Missed by GenomeX:")
            lines.append("")
            lines.append("| genome | CheckM2 contamination |")
            lines.append("|---|---|")
            for g, c in sorted(missed, key=lambda r: -r[1])[:20]:
                lines.append(f"| {g} | {c:.2f}% |")
        lines.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"genomes compared: {len(shared)}")
    for stat in STAT_TOLERANCE:
        n = len(stat_failures[stat])
        print(f"  {stat:<16} {'OK' if not n else f'{n} exceed tolerance'}")
    if diffs:
        deltas = [d for _, _, _, d in diffs]
        print(f"  completeness     mean abs diff {sum(abs(d) for d in deltas) / len(deltas):.2f} points")
    print(f"  verdicts: {dict(Counter(gx_all[g]['contamination']['verdict'] for g in shared))}")
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
