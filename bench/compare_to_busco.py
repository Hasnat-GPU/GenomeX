"""Marker-by-marker comparison of GenomeX against real BUSCO.

GenomeX calls its completeness scan "BUSCO-compatible". That is a claim, and this
script is the test of it: run both tools over the same genomes with the same
lineage, then compare not the headline percentages but the classification of each
individual marker.

Aggregate agreement is the weaker check and can hide real disagreement -- two
tools both reporting 100% complete may still disagree on which markers are
duplicated, which is precisely the number the contamination module consumes.

Comparing against BUSCO in two modes separates two different disagreements:

    genome    BUSCO predicts its own genes. End-to-end: what a user experiences.
    proteins  BUSCO scores GenomeX's gene calls. Isolates marker classification,
              so a difference here is about cutoffs, not about Prodigal.

Usage:
    python bench/compare_to_busco.py --genomex runs/demo-rhizobia \\
        --busco genome=bench/busco_out/genome \\
        --busco proteins=bench/busco_out/proteins \\
        --out docs/benchmark-busco.md
"""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

CLASSES = ["Complete", "Duplicated", "Fragmented", "Missing"]

MODE_BLURB = {
    "genome": (
        "BUSCO ran in `-m genome` mode and predicted its own genes. Every "
        "difference below is end-to-end: gene calling and marker classification "
        "combined, which is what a user choosing between the two tools gets."
    ),
    "proteins": (
        "BUSCO ran in `-m proteins` mode over GenomeX's own predicted proteome, "
        "so both tools saw identical genes. Differences here are attributable to "
        "marker classification alone -- score and length cutoffs, domain "
        "envelope handling, and the complete/fragmented boundary."
    ),
}


def load_genomex_markers(run_dir: Path) -> dict[str, dict[str, str]]:
    """genome -> {busco_id: status} from each genomes/<name>/markers.tsv."""
    out: dict[str, dict[str, str]] = {}
    for tsv in sorted(run_dir.glob("genomes/*/markers.tsv")):
        genome = tsv.parent.name
        status: dict[str, str] = {}
        for line in tsv.read_text(encoding="utf-8").splitlines()[1:]:
            if not line.strip():
                continue
            bid, st = line.split("\t")[:2]
            status[bid] = st
        out[genome] = status
    return out


def load_busco_markers(busco_dir: Path) -> dict[str, dict[str, str]]:
    """genome -> {busco_id: status} from each BUSCO run's full_table.tsv.

    BUSCO writes one row per copy for duplicated markers, so a marker seen
    Complete twice is Duplicated regardless of which row is read last.
    """
    out: dict[str, dict[str, str]] = {}
    for table in sorted(busco_dir.glob("*/run_*/full_table.tsv")):
        genome = table.parent.parent.name
        status: dict[str, str] = {}
        for line in table.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                continue
            bid, st = fields[0], fields[1]
            if bid in status and status[bid] == "Complete" and st == "Complete":
                st = "Duplicated"
            status[bid] = st
        out[genome] = status
    return out


def load_busco_summary(busco_dir: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    for js in sorted(busco_dir.glob("*/short_summary.*.json")):
        genome = js.parent.name
        data = json.loads(js.read_text(encoding="utf-8"))
        out[genome] = data.get("results", data)
    return out


def load_genomex_summary(run_dir: Path) -> dict[str, dict]:
    data = json.loads((run_dir / "genomex.json").read_text(encoding="utf-8"))
    return {g["genome"]: g["markers"] for g in data["genomes"]}


def confusion(gx: dict[str, str], bu: dict[str, str]) -> tuple[Counter, list[tuple]]:
    """Counter of (genomex_class, busco_class) plus the list of disagreements."""
    pairs: Counter = Counter()
    diffs: list[tuple] = []
    for bid in sorted(set(gx) | set(bu)):
        a, b = gx.get(bid, "Absent"), bu.get(bid, "Absent")
        pairs[(a, b)] += 1
        if a != b:
            diffs.append((bid, a, b))
    return pairs, diffs


def _pct(n: int, d: int) -> str:
    return f"{100.0 * n / d:.1f}%" if d else "-"


def _headline_rows(
    genomes: list[str], gx_sum: dict[str, dict], bu_sum: dict[str, dict], label: str
) -> list[str]:
    lines = ["| genome | tool | C | S | D | F | M | completeness |", "|---|---|---|---|---|---|---|---|"]
    for g in genomes:
        gs, bs = gx_sum.get(g, {}), bu_sum.get(g, {})
        n = gs.get("markers_total", 124)
        lines.append(
            f"| {g} | GenomeX | {gs.get('complete', '-')} | {gs.get('single_copy', '-')} | "
            f"{gs.get('duplicated', '-')} | {gs.get('fragmented', '-')} | {gs.get('missing', '-')} | "
            f"{gs.get('completeness_percent', '-')}% |"
        )
        c = bs.get("Complete BUSCOs")
        pct = f"{bs.get('Complete percentage')}%" if c is not None else "-"
        lines.append(
            f"| {g} | BUSCO ({label}) | {c} | {bs.get('Single copy BUSCOs')} | "
            f"{bs.get('Multi copy BUSCOs')} | {bs.get('Fragmented BUSCOs')} | "
            f"{bs.get('Missing BUSCOs')} | {pct} |"
        )
    return lines


def compare_mode(
    label: str,
    gx_markers: dict[str, dict[str, str]],
    gx_sum: dict[str, dict],
    bu_markers: dict[str, dict[str, str]],
    bu_sum: dict[str, dict],
) -> tuple[list[str], dict]:
    genomes = sorted(set(gx_markers) & set(bu_markers))
    lines: list[str] = [f"## BUSCO `-m {label}`", ""]
    lines.append(MODE_BLURB.get(label, ""))
    lines.append("")

    if not genomes:
        lines.append("_No genomes in common._")
        return lines, {}

    lines.append("### Headline numbers")
    lines.append("")
    lines.extend(_headline_rows(genomes, gx_sum, bu_sum, label))
    lines.append("")
    lines.append("Side by side, in BUSCO's own summary notation:")
    lines.append("")
    lines.append("| genome | GenomeX | BUSCO |")
    lines.append("|---|---|---|")
    for g in genomes:
        gx_str = gx_sum.get(g, {}).get("busco_style_string", "-")
        bu_str = bu_sum.get(g, {}).get("one_line_summary", "-")
        lines.append(f"| {g} | `{gx_str}` | `{bu_str}` |")
    lines.append("")

    total_pairs: Counter = Counter()
    all_diffs: dict[str, list[tuple]] = {}
    lines.append("### Per-marker agreement")
    lines.append("")
    lines.append("| genome | markers | identical calls | agreement |")
    lines.append("|---|---|---|---|")
    for g in genomes:
        pairs, diffs = confusion(gx_markers[g], bu_markers[g])
        total_pairs += pairs
        all_diffs[g] = diffs
        n = sum(pairs.values())
        same = sum(v for (a, b), v in pairs.items() if a == b)
        lines.append(f"| {g} | {n} | {same} | {_pct(same, n)} |")
    n_all = sum(total_pairs.values())
    same_all = sum(v for (a, b), v in total_pairs.items() if a == b)
    lines.append(f"| **all** | **{n_all}** | **{same_all}** | **{_pct(same_all, n_all)}** |")
    lines.append("")

    observed = sorted({a for a, _ in total_pairs} | {b for _, b in total_pairs})
    order = [c for c in CLASSES if c in observed] + [c for c in observed if c not in CLASSES]
    lines.append("### Confusion matrix, pooled")
    lines.append("")
    lines.append("Rows are GenomeX, columns BUSCO.")
    lines.append("")
    lines.append("| GenomeX \\ BUSCO | " + " | ".join(order) + " |")
    lines.append("|---" * (len(order) + 1) + "|")
    for a in order:
        lines.append(f"| **{a}** | " + " | ".join(str(total_pairs.get((a, b), 0)) for b in order) + " |")
    lines.append("")

    kinds = Counter((a, b) for diffs in all_diffs.values() for _, a, b in diffs)
    n_diff = sum(kinds.values())
    lines.append("### Where they disagree")
    lines.append("")
    if not n_diff:
        lines.append("No disagreements on any marker in any genome.")
    else:
        lines.append("| GenomeX called | BUSCO called | markers |")
        lines.append("|---|---|---|")
        for (a, b), n in kinds.most_common():
            lines.append(f"| {a} | {b} | {n} |")
        lines.append("")
        lines.append("<details><summary>Every disagreeing marker</summary>")
        lines.append("")
        lines.append("| genome | marker | GenomeX | BUSCO |")
        lines.append("|---|---|---|---|")
        for g in genomes:
            for bid, a, b in all_diffs[g]:
                lines.append(f"| {g} | `{bid}` | {a} | {b} |")
        lines.append("")
        lines.append("</details>")
    lines.append("")

    stats = {
        "genomes": len(genomes),
        "markers": n_all,
        "identical": same_all,
        "agreement": _pct(same_all, n_all),
        "disagreements": n_diff,
        "kinds": kinds,
    }
    return lines, stats


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--genomex", required=True, type=Path, help="a GenomeX --outdir")
    ap.add_argument(
        "--busco",
        required=True,
        action="append",
        metavar="LABEL=DIR",
        help="labelled BUSCO output directory, repeatable "
             "(genome=bench/busco_out/genome proteins=bench/busco_out/proteins)",
    )
    ap.add_argument("--out", type=Path, default=Path("docs/benchmark-busco.md"))
    ap.add_argument("--lineage", default="bacteria_odb10")
    args = ap.parse_args()

    modes: dict[str, Path] = {}
    for spec in args.busco:
        label, _, path = spec.partition("=")
        if not path:
            label, path = "busco", spec
        modes[label] = Path(path)

    gx_markers = load_genomex_markers(args.genomex)
    gx_sum = load_genomex_summary(args.genomex)
    if not gx_markers:
        print(f"no markers.tsv under {args.genomex}/genomes -- rerun the pipeline")
        return 1

    header: list[str] = ["# GenomeX vs BUSCO", ""]
    header.append(
        "GenomeX describes its completeness scan as *BUSCO-compatible*: HMMER run "
        "directly against the OrthoDB v10 profiles BUSCO ships, with BUSCO's own "
        "score and length cutoffs applied. This page tests that description against "
        f"BUSCO itself on `{args.lineage}`, marker by marker. Both tools read the "
        "same local lineage directory, so no disagreement here comes from profile "
        "or cutoff versions."
    )
    header.append("")
    header.append("Generated by `bench/compare_to_busco.py`; regenerate with `bench/run_busco.sh`.")
    header.append("")

    body: list[str] = []
    summary: dict[str, dict] = {}
    for label in sorted(modes, key=lambda x: (x != "proteins", x)):
        bm, bs = load_busco_markers(modes[label]), load_busco_summary(modes[label])
        if not bm:
            print(f"warning: no BUSCO output under {modes[label]} for mode {label!r}, skipping")
            continue
        lines, stats = compare_mode(label, gx_markers, gx_sum, bm, bs)
        body.extend(lines)
        if stats:
            summary[label] = stats

    if not summary:
        print("no BUSCO output found in any given directory")
        return 1

    if len(summary) > 1:
        header.append("## Summary")
        header.append("")
        header.append("| BUSCO mode | markers compared | identical calls | agreement |")
        header.append("|---|---|---|---|")
        for label, s in summary.items():
            header.append(f"| `{label}` | {s['markers']} | {s['identical']} | {s['agreement']} |")
        header.append("")

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text("\n".join(header + body) + "\n", encoding="utf-8")

    for label, s in summary.items():
        print(f"[{label}] {s['genomes']} genomes, {s['markers']} markers, "
              f"{s['identical']} identical ({s['agreement']}), {s['disagreements']} differ")
        for (a, b), n in s["kinds"].most_common(6):
            print(f"    GenomeX {a:<11} vs BUSCO {b:<11} {n}")
    print(f"written: {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
