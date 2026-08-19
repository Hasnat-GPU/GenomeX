"""Parsers for the CheckM2 comparison.

The deterministic-quantity check is only meaningful if the two tools' columns
are lined up correctly. A unit mismatch (GC as a fraction vs a percentage) would
silently report a 99% gap on every genome and look like a catastrophic bug in
GenomeX, so it is pinned here.
"""

import importlib.util
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bench_checkm2", ROOT / "bench" / "compare_to_checkm2.py")
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)

HEADER = (
    "Name\tCompleteness\tContamination\tCoding_Density\tContig_N50\tGenome_Size\t"
    "GC_Content\tTotal_Coding_Sequences\tTotal_Contigs\n"
)

GX = {
    "assembly": {"total_bp": 7_650_000, "n_contigs": 3, "gc_percent": 63.5, "n50": 4_152_217},
    "genes": {"n_genes": 6771, "coding_density": 0.88},
    "markers": {"completeness_percent": 100.0},
    "contamination": {"verdict": "clean"},
}


def _report(tmp_path, rows: list[str]) -> Path:
    p = tmp_path / "quality_report.tsv"
    p.write_text(HEADER + "\n".join(rows) + "\n")
    return p


def test_quality_report_is_keyed_by_genome_name(tmp_path):
    p = _report(
        tmp_path,
        [
            "g1\t100.0\t0.51\t0.88\t4152217\t7650000\t0.635\t6771\t3",
            "g2\t98.2\t7.30\t0.87\t100282\t9030000\t0.625\t8100\t320",
        ],
    )
    got = bench.load_checkm2(p)
    assert sorted(got) == ["g1", "g2"]
    assert got["g2"]["Contamination"] == "7.30"


def test_gc_expressed_as_a_fraction_is_converted_to_percent(tmp_path):
    p = _report(tmp_path, ["g1\t100.0\t0.5\t0.88\t4152217\t7650000\t0.635\t6771\t3"])
    cm = bench.load_checkm2(p)["g1"]
    _, checkm2_gc = bench.stat_pairs(GX, cm)["gc_percent"]
    assert checkm2_gc == 63.5


def test_gc_already_in_percent_is_left_alone(tmp_path):
    p = _report(tmp_path, ["g1\t100.0\t0.5\t0.88\t4152217\t7650000\t63.5\t6771\t3"])
    cm = bench.load_checkm2(p)["g1"]
    _, checkm2_gc = bench.stat_pairs(GX, cm)["gc_percent"]
    assert checkm2_gc == 63.5


def test_matching_deterministic_stats_show_no_gap(tmp_path):
    p = _report(tmp_path, ["g1\t100.0\t0.5\t0.88\t4152217\t7650000\t0.635\t6771\t3"])
    cm = bench.load_checkm2(p)["g1"]
    for stat, (x, y) in bench.stat_pairs(GX, cm).items():
        assert bench.relative_gap(x, y) == 0.0, f"{stat} should agree exactly"


def test_relative_gap_is_symmetric_and_scale_free():
    assert bench.relative_gap(100, 100) == 0.0
    assert bench.relative_gap(100, 90) == bench.relative_gap(90, 100)
    assert bench.relative_gap(100, 90) == 0.1
    assert math.isnan(bench.relative_gap(None, 5))


def test_missing_columns_do_not_crash_the_comparison(tmp_path):
    p = tmp_path / "quality_report.tsv"
    p.write_text("Name\tCompleteness\tContamination\ng1\t99.0\t1.2\n")
    cm = bench.load_checkm2(p)["g1"]
    pairs = bench.stat_pairs(GX, cm)
    assert pairs["genome_size"][1] is None
    assert math.isnan(bench.relative_gap(*pairs["genome_size"]))


# --------------------------------------------------- end-to-end page build

def _genomex_json(tmp_path, genomes: list[dict]) -> Path:
    import json
    d = tmp_path / "gx"
    d.mkdir()
    (d / "genomex.json").write_text(json.dumps({"genomes": genomes}), encoding="utf-8")
    return d


def _gx_entry(name: str, lineage: str | None = None, n_markers: int | None = None) -> dict:
    import copy
    g = copy.deepcopy(GX)
    g["genome"] = name
    if lineage:
        g["markers"]["lineage"] = lineage
        g["markers"]["markers_total"] = n_markers
    return g


def test_report_names_the_marker_set_it_used(tmp_path):
    """A contamination result quoted without its lineage cannot be compared to
    another one: 124 markers quantise at 0.81%, 688 at 0.15%."""
    gx = _genomex_json(tmp_path, [_gx_entry("g1", "burkholderiales_odb10", 688)])
    cm = _report(tmp_path, ["g1\t100.0\t0.5\t0.88\t4152217\t7650000\t63.5\t6771\t3"])
    out = tmp_path / "page.md"

    assert bench.main(["--genomex", str(gx), "--checkm2", str(cm), "--out", str(out)]) == 0
    text = out.read_text(encoding="utf-8")
    assert "burkholderiales_odb10" in text
    assert "688 markers" in text


def test_report_builds_when_the_lineage_is_absent(tmp_path):
    """Older run directories predate the field. Omitting the claim is right;
    crashing after a 90-minute sweep is not."""
    gx = _genomex_json(tmp_path, [_gx_entry("g1")])
    cm = _report(tmp_path, ["g1\t100.0\t0.5\t0.88\t4152217\t7650000\t63.5\t6771\t3"])
    out = tmp_path / "page.md"

    assert bench.main(["--genomex", str(gx), "--checkm2", str(cm), "--out", str(out)]) == 0
    assert "1 genomes compared." in out.read_text(encoding="utf-8")
