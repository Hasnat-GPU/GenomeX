"""The benchmark script produces a public claim, so its parsers are tested too.

A silent parsing bug here would make GenomeX look either better or worse than it
is against BUSCO, which is worse than having no benchmark at all.
"""

import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("bench_compare", ROOT / "bench" / "compare_to_busco.py")
bench = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bench)


FULL_TABLE_HEADER = (
    "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
)


def _busco_run(tmp_path, genome: str, rows: list[str]) -> Path:
    d = tmp_path / genome / "run_bacteria_odb10"
    d.mkdir(parents=True)
    (d / "full_table.tsv").write_text(FULL_TABLE_HEADER + "\n".join(rows) + "\n")
    return tmp_path


def test_busco_duplicated_markers_are_recognised_from_repeated_rows(tmp_path):
    """BUSCO writes one row per copy; two Complete rows mean Duplicated."""
    _busco_run(
        tmp_path,
        "g1",
        [
            "m1\tComplete\tctg1\t1\t900\t+\t250.0\t300",
            "m2\tComplete\tctg1\t1\t900\t+\t250.0\t300",
            "m2\tComplete\tctg7\t1\t900\t+\t249.0\t300",
            "m3\tFragmented\tctg1\t1\t400\t+\t80.0\t130",
            "m4\tMissing",
        ],
    )
    got = bench.load_busco_markers(tmp_path)["g1"]
    assert got == {
        "m1": "Complete",
        "m2": "Duplicated",
        "m3": "Fragmented",
        "m4": "Missing",
    }


def test_busco_native_duplicated_label_is_preserved(tmp_path):
    _busco_run(
        tmp_path,
        "g1",
        [
            "m1\tDuplicated\tctg1\t1\t900\t+\t250.0\t300",
            "m1\tDuplicated\tctg2\t1\t900\t+\t250.0\t300",
        ],
    )
    assert bench.load_busco_markers(tmp_path)["g1"] == {"m1": "Duplicated"}


def test_genomex_marker_table_is_read_back(tmp_path):
    d = tmp_path / "genomes" / "g1"
    d.mkdir(parents=True)
    (d / "markers.tsv").write_text(
        "busco_id\tstatus\tprotein_id\tcontig\tscore\tmatched_aa\n"
        "m1\tComplete\tctg1_5\tctg1\t250.0\t300\n"
        "m2\tDuplicated\tctg1_9\tctg1\t250.0\t300\n"
        "m2\tDuplicated\tctg7_2\tctg7\t249.0\t300\n"
        "m3\tMissing\t-\t-\t-\t-\n"
    )
    got = bench.load_genomex_markers(tmp_path)["g1"]
    assert got == {"m1": "Complete", "m2": "Duplicated", "m3": "Missing"}


def test_confusion_counts_pairs_and_lists_only_disagreements():
    gx = {"m1": "Complete", "m2": "Duplicated", "m3": "Missing", "m4": "Complete"}
    bu = {"m1": "Complete", "m2": "Complete", "m3": "Missing", "m4": "Fragmented"}
    pairs, diffs = bench.confusion(gx, bu)

    assert pairs[("Complete", "Complete")] == 1
    assert pairs[("Duplicated", "Complete")] == 1
    assert pairs[("Missing", "Missing")] == 1
    assert sorted(d[0] for d in diffs) == ["m2", "m4"]
    assert sum(pairs.values()) == 4


def test_marker_present_in_only_one_tool_is_not_silently_dropped():
    pairs, diffs = bench.confusion({"m1": "Complete"}, {"m2": "Complete"})
    assert pairs[("Complete", "Absent")] == 1
    assert pairs[("Absent", "Complete")] == 1
    assert len(diffs) == 2, "a marker missing from one side must count as a disagreement"


def test_busco_summary_json_keys_are_the_ones_busco_writes(tmp_path):
    """Guards against silently reading the wrong keys and printing empty cells.

    BUSCO writes 'Complete BUSCOs', not 'Complete'; the first version of this
    script read the latter and produced a table full of None.
    """
    d = tmp_path / "g1"
    d.mkdir()
    (d / "short_summary.specific.bacteria_odb10.g1.json").write_text(
        '{"results": {"one_line_summary": "C:99.2%[S:97.6%,D:1.6%],F:0.8%,M:0.0%,n:124",'
        ' "Complete BUSCOs": 123, "Complete percentage": 99.2, "Single copy BUSCOs": 121,'
        ' "Multi copy BUSCOs": 2, "Fragmented BUSCOs": 1, "Missing BUSCOs": 0, "n_markers": 124}}'
    )
    got = bench.load_busco_summary(tmp_path)["g1"]
    assert got["Complete BUSCOs"] == 123

    rows = bench._headline_rows(["g1"], {"g1": {"markers_total": 124}}, {"g1": got}, "proteins")
    busco_row = [r for r in rows if "BUSCO" in r][0]
    assert "123" in busco_row and "121" in busco_row
    assert "None" not in busco_row, "empty cells mean the key mapping drifted"
