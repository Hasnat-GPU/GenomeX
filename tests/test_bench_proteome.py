"""The proteome comparison is a measuring instrument, so it is tested like one.

Its output is an agreement count, and the failure that matters is silent: a
parser that reads the wrong column, or drops rows, reports *perfect* agreement
rather than an error. `compare_to_busco.py` once printed a whole table of `None`
for exactly that reason. These pin the two parsers against hand-built tables
whose answer is known.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))

from compare_proteome_to_busco import load_busco_full_table, load_genomex_table

BUSCO_TABLE = """\
# BUSCO version is: 5.8.3
# Busco id\tStatus\tSequence\tScore\tLength\tOrthoDB url\tDescription
m1\tComplete\tNP_1.1\t520.2\t726\thttp://x\tsomething
m2\tDuplicated\tNP_2.1\t485.6\t383\thttp://x\tsomething
m2\tDuplicated\tNP_3.1\t401.2\t377\thttp://x\tsomething
m3\tFragmented\tNP_4.1\t96.9\t120\thttp://x\tsomething
m4\tMissing\t-\t-\t-\thttp://x\tsomething
"""

GX_TABLE = """\
busco_id\tstatus\tprotein_id\tcontig\tscore\tmatched_aa
m1\tComplete\tNP_1.1\tNP\t520.2\t726
m2\tDuplicated\tNP_2.1\tNP\t485.6\t383
m2\tDuplicated\tNP_3.1\tNP\t401.2\t377
m3\tFragmented\tNP_4.1\tNP\t96.9\t120
m4\tMissing\t-\t-\t-\t-
"""


def test_busco_table_parses_status_and_every_sequence(tmp_path):
    run = tmp_path / "run_fungi_odb10"
    run.mkdir()
    (run / "full_table.tsv").write_text(BUSCO_TABLE, encoding="utf-8")
    status, seqs = load_busco_full_table(tmp_path)
    assert status == {"m1": "Complete", "m2": "Duplicated",
                      "m3": "Fragmented", "m4": "Missing"}
    # Both copies of the duplicated marker, or duplication is silently halved.
    assert sorted(seqs["m2"]) == ["NP_2.1", "NP_3.1"]
    assert "m4" not in seqs, "a Missing marker must name no sequence"


def test_genomex_table_parses_the_same_way(tmp_path):
    p = tmp_path / "markers.tsv"
    p.write_text(GX_TABLE, encoding="utf-8")
    status, seqs = load_genomex_table(p)
    assert status == {"m1": "Complete", "m2": "Duplicated",
                      "m3": "Fragmented", "m4": "Missing"}
    assert sorted(seqs["m2"]) == ["NP_2.1", "NP_3.1"]
    assert "m4" not in seqs


def test_the_two_parsers_agree_on_identical_content(tmp_path):
    """The comparison is only meaningful if the vocabularies really do line up."""
    run = tmp_path / "run_x"
    run.mkdir()
    (run / "full_table.tsv").write_text(BUSCO_TABLE, encoding="utf-8")
    gx = tmp_path / "markers.tsv"
    gx.write_text(GX_TABLE, encoding="utf-8")
    b_status, b_seqs = load_busco_full_table(tmp_path)
    g_status, g_seqs = load_genomex_table(gx)
    assert b_status == g_status
    assert {k: sorted(v) for k, v in b_seqs.items()} == {k: sorted(v) for k, v in g_seqs.items()}


def test_a_missing_full_table_is_an_error_not_an_empty_agreement(tmp_path):
    """Zero markers on both sides would otherwise read as 100% agreement."""
    import pytest
    with pytest.raises(SystemExit):
        load_busco_full_table(tmp_path)
