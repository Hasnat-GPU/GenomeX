"""Finding a marker set, and refusing to guess which one.

The marker set is the flag with the largest effect on what GenomeX reports, and
it was reachable only by typing an absolute path. These tests pin the two halves
of making it reachable: a bare name resolves, and *nothing else* does. A helper
that quietly turned `bacteria` into `bacteria_odb10` would be inferring a
lineage, which is the one thing this tool must never do -- a set from the wrong
clade reports core genes missing that the organism never had.
"""

from __future__ import annotations

import pytest

from genomex.markers import (
    DEFAULT_LINEAGE_NAME,
    Lineage,
    available_lineages,
    describe_lineage,
    resolve_lineage,
)
from genomex.pipeline import lineage_notice


def make_lineage_dir(root, name: str, n_markers: int, *, cfg_count: int | None = None):
    """A directory shaped like a BUSCO odb10 set, with `n_markers` profiles."""
    d = root / name
    (d / "hmms").mkdir(parents=True)
    for i in range(n_markers):
        (d / "hmms" / f"{i}at2.hmm").write_text("HMMER3/f\n", encoding="utf-8")
    count = n_markers if cfg_count is None else cfg_count
    (d / "dataset.cfg").write_text(
        f"name={name}\nspecies=fake\nnumber_of_BUSCOs={count}\n", encoding="utf-8"
    )
    return d


def _lineage(name: str, n: int) -> Lineage:
    from pathlib import Path
    return Lineage(path=Path("/db") / name, name=name, n_markers=n, scores={}, lengths={})


# ------------------------------------------------------------------ discovery

def test_available_lineages_orders_by_marker_count(tmp_path):
    make_lineage_dir(tmp_path, "fungi_odb10", 7)
    make_lineage_dir(tmp_path, "bacteria_odb10", 3)
    make_lineage_dir(tmp_path, "burkholderiales_odb10", 5)

    found = available_lineages(tmp_path)
    assert [i.name for i in found] == [
        "bacteria_odb10", "burkholderiales_odb10", "fungi_odb10"
    ]
    assert [i.n_markers for i in found] == [3, 5, 7]


def test_discovery_ignores_directories_that_are_not_lineage_sets(tmp_path):
    """`~/genomex-work/db` holds other things -- CheckM2's database lives beside
    the lineages. A directory with no `hmms/` is not a marker set."""
    make_lineage_dir(tmp_path, "bacteria_odb10", 3)
    (tmp_path / "checkm2").mkdir()
    (tmp_path / "notes.txt").write_text("x", encoding="utf-8")

    assert [i.name for i in available_lineages(tmp_path)] == ["bacteria_odb10"]


def test_a_missing_db_directory_lists_nothing_rather_than_raising(tmp_path):
    assert available_lineages(tmp_path / "nope") == []


def test_the_declared_marker_count_wins_over_the_file_count(tmp_path):
    """A set whose `dataset.cfg` disagrees with its `hmms/` directory is a partial
    download. Reporting the declared number makes that visible; counting the files
    would quietly report a complete-looking scan against a truncated set."""
    d = make_lineage_dir(tmp_path, "bacteria_odb10", 3, cfg_count=124)
    assert describe_lineage(d).n_markers == 124


# ----------------------------------------------------------------- resolution

def test_a_bare_name_resolves_against_the_db(tmp_path):
    made = make_lineage_dir(tmp_path, "fungi_odb10", 4)
    assert resolve_lineage("fungi_odb10", tmp_path) == made


def test_an_explicit_path_still_resolves(tmp_path):
    """The old spelling has to keep working: every script and doc that predates
    the name form passes a directory."""
    made = make_lineage_dir(tmp_path / "elsewhere", "fungi_odb10", 4)
    assert resolve_lineage(made, tmp_path) == made
    assert resolve_lineage(str(made), tmp_path) == made


def test_a_near_miss_name_is_refused_not_guessed(tmp_path):
    """`bacteria` is not `bacteria_odb10`, and prefix-matching it would be
    inference. The wrong clade reports fake incompleteness, so a name that does
    not exist is an error, never a best guess."""
    make_lineage_dir(tmp_path, "bacteria_odb10", 3)
    # Case is left out deliberately: whether `BACTERIA_ODB10` resolves is a
    # property of the filesystem, not of this code, and asserting it would fail
    # on Windows for a reason that is not a bug.
    for near_miss in ("bacteria", "bacteria_odb11", "odb10", "bacteria_odb10_v2"):
        with pytest.raises(FileNotFoundError):
            resolve_lineage(near_miss, tmp_path)


def test_the_refusal_names_what_is_actually_installed(tmp_path):
    """An error that only says "not found" leaves the user where they started:
    not knowing the flag takes a name or what names exist."""
    make_lineage_dir(tmp_path, "bacteria_odb10", 3)
    make_lineage_dir(tmp_path, "fungi_odb10", 7)

    with pytest.raises(FileNotFoundError) as exc:
        resolve_lineage("archaea_odb10", tmp_path)
    msg = str(exc.value)
    assert "archaea_odb10" in msg
    assert "bacteria_odb10" in msg and "fungi_odb10" in msg
    assert "7 markers" in msg


def test_an_empty_db_points_at_where_to_get_a_set(tmp_path):
    with pytest.raises(FileNotFoundError) as exc:
        resolve_lineage("bacteria_odb10", tmp_path)
    msg = str(exc.value)
    assert "busco-data.ezlab.org" in msg
    assert "GENOMEX_DB" in msg


# --------------------------------------------------------------- the run notice

def test_the_default_set_announces_what_it_costs():
    """A user who never finds `--lineage` still gets a verdict, and it looks the
    same as a good one. This is where they are told otherwise."""
    lines = lineage_notice(_lineage(DEFAULT_LINEAGE_NAME, 124))
    text = " ".join(lines)
    assert "124 markers" in text
    assert "0 of 4" in text and "4 of 4" in text
    assert "genomex lineages" in text


def test_a_chosen_set_is_named_and_left_alone():
    """The notice is for the default, which is a choice nobody made. Someone who
    passed `--lineage` has already made it, and does not need the lecture."""
    lines = lineage_notice(_lineage("burkholderiales_odb10", 688))
    assert lines == ["lineage: burkholderiales_odb10 (688 markers)"]
