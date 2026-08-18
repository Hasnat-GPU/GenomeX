"""Marker classification: BUSCO's cutoffs applied to a HMMER table.

The domtblout rows here are hand-built so the expected answer is known exactly:
which markers are single, duplicated, fragmented, and missing.
"""

from genomex.markers import Lineage, MarkerResult, parse_domtbl, sum_hmm_len

# Two markers. m1 expects a 100 aa protein (sigma 5 -> complete above 90 aa),
# m2 expects 200 aa (sigma 10 -> complete above 180 aa).
LINEAGE = Lineage(
    path=__import__("pathlib").Path("."),
    name="test_odb10",
    n_markers=3,
    scores={"m1": 50.0, "m2": 50.0, "m3": 50.0},
    lengths={"m1": (5.0, 100.0), "m2": (10.0, 200.0), "m3": (5.0, 100.0)},
)


def _row(protein: str, busco: str, score: float, hmm_from: int, hmm_to: int) -> str:
    """One domtblout row. Coordinates are HMM-profile coords, which is what
    BUSCO measures coverage on -- fields 16 and 17, not the sequence envelope."""
    f = ["-"] * 23
    f[0], f[3] = protein, busco
    f[2], f[5] = "300", "300"
    f[7] = str(score)
    f[15], f[16] = str(hmm_from), str(hmm_to)
    f[19], f[20] = str(hmm_from), str(hmm_to)
    return " ".join(f)


def _write(tmp_path, rows):
    p = tmp_path / "markers.domtbl"
    p.write_text("# hmmsearch domtblout\n" + "\n".join(rows) + "\n")
    return p


def test_hmm_coverage_merges_overlapping_domains():
    assert sum_hmm_len([(1, 50), (40, 80)]) == 80
    assert sum_hmm_len([(1, 50), (60, 80)]) == 71
    assert sum_hmm_len([]) == 0


def test_single_copy_marker(tmp_path):
    p = _write(tmp_path, [_row("c1_1", "m1", 120.0, 1, 99)])
    res = parse_domtbl(p, LINEAGE, {"c1_1": "c1"})
    assert res.single == ["m1"]
    assert res.duplicated == []
    assert set(res.missing) == {"m2", "m3"}
    assert res.completeness == round(100 / 3, 2)


def test_two_complete_copies_are_duplicated_and_traceable_to_contigs(tmp_path):
    p = _write(
        tmp_path,
        [_row("c1_1", "m1", 120.0, 1, 99), _row("c7_4", "m1", 118.0, 1, 99)],
    )
    res = parse_domtbl(p, LINEAGE, {"c1_1": "c1", "c7_4": "c7"})
    assert res.duplicated == ["m1"]
    assert sorted(res.duplicated_marker_contigs()["m1"]) == ["c1", "c7"]
    assert res.contig_marker_counts() == {"c1": 1, "c7": 1}


def test_short_alignment_is_fragmented_not_complete(tmp_path):
    # m2 expects 200 aa with sigma 10; 120 aa is below 200 - 2*10 = 180.
    p = _write(tmp_path, [_row("c1_2", "m2", 300.0, 1, 120)])
    res = parse_domtbl(p, LINEAGE, {"c1_2": "c1"})
    assert res.fragmented == ["m2"]
    assert res.single == [] and res.duplicated == []
    assert res.completeness == 0.0


def test_hit_below_score_cutoff_is_ignored(tmp_path):
    p = _write(tmp_path, [_row("c1_3", "m1", 49.9, 1, 99)])
    res = parse_domtbl(p, LINEAGE, {"c1_3": "c1"})
    assert set(res.missing) == {"m1", "m2", "m3"}
    assert res.hits == []


def test_multiple_domains_of_one_protein_sum_toward_completeness(tmp_path):
    # Two domains covering profile positions 1-60 and 55-100: 100 covered,
    # above the 90 position threshold.
    p = _write(
        tmp_path,
        [_row("c1_4", "m1", 120.0, 1, 60), _row("c1_4", "m1", 120.0, 55, 100)],
    )
    res = parse_domtbl(p, LINEAGE, {"c1_4": "c1"})
    assert res.single == ["m1"]
    assert res.hits[0].matched_aa == 100


def test_summary_reports_busco_style_string():
    res = MarkerResult(
        lineage="test_odb10", n_markers=100,
        single=["s"] * 90, duplicated=["d"] * 5, fragmented=["f"] * 3, missing=["m"] * 2,
    )
    s = res.summary()
    assert s["completeness_percent"] == 95.0
    assert s["duplication_percent"] == 5.0
    assert s["busco_style_string"] == "C:95.0%[S:90.0%,D:5.0%],F:3.0%,M:2.0%,n:100"
    assert "not BUSCO" in s["method"]


def test_empty_table_is_all_missing(tmp_path):
    p = _write(tmp_path, [])
    res = parse_domtbl(p, LINEAGE, {})
    assert len(res.missing) == 3
    assert res.completeness == 0.0


def test_marker_table_lists_every_marker_with_its_class(tmp_path):
    from genomex.markers import write_marker_table

    p = _write(
        tmp_path,
        [
            _row("c1_1", "m1", 120.0, 1, 99),    # single -> Complete
            _row("c1_2", "m2", 300.0, 1, 120),   # short  -> Fragmented
            _row("c9_5", "m3", 120.0, 1, 99),
            _row("c4_2", "m3", 119.0, 1, 99),    # two copies -> Duplicated
        ],
    )
    res = parse_domtbl(p, LINEAGE, {"c1_1": "c1", "c1_2": "c1", "c9_5": "c9", "c4_2": "c4"})
    out = write_marker_table(res, tmp_path / "markers.tsv")

    lines = out.read_text().strip().splitlines()
    header, rows = lines[0], lines[1:]
    assert header.split("\t")[:2] == ["busco_id", "status"]

    status = {r.split("\t")[0]: r.split("\t")[1] for r in rows}
    assert status == {"m1": "Complete", "m2": "Fragmented", "m3": "Duplicated"}
    # a duplicated marker contributes one row per copy, so its contigs stay visible
    m3_contigs = sorted(r.split("\t")[3] for r in rows if r.startswith("m3\t"))
    assert m3_contigs == ["c4", "c9"]


def test_marker_table_records_missing_markers_too(tmp_path):
    from genomex.markers import write_marker_table

    p = _write(tmp_path, [_row("c1_1", "m1", 120.0, 1, 99)])
    res = parse_domtbl(p, LINEAGE, {"c1_1": "c1"})
    out = write_marker_table(res, tmp_path / "markers.tsv")
    rows = out.read_text().strip().splitlines()[1:]
    missing = [r for r in rows if r.split("\t")[1] == "Missing"]
    assert len(missing) == 2
    assert all(r.split("\t")[2] == "-" for r in missing)


def test_weak_paralog_hits_do_not_count_as_copies(tmp_path):
    """A hit far below the marker's best score is a distant paralog, not a copy.

    Real case from Paraburkholderia: marker 1074831at2 had one hit at 296.8 and
    six between 17 and 30, all above a permissive score cutoff. Counting them as
    copies reported the marker duplicated and inflated the genome's duplication
    rate -- the number the contamination module then consumes.
    """
    rows = [_row("c1_1", "m1", 296.8, 1, 99)]
    rows += [
        _row(f"c1_{i + 2}", "m1", score, 1, 99)
        for i, score in enumerate([29.7, 24.1, 20.9, 19.7, 19.0, 17.4])
    ]
    res = parse_domtbl(_write(tmp_path, rows), LINEAGE, {})

    assert res.single == ["m1"], "one strong hit plus weak paralogs is a single copy"
    assert res.duplicated == []
    assert len(res.hits) == 1


def test_a_genuine_second_copy_still_counts(tmp_path):
    """The retention rule must not erase real duplication -- the contamination signal."""
    rows = [_row("c1_1", "m1", 300.0, 1, 99), _row("c9_4", "m1", 291.0, 1, 99)]
    res = parse_domtbl(_write(tmp_path, rows), LINEAGE, {"c1_1": "c1", "c9_4": "c9"})

    assert res.duplicated == ["m1"]
    assert sorted(res.duplicated_marker_contigs()["m1"]) == ["c1", "c9"]


def test_retention_threshold_boundary(tmp_path):
    """Exactly 85% of the best score is retained; just below it is not."""
    kept = parse_domtbl(
        _write(tmp_path, [_row("c1_1", "m1", 100.0, 1, 99), _row("c2_1", "m1", 85.0, 1, 99)]),
        LINEAGE, {},
    )
    dropped = parse_domtbl(
        _write(tmp_path, [_row("c1_1", "m1", 100.0, 1, 99), _row("c2_1", "m1", 84.9, 1, 99)]),
        LINEAGE, {},
    )
    assert kept.duplicated == ["m1"]
    assert dropped.single == ["m1"]


def test_hmm_coverage_not_sequence_envelope_decides_completeness(tmp_path):
    """A wide envelope over a short profile match must not read as complete.

    m2 expects 200 profile positions with sigma 10, so complete needs >= 180.
    This hit covers 120 profile positions inside a 400-residue envelope.
    """
    f = ["-"] * 23
    f[0], f[3] = "c1_9", "m2"
    f[7] = "300.0"
    f[15], f[16] = "1", "120"      # HMM coords: 120 positions
    f[19], f[20] = "1", "400"      # envelope: much wider
    res = parse_domtbl(_write(tmp_path, [" ".join(f)]), LINEAGE, {})
    assert res.fragmented == ["m2"]
    assert res.hits[0].matched_aa == 120
