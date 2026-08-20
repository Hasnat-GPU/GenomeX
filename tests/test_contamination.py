"""Ground-truth validation of the contamination detector.

Every test here plants a known answer and checks the detector recovers it.
The hard case (a contaminant at the *same* GC as its host) is included
deliberately: it is where GC-only screening fails and TNF has to carry the call,
and it is the honest measure of what this module can do.
"""

import json

import numpy as np
import pytest

from fixtures.synthetic import clean_genome, contaminated_genome
from genomex.contamination import (
    N_CANON,
    detect_contamination,
    tetranucleotide_freq,
)
from genomex.fasta import Assembly, write_fasta


def _load(tmp_path, genome):
    return Assembly.load(write_fasta(tmp_path / f"{genome.name}.fa", genome.records))


def _revcomp(s: str) -> str:
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


# ---------------------------------------------------------------- TNF vector

def test_tnf_vector_is_canonical_and_normalised():
    seq = "ACGTTGCAAGGCTTACGATCGGATCCAGTTACG" * 20
    v = tetranucleotide_freq(seq)
    assert v.shape == (N_CANON,) == (136,)
    assert v.sum() == pytest.approx(1.0)


def test_tnf_is_strand_independent():
    seq = "ACGTTGCAAGGCTTACGATCGGATCCAGTTACG" * 20
    assert np.allclose(tetranucleotide_freq(seq), tetranucleotide_freq(_revcomp(seq)))


def test_tnf_ignores_ambiguous_windows():
    clean = tetranucleotide_freq("ACGT" * 50)
    withn = tetranucleotide_freq("ACGT" * 25 + "N" * 10 + "ACGT" * 25)
    assert np.allclose(clean, withn, atol=0.02)


def test_tnf_of_empty_or_tiny_sequence_is_zero():
    assert tetranucleotide_freq("").sum() == 0
    assert tetranucleotide_freq("AC").sum() == 0


# ------------------------------------------------------- end-to-end verdicts

def test_single_organism_is_called_clean_with_no_false_positives(tmp_path):
    genome = clean_genome(seed=3)
    result = detect_contamination(_load(tmp_path, genome))
    assert result.verdict == "clean"
    assert result.n_suspect_contigs == 0, f"false positives: {result.suspect_contig_names()}"


def test_planted_contaminant_is_found_and_named(tmp_path):
    genome = contaminated_genome(n_host=12, n_foreign=3, seed=3)
    result = detect_contamination(_load(tmp_path, genome))

    assert result.verdict in ("likely", "possible")
    flagged = result.suspect_contig_names()
    truth = genome.foreign_contigs
    recall = len(flagged & truth) / len(truth)
    false_positives = flagged - truth
    assert recall == 1.0, f"missed planted contigs: {truth - flagged}"
    assert not false_positives, f"flagged host contigs: {false_positives}"
    assert result.bins["bimodal"] is True


def test_contaminant_at_identical_gc_is_still_detected(tmp_path):
    """The hard case: only k-mer composition separates host from contaminant."""
    genome = contaminated_genome(
        name="isogc", n_host=14, n_foreign=3, host_gc=0.62, foreign_gc=0.62, seed=7
    )
    result = detect_contamination(_load(tmp_path, genome))
    flagged = result.suspect_contig_names()
    truth = genome.foreign_contigs
    # GC alone cannot separate these, so any recovery at all is TNF doing the work.
    assert len(flagged & truth) >= 2, (
        f"TNF failed on the iso-GC case: flagged={flagged} truth={truth} "
        f"verdict={result.verdict}"
    )


def test_duplicated_markers_across_contigs_raise_the_verdict(tmp_path):
    """Composition alone says clean; duplicated single-copy markers must not be ignored."""
    genome = clean_genome(seed=5)
    asm = _load(tmp_path, genome)
    names = [c.name for c in asm.contigs]
    dup = {f"busco{i}": [names[0], names[1]] for i in range(5)}

    baseline = detect_contamination(asm)
    with_dups = detect_contamination(asm, duplicated_marker_contigs=dup)

    assert baseline.verdict == "clean"
    assert with_dups.verdict == "likely"
    assert any("duplicated" in r for r in with_dups.reasons)
    assert {names[0], names[1]} <= with_dups.flagged_contig_names()


def test_marker_duplication_alone_does_not_condemn_a_contig(tmp_path):
    """A chromosome sharing paralogs with a chromid is not a foreign contig.

    The genome-level verdict rises, but the contigs stay out of the contamination
    fraction -- otherwise every multipartite genome reports its own chromosome as
    a contaminant, and the comparative module then discards real genes as artefacts.
    """
    genome = clean_genome(seed=5)
    asm = _load(tmp_path, genome)
    names = [c.name for c in asm.contigs]
    dup = {f"busco{i}": [names[0], names[1]] for i in range(5)}

    result = detect_contamination(asm, duplicated_marker_contigs=dup)
    calls = {c.name: c.call for c in result.contigs}

    assert calls[names[0]] == "marker_conflict"
    assert calls[names[1]] == "marker_conflict"
    assert result.suspect_contig_names() == set()
    assert result.suspect_bp == 0
    assert set(result.summary()["marker_conflict_contigs"]) == {names[0], names[1]}


def test_distinct_sequence_without_core_markers_is_called_a_replicon(tmp_path):
    """The megaplasmid case: distinct composition, and positive evidence it is not
    a second organism -- the rest of the assembly carries the single-copy core and
    this sequence carries none of it. Plasmids and chromids do not carry the
    universal core; a foreign genome fragment does."""
    genome = contaminated_genome(
        name="plasmid", n_host=12, n_foreign=1, contig_len=60_000,
        host_gc=0.63, foreign_gc=0.52, seed=13,
    )
    asm = _load(tmp_path, genome)
    planted = next(iter(genome.foreign_contigs))
    # Core markers spread over the host contigs, none on the distinct sequence.
    host_contigs = [c.name for c in asm.contigs if c.name != planted]
    marker_counts = {name: 10 for name in host_contigs}

    result = detect_contamination(asm, contig_marker_counts=marker_counts)
    call = {c.name: c.call for c in result.contigs}[planted]

    assert call == "replicon_candidate"
    assert planted not in result.suspect_contig_names()
    assert planted in result.flagged_contig_names()
    assert any("plasmid" in r for r in result.reasons)


def test_without_a_marker_scan_no_plasmid_claim_is_made(tmp_path):
    """Absence of evidence is not evidence of absence.

    The same sequence, with no marker information supplied, must be reported as a
    contaminant candidate rather than explained away as a plasmid -- the plasmid
    reading rests on knowing the core markers are elsewhere, which an unscanned
    assembly cannot establish."""
    genome = contaminated_genome(
        name="plasmid", n_host=12, n_foreign=1, contig_len=60_000,
        host_gc=0.63, foreign_gc=0.52, seed=13,
    )
    asm = _load(tmp_path, genome)
    planted = next(iter(genome.foreign_contigs))

    result = detect_contamination(asm)
    assert {c.name: c.call for c in result.contigs}[planted] == "contaminant_candidate"
    assert planted in result.suspect_contig_names()


def test_foreign_sequence_carrying_core_markers_is_a_second_organism(tmp_path):
    """A compositionally distinct group carrying *duplicate* copies of the core is
    chromosomal sequence from another organism, and mass must not buy it a plasmid
    reading however large it is.

    Duplication is the load-bearing word. A second organism arrives with its own
    copy of a universal single-copy gene, so the assembly ends up holding two --
    one on host sequence, one on foreign. That is what is planted here."""
    genome = contaminated_genome(
        name="organism", n_host=12, n_foreign=3, contig_len=60_000,
        host_gc=0.63, foreign_gc=0.52, seed=13,
    )
    asm = _load(tmp_path, genome)
    marker_counts = {c.name: 10 for c in asm.contigs}  # core spread everywhere
    host = sorted(c.name for c in asm.contigs if c.name not in genome.foreign_contigs)
    # Each foreign contig holds a second copy of markers the host already has.
    duplicated = {
        f"BUSCO{k}": [host[k % len(host)], foreign]
        for foreign in sorted(genome.foreign_contigs)
        for k in range(4)
    }

    result = detect_contamination(
        asm, contig_marker_counts=marker_counts, duplicated_marker_contigs=duplicated,
    )
    calls = {c.name: c.call for c in result.contigs}
    for planted in genome.foreign_contigs:
        assert calls[planted] == "contaminant_candidate", planted
    assert genome.foreign_contigs <= result.suspect_contig_names()


def test_sole_copies_of_the_core_make_a_region_the_host_not_a_contaminant(tmp_path):
    """The same geometry with the duplication removed must not be contamination.

    A group whose core markers are the assembly's *only* copies cannot be a second
    organism: a second organism brings its own copies of a universal single-copy
    set, and they show up as duplicates. Sole copies mean this is the one genome
    present, and the odd composition is a prophage, an island, or the tail of the
    null -- so it is reported, and excluded from the contamination fraction.

    Calling it contamination is what let the family-wise error rate leak straight
    into the genome verdict: `docs/benchmark-fragmentation.md` measures the cost."""
    genome = contaminated_genome(
        name="island", n_host=12, n_foreign=3, contig_len=60_000,
        host_gc=0.63, foreign_gc=0.52, seed=13,
    )
    asm = _load(tmp_path, genome)
    marker_counts = {c.name: 10 for c in asm.contigs}

    result = detect_contamination(asm, contig_marker_counts=marker_counts)
    calls = {c.name: c.call for c in result.contigs}
    for planted in genome.foreign_contigs:
        assert calls[planted] == "atypical_host_region", planted
    assert result.suspect_contig_names() == set()
    assert result.n_suspect_contigs == 0
    # Still reported, with the count that justifies the call.
    assert genome.foreign_contigs <= result.flagged_contig_names()
    assert {r["contig"] for r in result.summary()["atypical_host_regions"]} == (
        genome.foreign_contigs
    )
    flags = {c.name: c.flags for c in result.contigs}
    assert any("sole_copy_core_markers" in f for f in flags[sorted(genome.foreign_contigs)[0]])


def test_too_few_contigs_returns_undetermined_not_a_guess(tmp_path):
    p = write_fasta(tmp_path / "tiny.fa", [("c1", "ACGT" * 2000), ("c2", "ACGT" * 2000)])
    result = detect_contamination(Assembly.load(p))
    assert result.verdict == "undetermined"
    assert result.assessed is False
    # Not 0. This assertion used to read `== 0`, which is the bug written down
    # as a test: zero suspect contigs is what a clean genome measures.
    assert result.n_suspect_contigs is None
    assert result.suspect_bp is None
    assert result.suspect_fraction is None


def test_an_abstention_serialises_as_null_and_never_as_a_clean_result(tmp_path):
    """The failure this pins: a machine consumer reads `suspect_fraction_percent`
    without first checking `verdict`, and a refusal to measure comes back as
    0.0% contamination -- a better result than most real genomes get.

    Every quantitative field goes to `null`, the lists included. An empty
    `replicon_candidates` says "we looked and found none"; on this path nothing
    was looked at."""
    p = write_fasta(tmp_path / "tiny.fa", [("c1", "ACGT" * 2000), ("c2", "ACGT" * 2000)])
    s = detect_contamination(Assembly.load(p)).summary()

    assert s["verdict"] == "undetermined"
    for field_name in ("suspect_contigs", "suspect_bp", "suspect_fraction_percent",
                       "replicon_candidates", "atypical_host_regions",
                       "marker_conflict_contigs"):
        assert s[field_name] is None, f"{field_name} = {s[field_name]!r}, should be null"

    # The refusal still ships its reason -- an abstention with no evidence
    # attached would trade one silent failure for another.
    assert s["reasons"] and "more sequence" in s["reasons"][0]
    assert json.loads(json.dumps(s))["suspect_fraction_percent"] is None


def test_a_clean_genome_still_reports_measured_zeros(tmp_path):
    """The other half: `null` must mean unmeasured, so a genome that really was
    scored and really has nothing foreign has to keep reporting 0 and `[]`. If
    both cases returned null the distinction would be worthless."""
    genome = clean_genome(n_contigs=10, seed=11)
    s = detect_contamination(_load(tmp_path, genome)).summary()

    assert s["verdict"] == "clean"
    assert s["suspect_contigs"] == 0
    assert s["suspect_fraction_percent"] == 0.0
    assert s["replicon_candidates"] == []


def test_short_contigs_are_reported_as_skipped(tmp_path):
    genome = clean_genome(n_contigs=6, contig_len=20_000, seed=9)
    records = genome.records + [(f"short{i}", "ACGT" * 100) for i in range(4)]
    p = write_fasta(tmp_path / "mix.fa", records)
    result = detect_contamination(Assembly.load(p), min_contig_length=3000)
    assert result.params["contigs_skipped_short"] == 4
    assert result.params["contigs_scored"] == 6


def test_verdict_is_deterministic_across_runs(tmp_path):
    genome = contaminated_genome(seed=11)
    asm = _load(tmp_path, genome)
    a = detect_contamination(asm)
    b = detect_contamination(asm)
    assert a.verdict == b.verdict
    assert [c.suspicion for c in a.contigs] == [c.suspicion for c in b.contigs]
    assert a.bins == b.bins


# --------------------------------------------- marker-set scale invariance

def test_duplication_verdict_is_a_rate_not_a_count(tmp_path):
    """The same genome scanned with a bigger marker set must not get a worse verdict.

    Swapping bacteria_odb10 (124 markers) for burkholderiales_odb10 (688)
    multiplies every duplication count by 5.5 without adding one base of
    contamination. Under the old absolute cutoffs that turned finished reference
    genomes -- carrying an ordinary ~1% of cross-replicon paralogs -- into
    `likely` contamination.
    """
    genome = clean_genome(seed=5)
    asm = _load(tmp_path, genome)
    names = [c.name for c in asm.contigs]

    # 1% of each set duplicated across two contigs: the same biology, twice.
    small = {f"busco{i}": [names[0], names[1]] for i in range(1)}
    large = {f"busco{i}": [names[0], names[1]] for i in range(7)}

    assert detect_contamination(asm, duplicated_marker_contigs=small,
                                total_markers=124).verdict == "clean"
    assert detect_contamination(asm, duplicated_marker_contigs=large,
                                total_markers=688).verdict == "clean"


def test_duplication_rate_above_threshold_still_fires(tmp_path):
    """Scale invariance must not be bought by going blind: 5% of a large set is
    contamination just as 5% of a small one is."""
    genome = clean_genome(seed=5)
    asm = _load(tmp_path, genome)
    names = [c.name for c in asm.contigs]
    dup = {f"busco{i}": [names[0], names[1]] for i in range(40)}  # 5.8% of 688

    result = detect_contamination(asm, duplicated_marker_contigs=dup, total_markers=688)
    assert result.verdict == "likely"
    assert result.params["cross_contig_duplication_percent"] == pytest.approx(5.81, abs=0.01)


def test_unknown_marker_set_falls_back_to_counts(tmp_path):
    """No denominator means no rate. Falling back to the previous absolute rule is
    honest; inventing a marker-set size would not be."""
    genome = clean_genome(seed=5)
    asm = _load(tmp_path, genome)
    names = [c.name for c in asm.contigs]
    dup = {f"busco{i}": [names[0], names[1]] for i in range(5)}

    result = detect_contamination(asm, duplicated_marker_contigs=dup)
    assert result.verdict == "likely"
    assert result.params["cross_contig_duplication_percent"] is None
    assert result.params["marker_set_size"] is None
