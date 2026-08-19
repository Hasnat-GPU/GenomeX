"""Ground-truth validation of the contamination detector.

Every test here plants a known answer and checks the detector recovers it.
The hard case (a contaminant at the *same* GC as its host) is included
deliberately: it is where GC-only screening fails and TNF has to carry the call,
and it is the honest measure of what this module can do.
"""

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
    """A compositionally distinct group that carries the single-copy core at the
    genome-wide rate is chromosomal sequence from another organism, and mass must
    not buy it a plasmid reading however large it is."""
    genome = contaminated_genome(
        name="organism", n_host=12, n_foreign=3, contig_len=60_000,
        host_gc=0.63, foreign_gc=0.52, seed=13,
    )
    asm = _load(tmp_path, genome)
    marker_counts = {c.name: 10 for c in asm.contigs}  # core spread everywhere

    result = detect_contamination(asm, contig_marker_counts=marker_counts)
    calls = {c.name: c.call for c in result.contigs}
    for planted in genome.foreign_contigs:
        assert calls[planted] == "contaminant_candidate", planted
    assert genome.foreign_contigs <= result.suspect_contig_names()


def test_too_few_contigs_returns_undetermined_not_a_guess(tmp_path):
    p = write_fasta(tmp_path / "tiny.fa", [("c1", "ACGT" * 2000), ("c2", "ACGT" * 2000)])
    result = detect_contamination(Assembly.load(p))
    assert result.verdict == "undetermined"
    assert result.n_suspect_contigs == 0


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
