"""Proteome input, and the fabrications it would otherwise produce.

Every test here pins a number that a wrong implementation returns *plausibly*.
That is the whole risk of this input path: nothing crashes when a proteome is
treated as an assembly. `assembly_stats` reports a GC content, `Assembly.load`
reports a contig count, and `parse_domtbl` reports a contig name. All three are
in range, none is a measurement, and a reader cannot tell.
"""

import pytest

from genomex.fasta import Assembly, assembly_stats, iter_fasta
from genomex.markers import UNKNOWN_CONTIG, Lineage, MarkerHit, MarkerResult, parse_domtbl
from genomex.proteome import (
    NotAnAssembly,
    NotAProteome,
    Proteome,
    assert_nucleotide,
    protein_only_fraction,
)

# Real amino-acid composition matters: the discriminator counts residues that
# have no nucleotide meaning, so a toy protein of only ACGT would be
# indistinguishable from DNA and *should* be.
PROTEINS = """\
>NP_009332.1 hypothetical protein [S. cerevisiae]
MSEQNNTEMTFQIQRIYTKDISFEAPNAPHVFQKDWQPEVKLDLDTASSQLADDVYEVVLR
VTVTASLGEETAFLCEVQQGGIFSIAGIEGTQMAHCLGAYCPNILFPYARECITSMVSRGT
>NP_010001.1 another protein
MGGRDAAWQPILEFKNAVTLGDKQWMSLPEQFKEYLQNDPHCFLIFEEGNKQWMDTVIPKQ
>NP_010002.1
MKKTTLLEQFPNVISQAGDKVWQEFLDNPHYALIFEEQGKQWMHTVIPK
"""

# Long enough to clear MIN_RESIDUES_FOR_ALPHABET_TEST -- below it the alphabet
# check abstains by design, so a short fixture would test nothing.
_CODING = ("ATGAGCGAACAAAACAACACCGAAATGACCTTTCAGATTCAGCGCATTTATACCAAAGATA"
           "TTAGCTTTGAAGCACCGAACGCACCGCATGTGTTTCAGAAAGATTGGCAGCCGGAAGTGAA")
NUCLEOTIDE = f">contig_1\n{_CODING}\n>contig_2\n{_CODING}\n"


LINEAGE = Lineage(
    path=__import__("pathlib").Path("."),
    name="test_odb10",
    n_markers=2,
    scores={"m1": 50.0, "m2": 50.0},
    lengths={"m1": (5.0, 100.0), "m2": (5.0, 100.0)},
)


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _records(path):
    return iter_fasta(path)


# --------------------------------------------------------------------------
# The discriminator


def test_protein_only_fraction_separates_the_two_alphabets(tmp_path):
    """A, C, G and T are also alanine, cysteine, glycine and threonine, so an
    ACGT-counting test cannot tell the alphabets apart. E/F/I/L/P/Q have no
    nucleotide reading at all."""
    prot, _ = protein_only_fraction(
        s for _h, s in _records(_write(tmp_path, "p.faa", PROTEINS))
    )
    nuc, _ = protein_only_fraction(
        s for _h, s in _records(_write(tmp_path, "n.fna", NUCLEOTIDE))
    )
    assert prot > 0.15, prot
    assert nuc == 0.0, nuc


def test_an_n_gapped_scaffold_is_not_mistaken_for_protein():
    """The reason the discriminator is not "how much of this is ACGTU".

    A scaffold padded with N scores 0.80 on plain ACGTU -- below a 0.90 bar, so
    a naive check would reject a perfectly ordinary assembly.
    """
    scaffold = "ACGT" * 20_000 + "N" * 20_000
    frac, total = protein_only_fraction([scaffold])
    assert total == 100_000
    assert frac == 0.0
    assert_nucleotide([scaffold], "scaffold.fna")  # must not raise


# --------------------------------------------------------------------------
# Refusals in both directions


def test_a_proteome_on_the_assembly_path_is_refused(tmp_path):
    """The refusal exists because the alternative is not a crash but a number.

    Handed the S. cerevisiae proteome, `assembly_stats` reports GC 35.36%
    against that genome's real 38.15% -- in range, and unfalsifiable by eye.
    """
    p = _write(tmp_path, "prot.faa", PROTEINS)
    asm = Assembly.load(p)
    with pytest.raises(NotAnAssembly, match="looks like protein"):
        assert_nucleotide((c.seq for c in asm.contigs), p)


def test_the_fabricated_stats_this_refusal_prevents(tmp_path):
    """Pins the failure mode itself, so it stays visible if the guard is moved.

    Without the guard these are what the report renders: a GC percentage from
    glycine and cysteine, a contig count that is a protein count, and a genome
    size that is a residue count.
    """
    asm = Assembly.load(_write(tmp_path, "prot.faa", PROTEINS))
    s = assembly_stats(asm)
    assert s["n_contigs"] == 3, "one 'contig' per protein"
    assert 20 < s["gc_percent"] < 60, (
        f"gc_percent {s['gc_percent']} is inside the plausible bacterial range, "
        "which is precisely why this cannot be left to a reader to catch"
    )


def test_a_nucleotide_file_on_the_proteome_path_is_refused(tmp_path):
    p = _write(tmp_path, "genome.faa", NUCLEOTIDE)
    with pytest.raises(NotAProteome, match="looks like nucleotide"):
        Proteome.load(p)


def test_a_short_peptide_of_nucleotide_legal_letters_is_still_accepted(tmp_path):
    """Sixteen letters are legal nucleotide codes, so a short peptide can be
    100% nucleotide-legal by accident.

    `MKVAAGCGTS` is entirely IUPAC-nucleotide: M and V are ambiguity codes as
    well as methionine and valine, K as well as lysine, S as well as serine.
    The alphabet test therefore cannot run on 20 residues, and must not pretend
    it can -- this is the same abstention as `min_contigs_for_zscores`.
    """
    p = _write(tmp_path, "tiny.faa", ">pep\nMKVAAGCGTSMKVAAGCGTS\n")
    assert len(Proteome.load(p)) == 1


def test_a_full_length_nucleotide_file_is_refused_where_a_peptide_is_not(tmp_path):
    """The gate is sample size, not leniency: past it, the refusal still fires."""
    long_dna = ">c1\n" + ("ACGTACGTAC" * 40) + "\n"   # 400 bases, over the gate
    with pytest.raises(NotAProteome, match="looks like nucleotide"):
        Proteome.load(_write(tmp_path, "big.faa", long_dna))


def test_an_empty_file_is_refused_rather_than_scored_as_zero(tmp_path):
    """Zero markers over zero proteins is 0% complete, which reads as a finding."""
    with pytest.raises(NotAProteome, match="no FASTA records"):
        Proteome.load(_write(tmp_path, "empty.faa", ""))


# --------------------------------------------------------------------------
# Honest stats


def test_stats_report_the_input_and_nothing_derived_from_a_genome(tmp_path):
    prot = Proteome.load(_write(tmp_path, "p.faa", PROTEINS))
    s = prot.stats()
    assert s["n_proteins"] == 3
    assert s["total_aa"] == sum(p.length_aa for p in prot.proteins)
    for forbidden in ("coding_density", "genes_per_mb", "gc_percent", "n50",
                      "genetic_code", "prodigal_procedure"):
        assert forbidden not in s, (
            f"{forbidden} has no value without a genome; emitting 0.0 would read "
            "as a measurement rather than an absence"
        )


def test_a_trailing_stop_codon_is_not_counted_as_a_residue(tmp_path):
    """Prodigal writes '*', NCBI does not. Counting it would make the same
    protein one residue longer depending on who called it."""
    a = Proteome.load(_write(tmp_path, "a.faa", ">x\nMKVLA*\n"))
    b = Proteome.load(_write(tmp_path, "b.faa", ">x\nMKVLA\n"))
    assert a.stats()["total_aa"] == b.stats()["total_aa"] == 5


# --------------------------------------------------------------------------
# The contig sentinel


def _domtbl(tmp_path, rows):
    """Minimal domtblout, same column layout as tests/test_markers.py."""
    lines = ["# hmmsearch domtblout"]
    for protein, busco, score, hmm_from, hmm_to in rows:
        f = ["-"] * 23
        f[0], f[3] = protein, busco
        f[2], f[5] = "300", "300"
        f[7] = str(score)
        f[15], f[16] = str(hmm_from), str(hmm_to)
        f[19], f[20] = str(hmm_from), str(hmm_to)
        lines.append(" ".join(f))
    p = tmp_path / "markers.domtbl"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return p


def test_marker_hits_without_a_contig_map_get_the_sentinel_not_a_guess(tmp_path):
    """The fabrication this replaced: `protein_id.rsplit("_", 1)[0]`.

    That heuristic is Prodigal's naming convention (`<contig>_<n>`). Applied to
    NCBI ids it turned 754 of 777 marker hits on the S. cerevisiae proteome into
    one fictitious contig called "NP".
    """
    dom = _domtbl(tmp_path, [("NP_009332.1", "m1", 120.0, 1, 99)])
    res = parse_domtbl(dom, LINEAGE, {})
    assert res.contigs_known is False
    assert all(h.contig == UNKNOWN_CONTIG for h in res.hits)
    assert not any(h.contig == "NP" for h in res.hits)


def test_contig_accessors_refuse_rather_than_aggregate_on_the_sentinel():
    """`{"m1": ["-", "-"]}` is not "both copies are on one contig".

    detect_contamination reads two copies on a single contig as affirmative
    evidence of paralogy and discards the marker, so a sentinel-keyed dict would
    launder missing information into a finding.
    """
    hits = [
        MarkerHit("m1", "NP_1.1", UNKNOWN_CONTIG, 500.0, 400, "complete"),
        MarkerHit("m1", "NP_2.1", UNKNOWN_CONTIG, 480.0, 390, "complete"),
    ]
    res = MarkerResult(lineage="x", n_markers=10, duplicated=["m1"], hits=hits,
                       contigs_known=False)
    assert res.contig_marker_counts() is None
    assert res.duplicated_marker_contigs() is None

    known = MarkerResult(lineage="x", n_markers=10, duplicated=["m1"],
                         hits=[MarkerHit("m1", "c1_4", "c1", 500.0, 400, "complete")])
    assert known.contig_marker_counts() == {"c1": 1}


def test_completeness_is_unaffected_by_contig_attribution():
    """Completeness must not depend on knowing where a marker sits, or the
    proteome path would measure something different from the assembly path."""
    hits = [MarkerHit("m1", "p1", UNKNOWN_CONTIG, 500.0, 400, "complete")]
    res = MarkerResult(lineage="x", n_markers=4, single=["m1"], hits=hits,
                       contigs_known=False)
    assert res.completeness == 25.0
