import gzip

from genomex.fasta import Assembly, assembly_stats, iter_fasta, write_fasta


def test_parses_multiline_records_and_descriptions(tmp_path):
    p = tmp_path / "x.fa"
    p.write_text(">c1 length=14 cov=3\nACGTACGT\nAAGG\nCC\n>c2\nGGGG\n")
    asm = Assembly.load(p)
    assert [c.name for c in asm.contigs] == ["c1", "c2"]
    assert asm.contigs[0].description == "length=14 cov=3"
    assert asm.contigs[0].length == 14
    assert asm.total_bp == 18


def test_reads_gzipped_input(tmp_path):
    p = tmp_path / "x.fa.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(">c1\nACGTACGTAC\n")
    assert [h for h, _ in iter_fasta(p)] == ["c1"]
    assert Assembly.load(p).total_bp == 10


def test_gc_ignores_ambiguous_bases():
    p = Assembly(path=__import__("pathlib").Path("."), name="t", contigs=[])
    from genomex.fasta import Contig

    c = Contig(name="c", seq="GGCCAANN")
    assert c.n_count == 2
    assert c.gc == 4 / 6  # 4 GC bases of 6 unambiguous ones; the two Ns do not count
    assert p.total_bp == 0


def test_n50_l50_against_hand_computed_values(tmp_path):
    # lengths 100, 60, 30, 10 -> total 200; half = 100 reached by the first contig
    records = [(f"c{i}", "A" * n) for i, n in enumerate([100, 60, 30, 10])]
    p = write_fasta(tmp_path / "y.fa", records)
    stats = assembly_stats(Assembly.load(p))
    assert stats["total_bp"] == 200
    assert stats["n_contigs"] == 4
    assert stats["n50"] == 100 and stats["l50"] == 1
    # 90% = 180 -> 100+60+30 = 190 passes at the third contig
    assert stats["n90"] == 30 and stats["l90"] == 3
    assert stats["longest_contig"] == 100


def test_gc_percent_is_computed_over_unambiguous_bases(tmp_path):
    p = write_fasta(tmp_path / "z.fa", [("c1", "GGCC" + "N" * 96)])
    stats = assembly_stats(Assembly.load(p))
    assert stats["gc_percent"] == 100.0
    assert stats["n_bases"] == 96


def test_roundtrip_write_then_read(tmp_path):
    records = [("a", "ACGT" * 30), ("b", "TTTT")]
    p = write_fasta(tmp_path / "rt.fa", records, width=17)
    assert list(iter_fasta(p)) == records
