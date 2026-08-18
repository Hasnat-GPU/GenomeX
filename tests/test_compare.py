"""Orthogroup partitioning and the explanations attached to strain-unique genes."""

from genomex.compare import Pangenome, explain_unique_genes
from genomex.genes import Gene


def _cluster(*members: str) -> list[str]:
    return list(members)


def test_core_accessory_unique_partition_matches_the_planted_answer():
    clusters = {
        "r1": _cluster("A|c1_1", "B|c1_1", "C|c1_1"),   # core
        "r2": _cluster("A|c1_2", "B|c1_2"),             # accessory
        "r3": _cluster("A|c1_3"),                       # unique to A
        "r4": _cluster("C|c1_9", "C|c1_10"),            # unique to C (a paralog pair)
    }
    pg = Pangenome(["A", "B", "C"], clusters)
    part = pg.partition()
    assert part["core"] == ["r1"]
    assert part["accessory"] == ["r2"]
    assert sorted(part["unique"]) == ["r3", "r4"]
    s = pg.summary()
    assert (s["core"], s["accessory"], s["strain_unique"]) == (1, 1, 2)
    assert s["orthogroups_total"] == 4


def test_genes_of_returns_only_that_genomes_members():
    clusters = {"r1": _cluster("A|c1_1", "B|c1_1"), "r2": _cluster("A|c1_2")}
    pg = Pangenome(["A", "B"], clusters)
    assert pg.genes_of("A", ["r1", "r2"]) == ["c1_1", "c1_2"]
    assert pg.genes_of("B", ["r2"]) == []


def _genes(contig: str, indices: list[int]) -> list[Gene]:
    return [
        Gene(f"{contig}_{i}", contig, i * 1000, i * 1000 + 900, 1, "00", 300)
        for i in indices
    ]


def test_consecutive_unique_genes_are_called_an_island():
    genes = _genes("c1", list(range(1, 21)))
    unique = ["c1_5", "c1_6", "c1_7", "c1_8"]
    out = explain_unique_genes(unique, genes)
    assert all(g.island_id is not None for g in out)
    assert {g.island_size for g in out} == {4}
    assert all("island" in g.explanation for g in out)


def test_scattered_unique_genes_are_not_an_island():
    genes = _genes("c1", list(range(1, 41)))
    unique = ["c1_3", "c1_15", "c1_31"]
    out = explain_unique_genes(unique, genes)
    assert all(g.island_id is None for g in out)
    assert all("isolated" in g.explanation for g in out)


def test_a_single_shared_gene_does_not_break_an_island():
    # island_max_gap=1 tolerates one interrupting shared gene.
    genes = _genes("c1", list(range(1, 21)))
    out = explain_unique_genes(["c1_5", "c1_6", "c1_8", "c1_9"], genes)
    assert len({g.island_id for g in out}) == 1
    assert out[0].island_size == 4


def test_genes_on_a_contaminated_contig_are_labelled_as_such_not_as_biology():
    genes = _genes("cX", [1, 2, 3, 4]) + _genes("c1", [1, 2])
    out = explain_unique_genes(
        ["cX_1", "cX_2", "cX_3", "cX_4", "c1_1"], genes, suspect_contigs={"cX"}
    )
    contaminated = [g for g in out if g.contig == "cX"]
    assert len(contaminated) == 4
    assert all(g.on_suspect_contig for g in contaminated)
    assert all(g.explanation == "contamination-suspect contig" for g in contaminated)
    # ...and the contamination label wins over the island label, which would
    # otherwise fire on this run of four consecutive genes.
    assert all("island" not in g.explanation for g in contaminated)


def test_islands_are_split_across_contigs():
    genes = _genes("c1", [1, 2, 3]) + _genes("c2", [1, 2, 3])
    out = explain_unique_genes(
        ["c1_1", "c1_2", "c1_3", "c2_1", "c2_2", "c2_3"], genes
    )
    islands = {g.island_id for g in out}
    assert len(islands) == 2, "genes on different contigs must not share an island"


def test_gc_deviation_marks_an_isolated_gene_as_an_acquisition_candidate(tmp_path):
    genes = _genes("c1", list(range(1, 21)))
    # Genome-wide GC ~62%, the unique gene sits at ~30%.
    records = [(f"c1_{i}", "GCGC" * 200) for i in range(1, 21)]
    records[6] = ("c1_7", "ATAT" * 190 + "GCGC" * 10)
    fna = tmp_path / "genes.fna"
    fna.write_text("".join(f">{h}\n{s}\n" for h, s in records))

    out = explain_unique_genes(["c1_7"], genes, genes_fna=fna)
    assert out[0].gc_deviation is not None
    assert out[0].gc_deviation < -5
    assert "atypical GC" in out[0].explanation


def test_island_enrichment_detects_real_clustering():
    """Genuinely clustered unique genes must beat the permutation null."""
    from genomex.compare import island_enrichment

    genes = _genes("c1", list(range(1, 101)))
    clustered = [f"c1_{i}" for i in range(10, 25)]  # one solid 15-gene block
    stats = island_enrichment(clustered, genes, trials=25)

    assert stats["genes_in_islands_observed"] == 15
    assert stats["enrichment"] > 1.5
    assert stats["informative"] is True


def test_island_enrichment_rejects_clustering_that_is_only_density():
    """When most genes are unique, runs are inevitable and must not be called islands."""
    from genomex.compare import island_enrichment

    genes = _genes("c1", list(range(1, 101)))
    # 80 of 100 genes unique, spread evenly: runs happen because there is no room not to.
    dense = [f"c1_{i}" for i in range(1, 101) if i % 5 != 0]
    stats = island_enrichment(dense, genes, trials=25)

    assert stats["enrichment"] is not None
    assert stats["enrichment"] < 1.5
    assert stats["informative"] is False


def test_island_enrichment_is_deterministic():
    from genomex.compare import island_enrichment

    genes = _genes("c1", list(range(1, 61)))
    unique = [f"c1_{i}" for i in (3, 4, 5, 6, 20, 40, 41, 42)]
    assert island_enrichment(unique, genes) == island_enrichment(unique, genes)
