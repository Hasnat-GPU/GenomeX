"""The fragmentation harness is a measuring instrument, so it is tested like one.

Its whole output is a false-positive *rate*, and a rate is only as good as the
independence of the draws behind it. Two properties carry that: different seeds
must give genuinely different shreddings, and one seed must always give the same
one. If seeds silently collapsed onto the same cut points, 320 draws would report
the confidence of 320 samples while holding the information of 40.

The marker mapping is pinned for the same reason. `marker_rate` is estimated from
markers located on contigs, and the sole-copy rule in `genomex.contamination`
compares a group's marker count against that rate -- so a harness that dropped or
double-counted markers when it cut would move a call without moving any code.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))

from fixtures.synthetic import clean_genome
from fragmentation_ladder import marker_counts_for, shred_with_spans
from genomex.fasta import Assembly, write_fasta


@pytest.fixture
def genome(tmp_path):
    g = clean_genome(name="host", n_contigs=4, contig_len=200_000, seed=5)
    return Assembly.load(write_fasta(tmp_path / "host.fa", g.records))


def _cut_points(asm, n, seed):
    _frag, spans = shred_with_spans(asm, n, seed=seed)
    return tuple(sorted((parent, lo, hi) for _c, parent, lo, hi in spans))


def test_the_same_seed_gives_the_same_shredding(genome):
    """Determinism, and the reason a `--seeds 8` run can be re-run and checked."""
    assert _cut_points(genome, 50, 0) == _cut_points(genome, 50, 0)


def test_different_seeds_give_different_shreddings(genome):
    """Otherwise the rate is a single draw wearing eight hats."""
    draws = {_cut_points(genome, 50, s) for s in range(8)}
    assert len(draws) == 8


def test_shredding_conserves_every_base(genome):
    frag, _spans = shred_with_spans(genome, 50, seed=3)
    assert frag.total_bp == genome.total_bp


def test_spans_tile_each_parent_without_gap_or_overlap(genome):
    """A marker's fragment is found by span containment, so the spans must be a
    partition -- a gap loses the marker, an overlap counts it twice."""
    _frag, spans = shred_with_spans(genome, 50, seed=3)
    by_parent: dict[str, list[tuple[int, int]]] = {}
    for _c, parent, lo, hi in spans:
        by_parent.setdefault(parent, []).append((lo, hi))
    lengths = {c.name: c.length for c in genome.contigs}
    for parent, pieces in by_parent.items():
        pieces.sort()
        assert pieces[0][0] == 0
        assert pieces[-1][1] == lengths[parent]
        for (_lo_a, hi_a), (lo_b, _hi_b) in zip(pieces, pieces[1:]):
            assert hi_a == lo_b, parent


def test_every_marker_lands_on_exactly_one_fragment(genome):
    """Total markers in must equal total markers out, at every rung."""
    markers = [(c.name, p) for c in genome.contigs
               for p in range(5_000, c.length, 20_000)]
    for n in (10, 50, 300):
        _frag, spans = shred_with_spans(genome, n, seed=1)
        counts = marker_counts_for(spans, markers)
        assert sum(counts.values()) == len(markers), n
