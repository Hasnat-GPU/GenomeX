"""The mixture harness is a measuring instrument, so it is tested like one.

A benchmark that silently builds the wrong thing produces numbers that look
fine and mean nothing -- `compare_to_busco.py` once read the wrong column and
printed a table of `None` before anyone noticed. These tests pin the properties
the mixture ladder's conclusions actually rest on: that the truth labels really
are the donor's sequence, that the requested contamination level is the level
delivered, and that the same inputs give the same mixture every time.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "bench"))

from fixtures.synthetic import clean_genome
from genomex.fasta import Assembly, write_fasta
from mixture_ladder import DONOR_PREFIX, MarkerHit, build_mixture


def _asm(tmp_path, genome):
    return Assembly.load(write_fasta(tmp_path / f"{genome.name}.fa", genome.records))


@pytest.fixture
def pair(tmp_path):
    host = _asm(tmp_path, clean_genome(name="host", n_contigs=6, contig_len=120_000, seed=1))
    donor = _asm(tmp_path, clean_genome(name="donor", n_contigs=6, contig_len=120_000, seed=2))
    return host, donor


def _markers(asm, per_contig=4):
    """Evenly spaced fake core markers, shared ids between host and donor."""
    out = []
    for c in asm.contigs:
        for i in range(per_contig):
            out.append(MarkerHit(f"m{i}_{c.name[-1]}", c.name,
                                 (i + 1) * c.length // (per_contig + 1)))
    return out


def test_zero_percent_mixture_adds_nothing(pair):
    host, donor = pair
    mix = build_mixture(host, donor, _markers(host), _markers(donor), 0.0, host_pieces=40)
    assert mix.donor_contigs == set()
    assert mix.donor_bp == 0
    assert mix.achieved_pct == 0.0
    assert mix.assembly.total_bp == host.total_bp


@pytest.mark.parametrize("pct", [2.0, 5.0, 10.0, 20.0])
def test_requested_level_is_the_level_delivered(pair, pct):
    """Within one fragment's worth. If the x-axis drifts, every recall number
    plotted against it is mislabelled."""
    host, donor = pair
    mix = build_mixture(host, donor, _markers(host), _markers(donor), pct, host_pieces=40)
    assert mix.achieved_pct == pytest.approx(pct, abs=1.0)


def test_donor_contigs_are_labelled_and_are_donor_sequence(pair):
    host, donor = pair
    mix = build_mixture(host, donor, _markers(host), _markers(donor), 10.0, host_pieces=40)
    by_name = {c.name: c for c in mix.assembly.contigs}

    assert mix.donor_contigs
    assert all(n.startswith(DONOR_PREFIX) for n in mix.donor_contigs)
    host_names = {c.name for c in mix.assembly.contigs} - mix.donor_contigs
    assert host_names and not any(n.startswith(DONOR_PREFIX) for n in host_names)

    donor_seqs = "".join(c.seq for c in donor.contigs)
    for name in mix.donor_contigs:
        assert by_name[name].seq in donor_seqs


def test_donor_fragments_are_not_longer_than_host_fragments(pair):
    """Length must not be a shortcut. A donor replicon left whole would be both
    an outlier in length and an overshoot of the target mass."""
    host, donor = pair
    mix = build_mixture(host, donor, _markers(host), _markers(donor), 10.0, host_pieces=40)
    by_name = {c.name: c for c in mix.assembly.contigs}
    host_max = max(len(by_name[c.name].seq)
                   for c in mix.assembly.contigs if c.name not in mix.donor_contigs)
    assert all(len(by_name[n].seq) <= host_max for n in mix.donor_contigs)


def test_donor_fragments_carry_markers_and_create_duplication(pair):
    """A real contaminant displaces part of the single-copy core. Withholding
    that would test a detector nobody runs."""
    host, donor = pair
    mix = build_mixture(host, donor, _markers(host), _markers(donor), 20.0, host_pieces=40)

    on_donor = sum(n for c, n in mix.marker_counts.items() if c in mix.donor_contigs)
    assert on_donor > 0
    cross = [b for b, cs in mix.duplicated.items() if len(set(cs)) > 1]
    assert cross, "no marker landed on both host and donor sequence"


def test_mixture_is_deterministic(pair):
    host, donor = pair
    hm, dm = _markers(host), _markers(donor)
    a = build_mixture(host, donor, hm, dm, 10.0, seed=3, host_pieces=40)
    b = build_mixture(host, donor, hm, dm, 10.0, seed=3, host_pieces=40)
    assert a.donor_contigs == b.donor_contigs
    assert a.donor_bp == b.donor_bp
    assert [c.name for c in a.assembly.contigs] == [c.name for c in b.assembly.contigs]
