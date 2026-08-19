"""Length-conditioned composition scoring.

A tetranucleotide frequency estimated from 3 kb of sequence is far noisier than
one estimated from 300 kb, because it is a multinomial sample: the expected
chi-square divergence from the true composition falls as 1/length. Scoring every
contig against one robust-z threshold therefore flags the short tail of any
fragmented assembly, and the flagging rate becomes a function of contig count
rather than of contamination.

That is not hypothetical. Shredding one finished genome into 10, 50, 100, 300 and
1000 pieces -- same DNA, same organism -- moved the old detector's suspect
fraction from 0.00% to 6.51% and its verdict from `possible` to `likely`
(`bench/fragmentation_ladder.py`).

The fix is to compare each contig against a null built from sequence known to be
native: windows tiled inside the assembly's own contigs. Windows of length W drawn
from one organism give the divergence distribution expected at length W when
nothing is foreign, including the effects theory cannot supply -- overlapping
k-mer correlation, coding structure, local isochore drift. A contig is then
flagged only when it diverges more than same-length native sequence does.

The theoretical anchor is E[D] ~ (C-1)/L for C categories, which the measured
curve should track; `null_curve` records both so the two can be compared.

Rejected alternatives, each measured rather than argued away:

* Column-wise permutation of the TNF matrix, block bootstrap over k-mer count
  blocks, and multinomial draws from the pooled composition all destroy the
  covariance between 4-mers that makes real composition coherent, and declare
  every clean genome contaminated at p = 0.005.
* Permuting principal-component scores independently has the opposite failure:
  p = 1.000 with no power against a planted contaminant.
* Hartigan's dip test needs a one-dimensional projection chosen from 136
  dimensions and weights a 3 kb contig like a 4 Mb chromosome. The gap statistic
  references a uniform box, which is meaningless on a simplex. Silhouette is
  pathological on the unbalanced splits that contamination actually produces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

K = 4
N_RAW = 4**K
_COMPLEMENT = {0: 3, 1: 2, 2: 1, 3: 0}

_LUT = np.full(256, 255, dtype=np.uint8)
for _base, _code in zip(b"ACGTacgt", [0, 1, 2, 3, 0, 1, 2, 3]):
    _LUT[_base] = _code


def _revcomp_index(idx: int) -> int:
    out = 0
    for i in range(K):
        out = (out << 2) | _COMPLEMENT[(idx >> (2 * i)) & 3]
    return out


_CANON_GROUP: dict[int, int] = {}
_canon_order: list[int] = []
for _i in range(N_RAW):
    _key = min(_i, _revcomp_index(_i))
    if _key not in _CANON_GROUP:
        _CANON_GROUP[_key] = len(_canon_order)
        _canon_order.append(_key)
N_CANON = len(_canon_order)  # 136
_RAW_TO_CANON = np.array(
    [_CANON_GROUP[min(i, _revcomp_index(i))] for i in range(N_RAW)], dtype=np.int32
)

#: GC fraction of each canonical 4-mer. GC count is invariant under reverse
#: complement, so this is well defined after canonicalisation -- which makes GC
#: an exact linear functional of the TNF vector rather than a separate statistic.
_GC_W = np.zeros(N_CANON, dtype=np.float64)
for _raw in range(N_RAW):
    _gc = sum(1 for i in range(K) if ((_raw >> (2 * i)) & 3) in (1, 2))
    _GC_W[_RAW_TO_CANON[_raw]] = _gc / K

#: Smallest window used to build the null, and the number of dyadic rungs above
#: it. 2 kb is the floor at which a 136-category multinomial is estimable at all.
BASE_WINDOW = 2000
MAX_RUNGS = 8
#: Minimum windows at a rung before its median/MAD is trusted over theory.
MIN_WINDOWS_PER_RUNG = 8


def kmer_codes(seq: str) -> np.ndarray:
    """Canonical 4-mer codes for every unambiguous window of `seq`."""
    codes = _LUT[np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)]
    if codes.size < K:
        return np.empty(0, dtype=np.int32)
    valid = codes != 255
    windows = np.lib.stride_tricks.sliding_window_view(codes, K)
    ok = np.lib.stride_tricks.sliding_window_view(valid, K).all(axis=1)
    if not ok.any():
        return np.empty(0, dtype=np.int32)
    w = windows[ok].astype(np.int32)
    raw = (w[:, 0] << 6) | (w[:, 1] << 4) | (w[:, 2] << 2) | w[:, 3]
    return _RAW_TO_CANON[raw]


def counts_from_codes(codes: np.ndarray) -> np.ndarray:
    return np.bincount(codes, minlength=N_CANON).astype(np.float64)


def frequencies(counts: np.ndarray) -> np.ndarray:
    total = counts.sum()
    return counts / total if total else counts


def gc_percent_from_freq(freq: np.ndarray) -> float:
    """GC as a linear functional of the canonical TNF vector."""
    return 100.0 * float(freq @ _GC_W)


def chi2_divergence(freq: np.ndarray, reference: np.ndarray, floor: float = 1e-9) -> float:
    """Pearson divergence of a composition from a reference composition.

    Scale-free in the categories and, under multinomial sampling of n k-mers from
    `reference`, has expectation (C-1)/n. That 1/n is exactly the length effect
    being corrected for.
    """
    ref = np.maximum(reference, floor)
    d = freq - ref
    return float(np.sum(d * d / ref))


def base_tiles(codes: np.ndarray, window: int = BASE_WINDOW) -> np.ndarray:
    """Counts for consecutive non-overlapping windows: shape (n_tiles, 136)."""
    n_tiles = codes.size // window
    if n_tiles == 0:
        return np.zeros((0, N_CANON), dtype=np.float64)
    trimmed = codes[: n_tiles * window].reshape(n_tiles, window)
    out = np.zeros((n_tiles, N_CANON), dtype=np.float64)
    for i in range(n_tiles):
        out[i] = np.bincount(trimmed[i], minlength=N_CANON)
    return out


def _aggregate(tiles: np.ndarray, factor: int) -> np.ndarray:
    """Merge `factor` adjacent base tiles, giving windows `factor` times longer.

    Dyadic aggregation means the whole ladder costs one pass over the sequence
    rather than one pass per rung.
    """
    n = tiles.shape[0] // factor
    if n == 0:
        return np.zeros((0, N_CANON), dtype=np.float64)
    return tiles[: n * factor].reshape(n, factor, N_CANON).sum(axis=1)


@dataclass
class NullCurve:
    """Divergence expected from native sequence, as a function of length."""

    log_lengths: np.ndarray = field(default_factory=lambda: np.empty(0))
    mu: np.ndarray = field(default_factory=lambda: np.empty(0))
    sigma: np.ndarray = field(default_factory=lambda: np.empty(0))
    gc_sigma: np.ndarray = field(default_factory=lambda: np.empty(0))
    n_windows: np.ndarray = field(default_factory=lambda: np.empty(0))
    theoretical: bool = False

    def __bool__(self) -> bool:
        return self.log_lengths.size > 0

    def at(self, length: int) -> tuple[float, float]:
        """(mu, sigma) of log divergence at this contig length.

        Flat extrapolation past the ends: beyond the observed rungs there is no
        evidence for a trend, and extrapolating one would invent precision.
        """
        if not self:
            return _theoretical_mu_sigma(length)
        x = math.log(max(length, 1))
        mu = float(np.interp(x, self.log_lengths, self.mu))
        sigma = float(np.interp(x, self.log_lengths, self.sigma))
        return mu, max(sigma, 1e-6)

    def gc_sigma_at(self, length: int) -> float:
        if not self or self.gc_sigma.size == 0:
            return 1.0
        x = math.log(max(length, 1))
        return max(float(np.interp(x, self.log_lengths, self.gc_sigma)), 0.05)

    def summary(self) -> dict:
        return {
            "rungs": int(self.log_lengths.size),
            "window_lengths": [int(round(math.exp(x))) for x in self.log_lengths],
            "windows_per_rung": [int(n) for n in self.n_windows],
            "log_divergence_median": [round(float(m), 3) for m in self.mu],
            "log_divergence_mad_sd": [round(float(s), 3) for s in self.sigma],
            "source": "theoretical" if self.theoretical else "within-genome windows",
        }


def _theoretical_mu_sigma(length: int) -> tuple[float, float]:
    """Fallback when a genome cannot supply enough native windows.

    Under multinomial sampling D has expectation (C-1)/n and relative standard
    deviation sqrt(2/(C-1)), so log D has an approximately constant spread.
    Anti-conservative on heterogeneous genomes, which is why callers record that
    it was used.
    """
    n = max(length - K + 1, 1)
    mu = math.log((N_CANON - 1) / n)
    sigma = math.sqrt(2.0 / (N_CANON - 1))
    return mu, sigma


def null_curve(
    tiles_by_contig: list[np.ndarray],
    totals_by_contig: list[np.ndarray],
    *,
    base_window: int = BASE_WINDOW,
    max_rungs: int = MAX_RUNGS,
    min_windows: int = MIN_WINDOWS_PER_RUNG,
) -> NullCurve:
    """Divergence of native windows from their own genome, rung by rung.

    Each window is scored against a centroid that excludes its parent contig, so
    a window is never compared against a reference it helped define -- the same
    leave-one-out convention the per-contig score uses. Without that, long contigs
    would score artificially low and the correction would reintroduce the very
    length bias it removes.
    """
    if not tiles_by_contig:
        return NullCurve(theoretical=True)

    grand = np.sum(totals_by_contig, axis=0)
    log_lengths, mus, sigmas, gc_sigmas, counts = [], [], [], [], []

    for rung in range(max_rungs):
        factor = 2**rung
        window = base_window * factor
        logs: list[float] = []
        gcs: list[float] = []
        for tiles, own in zip(tiles_by_contig, totals_by_contig):
            merged = _aggregate(tiles, factor)
            if merged.shape[0] < 2:
                continue
            ref_counts = grand - own
            if ref_counts.sum() <= 0:
                continue
            ref = frequencies(ref_counts)
            ref_gc = gc_percent_from_freq(ref)
            for row in merged:
                f = frequencies(row)
                d = chi2_divergence(f, ref)
                if d > 0:
                    logs.append(math.log(d))
                    gcs.append(gc_percent_from_freq(f) - ref_gc)
        if len(logs) < min_windows:
            continue
        arr = np.asarray(logs)
        med = float(np.median(arr))
        mad = float(np.median(np.abs(arr - med)))
        sd = 1.4826 * mad
        if sd <= 0:
            sd = float(arr.std()) or 1e-6
        g = np.asarray(gcs)
        gmed = float(np.median(g))
        gsd = 1.4826 * float(np.median(np.abs(g - gmed))) or float(g.std()) or 0.05
        log_lengths.append(math.log(window))
        mus.append(med)
        sigmas.append(sd)
        gc_sigmas.append(gsd)
        counts.append(len(logs))

    if not log_lengths:
        return NullCurve(theoretical=True)
    return NullCurve(
        log_lengths=np.asarray(log_lengths),
        mu=np.asarray(mus),
        sigma=np.asarray(sigmas),
        gc_sigma=np.asarray(gc_sigmas),
        n_windows=np.asarray(counts),
    )


def fwer_threshold(n_tests: int, alpha: float = 0.05) -> float:
    """Bonferroni-corrected z cutoff: the standard normal quantile of 1 - alpha/N.

    Replaces a fixed threshold of 4.0, which meant a different error rate on a
    3-contig genome than on a 1243-contig one -- so the same assembly, shredded
    finer, was held to a progressively laxer standard per contig and a stricter
    one overall. Stating alpha makes the trade explicit and makes the cutoff move
    with the number of tests, as it must.

    Deterministic bisection on math.erf; no scipy.
    """
    n = max(int(n_tests), 1)
    target = 1.0 - alpha / n
    lo, hi = 0.0, 40.0
    for _ in range(200):
        mid = (lo + hi) / 2
        cdf = 0.5 * (1.0 + math.erf(mid / math.sqrt(2.0)))
        if cdf < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2
