# Fragmentation invariance — specificity against constructed truth

Shred one finished genome into 10, 50, 100, 300 and 1000 contigs. Same DNA, same
organism, one composition; only the contig count changes. A contamination
detector must return `clean` at every rung, so every non-clean verdict here is a
false positive with no innocent reading available.

The truth costs nothing and, unlike the CheckM2 comparison, is not borrowed from
another tool. That makes this the one contamination benchmark in the repository
that may safely be *developed against*: the right answer follows from how the
input was built, not from a reference implementation that could itself be wrong.

```bash
python bench/fragmentation_ladder.py \
    --genomes ~/genomex-work/data/Bacteria/GCA_000069785.1.fna ... \
    --run-dir ~/genomex-work/runs/sweep \
    --seeds 8 \
    --out bench/fragmentation_seeds_after.tsv
```

## One draw is not a rate

The previous version of this measurement ran one shredding per rung on one
genome: 15 draws, of which 13 came back clean. That was reported as fragmentation
invariance restored, and the two rungs that flipped were written down as
"residual anomalies" to look at later.

They were not anomalies. They were the visible part of a rate.

A rung that fires on one shredding in five is indistinguishable from a rung that
never fires when the rung is only ever shredded once. Repeating each rung with
eight independent seeds, across the eight finished genomes in the collection with
five or fewer contigs, turns 15 draws into 320 and the anecdote into a number.

| | before | after |
|---|---|---|
| draws | 320 | 320 |
| false positives | **40 (12.5%)** | **0 (0.0%)** |
| 10 pieces | 1/64 (1.6%) | 0/64 |
| 50 pieces | 2/64 (3.1%) | 0/64 |
| 100 pieces | 6/64 (9.4%) | 0/64 |
| 300 pieces | 14/64 (21.9%) | 0/64 |
| 1000 pieces | 17/64 (26.6%) | 0/64 |
| genomes with at least one false positive | 7 of 8 | 0 of 8 |

The "before" column is not the old single-draw run re-read. It is a fresh 320-draw
run of `git archive HEAD` — the released v0.2.0 detector — driven by a
byte-identical copy of the harness, so exactly one thing differs between the two
columns. Per-draw records: [`../bench/fragmentation_seeds_before.tsv`](../bench/fragmentation_seeds_before.tsv)
and [`../bench/fragmentation_seeds_after.tsv`](../bench/fragmentation_seeds_after.tsv).

The rate rising from 1.6% to 26.6% down the ladder is the part that matters. A
flat 12.5% would be a calibration error. A rate that climbs with fragmentation is
the defect the ladder was built to catch, still present after the rebuild that was
supposed to remove it, and invisible to a single draw.

## Every one of the 40 came through the same branch

All 40 false-positive draws contained at least one contig called
`contaminant_candidate` by this arm of the call tree:

> the group is compositionally distinct **and** carries core single-copy markers
> at roughly the genome-wide rate, therefore it is chromosomal sequence,
> therefore — since plasmids do not carry the universal core — it is a second
> organism.

Every step of that is true except the conclusion. Chromosomal sequence in an
assembly of one organism is *that organism's* chromosome.

`GCA_000300095.1` shredded into 300 pieces, seed 0, is the whole failure in one
frame. Nineteen fragments totalling 625 kb were flagged as compositional
outliers, mostly pieces of the genuine 785 kb low-GC replicon that the detector
correctly reports as a `replicon_candidate` at every other rung. At this rung the
union-find merged one further fragment into the same group: `CP003863.1_p155`,
41 kb of ordinary chromosome carrying **34 core markers**. The group's marker
count jumped from near zero to 38 against an expectation of 10, the group cleared
"carries core markers at the genome-wide rate", and all nineteen contigs — the
megaplasmid pieces included — were relabelled contamination. 8.17% of the
assembly, verdict `likely`, on a finished reference genome.

Not one of those 38 markers was duplicated anywhere else in the assembly.

## The rule that replaced it

A second organism does not merely carry core genes. It carries its **own copies**
of a set defined to occur once per genome, and those copies duplicate the host's.
Duplication is the signal; presence is not.

So a compositionally distinct group is judged on which of the two it shows:

| group holds | reading | call |
|---|---|---|
| core markers duplicated elsewhere in the assembly | a second organism | `contaminant_candidate` |
| the assembly's **only** copies of core markers | this organism's own chromosome | `atypical_host_region` |
| no core markers, and ≥ 50 kb of mass | plasmid or chromid, unresolvable by composition | `replicon_candidate` |

The first row is the branch that already existed and that the mixture ladder has
ground truth for. The second is new. It is an abstention in the repository's usual
sense: the composition is odd, the tool says so, and it declines to call odd
composition a second genome when the marker evidence says the sequence is native.

This is a change of inference, not of threshold. No cutoff moved, and on any
assembly where the flagged group really does hold duplicate copies of the core the
behaviour is exactly what it was.

## The detector did not go quiet

Zero false positives would be easy to reach by flagging nothing. That is not what
happened. In the same 320 post-fix draws:

- 292 draws still find at least one compositional outlier group;
- 230 draws still report `replicon_candidate` contigs, 2099 of them in total.

The same 19 contigs of `GCA_000300095.1` at 300 pieces are still flagged, still
listed, still carry their z-scores and their marker counts. They are called
`atypical_host_region` instead of `contaminant_candidate`, and the genome verdict
is `clean`. Nothing was suppressed; one inference was corrected.

## What it cost on the recall side

Nothing that `bench/mixture_ladder.py` can measure.

Both mixture ladders were re-run against the fix — four host/donor pairs at 0, 1,
2, 5, 10 and 20% donor by mass, under `bacteria_odb10` and
`burkholderiales_odb10`. Every figure is **bit-for-bit identical** to the v0.2.0
run: same verdicts, same contig recall, same bp recall, same precision, same
host false positives, same duplicated-marker counts, in all 48 rows.

`donor_as_host_region` — donor contigs excused as the host's own atypical
chromosome, the exact failure this change could have introduced — is **0 in every
row of both files**. That is what the theory predicts: donor fragments cut from a
real genome carry core markers that duplicate the host's, so they meet the
duplication branch before the sole-copy branch is ever reached.

A cost that cannot be shown to exist is not the same as no cost. The honest
statement is that on 48 constructed mixtures spanning two marker sets and a
20-fold range of contamination, this change moved nothing, and the case it would
hurt — a contaminant so scarce or so divergent that it contributes core markers
without duplicating any of the host's — is not represented in this ladder.

## What this benchmark cannot tell you

- **It is silent on recall.** Nothing foreign is ever added. Detection is measured
  by [`benchmark-mixture.md`](benchmark-mixture.md), against sequence spliced in
  at a known share.
- **It is silent on real drafts.** A shredded finished genome has uniform
  coverage, no assembly error and no chimeric contigs. It bounds the false
  positives caused by fragmentation alone, not the false positives of a real
  assembly.
- **Eight genomes, one family.** All eight are Burkholderiaceae, because that is
  the collection this repository has. The 0/320 is a statement about these eight.
- **`atypical_host_region` is a candidate, not a finding.** The class says the
  markers on that sequence are the assembly's only copies. It does not establish
  that the region is a prophage or an acquired island; nothing here runs a null
  test on that, and the island null in the comparative step is a separate
  instrument answering a separate question.
