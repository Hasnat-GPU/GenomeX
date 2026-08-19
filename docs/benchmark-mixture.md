# Contamination recall against constructed mixtures

Every other contamination measurement in this repository either has no foreign
sequence in it or borrows its truth from another tool. `bench/fragmentation_ladder.py`
shreds a genome and checks the verdict holds — a specificity test, with nothing
foreign ever added. `docs/benchmark-contamination.md` compares against CheckM2,
but only four of those 72 genomes exceed 5% contamination, so recall there is an
average of four coin flips, and CheckM2 is a held-out test set that should be
spent once rather than developed against.

This page is the missing instrument: recall measured against truth the
repository constructs and therefore knows exactly.

## How the mixtures are built

Take a finished genome, shred it into 300 fragments (~29 kb each, an ordinary
draft contig), then splice in fragments cut from a **different** genome until
they account for a stated share of the total. Every contig is labelled by
construction, so recall and precision are exact rather than estimated.

Two choices decide whether the measurement means anything:

- **Donor fragments match host fragments in length.** Real contaminants often
  assemble shorter, because they are sequenced at lower depth. If donor
  fragments were shorter here, the length-conditioned null would have a shortcut
  to exploit and the result would flatter the detector. Matching the lengths
  removes it and measures composition alone.
- **Donor fragments carry their share of the single-copy core**, because they
  are cut from a real genome rather than simulated. A 5% mixture displaces part
  of the marker set onto foreign sequence. Withholding that would test a
  detector nobody runs.

The donor's share is a share of the *mixture*, not of the host — the donor mass
solves `d/(h+d) = pct/100`. Taking a percentage of the host instead delivers
16.7% when the row says 20%, which is a mistake this harness made until
`tests/test_bench_mixture.py` caught it.

Four host/donor pairs, chosen to span donor distance:

| host | donor | relationship | GC offset |
|---|---|---|---|
| GCF_000020045.1 *P. phymatum* | GCF_001449005.1 *P. caribensis* | same genus | 0.3 pt |
| GCF_000020045.1 *P. phymatum* | GCA_040965455.1 *Trinickia* sp. | same family | 0.7 pt |
| GCF_000020045.1 *P. phymatum* | GCA_001598055.1 *C. nantongensis* | same family | 4.4 pt |
| GCF_013366925.1 *P. youngii* | GCA_900249745.1 *C. phytorum* | same family | 4.5 pt |

## Result: genome-level detection

Pairs where the verdict is not `clean`, by donor share of the assembly:

| donor % | `bacteria_odb10` (124 markers) | `burkholderiales_odb10` (688 markers) |
|---|---|---|
| 0 (control) | 4/4 clean | 4/4 clean |
| 1 | 1/4 detected | 1/4 |
| 2 | 0/4 | **4/4** |
| 5 | 3/4 | **4/4** |
| 10 | 4/4 | 4/4 |
| 20 | 4/4 | 4/4 |

The controls matter as much as the detections: at 0% donor the same host, shredded
the same way, stays `clean` in every pair under both marker sets. The gain at 2%
and 5% is not bought with false positives.

Two changes produced that column, and they are separable:

1. **A lineage-specific marker set.** All 72 genomes in the collection are
   Burkholderiaceae, so `burkholderiales_odb10` applies and carries 688 markers
   against `bacteria_odb10`'s 124. One duplicated marker is then 0.15% of the
   core rather than 0.81%, which is what makes a 2% mixture expressible at all.
   The lineage is a CLI argument (`--lineage`); GenomeX does not guess it.
2. **A verdict rule stated as a rate.** Duplicated markers were previously
   *counted*: `cross_contig_dups >= 5` meant `likely`, a cutoff written against
   124 markers. Multiplying the marker set by 5.5 multiplies every count by 5.5
   without adding a base of contamination, and four finished reference genomes
   carrying an ordinary ~1% of cross-replicon paralogy jumped straight to
   `likely`. The thresholds are now the old cutoffs divided by 124 — a change of
   units, not a retuning. On `bacteria_odb10` the rule is unchanged.

## Result: localisation, and where it fails

Naming *which* contigs are foreign is a different question from detecting that
some are, and the answer depends entirely on donor distance.

| pair | GC offset | recall (bp) at 2% | at 5% | at 10% |
|---|---|---|---|---|
| *P. phymatum* ← *P. caribensis* | 0.3 pt | 0.00 | 0.00 | 0.00 |
| *P. phymatum* ← *Trinickia* sp. | 0.7 pt | 0.00 | 0.00 | 0.00 |
| *P. phymatum* ← *C. nantongensis* | 4.4 pt | 1.00 | 0.91 | 0.98 |
| *P. youngii* ← *C. phytorum* | 4.5 pt | 0.21 | 0.00 | 0.08 |

At a 4-point GC offset the donor is located almost perfectly and at precision
1.00 up to 5%. At a 0.3-point offset it is never located at all — not at 20%
donor, where the genome-level verdict is emphatically `likely` and 133 markers
are duplicated across contigs.

This is a structural limit, not a threshold that needs moving. Contig calls are
gated on compositional evidence, and a same-genus donor at matched GC and matched
fragment length has none to give: its tetranucleotide composition sits inside the
host's own spread. The duplicated markers prove *a* second organism is present;
they do not say which fragment carries it, because each duplicated marker puts
one copy on a host fragment and one on a donor fragment without labelling which
is which.

The fourth pair shows this is not purely about GC — a 4.5-point offset with a
different host localises poorly. Composition distance between these particular
genomes, not the GC summary, is what decides it.

## What this means for using the tool

- **The verdict is trustworthy well below where it used to be.** 2% mixtures are
  detected on 4/4 pairs with a lineage-specific marker set, against 0/4 before.
- **`suspect_contigs` is a lower bound, always.** When it is empty and the
  verdict is not `clean`, the evidence is marker duplication, which is real but
  does not localise. Do not read an empty suspect list as "nothing foreign".
- **Do not read `cross_contig_duplication_percent` as a contamination
  percentage.** A 5% mixture produced 4.8% to 11.6% of it across these pairs, and
  clean finished genomes sit near 1% from paralogy alone.
- **A lineage-specific set is worth using when the lineage is known.** It is a
  strict improvement here, but the caller must name it — guessing the lineage
  would violate the abstention rule the rest of the tool holds to.

## Reproducing

```bash
python -m genomex qc ~/genomex-work/data/Bacteria/*.fna \
    --lineage ~/genomex-work/db/burkholderiales_odb10 \
    --outdir  ~/genomex-work/runs/sweep_deep

python bench/mixture_ladder.py \
    --data-dir ~/genomex-work/data/Bacteria \
    --run-dir  ~/genomex-work/runs/sweep_deep --total-markers 688 \
    --pair GCF_000020045.1:GCF_001449005.1 \
    --pair GCF_000020045.1:GCA_040965455.1 \
    --pair GCF_000020045.1:GCA_001598055.1 \
    --pair GCF_013366925.1:GCA_900249745.1 \
    --out bench/mixture_ladder_burkholderiales_odb10.tsv
```

Full tables: [`bench/mixture_ladder_burkholderiales_odb10.tsv`](../bench/mixture_ladder_burkholderiales_odb10.tsv)
and [`bench/mixture_ladder_bacteria_odb10.tsv`](../bench/mixture_ladder_bacteria_odb10.tsv).
