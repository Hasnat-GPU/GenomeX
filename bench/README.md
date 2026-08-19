# Benchmarks

## GenomeX vs BUSCO

GenomeX calls its completeness scan *BUSCO-compatible*, not BUSCO. This directory
holds the test of that claim.

```bash
micromamba create -y -n busco -c conda-forge -c bioconda busco
bench/run_busco.sh                                    # real BUSCO, same lineage, offline
python bench/compare_to_busco.py \
    --genomex runs/demo-rhizobia \
    --busco   bench/busco_out \
    --out     docs/benchmark-busco.md
```

Both tools score against the same local `bacteria_odb10` directory, so any
disagreement comes from method rather than from profile or cutoff versions.

The comparison is marker by marker, not headline percentage. Two tools can both
report 100% complete while disagreeing about which markers are duplicated — and
duplication is the number the contamination module consumes, so that is exactly
where a silent disagreement would matter.

Results: [`docs/benchmark-busco.md`](../docs/benchmark-busco.md) — currently
496/496 markers identical in both modes, across four genomes.

That was not the first result. The initial run disagreed on 9 markers in one
genome *with identical gene calls*, which located two real defects in
`markers.py`: the missing 85%-of-best-hit retention rule, and completeness being
measured on the sequence envelope rather than the HMM profile. Both are fixed;
the benchmark is what found them.

BUSCO output (`busco_out/`) is not committed; it is regenerable and large.

## GenomeX vs CheckM2

Completeness now has a reference check; contamination does not. CheckM2 provides
one, over a whole collection rather than four genomes.

```bash
micromamba create -y -n checkm2 -c conda-forge -c bioconda checkm2
micromamba run -n checkm2 checkm2 database --download --path ~/genomex-work/db/checkm2
micromamba run -n checkm2 checkm2 predict --threads 8     --input ~/genomex-work/data/Bacteria --extension .fna     --output-directory ~/genomex-work/checkm2_out

python -m genomex qc ~/genomex-work/data/Bacteria/*.fna --outdir ~/genomex-work/runs/sweep
python bench/compare_to_checkm2.py     --genomex  ~/genomex-work/runs/sweep     --checkm2  ~/genomex-work/checkm2_out/quality_report.tsv     --out      docs/benchmark-contamination.md
```

The database download is about 2.9 GB.

This comparison deliberately separates two things. Genome size, contig count,
GC, N50 and CDS count are computed from the same FASTA by both tools and **must**
agree -- a gap there is a bug in GenomeX, and it is the only independent check on
`fasta.py` and `genes.py` in the repository. Completeness and contamination are
estimates from different methods on different bases; they should track each
other, and exact equality would be suspicious rather than reassuring.

Results: [`docs/benchmark-contamination.md`](../docs/benchmark-contamination.md).

## Fragmentation ladder — constructed specificity

```bash
python bench/fragmentation_ladder.py \
    --genomes ~/genomex-work/data/Bacteria/GCA_000069785.1.fna \
    --run-dir ~/genomex-work/runs/sweep \
    --out bench/fragmentation_after.tsv
```

Shred one finished genome into 10, 50, 100, 300 and 1000 contigs. Same DNA, same
organism, one composition — only the contig count changes, so the verdict must
not. This is ground truth that costs nothing and, crucially, is **independent of
CheckM2**: a threshold tuned after seeing the 72-genome comparison is fitted to
it, whereas this ladder can be developed against safely because the right answer
follows from how the input was built.

It is what exposed the length bias in the old detector: 2/15 rungs stable, and a
correlation of −0.43 between contig score and log length. After the rebuild,
13/15 and −0.09.

## Mixture ladder — constructed recall

```bash
python bench/mixture_ladder.py \
    --data-dir ~/genomex-work/data/Bacteria \
    --run-dir  ~/genomex-work/runs/sweep_deep --total-markers 688 \
    --pair GCF_000020045.1:GCF_001449005.1 \
    --out bench/mixture_ladder_burkholderiales_odb10.tsv
```

The fragmentation ladder measures specificity; nothing foreign is ever added, so
it cannot speak to recall. The CheckM2 comparison does measure recall, but only
four of the 72 genomes exceed 5% contamination — an average of four coin flips —
and CheckM2 is a held-out test set that should be spent once, not developed
against.

This harness supplies the missing instrument. Take a finished genome, shred it,
splice in fragments of a *different* genome until they account for a known share
of the total, and ask the detector to name them. Every contig is labelled by
construction, so recall and precision are exact.

Two choices keep it honest:

- **Donor fragments match host fragments in length.** Real contaminants often
  assemble shorter, and length would then be a shortcut the length-conditioned
  null could exploit. Matching removes it and measures composition alone.
- **Donor fragments carry their share of the single-copy core**, because they
  are cut from a real genome. Withholding that would test a detector nobody runs.

Donor distance is the second axis, and it is where the limit shows: a donor from
another genus at a 4-point GC offset is located contig by contig, while a
same-genus donor at a 0.3-point offset is detected at genome level and never
localised.

