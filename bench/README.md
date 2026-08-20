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
micromamba run -n checkm2 checkm2 predict --threads 8 \
    --input ~/genomex-work/data/Bacteria --extension .fna \
    --output-directory ~/genomex-work/checkm2_out

python -m genomex qc ~/genomex-work/data/Bacteria/*.fna --outdir ~/genomex-work/runs/sweep
python bench/compare_to_checkm2.py \
    --genomex  ~/genomex-work/runs/sweep \
    --checkm2  ~/genomex-work/checkm2_out/quality_report.tsv \
    --out      docs/benchmark-contamination.md
```

The database download is about 2.9 GB.

This comparison deliberately separates two things. Genome size, contig count,
GC, N50 and CDS count are computed from the same FASTA by both tools and **must**
agree -- a gap there is a bug in GenomeX, and it is the only independent check on
`fasta.py` and `genes.py` in the repository. Completeness and contamination are
estimates from different methods on different bases; they should track each
other, and exact equality would be suspicious rather than reassuring.

Results: [`docs/benchmark-contamination.md`](../docs/benchmark-contamination.md).

## Marker scan against BUSCO, on a proteome

```bash
micromamba run -n busco busco -i proteome.faa -l ~/genomex-work/db/fungi_odb10 \
    -o scer_prot --out_path ~/genomex-work/busco_fungi -m proteins --offline

python bench/compare_proteome_to_busco.py \
    --proteins ~/genomex-work/data/Fungi/GCF_000146045.2.faa \
    --lineage  ~/genomex-work/db/fungi_odb10 \
    --busco    ~/genomex-work/busco_fungi/scer_prot \
    --workdir  ~/genomex-work/fungal_probe
```

`compare_to_busco.py` above compares whole GenomeX runs, so it can only reach
lineages the pipeline runs end to end — today, bacteria. This one isolates the
marker scan: a protein FASTA and a lineage, scored both ways and diffed marker
by marker. That is what made it possible to validate `fungi_odb10` while GenomeX
still has no eukaryotic gene calling at all.

**758/758** on the *S. cerevisiae* S288C proteome, every marker naming the same
protein. It was 757/758 before: the one disagreement turned out to be a BUSCO
rule this code never had — a protein counts only under the marker it scores
highest against — which is inert on the bacterial collection and decisive on a
lineage full of paralogous families. See `docs/decisions.md` #15.

Exits non-zero on any disagreement, so it can be used as a check rather than
read by eye.

The `proteins` subcommand now reaches the same scan from the CLI, and reproduces
this result exactly — same statuses, same named proteins, 758/758 against BUSCO
5.8.3:

```bash
python -m genomex proteins ~/genomex-work/data/Fungi/GCF_000146045.2.faa \
    --lineage ~/genomex-work/db/fungi_odb10 \
    --outdir  ~/genomex-work/runs/scer_prot
```

The two differ in one column. This script let `parse_domtbl` fall back to
`protein_id.rsplit("_", 1)[0]` for the contig, which on NCBI identifiers put 754
of 777 hits on a fictitious contig named `NP`; the CLI path records `-` and sets
`contigs_known=False`. Keep both: the script isolates the scan for a lineage the
pipeline has no gene caller for, and it is the only place a new lineage can be
validated before anything else consumes it.

## Fragmentation ladder — constructed specificity

```bash
python bench/fragmentation_ladder.py \
    --genomes ~/genomex-work/data/Bacteria/GCA_000069785.1.fna ... \
    --run-dir ~/genomex-work/runs/sweep \
    --seeds 8 \
    --out bench/fragmentation_seeds_after.tsv
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

**Use `--seeds`.** Those 15-draw figures are one shredding per rung, and a rung
that fires on one shredding in five looks identical to a rung that never fires
when you only shred once. Eight seeds across the eight finished genomes in the
collection gives 320 draws, and 320 draws showed the post-rebuild detector still
producing a 12.5% false-positive rate that climbed from 1.6% at 10 pieces to
26.6% at 1000 — a defect the single-draw ladder had recorded as two anomalies.
Zero afterwards. Results and mechanism:
[`docs/benchmark-fragmentation.md`](../docs/benchmark-fragmentation.md).

`fragmentation_before.tsv` and `fragmentation_after.tsv` are the original
single-draw pair, kept as the record of the null rebuild.
`fragmentation_seeds_before.tsv` and `fragmentation_seeds_after.tsv` are the
320-draw pair, and their `before` half was produced from a clean checkout of the
v0.2.0 tag driven by the current harness, so the two halves differ in one thing.

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

