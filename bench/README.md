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
