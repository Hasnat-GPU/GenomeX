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

Results: [`docs/benchmark-busco.md`](../docs/benchmark-busco.md).

BUSCO output (`busco_out/`) is not committed; it is regenerable and large.
