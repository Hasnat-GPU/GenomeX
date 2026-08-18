# Contributing to GenomeX

The useful contributions here are mostly of one kind: **showing that a call is
wrong on real data.** GenomeX makes claims about genomes — this is complete, this
contig is foreign, these genes were acquired — and every one of them is falsifiable.

## What a good contribution looks like

**Best: a genome that breaks a call.** An assembly where the contamination
verdict, the replicon/contaminant split, or the island enrichment is
demonstrably wrong. Open an issue with the accession, the command, the output,
and what the right answer is and how you know. That is more valuable than a
feature, because the thresholds in `contamination.py` were set against four
genomes and a synthetic fixture — a sample of four.

**Also welcome:**

- A test that plants a known answer and shows the pipeline misses it. Add it to
  `tests/` even if you cannot fix it; a failing test with a clear ground truth is
  a complete contribution.
- Agreement or disagreement with a reference implementation — a genome where
  our completeness differs materially from real BUSCO, or our contamination call
  differs from CheckM2. Numbers, not impressions.
- A new lineage or organism class (fungal gene calling is the obvious gap; see
  the limits in the README).
- Documentation that makes a limit clearer.

## Development setup

```bash
micromamba create -y -n gx -c conda-forge -c bioconda \
    prodigal hmmer mmseqs2 mash fastani python numpy pytest
mkdir -p ~/genomex-work/db && cd ~/genomex-work/db
wget https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz
tar xzf bacteria_odb10.2024-01-08.tar.gz
```

```bash
micromamba run -n gx python -m pytest tests/ -q                  # everything
micromamba run -n gx python -m pytest tests/ -q -m "not tools"   # no external tools needed
```

Tests marked `tools` need Prodigal/HMMER/MMseqs2/fastANI on PATH and skip
automatically when they are absent. CI runs both sets.

## Standards this project holds itself to

1. **A verdict ships with its evidence.** Anything that returns a judgement
   returns the reasons alongside it. No bare booleans.
2. **State the null.** If a pattern could arise by chance, compare it to chance
   before naming it. `island_enrichment` exists because "1108 genes in islands"
   was meaningless without knowing that 600 were expected anyway.
3. **Abstain over guess.** Too few contigs, wrong lineage, unresolved ANI —
   return `undetermined` and say why. See the abstention list in the skill file.
4. **Name the method honestly.** It is a BUSCO-*compatible* scan, not BUSCO.
   Every report says so. Do not let a shorthand become a claim.
5. **Determinism.** No unseeded randomness anywhere in a scoring path: the
   2-means split seeds at the extremes, the permutation null takes a fixed seed.
   Two runs on one input must agree byte for byte.

## Pull requests

- Branch from `main`, one concern per PR.
- Add or update tests for behaviour changes; a threshold change needs a test
  that would fail under the old value.
- Run the full suite before pushing.
- Commit subjects in the imperative mood, under ~72 characters. Explain *why* in
  the body — the diff already shows what.
- If a change alters any number in the README or the skill file, update those in
  the same PR. Stale numbers are how a tool starts lying.
