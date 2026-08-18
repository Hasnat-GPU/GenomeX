# GenomeX

[![CI](https://github.com/Hasnat-GPU/GenomeX/actions/workflows/ci.yml/badge.svg)](https://github.com/Hasnat-GPU/GenomeX/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Status: alpha](https://img.shields.io/badge/status-alpha-orange.svg)](CHANGELOG.md)

Assemblies in, an evidence-linked answer out, for the three questions on the
whiteboard:

1. **How complete is this genome?** — single-copy ortholog recovery.
2. **Is there cross-contamination in it?** — composition + displaced markers.
3. **Why do two bacteria from the same environment carry different genes?** —
   ANI, orthogroups, and a cause attached to every strain-unique gene.

Real tools do the science: **Prodigal** (gene calling), **HMMER** against
**BUSCO odb10** profiles (completeness), **MMseqs2** (orthogroups), **fastANI**
(species boundary), **Mash** (optional triage). GenomeX is the layer that runs
them in one pass, joins their outputs, and refuses to overclaim.

## Layout

```
genomex/
  fasta.py          FASTA I/O, N50/L50/GC, sha256 of every input
  runtime.py        tool resolution + provenance (argv, exit code, wall time, version)
  genes.py          Prodigal wrapper; every protein traceable to its contig
  markers.py        BUSCO-compatible marker scan (HMMER + odb10 cutoffs)
  contamination.py  TNF + GC + displaced markers + bimodality  ← hero module 1
  compare.py        ANI, orthogroups, strain-unique genes with causes  ← hero module 2
  pipeline.py       orchestration; contamination feeds the comparison
  report.py         standalone HTML + JSON
tests/              40 tests: unit, ground-truth contamination, and end-to-end with real tools
.claude/skills/genomex/SKILL.md   the unified skill: how to run it and how to read it
```

## Setup

Science runs on Linux; Windows is the cockpit. Tools live in a micromamba
environment inside WSL:

```bash
micromamba create -y -n gx -c conda-forge -c bioconda prodigal hmmer mmseqs2 mash fastani
# marker set (bacteria; 124 profiles, ~2 MB)
mkdir -p ~/genomex-work/db && cd ~/genomex-work/db
wget https://busco-data.ezlab.org/v5/data/lineages/bacteria_odb10.2024-01-08.tar.gz
tar xzf bacteria_odb10.2024-01-08.tar.gz
```

No sudo, no Docker, no conda base install required.

## Run

```bash
python -m genomex run A.fna B.fna --outdir runs/demo --all-pairs
python -m genomex qc  A.fna       --outdir runs/qc          # per-genome only
```

From Windows PowerShell, `.\gx.ps1 run A.fna B.fna --outdir runs/demo` forwards
into WSL with the environment already set.

Outputs: `report.html`, `genomex.json`, `provenance.json`,
`genomes/<name>/contigs.tsv`, `pairs/<a>__vs__<b>/unique_to_*.tsv`.

Roughly one minute per 7 Mb genome. Re-running into the same `--outdir` reuses
gene calls and marker tables.

## Reading the output

`contigs.tsv` gives every contig a `call`:

| call | meaning |
|---|---|
| `core` | composition matches the rest of the assembly |
| `contaminant_candidate` | compositionally foreign **and** carrying single-copy markers duplicated elsewhere |
| `replicon_candidate` | large and compositionally distinct, but no displaced markers — a plasmid or second chromosome is at least as likely as a contaminant |

The `replicon_candidate` class exists because of what real data did to the first
version of this code: the finished *P. phenoliruptrix* BR3459a genome has a
785 kb replicon at 59.1% GC against a 63.5% GC chromosome, and the detector
called it contamination. It is a megaplasmid. Composition cannot tell those
apart, so the pipeline now says so instead of guessing.

Strain-unique genes are labelled `contamination-suspect contig`,
`clustered in genomic island (acquisition/HGT candidate)`,
`isolated gene with atypical GC`, or
`isolated gene, typical composition (divergence or gene loss)`.

Island calls are checked against a permutation null that holds each contig's gene
count and unique-gene count fixed and only reshuffles which genes are unique. On
the demo data:

| pair | ANI | island genes | expected by chance | enrichment |
|---|---|---|---|---|
| *P. phenoliruptrix* BR3459a vs *P. phenoliruptrix* | 98.2% | 1108 | 600 | **1.85x** |
| *P. phenoliruptrix* vs *P. phymatum* | 83.0% | 1772 | 1265 | 1.40x |
| *P. phenoliruptrix* vs *C. taiwanensis* | 77.9% | 4139 | 4135 | 1.00x |

Between genera, "islands" are an artefact of two thirds of the genes being unique.
The report says so instead of reporting 253 acquisition islands.

## Validated against BUSCO

GenomeX agrees with BUSCO 5.8.3 on **496 of 496 markers** across the four demo
genomes, in both BUSCO modes — scoring its own gene predictions and scoring
GenomeX's. Same lineage directory, same HMMER and Prodigal builds, so nothing in
that number comes from version drift.

| genome | GenomeX | BUSCO |
|---|---|---|
| *C. taiwanensis* | `C:99.19%[S:97.58%,D:1.61%],F:0.81%,M:0.0%` | `C:99.2%[S:97.6%,D:1.6%],F:0.8%,M:0.0%` |
| *P. phenoliruptrix* BR3459a | `C:100.0%[S:98.39%,D:1.61%],F:0.0%,M:0.0%` | `C:100.0%[S:98.4%,D:1.6%],F:0.0%,M:0.0%` |
| *P. phymatum* STM815 | `C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%` | `C:100.0%[S:100.0%,D:0.0%],F:0.0%,M:0.0%` |

It did not start there. The first benchmark showed GenomeX calling markers
*Duplicated* that BUSCO called *Complete* — 9 of 124 in one genome, even when
both tools scored identical proteins. Two causes, both now fixed and both
visible in [`docs/benchmark-busco.md`](docs/benchmark-busco.md):

1. **Bitscore retention.** BUSCO keeps only matches within 85% of a marker's
   best hit. One marker had a hit at 296.8 and six between 17 and 30, all above
   a permissive score cutoff; GenomeX counted all seven as copies.
2. **Coverage axis.** Completeness is measured on the HMM profile, not the
   sequence envelope. An envelope can run far past the modelled region, so
   fragments were passing the length test.

Duplication rates fell from 7–10% to 0–1.6%, and three of the four genomes moved
from *possible contamination* to *clean* — the earlier verdicts were driven by
markers that were never duplicated.

## Demo result

Four beta-rhizobia, `./demo.sh`, 23 s warm (2 min cold):

- three genomes clean; each finished one has its symbiosis megaplasmid found and
  called a `replicon_candidate` rather than contamination — CP003865.1 (785 kb,
  59.1% GC), NC_010627.1 (595 kb), CU633751.1 (557 kb, *pRalta*);
- the draft assembly is the only one called `likely` contaminated: 8 contigs,
  2.4% of its bases;
- ANI recovers the species boundary: 98.2% within *P. phenoliruptrix*, 83.0%
  across species, 77.9% across genera;
- 12,701 orthogroups over the four genomes; 1,900 core.

## Tests

```bash
micromamba run -n gx python -m pytest tests/ -q
```

The contamination tests plant a known answer and check it is recovered:

- a one-composition genome must come back **clean with zero false positives**;
- a genome with three planted foreign contigs must have **all three named**;
- a contaminant at the **same GC as its host** must still be found — the case
  where GC screening fails and only k-mer composition can work;
- five markers duplicated across two contigs must **raise a clean verdict to
  likely**, and name both contigs;
- fewer than three scorable contigs must return **undetermined**, not a guess.

## Limits

- Completeness is HMMER against BUSCO's profiles with BUSCO's cutoffs and
  filtering rules. It matches BUSCO on every marker tested so far, but it is a
  reimplementation: it does not re-predict genes per candidate region, so a
  genome outside the tested set may still diverge.
- Contamination detection is composition-based: blind to a close relative of the
  host, and unable to name a contaminant without a reference database.
- Contigs below `--min-contig-length` are not scored.
- No functional annotation. "Unique" means no sequence-cluster partner, not novel
  and not characterised.
- Genomic islands are runs of consecutive unique genes tested against a permutation
  null — candidates, not calls, and only where the null is beaten.
- Prodigal is prokaryotic; fungal genomes need a eukaryote-aware caller first.

See [`docs/example-report.html`](docs/example-report.html) for the full rendered
output of the demo run (download and open it — GitHub does not render HTML inline).

## Status

Alpha, but no longer uncalibrated. Completeness and duplication now match BUSCO
exactly on every marker of the demo set, and that comparison is reproducible
from `bench/`.

What remains unvalidated is the part with no reference implementation to check
against: the contamination thresholds, the replicon/contaminant split, and the
island null model were all set against four genomes and a synthetic fixture. The
algorithms are tested against planted ground truth and the plumbing runs
end-to-end in CI, but four genomes is not a validation set. Contamination has not
been compared against CheckM2.

If you have a genome where a call is wrong, that is the
[most useful contribution](CONTRIBUTING.md) you can make.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md). Wrong calls on real genomes are wanted
more than features. Security issues go through
[private reporting](SECURITY.md), not public issues.

## Citing

If GenomeX contributes to published work, cite the tools that do the science —
[Prodigal](https://doi.org/10.1186/1471-2105-11-119),
[HMMER](https://doi.org/10.1371/journal.pcbi.1002195),
[BUSCO/OrthoDB](https://doi.org/10.1093/molbev/msab199),
[MMseqs2](https://doi.org/10.1038/nbt.3988) and
[fastANI](https://doi.org/10.1038/s41467-018-07641-9) — alongside this
repository. See [CITATION.cff](CITATION.cff).

## License

[MIT](LICENSE). The bundled dataset metadata under `data/` comes from NCBI and
carries its own terms; the BUSCO `bacteria_odb10` marker set is downloaded at
setup time from [OrthoDB](https://busco-data.ezlab.org/) and is not redistributed
here.
