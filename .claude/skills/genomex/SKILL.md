---
name: genomex
description: Run and interpret the GenomeX genome-QC pipeline on bacterial (or fungal) assemblies - assembly stats, gene calling, BUSCO-style completeness, contamination detection, ANI, and core/accessory gene-content comparison between isolates. Use when asked to QC a genome, check completeness or contamination, compute ANI, compare two genomes' gene content, or ask why two isolates differ.
---

# GenomeX

One command takes assemblies from FASTA to an evidence-linked answer for three
questions: *how complete is this genome*, *is more than one organism in it*, and
*why do two isolates from the same environment carry different genes*.

## Before running

GenomeX shells out to Prodigal, HMMER, MMseqs2 and fastANI. On a Windows machine
those live in a micromamba environment inside WSL, not on Windows itself, so
invoke through the launchers rather than calling Python directly:

```bash
./gx.sh  run A.fna B.fna --outdir runs/NAME     # Git Bash, WSL, or Linux
.\gx.ps1 run A.fna B.fna --outdir runs/NAME     # PowerShell
```

Both resolve the repository path themselves and hop into WSL when needed.
Override `GENOMEX_WSL_DISTRO` (default `Ubuntu-24.04`) and `GENOMEX_ENV`
(default `gx`) if yours differ. Inside the environment the plain form works:
`micromamba run -n gx python -m genomex ...`.

Check these before promising a result:

| Precondition | Check | If it fails |
|---|---|---|
| Tools present | `micromamba run -n gx which prodigal hmmsearch mmseqs fastANI` | `micromamba create -y -n gx -c conda-forge -c bioconda prodigal hmmer mmseqs2 mash fastani` |
| Marker set present | `ls $HOME/genomex-work/db/bacteria_odb10/hmms \| wc -l` → 124 | download from `https://busco-data.ezlab.org/v5/data/lineages/` |
| Right lineage for the organism | bacteria → `bacteria_odb10`; fungi → `fungi_odb10` | a bacterial marker set scored against a fungal genome reports fake incompleteness — stop and say so |
| Genetic code | default 11 (bacteria); Mycoplasma/Spiroplasma need `--genetic-code 4` | wrong code inflates fragmented/missing counts |

## Running

```bash
# full pipeline: QC each genome, cluster proteomes, compare every pair
python -m genomex run A.fna B.fna --outdir runs/NAME --all-pairs

# QC only, no clustering or ANI
python -m genomex qc A.fna --outdir runs/NAME

# specific pairs out of a larger set
python -m genomex run *.fna --outdir runs/NAME --pair A B --pair A C
```

Useful flags: `--min-contig-length` (default 3000, floor for composition
analysis), `--min-seq-id` / `--coverage` (MMseqs2 clustering, default 0.5 / 0.8),
`--threads`, `--genetic-code`.

Outputs under `--outdir`: `report.html` (standalone), `genomex.json` (everything),
`provenance.json` (every tool invocation and version), `genomes/<name>/contigs.tsv`
(per-contig evidence), `pairs/<a>__vs__<b>/unique_to_*.tsv` (per-gene tables).

Runtime is roughly one minute per 7 Mb genome plus clustering. Re-running into the
same `--outdir` reuses existing gene calls and marker tables, so iterate freely.

## Interpreting the output

**Completeness** (`markers.completeness_percent`). Fraction of the lineage's
single-copy markers found complete. Above 95% is a good bacterial draft; below 90%
means genuinely missing core genes *or* the wrong lineage — check which before
reporting. This is a HMMER scan against BUSCO's odb10 profiles using BUSCO's own
cutoffs, plus its retention and coverage rules; call it "BUSCO-style", never
"BUSCO said".

**Duplication** (`markers.duplication_percent`). Markers found more than once,
after BUSCO's 85%-of-best-hit retention rule has removed distant paralogs. On
clean finished genomes expect 0–2%; anything above a few percent is worth
explaining. What matters for contamination is the *cross-contig* subset,
reported separately as `duplicated_markers_cross_contig`.

These numbers match BUSCO 5.8.3 on every marker of the demo set
(`docs/benchmark-busco.md`). Still say "BUSCO-style", not "BUSCO said" — it is a
reimplementation that agrees, not the same program.

**Contamination verdict**: `clean` / `possible` / `likely` / `undetermined`.
Read the `reasons` list, never the verdict alone. Two distinctions the pipeline
makes and you must preserve when reporting:

- `contaminant_candidate` — compositionally foreign *and* carrying single-copy
  markers that duplicate copies elsewhere. This is the real signal.
- `replicon_candidate` — large, compositionally distinct, but no displaced
  markers. In multipartite genomes (Paraburkholderia, Cupriavidus, Vibrio) this is
  usually a plasmid or a second chromosome. Composition cannot separate a
  megaplasmid from a contaminant; say so rather than picking one.
- `undetermined` — fewer than 8 scored contigs, so no outlier statistics are
  possible. A finished 3-replicon genome lands here. Report per-replicon
  composition instead of a verdict.

**ANI**. ≥95% is the operational species boundary. `null` means fastANI found too
little to align (below roughly 80% identity) — that means *distant*, not
identical, and the report says so.

**Why two genomes differ**. Each strain-unique gene carries one of four labels:

| Label | Means | How to report it |
|---|---|---|
| contamination-suspect contig | the gene sits on a contig flagged as foreign | not evidence of biology; exclude before counting |
| clustered in genomic island | part of a run of consecutive unique genes | acquisition candidate — prophage, ICE, symbiosis island |
| isolated gene with atypical GC | composition differs from the genome | acquisition candidate, weaker evidence |
| isolated gene, typical composition | neither | divergence, gene loss in the other isolate, or a clustering artefact |

**Check the island null test before repeating any island claim.** Each unique-gene
block carries `island_null_test`: observed island genes, the number expected when
the same unique genes are permuted across the same contigs, and the ratio. When a
large share of genes is strain-unique, runs of three arise by chance — on the demo
data the same-species pair enriches 1.85–2.43x (real), while cross-genus pairs sit
at 1.00–1.01x, meaning every "island" there is an artefact of density. If
`informative` is false, report the counts and say the clustering is
indistinguishable from chance; do not call them acquisitions.

A bare gene-count difference is not a finding. "Strain A carries 180 genes absent
from B, 140 of them in 6 islands at 2.4x the permutation null" is one.

## Abstain rather than guess

Say the analysis cannot answer, and why, when:

- the marker lineage does not match the organism's domain;
- fewer than 8 contigs pass `--min-contig-length` and a contamination verdict is
  requested (report composition per contig instead);
- completeness is below 90% and a comparative claim is requested — missing genes
  masquerade as strain-unique genes in the *other* genome;
- ANI is `null` and a same-species claim is requested;
- a *function* is requested for a unique gene. GenomeX does not annotate function.
  Unique means "no sequence cluster partner", not "novel" and not "characterised".

## What GenomeX does not do

No read assembly (input is assemblies), no functional annotation, no taxonomic
naming of a contaminant, no rRNA/tRNA calling, no plasmid-vs-chromosome
classification beyond the composition heuristic above, and no eukaryote-aware gene
calling — Prodigal is prokaryotic, so a fungal genome needs a different caller
before these numbers mean anything.
