# CLAUDE.md

Guidance for Claude Code working in this repository. Keep this file short — it
loads into every session. Detail belongs in the files it points to.

## What this is

GenomeX takes bacterial genome assemblies and answers three questions with the
evidence attached: how complete the genome is, whether it holds more than one
organism, and why two isolates from one environment carry different genes.

It orchestrates real tools — Prodigal, HMMER against BUSCO odb10 profiles,
MMseqs2, fastANI — and joins their output. It does not reimplement them.

## Running it

```bash
./gx.sh  run A.fna B.fna --outdir runs/NAME --all-pairs   # Git Bash, WSL, Linux
.\gx.ps1 run A.fna B.fna --outdir runs/NAME               # PowerShell
micromamba run -n gx python -m pytest tests/ -q           # inside the env
```

The science runs on Linux (WSL) in a micromamba environment; Windows is only a
cockpit. Local paths and environment names are in `.claude/state/ENVIRONMENT.md`
(untracked). `bench/README.md` covers the benchmarks.

## Standards this code holds itself to

These are not style preferences. They are why the tool can be trusted, and each
one exists because breaking it produced a wrong answer that looked right.

1. **A verdict ships with its evidence.** Anything returning a judgement returns
   the reasons beside it. No bare booleans, no unexplained categories.
2. **State the null.** If a pattern could arise by chance, compare it to chance
   before naming it. Genomic islands have a permutation null because "1108 genes
   in islands" was meaningless until we knew 600 were expected anyway. The
   contamination bimodality rule lacked one and produced an 86% false-positive
   rate — see `docs/benchmark-contamination.md`.
3. **Abstain over guess.** Too few contigs, wrong lineage, unresolved ANI →
   return `undetermined` and say why. The abstention list is in the skill file.
4. **Name the method honestly.** It is a BUSCO-*compatible* scan, not BUSCO.
   Every report says so. A shorthand must never harden into a claim.
5. **Determinism.** No unseeded randomness in any scoring path. Two runs on one
   input agree byte for byte.
6. **Never tune a threshold to a benchmark.** `docs/benchmark-*.md` are test
   sets. The moment one becomes a target it stops measuring anything. Fix the
   method, then re-measure and report whatever comes out.

## Current state of the claims

Read this before repeating any number to a user.

| Claim | Status |
|---|---|
| Assembly stats (size, contigs, GC, N50, CDS, coding density) | **verified** — exact agreement with CheckM2 on 72 genomes |
| Completeness / duplication | **verified** — 496/496 markers identical to BUSCO 5.8.3 |
| Completeness vs CheckM2 residual | explained: tracks fragmented count at r = −0.996 |
| ANI, orthogroups, core/accessory | unvalidated against a reference, but thin wrappers over fastANI/MMseqs2 |
| Genomic island calls | guarded by a permutation null; only trust when `informative` is true |
| **Contamination verdict** | **refuted** — flags 62/72 published assemblies CheckM2 calls clean. Do not rely on it; use `contigs.tsv` |
| replicon vs contaminant split | unvalidated, no reference implementation exists |

## Environment gotchas that have cost real time

- **Variables do not survive `wsl.exe -- bash -lc '...'`.** `D=/path; cmd $D`
  arrives with `$D` empty. Inline full paths instead.
- WSL here has **no passwordless sudo** and **no bzip2**. Install through
  micromamba, bootstrapped by extracting its tarball with Python's `bz2`.
- **Windows Python cannot open `\\wsl.localhost\...`** even though Git Bash `ls`
  can. Anything reading genome files runs Linux-side.
- A `nohup`-ed job launched by a one-shot `wsl.exe` call **dies when that call
  returns**. Keep the launching process alive instead.
- BUSCO writes `busco_downloads/` and `busco_<pid>.log` into its working
  directory. Both are gitignored.
- Heredocs through the Bash tool sometimes mangle long Python files; use the
  Write tool for those.

## Where to look next

- `.claude/state/PROGRESS.md` — what is in flight right now, and what is next.
  **Read it first when resuming work.** Untracked.
- `docs/decisions.md` — why the design is what it is, and what was rejected.
- `docs/benchmark-busco.md`, `docs/benchmark-contamination.md` — the evidence.
- `.claude/skills/genomex/SKILL.md` — how to run and interpret the pipeline,
  including when to refuse to answer.
