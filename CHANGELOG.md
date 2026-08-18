# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.1] - 2026-08-18

### Fixed

- **Marker duplication was overcounted, which propagated into contamination
  verdicts.** Benchmarking against BUSCO 5.8.3 showed GenomeX calling markers
  *Duplicated* that BUSCO called *Complete* — 9 of 124 in one genome even when
  both tools scored identical proteins. Two causes:
  - BUSCO retains only matches within 85% of a marker's best hit; GenomeX
    counted every hit above the score cutoff. One marker had a hit at 296.8 and
    six between 17 and 30, all counted as copies.
  - Completeness must be measured on the HMM profile axis, not the sequence
    envelope, which can extend well past the modelled region.

  Duplication rates on the demo set fell from 7–10% to 0–1.6%, and three of four
  genomes moved from *possible contamination* to *clean*. Agreement with BUSCO is
  now 496/496 markers in both `-m genome` and `-m proteins` mode.

### Added

- `bench/` — a reproducible BUSCO comparison: `run_busco.sh` runs BUSCO in both
  modes over the same local lineage, `compare_to_busco.py` diffs the calls marker
  by marker and writes `docs/benchmark-busco.md`. Its parsers are tested.
- `genomes/<name>/markers.tsv` — per-marker class, protein and contig, so
  disagreement can be located rather than merely totalled.

## [0.1.0] - 2026-08-18

First working version: assemblies in, an evidence-linked report out.

### Added

- Assembly statistics (N50/L50/N90, GC, contig distribution) with a SHA-256 of
  every input recorded in the run provenance.
- Gene prediction through Prodigal, with every protein traceable to its contig.
- BUSCO-compatible completeness: HMMER against OrthoDB v10 profiles using
  BUSCO's own score and length cutoffs, reported as a `C/S/D/F/M` string.
- Contamination detection from four independent signals — canonical
  tetranucleotide composition, GC, single-copy markers displaced across contigs,
  and a deterministic bimodality split — with per-contig evidence in `contigs.tsv`.
- Comparative genomics: fastANI against the 95% species boundary, MMseqs2
  orthogroups partitioned into core/accessory/strain-unique, and a cause
  attached to every strain-unique gene.
- `island_enrichment`: a permutation null for genomic-island calls, holding each
  contig's gene and unique-gene counts fixed.
- Standalone HTML report, `genomex.json`, and `provenance.json` recording every
  external invocation with argv, exit code, wall time and tool version.
- A Claude Code skill (`.claude/skills/genomex/`) covering preconditions,
  interpretation, and the conditions under which the pipeline should abstain.
- 40 tests: unit, ground-truth contamination validation, and end-to-end runs
  against the real toolchain.

### Notes on two corrections that real data forced

- Below 8 scored contigs, median/MAD z-scores are meaningless — a finished
  3-replicon genome produced z = 74. Absolute thresholds are used instead and
  the verdict is explicitly marked weak.
- A large, compositionally distinct contig carrying no displaced markers is
  called a `replicon_candidate`, not a contaminant: on the demo data those are
  the symbiosis megaplasmids of *P. phenoliruptrix*, *P. phymatum* and
  *C. taiwanensis*.

[Unreleased]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.1.1...HEAD
[0.1.1]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Hasnat-GPU/GenomeX/releases/tag/v0.1.0
