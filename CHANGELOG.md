# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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

[Unreleased]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/Hasnat-GPU/GenomeX/releases/tag/v0.1.0
