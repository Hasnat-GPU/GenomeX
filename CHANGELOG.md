# Changelog

All notable changes to this project are documented here. Format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); this project uses
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Nothing yet.

## [0.2.0] - 2026-08-20

### Added

- **`bench/mixture_ladder.py` — constructed ground truth for contamination
  recall.** Splices fragments of one genome into another at a known share of the
  total, so recall and precision are exact rather than borrowed. Donor fragments
  match host fragments in length (otherwise length is a shortcut the
  length-conditioned null could exploit) and carry their share of the single-copy
  core (otherwise it tests a detector nobody runs).

  Measured, four host/donor pairs, by donor share of the assembly: 2% detected in
  0/4 pairs with `bacteria_odb10` and 4/4 with `burkholderiales_odb10`; 5% in 3/4
  and 4/4; 0% controls clean in all pairs under both.

  Localisation is unchanged and is the standing limit: a donor a genus away at a
  4-point GC offset is named contig by contig at recall 0.91–1.00, while a
  same-genus donor at 0.3 points is never localised, even at 20% of the assembly.
  Evidence in `docs/benchmark-mixture.md`.
- `docs/benchmark-contamination-burkholderiales.md` — the 72-genome CheckM2
  comparison repeated against a 688-marker lineage-specific set. Agreement does
  not improve: recall at 5% stays 0.25 under both. The three genomes CheckM2 puts
  above 5% duplicate 0.15–0.58% of the core, at or below the ~1% finished
  reference genomes show from ordinary paralogy.
- The CheckM2 page names the marker set it was scanned against, and hand-written
  commentary below a marker now survives regeneration — a rerun once deleted 74
  lines explaining which numbers a reader must not trust.

### Fixed

- **Duplicated-marker thresholds are rates, not counts.** `cross_contig_dups >= 5`
  meant `likely`, a cutoff written against `bacteria_odb10`'s 124 markers with
  nothing recording the dependency. Scanning the same finished genome against a
  688-marker set multiplies every count by 5.5 without adding a base of
  contamination, and four reference genomes carrying ~1% cross-replicon paralogy
  were reported `likely` contaminated. The thresholds are the old cutoffs divided
  by 124, `>=` preserved, so on `bacteria_odb10` the rule is unchanged bit for bit.
- `GCA_022041995.1` carried the collection's highest duplication at 5.65% against
  CheckM2's 2.37%. Most of it was paralogy within single contigs, and against the
  lineage-specific set it falls to 2.18% — an artefact of scoring a
  Burkholderiaceae genome against a marker set defined across all bacteria.
- A `nan` recall no longer appears on a generated page when no genome exceeds the
  threshold being tested.

### Changed

- **Contamination scoring is now conditioned on contig length.** A
  tetranucleotide frequency is a multinomial sample whose divergence falls as
  1/length, so one threshold for a 3 kb contig and a 4 Mb chromosome flagged the
  short tail of every fragmented assembly. Each contig is scored against a null
  built from windows of its own length tiled inside the assembly's own contigs;
  the measured slope of log divergence against log window is −0.58 to −0.72,
  not the −1.0 sampling alone predicts, so the curve is measured rather than
  assumed.
- The fixed z of 4.0 is now a Bonferroni cutoff at α = 0.05, stated as an error
  rate and moving with the number of contigs tested.
- Bimodality no longer contributes to the verdict. It fired on 57 of 72 published
  genomes because a 2-means split always returns two bins. Retained as evidence.
- Flagged contigs are grouped into candidate replicons before judgement, so a
  plasmid is one object whether emitted as one 785 kb contig or thirty 27 kb ones.
- The replicon call rests on positive evidence: a group carrying the single-copy
  core at the genome-wide rate is another organism's chromosome; one carrying none
  is extrachromosomal. With no marker scan supplied, no plasmid claim is made.

  Measured on `bench/fragmentation_ladder.py`, which is constructed ground truth
  independent of CheckM2: 13 of 15 rungs clean, against 2 of 15 before, and
  corr(contig score, log length) from −0.43 to −0.09. Out of sample against
  CheckM2: 13 of 72 flagged where 62 were before, verdict medians monotone
  (0.93 / 1.23 / 2.13) where they had been flat and faintly inverted, and recall
  at CheckM2 ≥ 5% down from 1.00 to 0.25. The recall loss is reported alongside
  the precision gain because reporting either alone would misrepresent the change.

### Added

- `bench/fragmentation_ladder.py` — shreds a finished genome to 10/50/100/300/1000
  contigs and asserts the verdict does not move. Development ran against this
  rather than against the CheckM2 numbers, so no threshold is fitted to the
  benchmark it is later measured on.
- `genomex/composition.py` — length-conditioned scoring, the within-genome window
  null, and the recorded refutations of four alternative nulls and three classical
  unimodality statistics that were considered and measured.

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

[Unreleased]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.2.0...HEAD
[0.2.0]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.1.1...v0.2.0
[0.1.1]: https://github.com/Hasnat-GPU/GenomeX/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/Hasnat-GPU/GenomeX/releases/tag/v0.1.0
