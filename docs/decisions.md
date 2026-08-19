# Design decisions

Why the code is shaped this way, and what was rejected. Git history records what
changed; this records why, and what would have to be true to change it back.

Most entries exist because a first version was wrong on real data. That is the
useful part — the rejected option is usually the obvious one.

---

## 1. Verdicts carry their evidence

`ContaminationResult` returns `reasons`, `contigs.tsv` gives every contig its
flags, and `provenance.json` records every external invocation with argv, exit
code and version.

**Rejected:** returning a category alone. In genomics a wrong answer is
byte-identical in shape to a right one, so a reader who cannot audit the
reasoning cannot detect the failure. Every published verdict here has been wrong
at least once; the evidence is what let it be caught.

## 2. A pattern is not a finding until it beats chance

`compare.island_enrichment` permutes which genes are strain-unique, holding each
contig's gene and unique-gene counts fixed, and reports observed against
expected.

The first version reported 1108 unique genes in 73 genomic islands and called
them acquisitions. The null expects 600 by chance. Across genera the ratio is
1.00 — every island was an artefact of two-thirds of genes being unique.

**Consequence, learned late:** the contamination bimodality rule was written
without a null and flags 57 of 72 published assemblies. Same error, one module
over. See `benchmark-contamination.md`.

## 3. Composition cannot separate a plasmid from a contaminant

A large, compositionally distinct contig carrying no displaced single-copy
markers is a `replicon_candidate`, not a contaminant, and is excluded from the
contamination fraction.

**Why:** the finished *P. phenoliruptrix* BR3459a genome has a 785 kb replicon
at 59.1% GC against a 63.5% chromosome. The first version called it
contamination. It is the symbiosis megaplasmid. The same held for *P. phymatum*
(595 kb) and *C. taiwanensis* (*pRalta*, 557 kb).

**Rejected:** picking one interpretation. Composition genuinely cannot
distinguish them, so the class name says both and the report explains the
ambiguity.

## 4. Below eight contigs there are no outlier statistics

Under `min_contigs_for_zscores` the detector uses absolute thresholds and marks
its own verdict weak.

**Why:** a finished 3-replicon genome produced a robust z of 74 from a
near-zero MAD. Three points have no population to be an outlier against.

## 5. Marker duplication does not condemn a contig on its own

A contig with duplicated core markers but typical composition is
`marker_conflict`, not `contaminant_candidate`, and does not enter the
contamination fraction.

**Why:** in multipartite genomes the chromosome and the chromid legitimately
share paralogs of a few core markers. Treating that as foreign made every
*Paraburkholderia* report its own chromosome as a contaminant — and the
comparative module then discarded real genes as artefacts, because it consumes
the suspect-contig set.

The genome-level verdict still rises on cross-contig duplication. The
distinction is between genome-level evidence and per-contig blame.

## 6. BUSCO's filtering rules are part of BUSCO's method

`parse_domtbl` applies the score cutoff, then classifies on **HMM-profile**
coverage, then retains only hits within 85% of each marker's best bitscore.

**Why:** without retention, one marker with a hit at 296.8 and six between 17.4
and 29.7 counted as seven copies. Those are distant paralogs clearing a
permissive cutoff. And coverage measured on the sequence envelope — which runs
past the modelled region — let fragments pass the length test.

Agreement with BUSCO went from 9 disagreements in one genome to 496/496 across
four. Both rules were read out of BUSCO's `hmmer.py`, not inferred.

**Rejected:** implementing the cutoffs from the dataset files alone and assuming
the rest was detail. The filtering *is* the method.

## 7. Pairwise means pairwise

When a joint pangenome spans several genomes, "unique to A" for the pair (A, B)
means present in A and absent from B — not absent from every other genome.

**Why:** the latter silently answers a different question, and the count changes
depending on which unrelated genomes happen to be in the run.

## 8. Completeness counts complete markers only

Fragmented markers are reported separately and never folded into completeness.

**Why:** this is BUSCO's definition and the tool claims BUSCO compatibility.
CheckM2 effectively counts a gene broken across a contig boundary as present,
which is why the two disagree by up to 18.55 points — a gap that tracks the
fragmented count at r = −0.996. Neither is wrong; they answer different
questions. The fix is to report both, not to average them.

## 9. Windows is a cockpit; nothing scientific executes on it

The pipeline runs in a Linux environment; `gx.ps1` and `gx.sh` hop into it.

**Why:** Bioconda has no `win-64`. Building Windows-native would not add
friction, it would make the development environment structurally incapable of
matching the environment that produces publishable results.

## 10. Determinism over marginally better clustering

The 2-means split seeds at the extremes of the first principal component rather
than by k-means++; the permutation null takes a fixed seed.

**Why:** a reproducibility tool whose own output moves between runs is not one.
The marginal quality of a random restart is worth less than a byte-identical
rerun.

## 11. Benchmarks are test sets, not targets

`benchmark-busco.md` and `benchmark-contamination.md` are measurements. Fitting
thresholds until a confusion matrix improves would be training on the test set,
and the resulting number would mean nothing on any other collection.

When a benchmark fails, the method gets fixed and re-measured. If it still
fails, that is the finding, and it gets published — which is why the README
currently tells users not to rely on the contamination verdict.
