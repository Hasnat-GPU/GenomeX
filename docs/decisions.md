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

## 12. Thresholds on a marker set are rates, not counts

The contamination verdict tests `cross_contig_dups / n_markers` against 4.03%
and 1.61%, not the raw count against 5 and 2.

**Why:** those cutoffs were written against `bacteria_odb10` and its 124
markers, and nothing recorded that dependency. Scanning the same finished genome
against `burkholderiales_odb10` — 688 markers, more appropriate for this
collection, not one base of extra contamination — multiplies every count by 5.5,
and four reference genomes carrying an ordinary ~1% of cross-replicon paralogy
were reported as `likely` contaminated.

A count is not a property of the genome; it is a property of the genome and the
marker set together. The rate is the invariant: a mixture holding a foreign
genome at mass fraction *f* displaces roughly *f × N* of *N* single-copy
markers, so the rate estimates *f* whatever *N* is.

The new thresholds are the old ones divided by 124, `>=` preserved, so on
`bacteria_odb10` the rule is unchanged bit for bit. **Rejected:** recalibrating
the thresholds against the collection. That would have been fitting to the test
set, and there was no need — the defect was units, not calibration.

## 13. Recall needs constructed truth, not a borrowed number

`bench/mixture_ladder.py` splices fragments of one genome into another at a
known share of the total and asks the detector to name them.

**Why:** the only recall figure available before came from CheckM2, and just
four of those 72 genomes exceed 5% contamination — an average of four coin
flips, on a test set that should be spent once rather than developed against.
The fragmentation ladder is constructed truth but adds nothing foreign, so it
measures specificity and is silent on recall.

Two construction choices carry the result. Donor fragments match host fragments
in length, because real contaminants assemble shorter and length would otherwise
be a shortcut the length-conditioned null could exploit. Donor fragments carry
their share of the single-copy core, because they are cut from a real genome and
withholding that would test a detector nobody runs.

**Rejected:** simulating a contaminant from a composition model. The whole
question is whether real between-genome composition differences are separable,
and a simulated donor would answer a question about the simulator instead.

## 14. Presence of core markers is not evidence of a second organism

A compositionally distinct group that carries core single-copy markers is called
a `contaminant_candidate` only when those markers are **duplicated** elsewhere in
the assembly. When they are the assembly's only copies it is called an
`atypical_host_region` and excluded from the contamination fraction.

**Why:** the old rule reasoned that plasmids do not carry the universal core, so
a distinct group that does carry it must be chromosomal, and chromosomal sequence
with foreign composition must be a second organism. Every step holds except the
last. Chromosomal sequence in an assembly of one organism is *that organism's*
chromosome, and a prophage, a genomic island or a GC-skewed region is exactly
that: native sequence with atypical composition.

The discriminator was there all along and was being discarded. A second organism
does not merely carry core genes, it carries its **own copies** of a set defined
to occur once per genome, and those duplicate the host's. Presence is not the
signal; duplication is.

The cost of the old reading was measurable and large. The family-wise threshold
admits a compositional outlier on some share of clean assemblies by construction,
and this arm converted every such outlier into a genome-level verdict. Over 320
shreddings of eight finished genomes: 40 false positives, 12.5%, rising from 1.6%
at 10 pieces to 26.6% at 1000. Afterwards, zero — with 292 of the 320 draws still
flagging outlier groups, so nothing was silenced. Both mixture ladders came back
bit-for-bit identical, and no donor contig in 48 constructed rows was excused as a
host region. `docs/benchmark-fragmentation.md`.

It also cost biology. Genes on `contaminant_candidate` contigs are dropped by the
comparative step as artefacts, so a prophage — the textbook acquisition island —
was being deleted from the very analysis meant to find it.

**Rejected:** raising the outlier threshold until the ladder came back clean. The
outliers are real; the assembly does contain compositionally odd sequence, and
suppressing the flag would hide a true observation to fix a false conclusion. The
defect was the inference drawn from the flag, not the flag.

**Rejected:** dropping the arm entirely and letting such groups fall through to
`replicon_candidate`. That label means "plasmid or contaminant, unresolvable by
composition", and here the marker evidence *does* resolve it. Collapsing a
resolved case into an unresolved one throws away an answer the tool has.

## 15. A protein belongs to one marker

A protein that clears the score cutoff of several markers is counted only under
the one it scores highest against, and a protein already complete under one
marker is not also a fragment of another. Both run before the 85%-retention
step, which is where BUSCO runs them.

**Why:** the retention rule looks *within* a marker and asks whether a second
hit is a real copy or a distant paralog. It cannot see the other direction — a
hit that is a perfectly good protein which simply belongs to a different marker.
In `fungi_odb10`, three AAA-ATPases each clear the cutoff of the Pex1 marker as
well as of the AAA-ATPase marker, and score higher against the latter. Without
this rule Pex1 is reported duplicated when it is single.

That was the only marker of 758 on which GenomeX disagreed with BUSCO 5.8.3
scoring the same *S. cerevisiae* proteome. With the rule it is 758/758.

The rule was missing because the bacterial benchmark could not see it. Bacterial
core markers are largely ribosomal and rarely share a protein between families;
re-deriving all 72 genomes changes 1 marker call of 9,002 under
`bacteria_odb10` and 5 of 49,963 under `burkholderiales_odb10`, and the demo
set's 496/496 agreement is byte-identical either way. A defect invisible in one
lineage and decisive in another is an argument for validating against more than
one, not for assuming the first was representative.

**Rejected:** applying it after the retention step, which is where it would
naturally have gone in this code. BUSCO's order is the other way round, and
order changes the answer: a hit removed by the cross-marker rule must not first
be allowed to set a marker's best score, or the 85% floor is computed from a hit
that was never this marker's to begin with.

**Rejected:** treating the disagreement as a fungal quirk and special-casing
eukaryotic lineages. The rule is not about eukaryotes. It is about marker sets
containing paralogous families, which any sufficiently large set will.

## 16. A proteome is a different input, not a genome with holes

Protein FASTA input gets its own result type, its own JSON key (`proteomes`),
and its own report section. It is not a `GenomeResult` with the assembly fields
left empty.

**Why:** every downstream consumer of `genomes` — the report, the comparative
step, the benchmark harnesses — assumes an assembly behind each entry. An entry
that satisfies that assumption syntactically and violates it semantically is the
worst of the three options, because nothing fails. A proteome fed through
`Assembly.load` produces one "contig" per protein, a genome "size" that is a
residue count, and a GC content computed from glycine and cysteine: on the
*S. cerevisiae* proteome, **35.36% against that genome's real 38.15%**. In
range, plausible, and wrong.

The same reasoning drives what the report *says*. Contamination, ANI and islands
are listed with the reason each cannot be computed, rather than omitted. An
absent contamination section reads as "we looked and found nothing", which is a
stronger claim than any this input supports.

The `not_measured` block separates `impossible` from `unimplemented`, and the
distinction is load-bearing. Clustering proteomes into orthogroups needs no
nucleotides — MMseqs2 would do it today. What stops it is that the comparative
step attributes every unique gene to a contig and a coordinate. Calling that
"impossible" would dress a fact about this code as a fact about biology.

**Rejected:** auto-detecting the input type from the file and routing silently.
The two paths answer different questions and carry different guarantees, so
which one ran must be visible in the command the user typed, not inferred from
the file they passed. Alphabet detection is still used, but only to *refuse* —
never to switch paths.

**Rejected:** a shared `usable_for_comparative_analysis` boolean across both
result types. Boolean-shaped code is how the unassessed case came out green on
the assembly path for as long as it did: `undetermined` is not `likely`, so
`verdict != "likely"` passed it. Three states, or the absence of a check reads
as the success of one.

## 17. Refuse the wrong file type; never guess which one it is

Both input paths test the residue alphabet and refuse the other's file. The
discriminator is the share of residues that have no nucleotide meaning under any
IUPAC reading — E, F, I, L, P, Q, X, Z, J, O.

**Why:** the obvious test, "what fraction of this is ACGTU", fails on ordinary
input. A scaffolded assembly padded with N scores 0.80, so any threshold a
proteome clears would reject a real genome. Measured on the protein-only
alphabet instead: *S. cerevisiae* proteome 0.3535, its genome 0.0000, a
bacterial assembly 0.0000, a 20%-N scaffold 0.0000. The populations are 0.35
apart, so the exact cutoff is not load-bearing.

Below 200 residues the test abstains rather than answers. Sixteen of the
twenty-six letters are legal nucleotide ambiguity codes — M, K, V and S are
methionine, lysine, valine and serine *and* A/C, G/T, A/C/G, G/C — so a short
peptide can be entirely nucleotide-legal by chance. `MKVAAGCGTS` is. This is the
same abstention as `min_contigs_for_zscores`: below the sample size where a
statistic discriminates, do not claim it does.

**Rejected:** refusing on the file extension. `.fasta` is used for both, and an
extension is a claim by the person who named the file rather than a property of
its contents. Extension is a hint; alphabet is the fact.
