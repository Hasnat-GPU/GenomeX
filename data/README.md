# Dataset

`Bacteria.tar.gz` from the Box share (368 MB, 76 assemblies + metadata), extracted
to WSL ext4 for I/O speed:

    $GENOMEX_DATA/ (default: $HOME/genomex-work/data/Bacteria)     (772 MB extracted)

Only the metadata is copied into the repo (`mmc2_meta.tsv`, `acc.txt`); the FASTA
files stay on the Linux side because that is where the tools read them.

## What is in it

72 genomes with metadata, all beta-rhizobia -- legume nodule symbionts:

| Genus | Genomes |
|---|---|
| *Paraburkholderia* | 55 |
| *Cupriavidus* | 12 |
| *Burkholderia* | 2 |
| *Trinickia* | 1 |

Assembly level: 57 contig, 7 scaffold, 4 chromosome, 4 complete. Sizes 6.5–10.4 Mb,
3–1243 contigs, GC 62–67%.

This collection is a good fit for the third whiteboard question. Beta-rhizobia
carry their symbiosis genes (*nod*, *nif*, *fix*) on symbiosis islands and
megaplasmids that move by horizontal transfer, so two isolates from the same
nodule environment genuinely differ in gene content — and the difference has a
known mechanism to check the pipeline's explanations against.

## Demo selection

| Accession | Organism | Why |
|---|---|---|
| GCA_000300095.1 | *P. phenoliruptrix* BR3459a | finished, 3 replicons — exercises the few-contig and megaplasmid paths |
| GCA_040948545.1 | *P. phenoliruptrix* | draft, 77 contigs — same species, ANI ≈ 98% |
| GCF_000020045.1 | *P. phymatum* STM815 | finished, different species, same genus |
| GCA_000069785.1 | *C. taiwanensis* LMG 19424 | different genus, same environment — ANI below the species line |
