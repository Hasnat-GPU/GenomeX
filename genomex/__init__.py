"""GenomeX -- a unified genome QC and comparative-genomics pipeline.

Modules map one-to-one onto the whiteboard workflow:

    fasta          contigs in, assembly statistics out
    genes          ORF/CDS prediction (Prodigal)
    markers        single-copy ortholog recovery (BUSCO odb10 HMMs + HMMER)
    contamination  is this one organism, or more than one?
    compare        why do two isolates from one environment carry different genes?
    report         JSON + standalone HTML
"""

__version__ = "0.2.0"
