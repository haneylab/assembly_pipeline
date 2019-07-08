# Andrew Wilson
# andrew.wilson@alumni.ubc.ca
# Haney Lab - UBC

import os
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import math
# import snakemake


# o = open(snakemake.output[0], "w")
GC_content = {}
cov = {}
length = {}
seq_len = 0
sequence = {}

contigs = {}
for seq in SeqIO.parse(open("data/CH261/assembly/contigs.fasta", "r"), 'fasta'):
    info = seq.id.split("_")
    contigname = f"{info[0]}_{info[1]}"
    



#def is_phix()

