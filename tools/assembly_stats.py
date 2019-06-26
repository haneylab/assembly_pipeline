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
GC_skew = {}
cov = {}
length = {}
seq_len = 0
sequence = {}

contigs = []
with open("data/contigs.fasta", "r") as c:
    line = c.readline()
    while line:
        if line.startswith('>'):
            contigs.append(line.replace('>', '').rstrip())
            line = c.readline()
        else:
            seq_len += len(line.rstrip())
            line = c.readline()

for c in contigs:
    info = c.split('_')
    contigname = f"{info[0]}_{info[1]}"
    length[contigname] = info[3]
    cov[contigname] = info[5]
    if length[contigname] > 1000:



print(seq_len)


def gc_skew(seq:str, k=5):
    skew = {}
    for i in range(0, len(seq) - k):
        kmer = seq[i:i+k]
        gua = kmer.count("G")
        cyt = kmer.count("C")
        skew[i] = (gua - cyt) / (gua + cyt)
    return skew
