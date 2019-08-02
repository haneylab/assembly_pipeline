# Andrew Wilson
# andrew.wilson@alumni.ubc.ca
# Haney Lab - UBC

import os
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd


o = open(snakemake.output[0], "w")
cov = {}
lengths = {}
contigs = {}



def is_phix(file) -> list:
    phix_contigs = []
    for line in open(file, "r"):
        if line.startswith("#"):\
            continue
        else:
            hit = line.rstrip().split()
            if float(hit[13]) > 100:
                contig = hit[0].split("_")
                phix_contigs.append(f"{contig[0]}_{contig[1]}")
    return phix_contigs


for seq in SeqIO.parse(open(snakemake.input[0], "r"), 'fasta'):
    info = seq.id.split("_")
    contigname = f"{info[0]}_{info[1]}"
    contigs[contigname] = seq.seq
    lengths[contigname] = int(info[3])
    cov[contigname] = float(info[5])

phix = is_phix(snakemake.input[1])

filt_contigs = {}
for contig in contigs.keys():
    if lengths[contig] <= 1000:
        pass
    else:
        if contig not in phix:
            filt_contigs[contig] = contigs[contig]


o.write(f"Assembly statistics\n")
o.write(f"Total sequence length:\t{str(sum(lengths.values()))}\n")
o.write(f"Total number of contigs:\t{str(len(lengths.values()))}\n")
o.write(f"Number of contigs > 1000 bp:\t{str(len(filt_contigs.keys()))}\n")

n50 = 0
for i in sorted(lengths.values(),reverse=True):
        n50 += i
        if n50 > sum(lengths.values())/2:
            o.write(f"N50\t{str(i)}\n")
            break
o.close()

fasta = ""
for c in filt_contigs.keys():
    fasta += f">{c}\n{filt_contigs[c]}\n"

ff = open(snakemake.output[1], "w")
ff.write(fasta)
ff.close()


