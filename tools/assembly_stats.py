# Andrew Wilson
# andrew.wilson@alumni.ubc.ca
# Haney Lab - UBC

import os
from Bio import SeqIO
from Bio import SeqUtils
import pandas as pd
import math
import snakemake

with open(snakemake.output[0], "w") as o:
    GC_content = {}
    GC_skew = {}
    cov = {}
    length = {}
    contigs = []
    with open(snakemake.input[1], "r") as gff:
        line = gff.readline()
        while line:
            if line.startswith('##'):
                continue
            elif line.startswith(">"):
                break


