# A pipeline for bacterial genome assembly
# Andrew Wilson
# Haney Lab, UBC

SAMPLES = # Enter samples here ex. ["CH261", "CH253", "DC105", "PB126", "GXM4"]

temps = ["data/*/*.trimmed.*.fastq", "data/*/*.pear.*.fastq",
"data/*/*all.singles.fastq"]

rule all:
    input: expand("data/{samples}/{samples}.prokka/{samples}.fna", samples=SAMPLES)

rule clean:
    shell:
        "rm -rf data/*/assembly/ \
        rm -f {temps} \
        rm -rf data/*/bowtie2"

# Read trimming with trimmomatic
rule quality_trim:
    input:
        rf = "data/reads/{samples}_1.fastq",
        rr = "data/reads/{samples}_2.fastq"
    output:
        fp="data/{samples}/{samples}.trimmed.1.p.fastq",
        fu="data/{samples}/{samples}.trimmed.1.u.fastq",
        rp="data/{samples}/{samples}.trimmed.2.p.fastq",
        ru="data/{samples}/{samples}.trimmed.2.u.fastq",
    shell:
        "trimmomatic PE {input.rf} {input.rr} \
        {output.fp} {output.fu} {output.rp} {output.ru} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1:keepBothReads \
        SLIDINGWINDOW:4:15 MINLEN:50"

# Read merging in PEAR
rule pear:
    input:
        fw = "data/{samples}/{samples}.trimmed.1.p.fastq",
        rev = "data/{samples}/{samples}.trimmed.1.p.fastq"
    params:
        "data/{samples}/{samples}.pear"
    output:
        uf = "data/{samples}/{samples}.pear.unassembled.forward.fastq",
        ur = "data/{samples}/{samples}.pear.unassembled.reverse.fastq",
        a = "data/{samples}/{samples}.pear.assembled.fastq"
    shell:
        "pear -j 4 -f {input.fw} -r {input.rev} -o {params}"

# Putting all single-end reads into one file
rule file_prep:
    input:
        fw = "data/{samples}/{samples}.trimmed.1.u.fastq",
        rev = "data/{samples}/{samples}.trimmed.2.u.fastq",
        pear = "data/{samples}/{samples}.pear.assembled.fastq"
    output:
        "data/{samples}/{samples}.all.singles.fastq"
    shell:
        "cat {input.fw} {input.rev} > {output}"

# de novo Genome assembly with SPAdes
rule assemble:
    input:
        fw = "data/{samples}/{samples}.pear.unassembled.forward.fastq",
        rev = "data/{samples}/{samples}.pear.unassembled.reverse.fastq",
        singles = "data/{samples}/{samples}.all.singles.fastq"
    params:
        "data/{samples}/assembly"
    output:
        contigs = "data/{samples}/assembly/contigs.fasta"
    shell:
        "spades.py -m 16 -s {input.singles} -1 {input.fw} -2 {input.rev} \
        --careful --cov-cutoff auto -o {params}"

# Detection of PhiX contigs with nhmmer
rule PhiX_hits:
    input:
        "data/{samples}/assembly/contigs.fasta"
    output:
        "data/{samples}/assembly/{samples}.phiX.hits"
    shell:
        "nhmmer --tblout {output} tools/PhiX.fna {input}"

# Assembly statistics and contig filtering using a modified GetGenomeStats.py
rule assembly_stats:
    input:
        contigs = "data/{samples}/assembly/contigs.fasta",
        phix = "data/{samples}/assembly/{samples}.phiX.hits"
    output:
        "data/{samples}/assembly_stats.txt",
        "data/{samples}/{samples}.fasta"
    script:
        "tools/assembly_stats.py"

# Annotation of the filtered contigs
rule annotation:
    input:
        "data/{samples}/{samples}.fasta"
    params:
        outdir="data/{samples}/{samples}.prokka/",
        sample="{samples}"
    output:
        "data/{samples}/{samples}.prokka/{samples}.fna"
    shell:
        "prokka --force --outdir {params.outdir} --genus Pseudomonas --strain \
        {params.sample} --prefix {params.sample} --locustag {params.sample} \
        {input}"
