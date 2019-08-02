# A pipeline for bacterial genome assembly
# Andrew Wilson
# Haney Lab, UBC

sample = "{sample}"
d = f"data/{sample}"
temps = ["data/*/*.trimmed.*.fastq", "data/*/*.pear.*.fastq",
"data/*/*all.singles.fastq"]

rule all:
    input:
        f"{d}/assembly_stats.txt"
rule clean:
    shell:
        "rm -rf data/*/assembly/ \
        rm -f {temps} \
        rm -rf data/*/bowtie2"

# Read trimming with trimmomatic
rule quality_trim:
    input:
        rf = "data/reads/{sample}_1.fastq",
        rr = "data/reads/{sample}_2.fastq"
    output:
        fp=f"{d}/{sample}.trimmed.1.p.fastq",
        fu=f"{d}/{sample}.trimmed.1.u.fastq",
        rp=f"{d}/{sample}.trimmed.2.p.fastq",
        ru=f"{d}/{sample}.trimmed.2.u.fastq",
    log:
        f"{d}/logs/trim_log.txt"
    shell:
        "trimmomatic PE {input.rf} {input.rr} \
        {output.fp} {output.fu} {output.rp} {output.ru} \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:1:keepBothReads \
        SLIDINGWINDOW:4:15 MINLEN:50"

# Read merging in PEAR
rule pear:
    input:
        fw = f"{d}/{sample}.trimmed.1.p.fastq",
        rev = f"{d}/{sample}.trimmed.1.p.fastq"
    params:
        f"{d}/{sample}.pear"
    output:
        uf = f"{d}/{sample}.pear.unassembled.forward.fastq",
        ur = f"{d}/{sample}.pear.unassembled.reverse.fastq",
        a = f"{d}/{sample}.pear.assembled.fastq"
    log:
        f"{d}/logs/pear_log.txt"
    shell:
        "pear -j 4 -f {input.fw} -r {input.rev} -o {params}"

# Putting all single-end reads into one file
rule file_prep:
    input:
        fw = f"{d}/{sample}.trimmed.1.u.fastq",
        rev = f"{d}/{sample}.trimmed.2.u.fastq",
        pear = f"{d}/{sample}.pear.assembled.fastq"
    output:
        f"{d}/{sample}.all.singles.fastq"
    shell:
        "cat {input.fw} {input.rev} > {output}"

# de novo Genome assembly with SPAdes
rule assemble:
    input:
        fw = f"{d}/{sample}.pear.unassembled.forward.fastq",
        rev = f"{d}/{sample}.pear.unassembled.reverse.fastq",
        singles = f"{d}/{sample}.all.singles.fastq"
    params:
        f"{d}/assembly"
    output:
        contigs = f"{d}/assembly/contigs.fasta"
    log:
        f"{d}/logs/assembly_log.txt"
    shell:
        "spades.py -m 16 -s {input.singles} -1 {input.fw} -2 {input.rev} \
        --careful --cov-cutoff auto -o {params}"

# Detection of PhiX contigs with nhmmer
rule PhiX_hits:
    input:
        f"{d}/assembly/contigs.fasta"
    output:
        f"{d}/assembly/{sample}.phiX.hits"
    shell:
        "nhmmer --tblout {output} tools/PhiX.fna {input}"

# Generation of assembly statistics using a modified GetGenomeStats.py
rule assembly_stats:
    input:
        contigs = f"{d}/assembly/contigs.fasta",
        phix = f"{d}/assembly/{sample}.phiX.hits"
    output:
        f"{d}/assembly_stats.txt",
        f"{d}/{sample}.fasta"
    script:
        "tools/assembly_stats.py"

# Annotation of the filtered contigs
rule annotation:
    input:
        f"{d}/{sample}.fasta"
    params:
        outdir=f"{d}/{sample}.prokka/",
        sample=f"{sample}"
    output:
        f"{d}/{sample}.prokka/{sample}.fna"
    shell:
        "prokka --force --outdir {params.outdir} --genus Pseudomonas --strain \
        {params.sample} --prefix {params.sample} --locustag {params.sample} \
        {input}"
