# Bacterial Genome Assembly pipeline
This pipeline can be used to automatically assemble bacterial genomes sequenced on an Illumina sequencer.

## Installation
Ensure that you have a recent version of conda installed on your computer. I recommend the Python 3.7 version of Miniconda which you can install [here](https://docs.conda.io/en/latest/miniconda.html). Once you have conda, navigate to this repo in your command line and make a new conda environment with the command `conda env create -f dependencies.yml`. This will install all of the software needed to run the pipeline.

## Usage
1. Copy your read files into the `data/reads/` directory and rename them so that they fit the form **strain**_1.fastq and **strain**_2.fastq for forward and reverse reads, respectively.

2. Run the command `chmod u-w data/reads/*.fastq` so that you can't accidentally delete or rewrite your read files.

3. Open the `Snakefile` in your favourite text editor. Change the `SAMPLES` list on line 5 so that it contains your strains. Strain names need to be in quotations and separated by commas, and the entire list needs to be enclosed in square brackets.

4. Activate the conda enviroment by typing `conda activate assembly_pipeline`. In your terminal, type `snakemake` and wait for the pipeline to run.

5. To remove intermediate files, type `snakemake clean`. You can uninstall all of the software by deactivating the conda environment and deleting it with the commands:

```
conda deactivate assembly_pipeline
conda env remove -n assembly_pipeline
```

## Output
This pipeline produces genomes in the `data/` directory. Each genome is in a subdirectory named after the strain that it contains. There should be two files, `assembly_stats.txt` and a fasta file, as well as a directory which contains annotation information. 
