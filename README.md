# Nextflow RAD-Seq pipeline

## Overview

This is a RAD-Seq pipeline that generally follows

Data -> FastQC -> Trimmomatic -> FastQC -> process_radtags \
  -> ustacks -> cstacks -> sstacks -> tsv2bam -> gstacks \
  -> populations -> multiqc

## Running

1. Clone the repository

```bash
git clone --recursive https://github.com/bschiffthaler/nf-radseq.git
```

2. Place your raw data into `data/fastq`. The data can be raw fastq, or gzip compressed fastq (recommended).

3. Edit `data/metadata.csv`.
  3.1 Do not change the header row, it is used internally
  3.2 File paths should be either absolute or relative to the project root (this git repo)

4. Adapt `nextflow.config` if necessary

5. Start the pipeline: `nextflow run main.nf`

6. Generate template for materials and methods: `nextflow run materials.nf`

## Example metadata.csv

* *RF*: (string) Path to first PE read
* *RS*: (string) Path to second PE read
* *Id*: (string) A sample name. Avoid spaces and special characters. This should be unique
* *NumId*: (integer) A numeric sample id. This should be unique 
* *Population*: (string) the population this sample belongs to. If a sample belongs to multiple populations, enter all delimited by a `|`, e.g. `PopA|PopB|PopC`

```
RF,RS,Id,NumId,Population
./data/fastq/my_first_sample_1.fq.gz,./data/fastq/my_first_sample_2.fq.gz,first,1,default
./data/fastq/my_second_sample_1.fq.gz,./data/fastq/my_second_sample_2.fq.gz,second,2,default
```
