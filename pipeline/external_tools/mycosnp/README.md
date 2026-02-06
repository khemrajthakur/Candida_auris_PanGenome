# MycoSNP Variant Calling

Variant calling for Candida auris was performed using the MycoSNP
Nextflow pipeline developed by the CDC.

MycoSNP is executed as an external pipeline using Docker and Nextflow
and is not integrated directly into the main pan-genome workflow.

## Requirements
- Docker Desktop (running in background)
- WSL2 with Ubuntu
- Nextflow

## Pipeline repository
https://github.com/CDCgov/mycosnp-nf

## Basic usage

1. Create a working directory
mkdir -p mycosnp_pipeline
cd mycosnp_pipeline

2. Place FASTQ files in data/
mkdir data

3. Prepare samplesheet.csv

sample,fastq_1,fastq_2
SAMPLE1,/path/to/sample_1.fastq,/path/to/sample_2.fastq

4. Run MycoSNP

nextflow run CDCgov/mycosnp-nf \
  -profile docker \
  --input samplesheet.csv \
  --fasta candida_ref.fasta \
  --outdir results \
  --max_memory '7GB'

## Output
Variant calling results will be generated in the results/ directory.

