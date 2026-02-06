#!/bin/bash
set -e

# Quality control: FastQC + Trimmomatic

RAW_READS="../../data/raw_reads"
TRIMMED="../../data/trimmed_reads"
FASTQC_RAW="../../results/fastqc_raw"
FASTQC_TRIMMED="../../results/fastqc_trimmed"

THREADS=4
ADAPTERS="adapters.fa"

mkdir -p ${TRIMMED} ${FASTQC_RAW} ${FASTQC_TRIMMED}

echo "=============================="
echo "Running FastQC on raw reads"
echo "=============================="
fastqc ${RAW_READS}/*.fastq* -o ${FASTQC_RAW} -t ${THREADS}

echo "=============================="
echo "Running Trimmomatic"
echo "=============================="

for R1 in ${RAW_READS}/*_R1.fastq.gz
do
    SAMPLE=$(basename ${R1} _R1.fastq.gz)
    R2=${RAW_READS}/${SAMPLE}_R2.fastq.gz

    echo "Processing ${SAMPLE}"

    trimmomatic PE \
        -threads ${THREADS} \
        ${R1} ${R2} \
        ${TRIMMED}/${SAMPLE}_R1_paired.fastq.gz \
        ${TRIMMED}/${SAMPLE}_R1_unpaired.fastq.gz \
        ${TRIMMED}/${SAMPLE}_R2_paired.fastq.gz \
        ${TRIMMED}/${SAMPLE}_R2_unpaired.fastq.gz \
        ILLUMINACLIP:${ADAPTERS}:2:30:10 \
        LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36
done

echo "=============================="
echo "Running FastQC on trimmed reads"
echo "=============================="
fastqc ${TRIMMED}/*_paired.fastq.gz -o ${FASTQC_TRIMMED} -t ${THREADS}

echo "Quality control completed"

