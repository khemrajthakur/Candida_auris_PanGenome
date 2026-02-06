#!/bin/bash
set -e

# Pilon polishing

TRIMMED="../../data/trimmed_reads"
ASSEMBLY="../../data/assemblies"
PILON_OUT="../../data/pilon"
THREADS=4

mkdir -p ${PILON_OUT}

for DIR in ${ASSEMBLY}/*
do
    SAMPLE=$(basename ${DIR})
    CONTIGS=${DIR}/contigs.fasta

    echo "Polishing ${SAMPLE}"

    bwa index ${CONTIGS}

    bwa mem -t ${THREADS} \
        ${CONTIGS} \
        ${TRIMMED}/${SAMPLE}_R1_paired.fastq.gz \
        ${TRIMMED}/${SAMPLE}_R2_paired.fastq.gz | \
        samtools sort -o ${PILON_OUT}/${SAMPLE}.bam

    samtools index ${PILON_OUT}/${SAMPLE}.bam

    pilon \
        --genome ${CONTIGS} \
        --frags ${PILON_OUT}/${SAMPLE}.bam \
        --output ${SAMPLE}_pilon \
        --outdir ${PILON_OUT}
done

echo "Pilon completed"

