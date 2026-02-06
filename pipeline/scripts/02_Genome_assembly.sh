#!/bin/bash
set -e

# Assembly and assembly quality control
# SPAdes + QUAST
# Author: Milani Sharma

TRIMMED="../../data/trimmed_reads"
ASSEMBLY="../../data/assemblies"
QUAST_OUT="../../results/quast"

THREADS=4

mkdir -p ${ASSEMBLY} ${QUAST_OUT}

echo "=============================="
echo "Running SPAdes assemblies"
echo "=============================="

for R1 in ${TRIMMED}/*_R1_paired.fastq.gz
do
    SAMPLE=$(basename ${R1} _R1_paired.fastq.gz)
    R2=${TRIMMED}/${SAMPLE}_R2_paired.fastq.gz

    echo "Assembling ${SAMPLE}"

    spades.py \
        -1 ${R1} \
        -2 ${R2} \
        -o ${ASSEMBLY}/${SAMPLE} \
        -t ${THREADS}
done

echo "=============================="
echo "Running QUAST on assemblies"
echo "=============================="

quast.py ${ASSEMBLY}/*/contigs.fasta \
    -o ${QUAST_OUT} \
    -t ${THREADS}

echo "Assembly + QUAST completed"
