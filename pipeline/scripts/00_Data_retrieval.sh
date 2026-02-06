#!/bin/bash
set -e

# Download sequencing reads from NCBI SRA
# Requires: sra-tools

RAW_READS="../../data/raw_reads"
ACCESSIONS="../../pipeline/config/sra_accessions.txt"

mkdir -p ${RAW_READS}

echo "Starting read download..."

while read SRR
do
    echo "Downloading ${SRR}"
    prefetch ${SRR}
    fastq-dump --split-files --gzip ${SRR} -O ${RAW_READS}
done < ${ACCESSIONS}

echo "Read download completed"

