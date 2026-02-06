#!/bin/bash
set -e

# BUSCO analysis on Pilon-polished assemblies

PILON_DIR="../../data/pilon"
OUTDIR="../../results/busco"
LINEAGE="saccharomycetes_odb10"
THREADS=4

mkdir -p ${OUTDIR}

# Create summary file
SUMMARY=${OUTDIR}/busco_summary.csv
echo "Sample,Complete,Single_Copy,Duplicated,Fragmented,Missing" > ${SUMMARY}

echo "Running BUSCO on Pilon assemblies..."

for FASTA in ${PILON_DIR}/*_pilon.fasta
do
    SAMPLE=$(basename ${FASTA} _pilon.fasta)

    echo "BUSCO: ${SAMPLE}"

    busco \
        -i ${FASTA} \
        -o ${SAMPLE}_busco \
        -l ${LINEAGE} \
        -m genome \
        -c ${THREADS} \
        --out_path ${OUTDIR} \
        --quiet

    SUMMARY_FILE=${OUTDIR}/${SAMPLE}_busco/short_summary.specific.${LINEAGE}.${SAMPLE}_busco.txt

    COMPLETE=$(grep "Complete BUSCOs" ${SUMMARY_FILE} | awk '{print $1}')
    SINGLE=$(grep "Complete and single-copy BUSCOs" ${SUMMARY_FILE} | awk '{print $1}')
    DUP=$(grep "Complete and duplicated BUSCOs" ${SUMMARY_FILE} | awk '{print $1}')
    FRAG=$(grep "Fragmented BUSCOs" ${SUMMARY_FILE} | awk '{print $1}')
    MISS=$(grep "Missing BUSCOs" ${SUMMARY_FILE} | awk '{print $1}')

    echo "${SAMPLE},${COMPLETE},${SINGLE},${DUP},${FRAG},${MISS}" >> ${SUMMARY}
done

echo "BUSCO finished"
echo "Summary saved to: ${SUMMARY}"

