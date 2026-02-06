#!/bin/bash
set -e

# OrthoFinder pan-genome analysis
# Input: protein FASTA files from Funannotate
# Output: OrthoFinder results

ANNOTATION_DIR="../../data/annotations"
ORTHO_INPUT="../../data/orthofinder/proteins"
ORTHO_OUT="../../data/orthofinder/results"
THREADS=4

mkdir -p "${ORTHO_INPUT}"
mkdir -p "${ORTHO_OUT}"

echo "Collecting protein FASTA files..."

for dir in ${ANNOTATION_DIR}/*_cauris_annotation
do
    sample=$(basename "${dir}" _cauris_annotation)
    protein_fasta="${dir}/predict_results/"*.proteins.fa

    if ls ${protein_fasta} 1> /dev/null 2>&1; then
        cp ${protein_fasta} "${ORTHO_INPUT}/${sample}.faa"
        echo "Added ${sample}"
    fi
done

echo "Running OrthoFinder..."
orthofinder -f "${ORTHO_INPUT}" -t ${THREADS} -o "${ORTHO_OUT}"

echo "OrthoFinder completed"

