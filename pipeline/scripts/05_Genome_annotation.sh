#!/bin/bash
set -e

# Funannotate genome annotation for Candida auris
# Input: Pilon-polished FASTA files
# Output: Annotated genomes

WORK_DIR="../../data/pilon"
OUTPUT_BASE="../../data/annotations"
THREADS=8

# Funannotate database must be set before running
# export FUNANNOTATE_DB=/path/to/funannotate_db

mkdir -p "${OUTPUT_BASE}"
LOG_FILE="${OUTPUT_BASE}/funannotate.log"

echo "Starting Funannotate annotation"
echo "==============================="

cd "${WORK_DIR}"

for fasta in *_pilon.fasta; do

    base_name=$(basename "${fasta}" _pilon.fasta)
    output_dir="${OUTPUT_BASE}/${base_name}_cauris_annotation"

    echo "Processing ${base_name}"

    # Fix contig names
    fixed_fasta="${base_name}_fixed.fasta"
    awk '/^>/ {print ">NODE_" ++i; next} {print}' "${fasta}" > "${fixed_fasta}"

    # Clean assembly
    cleaned_fasta="${base_name}_cleaned.fasta"
    funannotate clean -i "${fixed_fasta}" -o "${cleaned_fasta}" --minlen 500

    # Mask repeats
    masked_fasta="${base_name}_masked.fasta"
    funannotate mask -i "${cleaned_fasta}" -o "${masked_fasta}" --cpus ${THREADS}

    # Gene prediction
    funannotate predict \
        -i "${masked_fasta}" \
        -o "${output_dir}" \
        --species "candida_albicans" \
        --protein_evidence "${FUNANNOTATE_DB}/trained_species/candida_auris/cauris_proteins.fa" \
        --optimize_augustus \
        --cpus ${THREADS}

    # Functional annotation
    funannotate annotate -i "${output_dir}" --cpus ${THREADS}

    # Clean temporary files
    rm -f "${fixed_fasta}" "${cleaned_fasta}" "${masked_fasta}"

    echo "Done ${base_name}"
done

echo "==============================="
echo "Funannotate completed"

