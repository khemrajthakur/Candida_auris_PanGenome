# Genome analysis pipeline – Candida auris pan-genome

This directory contains the stepwise Linux-based pipeline used for
genome assembly, polishing, annotation, pan-genome analysis, and
variant calling of *Candida auris* isolates.

All scripts are designed to be run sequentially and were tested on
Ubuntu (WSL) using Conda and Docker.

---

## Pipeline overview

The pipeline follows this order:

1. Data retrieval  
2. Quality control and read trimming  
3. De novo genome assembly  
4. Assembly polishing and completeness assessment  
5. Genome annotation  
6. Pan-genome analysis  
7. Variant calling and phylogeny

---

## Script order and description

### `00_Data_retrieval.sh`
- Downloads raw sequencing data from NCBI (SRA)
- Converts SRA files to paired-end FASTQ format

### `01_Quality_control.sh`
- Performs quality assessment using FastQC
- Trims adapters and low-quality bases using Trimmomatic

### `02_Genome_assembly.sh`
- De novo genome assembly using SPAdes
- Generates contig-level assemblies for each isolate
- Assembly quality assessment using QUAST

### `03_Assembly_polishing.sh`
- Aligns reads back to assemblies using BWA
- Polishes assemblies using Pilon

### `04_Assembly_completeness.sh`
- Assesses genome completeness using BUSCO
- BUSCO analysis is performed on Pilon-polished assemblies only

### `05_Genome_annotation.sh`
- Structural and functional genome annotation using Funannotate
- Generates predicted gene models and protein sequences

### `06_Pangenome_analysis.sh`
- Pan-genome inference using OrthoFinder
- Generates orthogroups and gene presence–absence matrices

---

## External tools

The following tools are required or used by the pipeline:

- FastQC  
- Trimmomatic  
- SPAdes  
- QUAST  
- BWA  
- Pilon  
- BUSCO  
- Funannotate  
- OrthoFinder  
- MycoSNP (Nextflow + Docker)

