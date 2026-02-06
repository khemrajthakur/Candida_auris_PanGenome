# Candida auris pan-genome analysis pipeline

This repository contains a complete end-to-end pipeline for
pan-genome analysis of *Candida auris*, including genome assembly,
annotation, population structure, phylogenetics, and variant analysis.

## Repository structure

- `pipeline/`  
  Scripts for genome assembly, polishing, annotation, pan-genome
  inference, and variant calling.

- `downstream_analysis/`  
  R scripts for visualization and statistical analysis, including
  pan-genome composition, presence–absence, phylogeny, DAPC, and
  variant burden heatmaps.

- `data/`  
  Input data directory (not tracked on GitHub).

- `results/`  
  Output directory generated during pipeline execution (not tracked).

## Requirements
- Linux (tested on Ubuntu/WSL)
- Conda
- Docker
- Nextflow
- R (≥4.2)

## How to run
See individual README files inside `pipeline/` and
`downstream_analysis/` for step-by-step instructions.

## Citation
If you use this pipeline, please cite:
[Your paper / preprint / thesis here]

