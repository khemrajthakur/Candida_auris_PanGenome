# Pan-genome analysis and visualization

This folder contains R scripts for visualizing pan-genome composition and
clade-wise gene family dynamics in *Candida auris* using OrthoFinder output.

## Input
- `Orthogroups.GeneCount.tsv` – gene family counts per genome (OrthoFinder)
- `Orthogroups.tsv` – gene presence/absence across genomes (OrthoFinder)

(Input files are not included in the GitHub repository.)

## Scripts
- `pangenome_pie_chart.R`  
  Generates a donut / pie chart showing the distribution of core,
  accessory, and singleton gene families.

- `presence_absence_plot.R`  
  Identifies clade-specific gene family gains and losses and generates
  presence–absence plots across clades (CL1–CL6).

## Output
- Publication-quality pan-genome figures (PNG)
- Summary tables of clade-specific gains and losses (TSV)

## Usage
Run the scripts from this directory:

