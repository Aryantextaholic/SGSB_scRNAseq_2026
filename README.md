Pipeline:
10x → QC → DoubletFinder → Integration → Harmony → Clustering → Marker discovery → Monocle3 pseudotime

Species:
_Nezara viridula_ midgut single-cell RNA-seq

Author:
Surjeet Kumar Arya, PhD
University of Kentucky

# SGSB midgut scRNA-seq analysis (Seurat + Harmony + DoubletFinder + Monocle3)

This repository contains the R workflow used for single-cell RNA-seq analysis of the southern green stink bug (SGSB), *Nezara viridula*, midgut.


## Methods summary

Single-cell RNA-seq data were processed using Seurat v5.
Quality control thresholds were applied on gene counts,
UMI counts, and mitochondrial gene percentages.

Doublets were removed using DoubletFinder.
Batch correction was performed using Harmony.
Clusters were identified using shared nearest neighbor graphs.
Pseudotime trajectories were reconstructed using Monocle3.


## Contents
- `SGSB_final_midgut_scRNA_2026.R` — end-to-end analysis script (QC → normalization/integration → clustering/UMAP → markers → pseudotime)

## Requirements
- R (>= 4.2 recommended)
- Key packages: Seurat, harmony, DoubletFinder, monocle3, SeuratWrappers, tidyverse, ggplot2

## Data
Large input objects (Cell Ranger matrices / Seurat `.rds`) are not included in this repository.
Update file paths inside the script to point to your local inputs.

## Outputs
The script generates:
- UMAP / cluster plots
- marker tables per cluster
- pseudotime trajectory plots (Monocle3)

## Citation
If you use this workflow, please cite:
Arya SK et al. (manuscript in preparation / under review).


## Data Availability

Raw sequencing data are not included in this repository due to file size and institutional policies.

Data are available from:
- NCBI SRA: [accession will be added after publication]
- Institutional storage: available upon reasonable request
- Processed objects (.rds) available on request for reproducibility
