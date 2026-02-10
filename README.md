# SGSB midgut scRNA-seq analysis (Seurat + Harmony + DoubletFinder + Monocle3)

This repository contains the R workflow used for single-cell RNA-seq analysis of the southern green stink bug (SGSB), *Nezara viridula*, midgut.

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
