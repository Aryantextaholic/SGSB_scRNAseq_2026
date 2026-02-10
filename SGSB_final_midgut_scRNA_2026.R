# =========================================================
# SGSB midgut scRNA-seq workflow (10x -> Seurat -> Harmony -> Monocle3)
# Author: Dr. Surjeet Kumar Arya
# Repo: https://github.com/Aryantextaholic/SGSB_scRNAseq_2026
# =========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(harmony)
  library(DoubletFinder)
  library(monocle3)
  library(SeuratWrappers)
})

# -----------------------
# 0) User settings
# -----------------------
BASE_DIR   <- "/Users/surjeetarya/Desktop/All_folders/SGSB/SGSB_re_analysis/sgsb_MT/"
REP1_DIR   <- file.path(BASE_DIR, "filtered_feature_bc_matrix_rep1_MT")
REP2_DIR   <- file.path(BASE_DIR, "filtered_feature_bc_matrix_rep2_merged_MT")

MT_GENE_CSV <- file.path(BASE_DIR, "sgsb_MT.csv")   # CSV with column GeneName
OUTDIR      <- "results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# QC thresholds (change as needed)
QC_REP1 <- list(feature_min=300, feature_max=1500, count_min=300,  count_max=6000,  percent_mt_max=10)
QC_REP2 <- list(feature_min=800, feature_max=1500, count_min=1000, count_max=5000,  percent_mt_max=10)

# Integration / clustering settings
DIMS_USE <- 1:20
RESOLUTION <- 0.5
HARMONY_THETA <- 1

# DoubletFinder settings
DF_PCS <- 1:20
DF_pN  <- 0.25
DF_rate <- 0.10  # expected doublet rate (10%)

# -----------------------
# 1) Input checks
# -----------------------
stopifnot(dir.exists(REP1_DIR), dir.exists(REP2_DIR), file.exists(MT_GENE_CSV))

mt_genes <- read.csv(MT_GENE_CSV)$GeneName
if (length(mt_genes) < 10) warning("MT gene list seems small; check MT_GENE_CSV column GeneName.")

# -----------------------
# 2) Helper functions
# -----------------------
read_and_preprocess <- function(data_dir, project_name, mt_genes) {
  counts <- Read10X(data.dir = data_dir)
  seu <- CreateSeuratObject(counts = counts, project = project_name)
  seu <- PercentageFeatureSet(seu, features = mt_genes, col.name = "percent.mt")
  return(seu)
}

quality_control <- function(seu, qc) {
  subset(
    seu,
    subset =
      nFeature_RNA >= qc$feature_min &
      nFeature_RNA <= qc$feature_max &
      nCount_RNA   >= qc$count_min &
      nCount_RNA   <= qc$count_max &
      percent.mt   <= qc$percent_mt_max
  )
}

run_df_single <- function(seu, pcs=1:20, pN=0.25, expected_rate=0.10) {
  # standard log-normalize + HVG + PCA + clustering to get seurat_clusters
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu)
  seu <- ScaleData(seu, vars.to.regress = c("nCount_RNA","percent.mt"))
  seu <- RunPCA(seu, verbose = FALSE)
  seu <- FindNeighbors(seu, dims = pcs)
  seu <- FindClusters(seu, resolution = 0.5)
  seu <- RunUMAP(seu, dims = pcs)

  sweep.res.list <- paramSweep_v3(seu, PCs = pcs, sct = FALSE)
  sweep.stats    <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn          <- find.pK(sweep.stats)

  best_pK <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::pull(pK) %>%
    as.character() %>%
    as.numeric()

  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)
  nExp_poi <- round(expected_rate * nrow(seu@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

  seu <- doubletFinder_v3(
    seu,
    PCs = pcs,
    pN = pN,
    pK = best_pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  # find the DF classification column that was created
  df_class_cols <- grep("^DF.classifications", colnames(seu@meta.data), value = TRUE)
  if (length(df_class_cols) == 0) stop("DoubletFinder did not create DF.classifications column.")
  df_col <- df_class_cols[length(df_class_cols)]  # take the newest
  seu$DF_class <- seu@meta.data[[df_col]]

  # keep singlets
  seu_clean <- subset(seu, subset = DF_class == "Singlet")
  return(list(seu_raw = seu, seu_clean = seu_clean, best_pK = best_pK, df_col = df_col))
}

# -----------------------
# 3) Read + QC each replicate
# -----------------------
rep1 <- read_and_preprocess(REP1_DIR, "SGSB_1", mt_genes)
rep2 <- read_and_preprocess(REP2_DIR, "SGSB_2", mt_genes)

rep1_qc <- quality_control(rep1, QC_REP1)
rep2_qc <- quality_control(rep2, QC_REP2)

saveRDS(rep1_qc, file.path(OUTDIR, "SGSB_rep1_QC.rds"))
saveRDS(rep2_qc, file.path(OUTDIR, "SGSB_rep2_QC.rds"))

# -----------------------
# 4) DoubletFinder per replicate (recommended)
# -----------------------
df1 <- run_df_single(rep1_qc, pcs=DF_PCS, pN=DF_pN, expected_rate=DF_rate)
df2 <- run_df_single(rep2_qc, pcs=DF_PCS, pN=DF_pN, expected_rate=DF_rate)

saveRDS(df1$seu_raw,   file.path(OUTDIR, "SGSB_rep1_after_DF_raw.rds"))
saveRDS(df1$seu_clean, file.path(OUTDIR, "SGSB_rep1_after_DF_singlets.rds"))
saveRDS(df2$seu_raw,   file.path(OUTDIR, "SGSB_rep2_after_DF_raw.rds"))
saveRDS(df2$seu_clean, file.path(OUTDIR, "SGSB_rep2_after_DF_singlets.rds"))

message("Rep1 best pK: ", df1$best_pK, " | DF column: ", df1$df_col)
message("Rep2 best pK: ", df2$best_pK, " | DF column: ", df2$df_col)

# -----------------------
# 5) Merge replicates + add metadata
# -----------------------
merged <- merge(df1$seu_clean, y = df2$seu_clean, add.cell.ids = c("Rep1","Rep2"), project = "SGSB_merged")

merged$sample <- rownames(merged@meta.data)
merged@meta.data <- separate(merged@meta.data, col = "sample", into = c("Sampletype","Barcode"), sep = "_", remove = FALSE)

saveRDS(merged, file.path(OUTDIR, "SGSB_merged_singlets.rds"))

# -----------------------
# 6) Integration (CCA) + optional Harmony
# -----------------------
obj.list <- SplitObject(merged, split.by = "Sampletype")
obj.list <- lapply(obj.list, NormalizeData)
obj.list <- lapply(obj.list, FindVariableFeatures)

features <- SelectIntegrationFeatures(object.list = obj.list)
anchors  <- FindIntegrationAnchors(object.list = obj.list, anchor.features = features)
integrated <- IntegrateData(anchors)

DefaultAssay(integrated) <- "integrated"
integrated <- ScaleData(integrated, verbose = FALSE)
integrated <- RunPCA(integrated, verbose = FALSE)

# Option A: cluster using integrated PCA
integrated <- FindNeighbors(integrated, dims = DIMS_USE)
integrated <- FindClusters(integrated, resolution = RESOLUTION)
integrated <- RunUMAP(integrated, dims = DIMS_USE)

# Option B: Harmony on PCA (batch correction)
integrated <- RunHarmony(
  object = integrated,
  group.by.vars = "Sampletype",
  reduction = "pca",
  dims.use = DIMS_USE,
  theta = HARMONY_THETA
)

integrated <- RunUMAP(integrated, reduction = "harmony", dims = DIMS_USE, reduction.name = "umap_harmony")
integrated <- FindNeighbors(integrated, reduction = "harmony", dims = DIMS_USE)
integrated <- FindClusters(integrated, resolution = RESOLUTION)

saveRDS(integrated, file.path(OUTDIR, "SGSB_integrated_harmony.rds"))

p_umap <- DimPlot(integrated, reduction = "umap_harmony", label = TRUE, repel = TRUE) + NoLegend()
ggsave(file.path(OUTDIR, "UMAP_harmony.png"), p_umap, width = 7, height = 6, dpi = 300)

# -----------------------
# 7) Markers
# -----------------------
DefaultAssay(integrated) <- "RNA"
integrated <- SCTransform(integrated, verbose = FALSE)
integrated <- PrepSCTFindMarkers(integrated)

# example markers per cluster
markers <- FindAllMarkers(integrated, only.pos = TRUE)
write.csv(markers, file.path(OUTDIR, "markers_FindAllMarkers.csv"), row.names = FALSE)

# -----------------------
# 8) Monocle3 pseudotime
# -----------------------
# Use harmony UMAP coordinates if you want:
# Embeddings(integrated, "umap_harmony")
cds <- as.cell_data_set(integrated)
reducedDims(cds)$UMAP <- Embeddings(integrated, "umap_harmony")
colData(cds)$ident <- as.character(Idents(integrated))

cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds, use_partition = FALSE)

# root example: use SC/EB if present
if ("SC/EB" %in% unique(colData(cds)$ident)) {
  root_cells <- rownames(subset(as.data.frame(colData(cds)), ident == "SC/EB"))
  cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)
} else {
  cds <- order_cells(cds, reduction_method = "UMAP")
}

p_pt <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE)
ggsave(file.path(OUTDIR, "Monocle3_pseudotime.png"), p_pt, width = 7, height = 6, dpi = 300)

deg_pseudo <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
write.csv(deg_pseudo, file.path(OUTDIR, "Monocle3_graph_test.csv"), row.names = FALSE)

saveRDS(cds, file.path(OUTDIR, "SGSB_monocle3_cds.rds"))

message("DONE. Checkpoints saved as .rds in: ", OUTDIR)
