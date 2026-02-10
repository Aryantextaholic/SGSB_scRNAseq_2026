
library(Seurat)
library(patchwork)
library(magrittr)
library(gridExtra)
library(harmony)
library(dplyr)
library(tidyr)
library(DoubletFinder)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(tidyverse)
library(ggplot2)
library(cowplot)

SGSB_integ <- readRDS("/Users/surjeetarya/Desktop/All_files/All_single_data/Insects_data/SGSB/Paper_revised/Final_SGSB_Revissed_2024/SGSB_New/SGSB_Final_tutorial_integ_prep_2.rds")

dim1 <- DimPlot(SGSB_integ, reduction = "UMAP", group.by = "seurat_clusters")

dim1

## Run Harmony

SGSB_integ <- RunHarmony(SGSB_integ, features = VariableFeatures(integrated),group.by.vars = "Sampletype",  theta = 1)

SGSB_integ <- SCTransform(SGSB_integ, verbose = FALSE)

# Prepare the data for finding markers
SGSB_integ <- PrepSCTFindMarkers(SGSB_integ)


## Assuming you have loaded Seurat and ggplot2 libraries
library(Seurat)
library(ggplot2)

new.cluster.ids <- c(" IEC", "EC like 1", "EC like 2", "pEC1", "pEC2", "mEC1", "SC/EB", "mEC2", "GC", "ECM-EpC", "EE", "VM")
names(new.cluster.ids) <- levels(SGSB_integ)
SGSB_integ <- RenameIdents(SGSB_integ, new.cluster.ids)

# Create the DimPlot
dim_plot <- DimPlot(SGSB_integ, reduction = "UMAP", label = TRUE, pt.size = 0.5) + NoLegend()

# Modify the theme to make the cluster IDs bold
dim_plot_sgsb <- dim_plot + theme(strip.text = element_text(face = "bold"))

dim_plot_sgsb


######

dim_plot_sgsb <- DimPlot(
  SGSB_integ,
  reduction = "UMAP",
  label = TRUE,
  repel = TRUE,      # avoids label overlap
  label.size = 7,    # big labels
  pt.size = 0.4
) +
  NoLegend() +
  theme(
    text = element_text(face = "bold"),
    axis.title = element_text(face = "bold", size = 16),
    axis.text  = element_text(face = "bold", size = 16)
  )

dim_plot_sgsb




######

saveRDS(SGSB_integ, file = "/Users/surjeetarya/Desktop/All_single_data/SGSB/Code_for_analysis/SGSB_Final_tutorial_Idents.rds")


DimPlot(SGSB_integ, split.by = 'Sampletype')

DimPlot(SGSB_integ, group.by = c("orig.ident", "ident"))

DimPlot(SGSB_integ, group.by =  "ident")


# Prepare the data for finding markers
SGSB_integ <- PrepSCTFindMarkers(SGSB_integ)

Violin_before <- VlnPlot(SGSB_integ, features =c('nFeature_RNA','nCount_RNA','percent.mt'))
Violin_before

# -----------------
# Additional Plots
# -----------------

# Violin plot for specific gene
VlnPlot(SGSB_integ, features = c("NEZAVI-LOCUS3023"),  log = TRUE)



#####

VlnPlot(
  SGSB_integ,
  features = "NEZAVI-LOCUS7980",
  log = TRUE
) +
  theme(
    axis.title.x = element_text(face = "bold", size = 18),
    axis.title.y = element_text(face = "bold", size = 18),
    axis.text.x  = element_text(face = "bold", size = 18),
    axis.text.y  = element_text(face = "bold", size = 18)
  )



#####


genelist2 <- read.csv("ALKP_genes.csv")

genelist2 <- read.csv("ATP_binding_cassettes.csv")

genelist2 <- read.csv("ABCA_transporter.csv")

gene_ID <- genelist2$GeneNames

#head(genelist2)

VlnPlot(SGSB_integ, features = gene_ID,  log = TRUE)

# Feature plot for specific gene
FeaturePlot(SGSB_integ, features = gene_ID)



#NEZAVI-LOCUS14563: Aminopeptidase N
#NEZAVI-LOCUS3639: DE-Cadherins
#NEZAVI-LOCUS15376: Alkaline phosphatase
#NEZAVI-LOCUS12249: ATP-binding cassettes A





######

# Finding markers for cluster 1
clustermEC2.markers <- FindMarkers(SGSB_integ, ident.1 = "mEC2", ident.2 = NULL, only.pos = FALSE)

# Writing cluster markers to a CSV file
write.csv(clustermEC2.markers, "/Users/surjeetarya/Library/CloudStorage/OneDrive-UniversityofKentucky/Desktop/Desktope_2026/Manuscript/SGSB/Manuscript_SGSB/Markers/cluster_mEC2.csv")




######





# Violin plot for specific gene
VlnPlot(SGSB_integ, features = c("NEZAVI-LOCUS280","NEZAVI-LOCUS3642","NEZAVI-LOCUS15255"), slot = "counts", log = TRUE)


VlnPlot()

# Feature plot for specific gene
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS2664")) # EC
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS10578")) # EB
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS1470")) # VM
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS13387")) # EE
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS6785"))  # GC
FeaturePlot(SGSB_integ, features = c("NEZAVI-LOCUS13292")) # SC

################################ Trajectory analysis ################################################

library(ggplot2)
library(patchwork)
library(SeuratWrappers)


install.packages("remotes")

remotes::install_github("satijalab/seurat-wrappers")

## install in terminal for mac: brew install hdf5 pkg-config


devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("devtools")
remotes::install_github("bnprks/BPCells/r")

install.packages("igraph", type = "source", dependencies = TRUE)



####################### Analysis ######

library(monocle3)
library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)

# 1) Convert Seurat -> Monocle3
cds <- as.cell_data_set(SGSB_integ)


Reductions(SGSB_integ)


# 2) Put Seurat UMAP coordinates into Monocle3 BEFORE clustering
reducedDims(cds)$UMAP <- Embeddings(SGSB_integ, "UMAP")

# 3) Bring Seurat cluster identities into colData(cds)
colData(cds)$cluster <- as.character(Idents(SGSB_integ))  # make sure identities are set in Seurat

# Optional: keep nCount_RNA if present
# View(colData(cds)$nCount_RNA)

# 4) Cluster cells in monocle3 space using UMAP
cds <- cluster_cells(cds, reduction_method = "UMAP")

# 5) Plots
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
p1 + p2

# 6) Learn trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

# 7) Choose root cells from a known cluster label (example: "SC")

#root_cells <- rownames(subset(colData(cds), cluster %in% c("SC (4)")))

root_cells <- rownames(subset(as.data.frame(colData(cds)), cluster == "SC/EB"))

cds <- order_cells(cds, reduction_method = "UMAP", root_cells = root_cells)

plot_cells(cds,
           color_cells_by = "pseudotime",
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)


p_pt <- plot_cells(cds,
                   color_cells_by = "pseudotime",
                   label_groups_by_cluster = FALSE,
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   label_leaves = FALSE)

p_cl <- plot_cells(cds,
                   color_cells_by = "cluster",
                   label_groups_by_cluster = TRUE,
                   label_branch_points = FALSE,
                   label_roots = FALSE,
                   label_leaves = FALSE,
                   group_label_size = 5)

p_pt | p_cl



colnames(colData(cds))



####

library(patchwork)

p_pt <- plot_cells(
  cds,
  color_cells_by = "monocle3_pseudotime",   # or "pseudotime" if that works in your cds
  label_cell_groups = FALSE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE
)

p_pt <- p_pt +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )



p_pt

p_id <- plot_cells(
  cds,
  color_cells_by = "ident",                 # <-- your cluster identity labels
  label_cell_groups = TRUE,                 # <-- this labels groups from colData
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  group_label_size = 6
)

p_pt | p_id

cds_epi <- cds[, !colData(cds)$ident %in% c("VM")]

plot_cells(
  cds_epi,
  color_cells_by = "ident",
  label_cell_groups = TRUE
)

plot_cells(
  cds,
  color_cells_by = "ident",
  label_cell_groups = TRUE,
  label_branch_points = FALSE,
  label_roots = FALSE,
  label_leaves = FALSE,
  group_label_size = 6
)












# 8) Save
saveRDS(cds, file = "/Users/surjeetarya/Desktop/All_single_data/SGSB/Code_for_analysis/SGSB_Final_tutorial_Idents_Pseudotime_2.rds")

# 9) Pseudotime in metadata
colData(cds)$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

ggplot(data.pseudo,
       aes(monocle3_pseudotime,
           reorder(cluster, monocle3_pseudotime, median),
           fill = cluster)) +
  geom_boxplot()

# 10) Genes changing with pseudotime
deg_bcells <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)

deg_bcells %>%
  arrange(q_value) %>%
  filter(status == "OK") %>%
  head()





###################### End ##############







#############################

cds <- as.cell_data_set(SGSB_integ)
cds <- cluster_cells(cds)
p1 <- plot_cells(cds, show_trajectory_graph = FALSE)
p2 <- plot_cells(cds, color_cells_by = "partition", show_trajectory_graph = FALSE)
wrap_plots(p1, p2)

View(cds@colData$nCount_RNA)


# assign paritions
reacreate.partition <- c(rep(1,length(cds@colData@rownames)))
names(reacreate.partition) <- cds@colData@rownames
reacreate.partition <- as.factor(reacreate.partition)

cds@clusters$UMAP$partitions <- reacreate.partition

# Assign the cluster info 

SGSB_integ@active.ident

list_cluster <- SGSB_integ@active.ident
cds@clusters$UMAP$clusters <- list_cluster


# Assign UMAP coordinate - cell embeddings

cds@int_colData@listData$reducedDims$UMAP <- SGSB_integ@reductions$UMAP@cell.embeddings


# plot

cluster.before.trajectory <- plot_cells(cds,
                                        color_cells_by = 'cluster',
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster.names <- plot_cells(cds,
                            color_cells_by = "cluster",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  scale_color_manual(values = c('red', 'blue', 'green', 'maroon', 'yellow','pink','black','orange','darkgreen','violet','brown','lightgreen')) +
  theme(legend.position = "right")

cluster.before.trajectory | cluster.names



# ...3. Learn trajectory graph ------------------------
cds <- learn_graph(cds, use_partition = FALSE)

v1 <- plot_cells(cds,
                 color_cells_by = 'cluster',
                 label_groups_by_cluster = FALSE,
                 label_branch_points = FALSE,
                 label_roots = FALSE,
                 label_leaves = FALSE,
                 group_label_size = 5)

v1


# ...4. Order the cells in pseudotime -------------------

cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds) == 'SC']))

table(clusters(cds))

plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = FALSE,
           label_branch_points = FALSE,
           label_roots = FALSE,
           label_leaves = FALSE)



saveRDS(cds, file = "/Users/surjeetarya/Desktop/All_single_data/SGSB/Code_for_analysis/SGSB_Final_tutorial_Idents_Pseudotime_2.rds")


# cells ordered by monocle3 pseudotime

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))

#ggplot(data.pseudo, aes(monocle3_pseudotime, ident, fill = ident)) +
geom_boxplot()

ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(ident, monocle3_pseudotime, median), fill = ident)) +
  geom_boxplot()


# ...5. Finding genes that change as a function of pseudotime --------------------
deg_bcells <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 4)

class(cds)
dim(cds)


results_1 <- deg_bcells

results_1

deg_bcells %>% 
  arrange(q_value) %>% 
  filter(status == 'OK') %>% 
  head()

