####################################################################################################


# Tutorial: https://satijalab.org/seurat/v3.0/immune_alignment.html

cd /home/muhammad/Downloads/Projects/CS_mm_PD_SNCA/contrast_processing
# scp -r iris-cluster:/scratch/users/mali/Tasks/CS_mm_PD_SNCA/alevin_salmon_processing/MAE*/ ./

library(Seurat)
library(tximport)
library(dplyr)
library(cowplot)
library(ggplot2)

# read gene name mapping file
gene_mapping <- read.table("/home/muhammad/Downloads/Projects/CS_mm_PD_SNCA/EnsT_EnsG_GeneName.tsv", sep="\t", header=T, stringsAsFactors=F)
gene_mapping <- gene_mapping[,c(2,3)]

# mapping table with "key-value" pairs
lookUp1 <- setNames(as.character(gene_mapping$gene_name), gene_mapping$gene_id)

MAE1 <- file.path("MAE1/alevin/quants_mat.gz")
file.exists(MAE1) # TRUE
MAE2 <- file.path("MAE2/alevin/quants_mat.gz")
file.exists(MAE2) # TRUE
MAE3 <- file.path("MAE3/alevin/quants_mat.gz")
file.exists(MAE3) # TRUE
MAE4 <- file.path("MAE4/alevin/quants_mat.gz")
file.exists(MAE4) # TRUE

# import count matrix
MAE1_txi <- tximport(MAE1, type="alevin")
MAE2_txi <- tximport(MAE2, type="alevin")
MAE3_txi <- tximport(MAE3, type="alevin")
MAE4_txi <- tximport(MAE4, type="alevin")

# Convert Ensembl IDs to gene symbol in the row.names of count matrix:
res <- lapply(row.names(MAE1_txi$counts), function(i) lookUp1[i])
gene_symbol <- unlist(res, use.names=FALSE)
row.names(MAE1_txi$counts) <- gene_symbol
MAE1_WT_F <- MAE1_txi$counts
rm(MAE1_txi)

res <- lapply(row.names(MAE2_txi$counts), function(i) lookUp1[i])
gene_symbol <- unlist(res, use.names=FALSE)
row.names(MAE2_txi$counts) <- gene_symbol
MAE2_TG_F <- MAE2_txi$counts
rm(MAE2_txi)

res <- lapply(row.names(MAE3_txi$counts), function(i) lookUp1[i])
gene_symbol <- unlist(res, use.names=FALSE)
row.names(MAE3_txi$counts) <- gene_symbol
MAE3_WT_M <- MAE3_txi$counts
rm(MAE3_txi)

res <- lapply(row.names(MAE4_txi$counts), function(i) lookUp1[i])
gene_symbol <- unlist(res, use.names=FALSE)
row.names(MAE4_txi$counts) <- gene_symbol
MAE4_TG_M <- MAE4_txi$counts
rm(MAE4_txi)

dim(MAE1_WT_F);dim(MAE2_TG_F);dim(MAE3_WT_M);dim(MAE4_TG_M);
# [1] 54281  1378
# [1] 54281  2457
# [1] 54281  1345
# [1] 54281   754
save(MAE1_WT_F, MAE2_TG_F, MAE3_WT_M, MAE4_TG_M, file="count_matrix_all.RData")

# load("count_matrix_all.RData")
# Set up control object
ctrl.F <- CreateSeuratObject(counts = MAE1_WT_F, project = "MAE1.CTRL", min.cells = 5)
ctrl.F$stim <- "CTRL.F"
ctrl.F <- subset(ctrl.F, subset = nFeature_RNA > 500)
ctrl.F <- NormalizeData(ctrl.F, verbose = FALSE)
ctrl.F <- FindVariableFeatures(ctrl.F, selection.method = "vst", nfeatures = 2000)

ctrl.M <- CreateSeuratObject(counts = MAE3_WT_M, project = "MAE3.CTRL", min.cells = 5)
ctrl.M$stim <- "CTRL.M"
ctrl.M <- subset(ctrl.M, subset = nFeature_RNA > 500)
ctrl.M <- NormalizeData(ctrl.M, verbose = FALSE)
ctrl.M <- FindVariableFeatures(ctrl.M, selection.method = "vst", nfeatures = 2000)

# Set up transgenic object
stim.F <- CreateSeuratObject(counts = MAE2_TG_F, project = "MAE2.TG", min.cells = 5)
stim.F$stim <- "TG.F"
stim.F <- subset(stim.F, subset = nFeature_RNA > 500)
stim.F <- NormalizeData(stim.F, verbose = FALSE)
stim.F <- FindVariableFeatures(stim.F, selection.method = "vst", nfeatures = 2000)

stim.M <- CreateSeuratObject(counts = MAE4_TG_M, project = "MAE4.TG", min.cells = 5)
stim.M$stim <- "TG.M"
stim.M <- subset(stim.M, subset = nFeature_RNA > 500)
stim.M <- NormalizeData(stim.M, verbose = FALSE)
stim.M <- FindVariableFeatures(stim.M, selection.method = "vst", nfeatures = 2000)

save(ctrl.F, ctrl.M, stim.F, stim.M, file="seurat_objects_all.RData")

#### Perform integration #####

# load("seurat_objects_all.RData")
scNorm.anchors <- FindIntegrationAnchors(object.list = list(ctrl.F, ctrl.M, stim.F, stim.M), dims = 1:20)
scNorm <- IntegrateData(anchorset = scNorm.anchors, dims = 1:20)
save.image("mm_Cortex.RData")

table(scNorm$orig.ident) # these are project IDs

# MAE1.CTRL   MAE2.TG       MAE3.CTRL   MAE4.TG 
#              909          1059                    799             733

table(scNorm$stim) # these are stim IDs

# CTRL.F     CTRL.M        TG.F     TG.M 
#       909            799        1059        733


###################################
#### Perform an integrated analysis ####
###################################


DefaultAssay(scNorm) <- "integrated"

# Run the standard workflow for visualization and clustering
scNorm <- ScaleData(scNorm, verbose = FALSE)
scNorm <- RunPCA(scNorm, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
scNorm <- RunUMAP(scNorm, reduction = "pca", dims = 1:20)
scNorm <- FindNeighbors(scNorm, reduction = "pca", dims = 1:20)
scNorm <- FindClusters(scNorm, resolution = 0.5)
# Number of nodes: 3500
# Number of edges: 164578
# Maximum modularity in 10 random starts: 0.9110
# Number of communities: 13

# Make silhoutte plots at different clustering resolutions (of seurat) to identify best cluster number/size

library(cluster)

dist.matrix <- dist(x = Embeddings(object = scNorm[[reduction]])[, dims])
reduction <- "pca"
dims <- 1:20

scNorm <- FindClusters(scNorm, resolution = 0.01)
clusters <- scNorm$integrated_snn_res.0.01
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_7_cluster_res0.01 # Average silhoutte width: 0.55

scNorm <- FindClusters(scNorm, resolution = 0.05)
clusters <- scNorm$integrated_snn_res.0.05
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_8_cluster_res0.05 # Average silhoutte width: 0.57
# manually take screenshot of the plot as when saved in a png file colors get affected

scNorm <- FindClusters(scNorm, resolution = 0.1)
clusters <- scNorm$integrated_snn_res.0.1
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_10_cluster_res0.1 # Average silhoutte width: 0.45

scNorm <- FindClusters(scNorm, resolution = 0.2)
clusters <- scNorm$integrated_snn_res.0.2
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_10_cluster_res0.2 # Average silhoutte width: 0.45

scNorm <- FindClusters(scNorm, resolution = 0.3)
clusters <- scNorm$integrated_snn_res.0.3
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_11_cluster_res0.3 # Average silhoutte width: 0.42

scNorm <- FindClusters(scNorm, resolution = 0.4)
clusters <- scNorm$integrated_snn_res.0.4
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_12_cluster_res0.4 # Average silhoutte width: 0.42

scNorm <- FindClusters(scNorm, resolution = 0.5)
clusters <- scNorm$integrated_snn_res.0.5
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_13_cluster_res0.5 # Average silhoutte width: 0.43

scNorm <- FindClusters(scNorm, resolution = 0.6)
clusters <- scNorm$integrated_snn_res.0.6
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_14_cluster_res0.6 # Average silhoutte width: 0.41

scNorm <- FindClusters(scNorm, resolution = 0.7)
clusters <- scNorm$integrated_snn_res.0.7
sil <- silhouette(x = as.numeric(x = as.factor(x = clusters)), dist = dist.matrix)
plot(sil) # silhoutte_17_cluster_res0.7 # Average silhoutte width: 0.31


scNorm <- FindClusters(scNorm, resolution = 0.05)
# Number of nodes: 3500
# Number of edges: 164578
# Maximum modularity in 10 random starts: 0.9862
# Number of communities: 8

table(scNorm$seurat_clusters)

#       0        1        2       3        4      5       6       7 
#  1442    705     675   304     142    87     74      71

p1 <- DimPlot(scNorm, reduction = "umap", group.by = "stim")
p2 <- DimPlot(scNorm, reduction = "umap", label = TRUE)
png("DimReductionPlot_mmCortex.png", width=12, height=6, units="in", res=300)
plot_grid(p1, p2)
dev.off()

# To visualize the two conditions side-by-side, we can use the split.by argument to show each condition colored by cluster.
png("DimReductionPlot_conditionSpecific_mmCortex.png", width=14, height=6, units="in", res=300)
DimPlot(scNorm, reduction = "umap", split.by = "stim")
dev.off()


########################################
#### Identify conserved cell type markers ####
########################################


# identify markers for each cluster / cell type
DefaultAssay(scNorm) <- "RNA"
cortex_markers_0 <- FindConservedMarkers(scNorm, ident.1 = 0, grouping.var = "stim", verbose = FALSE)
cortex_markers_1 <- FindConservedMarkers(scNorm, ident.1 = 1, grouping.var = "stim", verbose = FALSE)
cortex_markers_2 <- FindConservedMarkers(scNorm, ident.1 = 2, grouping.var = "stim", verbose = FALSE)
cortex_markers_3 <- FindConservedMarkers(scNorm, ident.1 = 3, grouping.var = "stim", verbose = FALSE)
cortex_markers_4 <- FindConservedMarkers(scNorm, ident.1 = 4, grouping.var = "stim", verbose = FALSE)
cortex_markers_5 <- FindConservedMarkers(scNorm, ident.1 = 5, grouping.var = "stim", verbose = FALSE)
cortex_markers_6 <- FindConservedMarkers(scNorm, ident.1 = 6, grouping.var = "stim", verbose = FALSE)
cortex_markers_7 <- FindConservedMarkers(scNorm, ident.1 = 7, grouping.var = "stim", verbose = FALSE)
cortex_markers_0[1:10,1:4]
cortex_markers_1[1:10,1:4]
cortex_markers_2[1:10,1:4]
cortex_markers_3[1:10,1:4]
cortex_markers_4[1:10,1:4]
cortex_markers_5[1:10,1:4]
cortex_markers_6[1:10,1:4]
cortex_markers_7[1:10,1:4]

save.image("mm_Cortex_v2.RData")

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
# load("mm_Cortex_v2.RData")

# Aim: Manually go through marker genes and find in which cell they are uniquely expressed and then label the cell types accordingly.

# http://bio-bigdata.hrbmu.edu.cn/CellMarker/index.jsp # https://academic.oup.com/nar/article/47/D1/D721/5115823 # https://panglaodb.se/markers.html
cluster 0 = P2ry12 (also Hexb, Ctss, Cx3cr1, C1qa, C1qb) = microglial cells
cluster 1 = Flt1 (also Cldn5,  Ly6c1, Ly6a, Adgrf5, Abcb1a) = endothelial cell # https://www.novusbio.com/products/vegfr1-flt-1-antibody_nb100-527 # http://bio-bigdata.hrbmu.edu.cn/CellMarker/search.jsp?quickSearchInfo=Flt1
cluster 2 = Clu (also Aldoc, Gja1) = quiescent NSC # http://bio-bigdata.hrbmu.edu.cn/CellMarker/search.jsp?quickSearchInfo=clu
cluster 2 = Pla2g7 (also Plpp3) = Bergmann glial cell # http://bio-bigdata.hrbmu.edu.cn/CellMarker/
cluster 3 = Cldn11 (also Ermn, Gm21984, Gjc3, Mog, Mag) = Oligodendrocyte
cluster 4 = Snap25 (neuron), Scn2b, Tnr (type 1C- 2- spiral ganglion neuron)  = ganglion neuron type2
cluster 5 = Dlx1,Evf1 or Dlx6os1 (neuroblast), Cd24a, Tubb3 (type 1 spiral gang Neuron), Sox11 (early precursor cell) = neuroblast
cluster 6 = Fabp7 (NSC, Olfactory ensheathing glia), Apod (Oligodendrocyte), Ncald (type 1C spiral gang neuron), Nnat (type 2 spiral gang neuron), Slc35f1, Matn4 (oligo dendro precursor cell) = oligo dendro precursor
cluster 7 = Igf2 (Mural cell), Htr2c (type 1B spiral gang N), Cox8b, Atp2b3 (type 1 spiral gang N), Otx2 (Ventral otocyst) = ganglion neuron type1

# annotate the clusters as specific cell types
scNorm <- RenameIdents(scNorm, `0` = "microglial cells", `1` = "endothelial cell", `2` = "quiescent NSC", `3` = "oligodendrocyte", 
 `4` = "ganglion neuron type2", `5` = "neuroblast", `6` = "oligodendrocyte precursor", `7` = "ganglion neuron type1")

# we can explore these marker genes for each cluster
png("Clusters_marker_mmCortex.png", width=10, height=6, units="in", res=300)
FeaturePlot(scNorm, features = c("P2ry12", "Flt1", "Clu", "Pla2g7", "Cldn11", "Tnr", "Dlx1", "Slc35f1", "Atp2b3"), min.cutoff = "q9")
dev.off()

# To visualize the cell type labelled clusters.
png("Labelled_clusters_mmCortex.png", width=10, height=6, units="in", res=300)
DimPlot(scNorm, label = TRUE)
dev.off()


Idents(scNorm) <- factor(Idents(scNorm), levels = c("microglial cells", "endothelial cell", "quiescent NSC", "oligodendrocyte",
 "ganglion neuron type2", "neuroblast", "oligodendrocyte precursor", "ganglion neuron type1"))
# markers.to.plot <- c("P2ry12", "Hexb", "Ctss", "Flt1", "Cldn5", "Ly6c1", "Clu", "Aldoc", "Gja1", "Cldn11", "Ermn", "Gm21984",
# "Snap25", "Scn2b", "Tnr", "Dlx1", "Dlx6os1", "Slc35f1", "Matn4", "Apod", "Cox8b", "Atp2b3", "Htr2c")
markers.to.plot <- c("P2ry12", "Hexb", "Ctss", "Flt1", "Cldn5", "Ly6c1", "Pla2g7", "Plpp3", "Slc1a2", "Cldn11", "Ermn", "Gm21984",
 "Snhg11",  "Kcnq1ot1", "Meg3", "Dlx1", "Dlx6os1", "Tubb3", "Fabp7", "Igfbp5", "Ift57", "Cox8b", "Atp2b3", "Htr2c")
 
png("Cluster_marker_acrossCondition_mmCortex.png", width=12, height=10, units="in", res=300)
DotPlot(scNorm, features = rev(markers.to.plot), cols = c("blue", "darkred", "darkgreen", "orange"), dot.scale = 8, split.by = "stim") + RotatedAxis()
dev.off()
 
idents_later_use <- Idents(scNorm)
 
####################################################### 
#### Identify differential expressed genes across conditions ####
####################################################### 


# cluster 0 = P2ry12 (also Hexb, Ctss, Cx3cr1, C1qa, C1qb) = microglial cells
# cluster 1 = Flt1 (also Cldn5,  Ly6c1, Ly6a, Adgrf5, Abcb1a) = endothelial cell

microglial.cells <- subset(scNorm, idents = "microglial cells")
Idents(microglial.cells) <- "stim"
avg.microglial.cells <- log1p(AverageExpression(microglial.cells, verbose = FALSE)$RNA)
avg.microglial.cells$gene <- rownames(avg.microglial.cells)
# avg.microglial.cells$Avg <- rowMeans(avg.microglial.cells)
# avg.microglial.cells <- avg.microglial.cells[order(-avg.microglial.cells$Avg),]

endothelial.cells <- subset(scNorm, idents = "endothelial cell")
Idents(endothelial.cells) <- "stim"
avg.endothelial.cells <- log1p(AverageExpression(endothelial.cells, verbose = FALSE)$RNA)
avg.endothelial.cells$gene <- rownames(avg.endothelial.cells)

gang_neuron_2 <- subset(scNorm, idents = "ganglion neuron type2")
Idents(gang_neuron_2) <- "stim"
avg_gang_neuron_2 <- log1p(AverageExpression(gang_neuron_2, verbose = FALSE)$RNA)
avg_gang_neuron_2$Avg <- rowMeans(avg_gang_neuron_2)
avg_gang_neuron_2 <- avg_gang_neuron_2[order(-avg_gang_neuron_2$Avg),] # some of the top ranked genes are CT-specific while most are house-keeping kind of genes.
avg_gang_neuron_2$gene <- rownames(avg_gang_neuron_2)


genes.to.label = c("P2ry12", "Hexb", "Ctss", "Cx3cr1", "C1qa", "Flt1", "Cldn5", "Ly6c1", "Adgrf5", "Abcb1a")
p1 <- ggplot(avg.microglial.cells, aes(CTRL.F, TG.F)) + geom_point() + ggtitle("Microglial Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.endothelial.cells, aes(CTRL.M, TG.M)) + geom_point() + ggtitle("Endothelial Cells")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p3 <- ggplot(avg.microglial.cells, aes(CTRL.M, TG.M)) + geom_point() + ggtitle("Microglial Cells")
p3 <- LabelPoints(plot = p3, points = genes.to.label, repel = TRUE)
p4 <- ggplot(avg.endothelial.cells, aes(CTRL.F, TG.F)) + geom_point() + ggtitle("Endothelial Cells")
p4 <- LabelPoints(plot = p4, points = genes.to.label, repel = TRUE)
png("Average_markerExpression_acrossCTs_mmCortex.png", width=14, height=12, units="in", res=300)
plot_grid(p1, p2, p3, p4)
dev.off() 

# what genes change in different conditions for cells of the same type
scNorm$celltype.stim <- paste(Idents(scNorm), scNorm$stim, sep = ".")
scNorm$celltype <- Idents(scNorm)
Idents(scNorm) <- "celltype.stim"


################################### DEGs based on within-CT comparison that are geneder-specific ################################


# see columns that can be used to make all cell type comparisons
# all CTs are divided into 4 groups (CTRL.F, CTRL.M, TG.F, TG.M)
table(scNorm$celltype.stim)

oligodendrocyte_F <- FindMarkers(scNorm, ident.1 = "oligodendrocyte.TG.F", ident.2 = "oligodendrocyte.CTRL.F", verbose = FALSE)
dim(oligodendrocyte_F) # 964   5
nrow(oligodendrocyte_F[oligodendrocyte_F$p_val_adj < 0.05,]) # 1
nrow(oligodendrocyte_F[oligodendrocyte_F$p_val < 0.05,]) # 251
oligodendrocyte_M <- FindMarkers(scNorm, ident.1 = "oligodendrocyte.TG.M", ident.2 = "oligodendrocyte.CTRL.M", verbose = FALSE)
dim(oligodendrocyte_M) # 1061    5
nrow(oligodendrocyte_M[oligodendrocyte_M$p_val_adj < 0.05,]) # 2
nrow(oligodendrocyte_M[oligodendrocyte_M$p_val < 0.05,]) # 308

endothelial_cell_F <- FindMarkers(scNorm, ident.1 = "endothelial cell.TG.F", ident.2 = "endothelial cell.CTRL.F", verbose = FALSE)
dim(endothelial_cell_F) # 538   5
nrow(endothelial_cell_F[endothelial_cell_F$p_val_adj < 0.05,]) # 7
nrow(endothelial_cell_F[endothelial_cell_F$p_val < 0.05,]) # 309
endothelial_cell_M <- FindMarkers(scNorm, ident.1 = "endothelial cell.TG.M", ident.2 = "endothelial cell.CTRL.M", verbose = FALSE)
dim(endothelial_cell_M) # 685    5
nrow(endothelial_cell_M[endothelial_cell_M$p_val_adj < 0.05,]) # 6
nrow(endothelial_cell_M[endothelial_cell_M$p_val < 0.05,]) # 414

microglial_cell_F <- FindMarkers(scNorm, ident.1 = "microglial cells.TG.F", ident.2 = "microglial cells.CTRL.F", verbose = FALSE)
dim(microglial_cell_F) # 199   5
nrow(microglial_cell_F[microglial_cell_F$p_val_adj < 0.05,]) # 29
nrow(microglial_cell_F[microglial_cell_F$p_val < 0.05,]) # 166
microglial_cell_M <- FindMarkers(scNorm, ident.1 = "microglial cells.TG.M", ident.2 = "microglial cells.CTRL.M", verbose = FALSE)
dim(microglial_cell_M) # 229    5
nrow(microglial_cell_M[microglial_cell_M$p_val_adj < 0.05,]) # 16
nrow(microglial_cell_M[microglial_cell_M$p_val < 0.05,]) # 186

neuroblast_F <- FindMarkers(scNorm, ident.1 = "neuroblast.TG.F", ident.2 = "neuroblast.CTRL.F", verbose = FALSE)
dim(neuroblast_F) # 2590   5
nrow(neuroblast_F[neuroblast_F$p_val_adj < 0.05,]) # 0
nrow(neuroblast_F[neuroblast_F$p_val < 0.05,]) # 347
neuroblast_M <- FindMarkers(scNorm, ident.1 = "neuroblast.TG.M", ident.2 = "neuroblast.CTRL.M", verbose = FALSE)
dim(neuroblast_M) # 3425    5
nrow(neuroblast_M[neuroblast_M$p_val_adj < 0.05,]) # 0
nrow(neuroblast_M[neuroblast_M$p_val < 0.05,]) # 319

oligodendrocyte_precursor_F <- FindMarkers(scNorm, ident.1 = "oligodendrocyte precursor.TG.F", ident.2 = "oligodendrocyte precursor.CTRL.F", verbose = FALSE)
dim(oligodendrocyte_precursor_F) # 3698   5
nrow(oligodendrocyte_precursor_F[oligodendrocyte_precursor_F$p_val_adj < 0.05,]) # 0
nrow(oligodendrocyte_precursor_F[oligodendrocyte_precursor_F$p_val < 0.05,]) # 673
oligodendrocyte_precursor_M <- FindMarkers(scNorm, ident.1 = "oligodendrocyte precursor.TG.M", ident.2 = "oligodendrocyte precursor.CTRL.M", verbose = FALSE)
dim(oligodendrocyte_precursor_M) # 4918    5
nrow(oligodendrocyte_precursor_M[oligodendrocyte_precursor_M$p_val_adj < 0.05,]) # 0
nrow(oligodendrocyte_precursor_M[oligodendrocyte_precursor_M$p_val < 0.05,]) # 707

quiescent_NSC_F <- FindMarkers(scNorm, ident.1 = "quiescent NSC.TG.F", ident.2 = "quiescent NSC.CTRL.F", verbose = FALSE)
dim(quiescent_NSC_F) # 365   5
nrow(quiescent_NSC_F[quiescent_NSC_F$p_val_adj < 0.05,]) # 6
nrow(quiescent_NSC_F[quiescent_NSC_F$p_val < 0.05,]) # 218
quiescent_NSC_M <- FindMarkers(scNorm, ident.1 = "quiescent NSC.TG.M", ident.2 = "quiescent NSC.CTRL.M", verbose = FALSE)
dim(quiescent_NSC_M) # 634    5
nrow(quiescent_NSC_M[quiescent_NSC_M$p_val_adj < 0.05,]) # 3
nrow(quiescent_NSC_M[quiescent_NSC_M$p_val < 0.05,]) # 301

ganglion_neuron_type1_F <- FindMarkers(scNorm, ident.1 = "ganglion neuron type1.TG.F", ident.2 = "ganglion neuron type1.CTRL.F", verbose = FALSE)
dim(ganglion_neuron_type1_F) # 3444   5
nrow(ganglion_neuron_type1_F[ganglion_neuron_type1_F$p_val_adj < 0.05,]) # 0
nrow(ganglion_neuron_type1_F[ganglion_neuron_type1_F$p_val < 0.05,]) # 181
ganglion_neuron_type1_M <- FindMarkers(scNorm, ident.1 = "ganglion neuron type1.TG.M", ident.2 = "ganglion neuron type1.CTRL.M", verbose = FALSE)
dim(ganglion_neuron_type1_M) # 2322    5
nrow(ganglion_neuron_type1_M[ganglion_neuron_type1_M$p_val_adj < 0.05,]) # 0
nrow(ganglion_neuron_type1_M[ganglion_neuron_type1_M$p_val < 0.05,]) # 242

ganglion_neuron_type2_F <- FindMarkers(scNorm, ident.1 = "ganglion neuron type2.TG.F", ident.2 = "ganglion neuron type2.CTRL.F", verbose = FALSE)
dim(ganglion_neuron_type2_F) # 2374   5
nrow(ganglion_neuron_type2_F[ganglion_neuron_type2_F$p_val_adj < 0.05,]) # 0
nrow(ganglion_neuron_type2_F[ganglion_neuron_type2_F$p_val < 0.05,]) # 692
ganglion_neuron_type2_M <- FindMarkers(scNorm, ident.1 = "ganglion neuron type2.TG.M", ident.2 = "ganglion neuron type2.CTRL.M", verbose = FALSE)
dim(ganglion_neuron_type2_M) # 2772    5
nrow(ganglion_neuron_type2_M[ganglion_neuron_type2_M$p_val_adj < 0.05,]) # 1
nrow(ganglion_neuron_type2_M[ganglion_neuron_type2_M$p_val < 0.05,]) # 586

# Make a table of differential expression result for comparison between all the cell types
# See for-loop code at the bottom of this document.


#################################### DEGs based on CT comparison with all remaining CTs (clusters) ###################################


# get the names of all cell types along with conditions and gender
identity_names <- names(table(scNorm$celltype.stim))

# get markers of every cell type / cluster
microglial_cells <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "microglial cells")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "microglial cells")], verbose = FALSE)
dim(microglial_cells) # 2158    5
nrow(microglial_cells[microglial_cells$p_val_adj < 0.05,]) # 2075
nrow(microglial_cells[microglial_cells$p_val < 0.05,]) # 2156

endothelial_cell <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "endothelial cell")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "endothelial cell")], verbose = FALSE)
dim(endothelial_cell) # 1824    5
nrow(endothelial_cell[endothelial_cell$p_val_adj < 0.05,]) # 1676
nrow(endothelial_cell[endothelial_cell$p_val < 0.05,]) # 1813

quiescent_NSC <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "quiescent NSC")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "quiescent NSC")], verbose = FALSE)
dim(quiescent_NSC) # 2132    5
nrow(quiescent_NSC[quiescent_NSC$p_val_adj < 0.05,]) # 1965
nrow(quiescent_NSC[quiescent_NSC$p_val < 0.05,]) # 2125

oligodendrocyte <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "oligodendrocyte")], verbose = FALSE)
dim(oligodendrocyte) # 1986    5
nrow(oligodendrocyte[oligodendrocyte$p_val_adj < 0.05,]) # 1488
nrow(oligodendrocyte[oligodendrocyte$p_val < 0.05,]) # 1936

ganglion_neuron_t2 <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type2")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "ganglion neuron type2")], verbose = FALSE)
dim(ganglion_neuron_t2) # 2394    5
nrow(ganglion_neuron_t2[ganglion_neuron_t2$p_val_adj < 0.05,]) # 1093
nrow(ganglion_neuron_t2[ganglion_neuron_t2$p_val < 0.05,]) # 2195

neuroblast <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "neuroblast")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "neuroblast")], verbose = FALSE)
dim(neuroblast) # 2553    5
nrow(neuroblast[neuroblast$p_val_adj < 0.05,]) # 1040
nrow(neuroblast[neuroblast$p_val < 0.05,]) # 2177

oligodendro_precursor <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte precursor")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "oligodendrocyte precursor")], verbose = FALSE)
dim(oligodendro_precursor) # 1839    5
nrow(oligodendro_precursor[oligodendro_precursor$p_val_adj < 0.05,]) # 439
nrow(oligodendro_precursor[oligodendro_precursor$p_val < 0.05,]) # 1310

ganglion_neuron_t1 <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type1")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) != "ganglion neuron type1")], verbose = FALSE)
dim(ganglion_neuron_t1) # 2502    5
nrow(ganglion_neuron_t1[ganglion_neuron_t1$p_val_adj < 0.05,]) # 845
nrow(ganglion_neuron_t1[ganglion_neuron_t1$p_val < 0.05,]) # 1904

# plot the markers
ct_markers <- c(row.names(microglial_cells)[1:5], row.names(endothelial_cell)[1:5], row.names(quiescent_NSC)[1:5], row.names(oligodendrocyte)[1:5],
  row.names(ganglion_neuron_t2)[1:5], row.names(neuroblast)[1:5], row.names(oligodendro_precursor)[1:5], row.names(ganglion_neuron_t1)[1:5])

# idents_later_use <- Idents(scNorm)
Idents(scNorm) <- idents_later_use
 
png("Cluster_marker_acrossCondition_mmCortex_v2.png", width=14, height=10, units="in", res=300)
DotPlot(scNorm, features = rev(ct_markers), cols = c("blue", "darkred", "darkgreen", "orange"), dot.scale = 8, split.by = "stim") + RotatedAxis()
dev.off()


################################### get gender-specific markers (DEGs) for every cell type / cluster (with-in CT) ###################################


# change back the Idents of scNorm
scNorm$celltype.stim <- paste(Idents(scNorm), scNorm$stim, sep = ".")
scNorm$celltype <- Idents(scNorm)
Idents(scNorm) <- "celltype.stim"


microglial_cells_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "microglial cells" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "microglial cells" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(microglial_cells_gender) # 347   5
nrow(microglial_cells_gender[microglial_cells_gender$p_val_adj < 0.05,]) # 157
nrow(microglial_cells_gender[microglial_cells_gender$p_val < 0.05,]) # 343

endothelial_cells_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "endothelial cell" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "endothelial cell" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(endothelial_cells_gender) # 432   5
nrow(endothelial_cells_gender[endothelial_cells_gender$p_val_adj < 0.05,]) # 77
nrow(endothelial_cells_gender[endothelial_cells_gender$p_val < 0.05,]) # 382

quiescent_NSC_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "quiescent NSC" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "quiescent NSC" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(quiescent_NSC_gender) # 524   5
nrow(quiescent_NSC_gender[quiescent_NSC_gender$p_val_adj < 0.05,]) # 80
nrow(quiescent_NSC_gender[quiescent_NSC_gender$p_val < 0.05,]) # 479

oligodendrocyte_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(oligodendrocyte_gender) # 779   5
nrow(oligodendrocyte_gender[oligodendrocyte_gender$p_val_adj < 0.05,]) # 31
nrow(oligodendrocyte_gender[oligodendrocyte_gender$p_val < 0.05,]) # 526

ganglion_neuron_t2_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type2" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type2" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(ganglion_neuron_t2_gender) # 2369   5
nrow(ganglion_neuron_t2_gender[ganglion_neuron_t2_gender$p_val_adj < 0.05,]) # 12
nrow(ganglion_neuron_t2_gender[ganglion_neuron_t2_gender$p_val < 0.05,]) # 1238

neuroblast_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "neuroblast" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "neuroblast" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(neuroblast_gender) # 2281   5
nrow(neuroblast_gender[neuroblast_gender$p_val_adj < 0.05,]) # 6
nrow(neuroblast_gender[neuroblast_gender$p_val < 0.05,]) # 721

oligodendro_precursor_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte precursor" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "oligodendrocyte precursor" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(oligodendro_precursor_gender) # 1968   5
nrow(oligodendro_precursor_gender[oligodendro_precursor_gender$p_val_adj < 0.05,]) # 3
nrow(oligodendro_precursor_gender[oligodendro_precursor_gender$p_val < 0.05,]) # 423

ganglion_neuron_t1_gender <- FindMarkers(scNorm, ident.1 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type1" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "M")], 
 ident.2 = identity_names[which(sapply(strsplit(as.character(identity_names), "\\."), "[", 1) == "ganglion neuron type1" & sapply(strsplit(as.character(identity_names), "\\."), "[", 3) == "F")], verbose = FALSE)
dim(ganglion_neuron_t1_gender) # 1737   5
nrow(ganglion_neuron_t1_gender[ganglion_neuron_t1_gender$p_val_adj < 0.05,]) # 1
nrow(ganglion_neuron_t1_gender[ganglion_neuron_t1_gender$p_val < 0.05,]) # 302


# plot the markers
ct_markers <- c(row.names(microglial_cells_gender)[1:5], row.names(endothelial_cells_gender)[1:5], row.names(quiescent_NSC_gender)[1:5], row.names(oligodendrocyte_gender)[1:5],
  row.names(ganglion_neuron_t2_gender)[1:5], row.names(neuroblast_gender)[1:5], row.names(oligodendro_precursor_gender)[1:5], row.names(ganglion_neuron_t1_gender)[1:5])
ct_markers <- unique(ct_markers)

# idents_later_use <- Idents(scNorm)
Idents(scNorm) <- idents_later_use
 
png("Cluster_marker_acrossCondition_mmCortex_Gender.png", width=14, height=10, units="in", res=300)
DotPlot(scNorm, features = rev(ct_markers), cols = c("blue", "darkred", "darkgreen", "orange"), dot.scale = 8, split.by = "stim") + RotatedAxis()
dev.off()


################################################# Save the datasets ############################################################


save.image("mm_Cortex_v3.RData")
save(microglial_cells, endothelial_cell, quiescent_NSC, oligodendrocyte, ganglion_neuron_t2, neuroblast, oligodendro_precursor, ganglion_neuron_t1,
 microglial_cells_gender, endothelial_cells_gender, quiescent_NSC_gender, oligodendrocyte_gender, ganglion_neuron_t2_gender, neuroblast_gender, 
 oligodendro_precursor_gender, ganglion_neuron_t1_gender, identity_names, 
 oligodendrocyte_F, oligodendrocyte_M, endothelial_cell_F, endothelial_cell_M, microglial_cell_F, microglial_cell_M, 
 neuroblast_F, neuroblast_M, oligodendrocyte_precursor_F, oligodendrocyte_precursor_M, quiescent_NSC_F, quiescent_NSC_M,
 ganglion_neuron_type1_F, ganglion_neuron_type1_M, ganglion_neuron_type2_F, ganglion_neuron_type2_M, file="DEGs_mmCortex.RData")

CT_specific_local <- list(ganglion_neuron_type1_F, ganglion_neuron_type1_M, ganglion_neuron_type2_F, ganglion_neuron_type2_M, microglial_cell_M, microglial_cell_F,
 neuroblast_F, neuroblast_M, oligodendrocyte_F, oligodendrocyte_M, quiescent_NSC_F, quiescent_NSC_M, endothelial_cell_M, endothelial_cell_F,
  oligodendrocyte_precursor_M, oligodendrocyte_precursor_F)
names(CT_specific_local) <- c("ganglion_neuron_type1_F", "ganglion_neuron_type1_M", "ganglion_neuron_type2_F", "ganglion_neuron_type2_M", "microglial_cell_M", "microglial_cell_F",
 "neuroblast_F", "neuroblast_M", "oligodendrocyte_F", "oligodendrocyte_M", "quiescent_NSC_F", "quiescent_NSC_M", "endothelial_cell_M", "endothelial_cell_F",
  "oligodendrocyte_precursor_M", "oligodendrocyte_precursor_F")
 
for (i in 1:length(CT_specific_local)){
  name <- names(CT_specific_local)[[i]]
  write.table(CT_specific_local[[i]], file=paste0(name, "_CT_specific_local_level1.txt", sep=""), sep="\t", quote=F)
}
 

CT_specific_global <- list(endothelial_cell, microglial_cells, neuroblast, oligodendro_precursor, oligodendrocyte, quiescent_NSC, ganglion_neuron_t1, ganglion_neuron_t2)
names(CT_specific_global) <- c("endothelial_cell", "microglial_cells", "neuroblast", "oligodendro_precursor", "oligodendrocyte", "quiescent_NSC", "ganglion_neuron_t1", "ganglion_neuron_t2")

for (i in 1:length(CT_specific_global)){
  name <- names(CT_specific_global)[[i]]
  write.table(CT_specific_global[[i]], file=paste0(name, "_CT_specific_global_level2.txt", sep=""), sep="\t", quote=F)
}


gender_specific <- list(microglial_cells_gender, endothelial_cells_gender, quiescent_NSC_gender, oligodendrocyte_gender, ganglion_neuron_t2_gender, neuroblast_gender, 
 oligodendro_precursor_gender, ganglion_neuron_t1_gender)
names(gender_specific) <- c("microglial_cells_gender", "endothelial_cells_gender", "quiescent_NSC_gender", "oligodendrocyte_gender", "ganglion_neuron_t2_gender", "neuroblast_gender", 
 "oligodendro_precursor_gender", "ganglion_neuron_t1_gender")
 
for (i in 1:length(gender_specific)){
  name <- names(gender_specific)[[i]]
  write.table(gender_specific[[i]], file=paste0(name, "_CT_specific_level3.txt", sep=""), sep="\t", quote=F)
}
 


################################################# Plot the intersection of DEGs ###################################################


library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
# load("mm_Cortex_v3.RData")

load("DEGs_mmCortex.RData")

# CT-specific DEGs w.r.t. All othe CTs (clusters) as a background set 
CT_specific <- c("endothelial_cell", "microglial_cells", "neuroblast", "oligodendro_precursor", "oligodendrocyte", "quiescent_NSC", "ganglion_neuron_t1", "ganglion_neuron_t2")

endothelial_cell <- endothelial_cell[endothelial_cell$p_val_adj < 0.05,]
endothelial_cell <- row.names(endothelial_cell)
length(endothelial_cell) # 1676
microglial_cells <- microglial_cells[microglial_cells$p_val_adj < 0.05,]
microglial_cells <- row.names(microglial_cells)
length(microglial_cells) # 2075
neuroblast <- neuroblast[neuroblast$p_val_adj < 0.05,]
neuroblast <- row.names(neuroblast)
length(neuroblast) # 1040
oligodendro_precursor <- oligodendro_precursor[oligodendro_precursor$p_val_adj < 0.05,]
oligodendro_precursor <- row.names(oligodendro_precursor)
length(oligodendro_precursor) # 439
oligodendrocyte <- oligodendrocyte[oligodendrocyte$p_val_adj < 0.05,]
oligodendrocyte <- row.names(oligodendrocyte)
length(oligodendrocyte) # 1488
quiescent_NSC <- quiescent_NSC[quiescent_NSC$p_val_adj < 0.05,]
quiescent_NSC <- row.names(quiescent_NSC)
length(quiescent_NSC) # 1965
ganglion_neuron_t1 <- ganglion_neuron_t1[ganglion_neuron_t1$p_val_adj < 0.05,]
ganglion_neuron_t1 <- row.names(ganglion_neuron_t1)
length(ganglion_neuron_t1) # 845
ganglion_neuron_t2 <- ganglion_neuron_t2[ganglion_neuron_t2$p_val_adj < 0.05,]
ganglion_neuron_t2 <- row.names(ganglion_neuron_t2)
length(ganglion_neuron_t2) # 1093

# overlap plot # common genes plot 
library(UpSetR)
x = list(endothelial = endothelial_cell, microgial = microglial_cells, neuroblasts = neuroblast, 
	oligo_dendro_Prec = oligodendro_precursor, oligo_dendro = oligodendrocyte, qNSC = quiescent_NSC, 
	gang_NeuronT1 = ganglion_neuron_t1, gang_NeuronT2 = ganglion_neuron_t2)
png("Overlap_DEGs_CTspecific.png", width=14, height=10, units="in", res=300)
upset(fromList(x), sets=c("endothelial", "microgial", "neuroblasts", "oligo_dendro_Prec", "oligo_dendro", "qNSC", "gang_NeuronT1", "gang_NeuronT2"), order.by = "degree")
dev.off()

# Gender-specific DEGs w.r.t. within-CT as a background set 
gender_specific <- c("microglial_cells_gender", "endothelial_cells_gender", "quiescent_NSC_gender", "oligodendrocyte_gender", "ganglion_neuron_t2_gender", "neuroblast_gender", 
 "oligodendro_precursor_gender", "ganglion_neuron_t1_gender")
 
endothelial_cells_gender <- endothelial_cells_gender[endothelial_cells_gender$p_val_adj < 0.05,]
endothelial_cells_gender <- row.names(endothelial_cells_gender)
length(endothelial_cells_gender) # 77
microglial_cells_gender <- microglial_cells_gender[microglial_cells_gender$p_val_adj < 0.05,]
microglial_cells_gender <- row.names(microglial_cells_gender)
length(microglial_cells_gender) # 157
neuroblast_gender <- neuroblast_gender[neuroblast_gender$p_val_adj < 0.05,]
neuroblast_gender <- row.names(neuroblast_gender)
length(neuroblast_gender) # 6
oligodendro_precursor_gender <- oligodendro_precursor_gender[oligodendro_precursor_gender$p_val_adj < 0.05,]
oligodendro_precursor_gender <- row.names(oligodendro_precursor_gender)
length(oligodendro_precursor_gender) # 3
oligodendrocyte_gender <- oligodendrocyte_gender[oligodendrocyte_gender$p_val_adj < 0.05,]
oligodendrocyte_gender <- row.names(oligodendrocyte_gender)
length(oligodendrocyte_gender) # 31
quiescent_NSC_gender <- quiescent_NSC_gender[quiescent_NSC_gender$p_val_adj < 0.05,]
quiescent_NSC_gender <- row.names(quiescent_NSC_gender)
length(quiescent_NSC_gender) # 80
ganglion_neuron_t1_gender <- ganglion_neuron_t1_gender[ganglion_neuron_t1_gender$p_val_adj < 0.05,]
ganglion_neuron_t1_gender <- row.names(ganglion_neuron_t1_gender)
length(ganglion_neuron_t1_gender) # 1
ganglion_neuron_t2_gender <- ganglion_neuron_t2_gender[ganglion_neuron_t2_gender$p_val_adj < 0.05,]
ganglion_neuron_t2_gender <- row.names(ganglion_neuron_t2_gender)
length(ganglion_neuron_t2_gender) # 12

# overlap plot # common genes plot 
library(UpSetR)
x = list(endothelial = endothelial_cells_gender, microgial = microglial_cells_gender, neuroblasts = neuroblast_gender, 
	oligo_dendro_Prec = oligodendro_precursor_gender, oligo_dendro = oligodendrocyte_gender, qNSC = quiescent_NSC_gender, 
	gang_NeuronT1 = ganglion_neuron_t1_gender, gang_NeuronT2 = ganglion_neuron_t2_gender)
png("Overlap_DEGs_genderSpecific.png", width=14, height=10, units="in", res=300)
upset(fromList(x), sets=c("endothelial", "microgial", "neuroblasts", "oligo_dendro_Prec", "oligo_dendro", "qNSC", "gang_NeuronT1", "gang_NeuronT2"), order.by = "degree")
dev.off()


###########################################################################################################################################################


#########################################################################################################################################
################################################### Pseudotime: Build single-cell trajectories ###################################################
#########################################################################################################################################


# http://cole-trapnell-lab.github.io/monocle-release/docs/#constructing-single-cell-trajectories

library(Seurat)
library(monocle)
library(DDRTree)
library(pheatmap)
library(reshape2)

setwd("/home/muhammad/Downloads/Projects/CS_mm_PD_SNCA/contrast_processing")

# Importing Seurat object
load("scNorm_for_monocle.RData")

data <- as(as.matrix(scNorm[["RNA"]]@counts), 'sparseMatrix') # get Seurat's count matrix
pd <- new('AnnotatedDataFrame', data = scNorm@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

#Construct monocle cds
monocle_cds <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily=negbinomial.size())
save(monocle_cds, file="monocle_cds_v1.RData")

# Estimate size factors and dispersions
monocle_cds <- estimateSizeFactors(monocle_cds)
monocle_cds <- estimateDispersions(monocle_cds) # Removing 32 outliers
save(monocle_cds, file="monocle_cds_v2.RData")

monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
print(head(fData(monocle_cds)))

# subset those genes that are expressed in more than 10 cells
expressed_genes <- row.names(subset(fData(monocle_cds),
    num_cells_expressed >= 10))
length(expressed_genes) # 15965

pData(monocle_cds)$Total_mRNAs <- Matrix::colSums(exprs(monocle_cds))

monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs < 1e6]
dim(monocle_cds)
# Features  Samples 
#    17988     3500

upper_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) +
            2*sd(log10(pData(monocle_cds)$Total_mRNAs)))
upper_bound # 4156.203
lower_bound <- 10^(mean(log10(pData(monocle_cds)$Total_mRNAs)) -
            2*sd(log10(pData(monocle_cds)$Total_mRNAs)))
lower_bound # 596.5125

# distribution of mRNA totals across the cells
png("mRNA_dist_across_samples_monocle.png", width=12, height=8, units="in", res=300)
qplot(Total_mRNAs, data = pData(monocle_cds), color = stim, geom ="density") +
geom_vline(xintercept = lower_bound) +
geom_vline(xintercept = upper_bound)
dev.off()

# remove cells with either very low mRNA recovery or far more mRNA than the typical cell.
monocle_cds <- monocle_cds[,pData(monocle_cds)$Total_mRNAs > lower_bound &
      pData(monocle_cds)$Total_mRNAs < upper_bound]
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
dim(monocle_cds)
# Features  Samples 
#    17988     3370

# Log-transform each value in the expression matrix.
L <- log(exprs(monocle_cds[expressed_genes,]))

# Standardize each gene, so that they are all on the same scale,
# Then melt the data with plyr so we can plot it easily
melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L)))) # plotting as it is gives out of memory error, therefore removing NAs.
dim(melted_dens_df) # 53802050        3
melted_dens_df <- na.omit(melted_dens_df)
dim(melted_dens_df) # 3370    3

# Plot the distribution of the standardized gene expression values.
png("standardized_mRNA_dist_monocle.png", width=12, height=8, units="in", res=300)
qplot(value, geom = "density", data = melted_dens_df) +
stat_function(fun = dnorm, size = 0.5, color = 'red') +
xlab("Standardized log(FPKM)") +
ylab("Density")
dev.off()

save.image(file="monocle_cds_v3.RData")


# Classifying cells by type


#Ctss_id <- row.names(subset(fData(monocle_cds), gene_short_name == "Ctss"))
C1qa_id <- row.names(subset(fData(monocle_cds), gene_short_name == "C1qa"))
Flt1_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Flt1"))
Plpp3_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Plpp3"))
Cldn11_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Cldn11"))
Snhg11_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Snhg11"))
Dlx1_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Dlx1"))
Npy_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Npy"))
Ecrg4_id <- row.names(subset(fData(monocle_cds),gene_short_name == "Ecrg4"))

# make CellTypeHierarchy object
cth <- newCellTypeHierarchy()
cth <- addCellType(cth, "microglial_cells", classify_func = function(x) 
 { x[C1qa_id,] >= 1 })
cth <- addCellType(cth, "endothelial_cells", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] > 1 })
cth <- addCellType(cth, "quiescent_NSC", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] > 1 })
cth <- addCellType(cth, "oligodendrocyte", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] < 1 & x[Cldn11_id,] > 1 })
cth <- addCellType(cth, "ganglion_neuronT2", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] < 1 & x[Cldn11_id,] < 1 & x[Snhg11_id,] > 1 })
cth <- addCellType(cth, "neuroblasts", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] < 1 & x[Cldn11_id,] < 1 & x[Snhg11_id,] < 1 & x[Dlx1_id,] > 1 })
cth <- addCellType(cth, "oligodendrocyte_precursors", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] < 1 & x[Cldn11_id,] < 1 & x[Snhg11_id,] < 1 & x[Dlx1_id,] < 1 & x[Npy_id,] > 1 })
cth <- addCellType(cth, "ganglion_neuronT1", classify_func = function(x)
 { x[C1qa_id,] < 1 & x[Flt1_id,] < 1 & x[Plpp3_id,] < 1 & x[Cldn11_id,] < 1 & x[Snhg11_id,] < 1 & x[Dlx1_id,] < 1 & x[Npy_id,] < 1 & x[Ecrg4_id,] > 1 })

# classify cells
monocle_cds <- classifyCells(monocle_cds, cth, 0.1)

table(pData(monocle_cds)$CellType)
#         endothelial_cells          ganglion_neuronT1 
#                       528                         40 
#         ganglion_neuronT2           microglial_cells 
#                        41                       1398 
#               neuroblasts            oligodendrocyte 
#                        35                        168 
#oligodendrocyte_precursors              quiescent_NSC 
#                        27                        556 
#                   Unknown 
#                       433 

save.image(file="monocle_cds_v4.RData")

load("monocle_cds_v3.RData")
library(monocle)

png("cell_classification_piechart_monocle.png", width=12, height=8, units="in", res=300)
pie <- ggplot(pData(monocle_cds),
aes(x = factor(1), fill = factor(CellType))) + geom_bar(width = 1)
pie + coord_polar(theta = "y") +
theme(axis.title.x = element_blank(), axis.title.y = element_blank())
dev.off()


# Clustering cells using marker genes


marker_diff <- markerDiffTable(monocle_cds[expressed_genes,], cth,
            residualModelFormulaStr = "~stim + num_genes_expressed", cores = 1)

candidate_clustering_genes <- row.names(subset(marker_diff, qval < 0.01))
length(candidate_clustering_genes) # 7838
marker_spec <-calculateMarkerSpecificity(monocle_cds[candidate_clustering_genes,], cth)
head(selectTopMarkers(marker_spec, 3))
dim(selectTopMarkers(marker_spec, 3)) # 103   3
dim(selectTopMarkers(marker_spec)) # 145   3
selectTopMarkers(marker_spec) # Top Markers for each cell type

# The "specificity" score is calculated using the metric described in Cabili et al and can range from zero to one. 
# The closer it is to one, the more restricted it is to the cell type in question.
clustering_markers <- selectTopMarkers(marker_spec)
write.table(clustering_markers, file="clustering_markers_monocle.txt", sep="\t", row.names=F, quote=F)

# To cluster the cells, we'll choose the top 500 markers for each of these cell types:
dim(selectTopMarkers(marker_spec, 500)) # 4000    3
semisup_clustering_genes <- unique(selectTopMarkers(marker_spec, 500)$gene_id)
monocle_cds <- setOrderingFilter(monocle_cds, semisup_clustering_genes)

png("gene_variability_monocle.png", width=12, height=8, units="in", res=300)
plot_ordering_genes(monocle_cds)
dev.off()
png("pca_variance_explained_monocle.png", width=12, height=8, units="in", res=300)
plot_pc_variance_explained(monocle_cds, return_all = F)
dev.off()

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, num_dim=3,
  norm_method = 'log', reduction_method = 'tSNE', perplexity = 21,
  residualModelFormulaStr = "~stim + num_genes_expressed", verbose = T)

monocle_cds <- clusterCells(monocle_cds, num_clusters = 8)
# Distance cutoff calculated to 4.53965

png("pca_clusters_monocle.png", width=12, height=8, units="in", res=300)
plot_cell_clusters(monocle_cds, 1, 2, color = "CellType")
dev.off()


# Trajectory step 1: Feature selection = choose genes that define a cell's progress

diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,], fullModelFormulaStr = "~stim")
dim(diff_test_res) # 15965     7
ordering_genes <- row.names(subset(diff_test_res, qval < 0.01))
length(ordering_genes) # 1501

monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
png("empirical_dispersion_mean_mRNA_monocle.png", width=12, height=8, units="in", res=300)
plot_ordering_genes(monocle_cds)
dev.off()


# Trajectory step 2: reduce data dimensionality 

monocle_cds <- reduceDimension(monocle_cds, max_components = 2, method = 'DDRTree')

# Trajectory step 3: order cells along the trajectory

monocle_cds <- orderCells(monocle_cds)

# Once the cells are ordered, we can visualize the trajectory in the reduced dimensional space. 
png("cell_trajectory_by_stimulation_monocle.png", width=12, height=8, units="in", res=300)
plot_cell_trajectory(monocle_cds, color_by = "stim")
dev.off()

# colnames(pData(monocle_cds))
png("cell_trajectory_by_CT_monocle.png", width=12, height=8, units="in", res=300)
plot_cell_trajectory(monocle_cds, color_by = "celltype")
dev.off()


save.image(file="monocle_cds_v5.RData")


GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$celltype, pData(cds)$stim)[,"TG.F"]
    return(as.numeric(which(T0_counts == max(T0_counts))))
  } else {
    return (1)
  }
}

# check root index for each cell type
cth
# find what is the GM state
# GM_state(monocle_cds)
# [1] "microglial cells"
# pass the index of GM_state from "cth" to the orderCells function
# monocle_cds <- orderCells(monocle_cds, root_state = "1")

# pass the index of GM_state from "cth" to the orderCells function
monocle_cds <- orderCells(monocle_cds, root_state = GM_state(monocle_cds))

png("cell_trajectory_by_Pseudotime.png", width=12, height=8, units="in", res=300)
plot_cell_trajectory(monocle_cds, color_by = "Pseudotime")
dev.off()

png("cell_trajectory_CTspecific.png", width=14, height=8, units="in", res=300)
plot_cell_trajectory(monocle_cds, color_by = "celltype") + facet_wrap(~celltype, nrow = 1)
dev.off()

save.image(file="monocle_cds_v6.RData")


###########################################################################################################################################################

# Make a table of differential expression result for comparison between all the cell types

all_CTs_names <- unique(scNorm$celltype.stim)
all_CTs_names <- unique(sapply(strsplit(as.character(all_CTs_names), "\\."), "[", 1))

ctrl_m <- ".CTRL.M"
ctrl_f <- ".CTRL.F"
tg_m <- ".TG.M"
tg_f <- ".TG.F"

table = data.frame()

for (i in 1:length(all_CTs_names)){
 print(i)
 name_f <- paste0(all_CTs_names[i], "_F", sep="")
 female <- FindMarkers(scNorm, ident.1 = paste0(all_CTs_names[i], tg_f, sep=""), ident.2 = paste0(all_CTs_names[i], ctrl_f, sep=""), verbose = FALSE)
 total_genes_f <- nrow(female)
 significant_genes_f <- nrow(female[female$p_val_adj < 0.05, ])
 nominal_significant_genes_f <- nrow(female[female$p_val < 0.05, ])

 female_vec <- cbind(name_f, total_genes_f, significant_genes_f, nominal_significant_genes_f)

 name_m <- paste0(all_CTs_names[i], "_M", sep="")
 male <- FindMarkers(scNorm, ident.1 = paste0(all_CTs_names[i], tg_m, sep=""), ident.2 = paste0(all_CTs_names[i], ctrl_m, sep=""), verbose = FALSE)
 total_genes_m <- nrow(male)
 significant_genes_m <- nrow(male[male$p_val_adj < 0.05, ])
 nominal_significant_genes_m <- nrow(male[male$p_val < 0.05, ])

 male_vec <- cbind(name_m, total_genes_m, significant_genes_m, nominal_significant_genes_m)

 complete_df <- rbind(female_vec, male_vec)

 ifelse(i == 1, table <- complete_df, table <- rbind(table, complete_df))
}
dim(table) # 16   4
colnames(table) <- c("cell_type", "total_DEGs", "p_val_adj", "p_val")
write.table(table, file="CT_specific_DE_Result.txt", sep="\t", row.names=F, quote=F)
 
####################################################################################################
