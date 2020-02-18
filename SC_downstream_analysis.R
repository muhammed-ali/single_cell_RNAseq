############################################################################
############################################################################
############################# Seurat Tutorial ##################################
############################################################################
############################################################################

# guide: https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html


################################################
########### Setup the Seurat Object ################
################################################


library(Seurat)
library(tximport)
library(dplyr)
files <- file.path("alevin_quants/alevin/quants_mat.gz")
file.exists(files) # TRUE

# read gene name mapping file
gene_mapping <- read.table("/home/muhammad/Downloads/Projects/CS_mm_PD_SNCA/EnsT_EnsG_GeneName.tsv", sep="\t", header=T, stringsAsFactors=F)
gene_mapping <- gene_mapping[,c(2,3)]

# mapping table with "key-value" pairs
lookUp1 <- setNames(as.character(gene_mapping$gene_name), gene_mapping$gene_id)

# import count matrix
txi <- tximport(files, type="alevin")

# Convert Ensembl IDs to gene symbol in the row.names of count matrix:
res <- lapply(row.names(txi$counts), function(i) lookUp1[i])
gene_symbol <- unlist(res, use.names=FALSE)
row.names(txi$counts) <- gene_symbol

# Create Seurat object
pbmc <- CreateSeuratObject(counts = txi$counts , min.cells = 3, min.features = 200, project = "dropSeq_Mm_cortex")
pbmc
# An object of class Seurat 
# 23087 features across 4188 samples within 1 assay 
# Active assay: RNA (23087 features)


#####################################
########### Filter cells ################
#####################################


# Aim: Filter cells with unique feature count > 2500, or less than 200, or the cells that have more than 5% mitochondrial counts

# We use the set of all genes starting with mt- as a set of mitochondrial genes
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^mt-")

# Visualize QC metrics as a violin plot
png("QC_Metrics_mm_Cortex.png", width=8, height=6, units="in", res=300)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# FeatureScatter is typically used to visualize feature-feature relationships
png("FeatureScatter_mm_Cortex.png", width=12, height=6, units="in", res=300)
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
dev.off()

# filter cells that have unique feature counts over 2,500 or less than 200
# filter cells that have >5% mitochondrial counts
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)


#######################
## Normalizing the data ###
#######################


# Aim: Employ a global-scaling normalization method “LogNormalize” that normalizes the feature expression measurements 
# for each cell by the total expression, multiplies this by a scale factor (10,000 by default), and log-transforms the result

# Normalized values are stored in pbmc[["RNA"]]@data
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)


#######################################################################
########### Identification of highly variable features (feature selection) ###########
#######################################################################


# Aim: Calculate a subset of features that exhibit high cell-to-cell variation in the dataset 
# (i.e, they are highly expressed in some cells, and lowly expressed in others)

# return 2,000 features per dataset. These will be used in downstream analysis, like PCA.
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(pbmc), 10)
top10
# [1] "Plp1"      "Ttr"       "Enpp2"     "Ptgds"     "Il1b"      "Hist1h2ap"
# [7] "Ermn"      "Apod"      "Ly6c1"     "Bsg" 

# plot variable features with and without labels
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
png("VariableFeatures_mm_Cortex.png", width=14, height=6, units="in", res=300)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#######################
#### Scaling the data ####
#######################


# Aim: Apply a linear transformation (‘scaling’) that is a standard pre-processing step prior to dimensional reduction techniques like PCA.
# Shifts the expression of each gene, so that the mean expression across cells is 0.
# Scales the expression of each gene, so that the variance across cells is 1.

# The results of this are stored in pbmc[["RNA"]]@scale.data
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc, features = all.genes)


#####################################################
########### Perform linear dimensional reduction ###########
#####################################################


# Aim: perform PCA on the scaled data for reducing the dimensions of data.
# DimHeatmap: Both cells and features are ordered according to their PCA scores. 

# perform linear dimensional reduction
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))

# Examine and visualize PCA results a few different ways
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)

# PC_ 1 
# Positive:  Ptn, Dbi, Gpm6b, Ttyh1, Gpm6a 
# Negative:  Nfkbia, Egr1, Atf3, Klf6, Gm43305 
# PC_ 2 
# Positive:  Stmn4, Syt11, Elavl3, Sox11, Fez1 
# Negative:  Cldn5, Igfbp7, Ly6c1, Itm2a, Adgrf5 
# PC_ 3 
# Positive:  Gja1, Id4, Pla2g7, Aldoc, Clu 
# Negative:  Cldn11, Sept4, Mal, Ermn, Ugt8a 
# PC_ 4 
# Positive:  Gja1, Htra1, Car2, Aldoc, Clu 
# Negative:  Dlx1, Sox11, Pbk, Top2a, Mki67 
# PC_ 5 
# Positive:  Calml4, Ecrg4, 2900040C04Rik, Clic6, Rsph1 
# Negative:  Sparcl1, Ptn, S1pr1, Cldn5, Plpp3 

# Seurat provides several useful ways of visualizing both cells and features that define the PCA, including VizDimReduction, DimPlot , and DimHeatmap

png("VizDimLoadings_mm_Cortex.png", width=14, height=6, units="in", res=300)
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
dev.off()

png("DimPlot_mm_Cortex.png", width=8, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "pca")
dev.off()

png("DimHeatmap_PC1_500cells_mm_Cortex.png", width=8, height=6, units="in", res=300)
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
dev.off()

png("DimHeatmap_PC15_500cells_mm_Cortex.png", width=12, height=12, units="in", res=300)
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)
dev.off()


##############################################
#### Determine the ‘dimensionality’ of the dataset ####
##############################################


# Aim: To overcome the extensive technical noise in any single feature for scRNA-seq data, 
# Seurat clusters cells based on their PCA scores, with each PC essentially representing a ‘metafeature’ 
# that combines information across a correlated feature set.
# JackStrawPlot: comparing the distribution of p-values for each PC with a uniform distribution (dashed line).
# ElbowPlot: a ranking of principle components based on the percentage of variance explained by each one.

# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
pbmc <- JackStraw(pbmc, num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)

png("JackStrawPlot_PC15_mm_Cortex.png", width=10, height=6, units="in", res=300)
JackStrawPlot(pbmc, dims = 1:15)
dev.off()
# Warning message: Removed 21701 rows containing missing values (geom_point).

# ‘Elbow plot’: a ranking of principle components based on the percentage of variance explained by each one
png("ElbowPlot_PC15_mm_Cortex.png", width=10, height=6, units="in", res=300)
ElbowPlot(pbmc)
dev.off()
# we can observe an ‘elbow’ around PC9-10, suggesting that the majority of true signal is captured in the first 10 PCs


#######################
#### Cluster the cells ####
#######################


# Aim: Seurat applies a graph-based clustering approach.
# Briefly, these methods embed cells in a graph structure - for example a K-nearest neighbor (KNN) graph, 
# with edges drawn between cells with similar feature expression patterns, and then attempt to partition 
# this graph into highly interconnected ‘quasi-cliques’ or ‘communities’.
# First construct a KNN graph based on the euclidean distance in PCA space, and refine the edge weights 
# between any two cells based on the shared overlap in their local neighborhoods (Jaccard similarity)

# WE CHOOSE PC-10 as a CUTOFF (based on Elbow and JackStrawPlot)
pbmc <- FindNeighbors(pbmc, dims = 1:10) # optimal: pbmc <- FindNeighbors(pbmc, dims = 1:9)

pbmc <- FindClusters(pbmc, resolution = 0.5) # optimal: pbmc <- FindClusters(pbmc, resolution = 0.4)

# Number of nodes: 930
# Number of edges: 31188
# Maximum modularity in 10 random starts: 0.7632
# Number of communities: 7

# Look at cluster IDs of the first 5 cells
head(Idents(pbmc), 5)

# CTGTCTCTTATA GAGCCCAGCCCC CCACGCGCGTTG CGAAGATCAGGG TTAACGCAGGCG 
#            1            1            4            3            6 
# Levels: 0 1 2 3 4 5 6

# Look at the number of cells in each cluster
table(Idents(pbmc))

#  microglia cells  microglial cells    Purkinje cells  Purkinje neurons 
#               364               353                62                53 
# endothelial cells  oligodendrocytes   ependymal cells 
#                38                36                24 



##################################################
### Run non-linear dimensional reduction (UMAP/tSNE) ###
##################################################


# Aim: Learn the underlying manifold of the data in order to  place similar cells together in low-dimension space.
# Non-linear dimensional reduction techniques, such as tSNE and UMAP, to visualize and explore these datasets.

# As input to the UMAP and tSNE, we suggest using the same PCs as input to the clustering analysis
pbmc <- RunUMAP(pbmc, dims = 1:10)

# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
png("DimPlot_umap_mm_Cortex.png", width=10, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "umap")
dev.off()

pbmc <- RunTSNE(pbmc, dims = 1:10)
png("DimPlot_tSNE_mm_Cortex.png", width=10, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "tsne")
dev.off()

pbmc <- RunPCA(pbmc, dims = 1:10)
png("DimPlot_pca_mm_Cortex.png", width=10, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "pca")
dev.off()

# save the object at this point so that it can easily be loaded back in without having to rerun the computationally intensive steps performed above
saveRDS(pbmc, file = "pbmc_mm_Cortex.rds")


##########################################################
### Finding diffierentially expressed features (cluster biomarkers) ###
##########################################################


Aim: Find the features (genes) that are uniquely expressed in one cluster (probable cell sub-population) w.r.t. rest of the cell populations (clusters).

# if you want to load the previously stored object
# pbmc <- readRDS("pbmc_mm_Cortex.rds")

# find all markers of cluster 1
cluster1.markers <- FindMarkers(pbmc, ident.1 = 1, min.pct = 0.25)
dim(cluster1.markers) # 286   5
head(cluster1.markers, n = 5)

#               p_val avg_logFC pct.1 pct.2    p_val_adj
# P2ry12 2.413380e-34 0.6592529 0.963 0.631 5.571771e-30
# Hexb   1.484819e-30 0.5431115 0.992 0.712 3.428001e-26
# Cx3cr1 1.876908e-30 0.5655268 0.969 0.669 4.333217e-26
# Ctss   5.979815e-30 0.5417822 0.989 0.688 1.380560e-25
# Ctsd   1.167179e-25 0.5200062 0.955 0.704 2.694665e-21

# find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)
dim(cluster5.markers) # 1001    5
head(cluster5.markers, n = 5)

# find all markers distinguishing cluster 4 from rest of all the clusters (0,1,2,3,5,6)
cluster4.markers <- FindMarkers(pbmc, ident.1 = 4, ident.2 = c(0,1,2,3,5,6), min.pct = 0.25)
dim(cluster4.markers) # 917   5
nrow(cluster5.markers[cluster5.markers$p_val < 0.05,]) # 747
nrow(cluster5.markers[cluster5.markers$p_val_adj < 0.05,]) # 291

# find markers for every cluster compared to all remaining cells, report only the positive ones
pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
dim(pbmc.markers) # 3045    7
table(pbmc.markers$cluster)
#    0     1    2     3    4     5    6 
# 208 195 384 638 447 449 724
nrow(pbmc.markers[pbmc.markers$p_val < 0.05,]) # 3045
nrow(pbmc.markers[pbmc.markers$p_val_adj < 0.05,]) # 1813
write.table(pbmc.markers, file="all_clusters_markers_mm_cortex.txt", sep="\t", quote=F)

# Top 2 genes from each cluster
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)

# ROC test returns the ‘classification power’ for any individual marker (ranging from 0 - random, to 1 - perfect)
# change "ident.1" for the cluster you want and get the results: (example, Cldn5 = cluster 4, Aldoc = cluster 2)
cluster1.markers <- FindMarkers(pbmc, ident.1 = 0, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
dim(cluster1.markers) # 383   5

# visualizing marker expression across different clusters
png("VlnPlot_markerCluster_expression_mm_Cortex.png", width=10, height=6, units="in", res=300)
VlnPlot(pbmc, features = c("Cldn5", "Aldoc"))
dev.off()

# you can plot raw counts as well
png("VlnPlot_markerCluster_rawCounts_mm_Cortex.png", width=10, height=6, units="in", res=300)
VlnPlot(pbmc, features = c("Cldn5", "Aldoc"), slot = "counts", log = TRUE)
dev.off()

png("FeaturePlot_mm_Cortex.png", width=12, height=12, units="in", res=300)
FeaturePlot(pbmc, features = c("Atf3", "Nfkbiz", "P2ry12", "Rnase4", "Aldoc", "Slc1a2", "Meg3", "Hist1h2ap", "Cldn5", "Ly6c1", "Ptgds", "Plp1"))
dev.off()

# DoHeatmap generates an expression heatmap for given cells and features. 
# In this case, we are plotting the top 20 markers (or all markers if less than 20) for each cluster.

png("DoHeatmap_mm_Cortex.png", width=12, height=10, units="in", res=300)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()
dev.off()


######################################################
######### Assigning cell type identity to clusters ############
######################################################


# Aim: Manually go through marker genes and find in which cell they are uniquely expressed and then label the cell types accordingly.

# https://science.sciencemag.org/content/347/6226/1138


# Aif1 = cluster 0 = microglia # https://science.sciencemag.org/content/347/6226/1138
# P2ry12 (also Hexb, Ctss) = cluster 1 = a specific marker for microglial cells # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5223388/ # https://www.cell.com/cell-reports/pdf/S2211-1247(17)31702-3.pdf
# Aldoc = cluster 2 = cerebellar Purkinje cells # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0086679
# Tuba1a = cluster 3 = is expressed in cells positive for neuronal markers # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5648057/
# Rtn1 = cluster 3 = mostly expressed by neurons, marker for dendrites of Purkinje neurons # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5522448/
# cldn5 = cluster 4 = Endothelial Cells in the Cerebral Cortex # https://www.ncbi.nlm.nih.gov/pubmed/27847469 # https://www.nature.com/articles/mp2017156.pdf?origin=ppub
# Mbp = cluster 5 = oligodendrocytes markers such as MBP # https://www.frontiersin.org/articles/10.3389/fnana.2018.00090/full
# Ttr (also Enpp2) = cluster 6 = ependymal cells # For example, Fol1, Igfbp2, Ttr, Enpp2 ... # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3510163/
#  Extra # Egr1, Junb, Fos = cluster 0 = Astroependymal Cells # https://www.sciencedirect.com/science/article/pii/S009286741830789X


new.cluster.ids <- c("microglia cells", "microglial cells", "Purkinje cells", "Purkinje neurons", "endothelial cells", "oligodendrocytes", "ependymal cells")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)

png("DimPlot_labelled_umap_mm_Cortex.png", width=8, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) 
dev.off()

png("DimPlot_labelled_tsne_mm_Cortex.png", width=8, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "tsne", label = TRUE, pt.size = 0.5) 
dev.off()

png("DimPlot_labelled_pca_mm_Cortex.png", width=8, height=6, units="in", res=300)
DimPlot(pbmc, reduction = "pca", label = TRUE, pt.size = 0.5) 
dev.off()

# visualizing marker expression across different clusters
png("VlnPlot_markerCluster_expression_mm_Cortex.png", width=10, height=6, units="in", res=300)
VlnPlot(pbmc, features = c("Cldn5", "Aif1"))
dev.off()

saveRDS(pbmc, "pbmc_mm_Cortex_final.rds")


####################################################################################################
