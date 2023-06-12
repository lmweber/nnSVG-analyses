#############################################
# Script for downstream clustering comparison
# Lukas Weber, updated Mar 2023
#############################################

# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(scater)
library(scran)


dir_clustering <- here("outputs", "downstream_clustering")


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
spe_list <- list(
  humanDLPFC_nnSVG = readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds")), 
  humanDLPFC_SPARKX = readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX_noFilt.rds")), 
  humanDLPFC_HVGs = readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds")), 
  humanDLPFC_MoransI = readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds"))
)

res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX_noFilt.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds"))), 
  humanDLPFC_MoransI = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["humanDLPFC_MoransI"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_MoransI"]]), "_MoransI")[-(1:2)]


# note filtering per method

table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)

table(res_list$humanDLPFC_MoransI$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_MoransI$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)


spe_out <- list()


# ----------------------------
# downstream clustering: nnSVG
# ----------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# nnSVG: top 1000 SVGs
ix <- which(res_list$humanDLPFC_nnSVG$rank_nnSVG <= 1000)
top <- res_list$humanDLPFC_nnSVG$gene_id[ix]

spe <- spe_list$humanDLPFC_nnSVG

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_nnSVG <- spe


# ------------------------------
# downstream clustering: SPARK-X
# ------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# SPARK-X: top 1000 SVGs
ix <- which(res_list$humanDLPFC_SPARKX$rank_SPARKX <= 1000)
top <- res_list$humanDLPFC_SPARKX$gene_id[ix]

spe <- spe_list$humanDLPFC_SPARKX

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(123456)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_SPARKX <- spe


# ---------------------------
# downstream clustering: HVGs
# ---------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# top 1000 HVGs
ix <- which(res_list$humanDLPFC_HVGs$rank_HVGs <= 1000)
top <- res_list$humanDLPFC_HVGs$gene_id[ix]

spe <- spe_list$humanDLPFC_HVGs

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(2)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_HVGs <- spe


# --------------------------------
# downstream clustering: Moran's I
# --------------------------------

# calculate downstream clustering on top 1000 SVGs or HVGs

# Moran's I: top 1000 SVGs
ix <- which(res_list$humanDLPFC_MoransI$rank_MoransI <= 1000)
top <- res_list$humanDLPFC_MoransI$gene_id[ix]

spe <- spe_list$humanDLPFC_MoransI

# dimensionality reduction

# note: selected random seeds to get equal number of clusters per method

# compute PCA
set.seed(1234567)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object
spe_out$spe_MoransI <- spe


# ------------
# save results
# ------------

# save colData and spatialCoords instead of full SPE objects to save space

coldata_out <- spatialcoords_out <- list()

coldata_out$nnSVG <- colData(spe_out$spe_nnSVG)
coldata_out$SPARKX <- colData(spe_out$spe_SPARKX)
coldata_out$HVGs <- colData(spe_out$spe_HVGs)
coldata_out$MoransI <- colData(spe_out$spe_MoransI)

spatialcoords_out$nnSVG <- spatialCoords(spe_out$spe_nnSVG)
spatialcoords_out$SPARKX <- spatialCoords(spe_out$spe_SPARKX)
spatialcoords_out$HVGs <- spatialCoords(spe_out$spe_HVGs)
spatialcoords_out$MoransI <- spatialCoords(spe_out$spe_MoransI)

res_out <- list(
  coldata_out = coldata_out, 
  spatialcoords_out = spatialcoords_out
)


fn <- here(dir_clustering, "res_downstream_clustering.rds")
saveRDS(res_out, file = fn)

