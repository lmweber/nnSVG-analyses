#############################################
# Script for downstream clustering comparison
# Lukas Weber, Jan 2023
#############################################

# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(scater)
library(scran)
library(ggplot2)


# directory to save plots
dir_plots <- here(file.path("plots", "clustering_comparison"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
spe_list <- list(
  humanDLPFC_nnSVG = readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds")), 
  humanDLPFC_SPARKX = readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX.rds")), 
  humanDLPFC_HVGs = readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds")), 
  humanDLPFC_MoransI = readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds"))
)

res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds"))), 
  humanDLPFC_MoransI = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["humanDLPFC_MoransI"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_MoransI"]]), "_MoransI")[-(1:2)]


# note filtering per method: same for both nnSVG and SPARK-X

table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_nnSVG$gene_id == res_list$humanDLPFC_SPARKX$gene_id)

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
# note: selected k and random seeds to get equal number of clusters per method
g <- buildSNNGraph(spe, k = 5, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object for plotting
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
# note: selected k and random seeds to get equal number of clusters per method
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object for plotting
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

# compute PCA
set.seed(12345)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
# note: selected k and random seeds to get equal number of clusters per method
g <- buildSNNGraph(spe, k = 4, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object for plotting
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

# compute PCA
set.seed(1)
spe <- runPCA(spe, subset_row = top)
# compute UMAP on top 50 PCs
set.seed(123)
spe <- runUMAP(spe, dimred = "PCA")
# update column names
colnames(reducedDim(spe, "UMAP")) <- paste0("UMAP", 1:2)

# clustering

# graph-based clustering
set.seed(123)
# note: selected k and random seeds to get equal number of clusters per method
g <- buildSNNGraph(spe, k = 10, use.dimred = "PCA")
g_walk <- igraph::cluster_walktrap(g)
clus <- g_walk$membership
colLabels(spe) <- factor(clus)


# store object for plotting
spe_out$spe_MoransI <- spe


# ---------------------
# evaluations and plots
# ---------------------


### to do: combined plot showing spot plots
### to do: calculate adjusted Rand index and plot summary

pal <- unname(palette.colors(12, "Polychrome 36"))

plotSpots(spe_out$spe_nnSVG, annotate = "label", palette = pal)
plotSpots(spe_out$spe_HVGs, annotate = "label", palette = pal)
plotSpots(spe_out$spe_SPARKX, annotate = "label", palette = pal)
plotSpots(spe_out$spe_MoransI, annotate = "label", palette = pal)

