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
library(mclust)
library(ggplot2)


# directory to save plots
dir_plots <- here(file.path("plots", "downstream_clustering"))


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


# --------------
# match clusters
# --------------

# match clusters to ground truth layers for each method

match_nnSVG <- c(5, 8, 1, 7, 2, 6, 3, 4)
colData(spe_out$spe_nnSVG)$label <- factor(
  colData(spe_out$spe_nnSVG)$label, levels = match_nnSVG)

match_SPARKX <- c(1, 6, 2, 8, 3, 4, 7, 5)
colData(spe_out$spe_SPARKX)$label <- factor(
  colData(spe_out$spe_SPARKX)$label, levels = match_SPARKX)

match_HVGs <- c(1, 8, 2, 4, 3, 6, 7, 5)
colData(spe_out$spe_HVGs)$label <- factor(
  colData(spe_out$spe_HVGs)$label, levels = match_HVGs)

match_MoransI <- c(5, 2, 8, 3, 4, 7, 1, 6)
colData(spe_out$spe_MoransI)$label <- factor(
  colData(spe_out$spe_MoransI)$label, levels = match_MoransI)


# ----------
# spot plots
# ----------

# plot clustering

# LIBD layer colors palette
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")

for (i in seq_along(spe_out)) {
  
  df <- as.data.frame(cbind(colData(spe_out[[i]]), spatialCoords(spe_out[[i]])))
  
  ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                 color = label)) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human DLPFC: ", gsub("spe_", "", names(spe_out[i])))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, 
                  paste0("clustering_humanDLPFC_", gsub("spe_", "", names(spe_out[i]))))
  ggsave(paste0(fn, ".pdf"), width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}


# -------------------
# adjusted Rand index
# -------------------

# calculate adjusted Rand index to evaluate clustering performance


# re-label factor levels combine clusters 7 and 8 into a single cluster
# representing white matter

colData(spe_out$spe_nnSVG)$label <- factor(
  colData(spe_out$spe_nnSVG)$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

colData(spe_out$spe_SPARKX)$label <- factor(
  colData(spe_out$spe_SPARKX)$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

colData(spe_out$spe_HVGs)$label <- factor(
  colData(spe_out$spe_HVGs)$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

colData(spe_out$spe_MoransI)$label <- factor(
  colData(spe_out$spe_MoransI)$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))


# calculate adjusted Rand index

ari_nnSVG <- adjustedRandIndex(colData(spe_out$spe_nnSVG)$label, 
                               colData(spe_out$spe_nnSVG)$ground_truth)

ari_SPARKX <- adjustedRandIndex(colData(spe_out$spe_SPARKX)$label, 
                                colData(spe_out$spe_SPARKX)$ground_truth)

ari_HVGs <- adjustedRandIndex(colData(spe_out$spe_HVGs)$label, 
                              colData(spe_out$spe_HVGs)$ground_truth)

ari_MoransI <- adjustedRandIndex(colData(spe_out$spe_MoransI)$label, 
                                 colData(spe_out$spe_MoransI)$ground_truth)


# plot adjusted Rand index

df <- data.frame(
  method = c("nnSVG", "SPARKX", "HVGs", "MoransI"), 
  ARI = c(ari_nnSVG, ari_SPARKX, ari_HVGs, ari_MoransI)
)
df$method <- factor(df$method, levels = df$method)


pal_methods <- c("blue3", "deepskyblue2", "darkorange", "firebrick3")


ggplot(df, aes(x = method, y = ARI, shape = method, color = method)) + 
  geom_point(stroke = 1.5, size = 2) + 
  scale_shape_manual(values = c(4, 3, 1, 2)) + 
  scale_color_manual(values = pal_methods) + 
  ylim(c(0, 1)) + 
  ggtitle("Downstream clustering performance") + 
  labs(y = "Adjusted Rand Index") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())

fn <- file.path(dir_plots, "summary_clustering_performance_ARI")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4.75, height = 3.5)

