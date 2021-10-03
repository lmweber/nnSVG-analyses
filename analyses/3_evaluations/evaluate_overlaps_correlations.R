########################################
# Script to calculate evaluation metrics
# Lukas Weber, Oct 2021
########################################

# to run in interactive session on cluster:
# module load conda_R/4.1.x
# Rscript filename.R


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ochRe)


# ------------
# load results
# ------------

# mOB dataset

spe_mOB_nnSVG <- readRDS(here("outputs", "results", "spe_mOB_nnSVG.rds"))
spe_mOB_nnSVG_clusters <- readRDS(here("outputs", "results", "spe_mOB_nnSVG_clusters.rds"))
spe_mOB_HVGs <- readRDS(here("outputs", "results", "spe_mOB_HVGs.rds"))

stopifnot(all(rowData(spe_mOB_HVGs)$gene_name == rowData(spe_mOB_nnSVG)$gene_name))
stopifnot(all(rowData(spe_mOB_HVGs)$gene_name == rowData(spe_mOB_nnSVG_clusters)$gene_name))

# create combined data frame of results

res_mOB_nnSVG <- 
  rowData(spe_mOB_nnSVG)[, c("gene_name", "LR_stat", "rank", "pval", "padj")]
colnames(res_mOB_nnSVG) <- 
  c("gene_name", "stat_nnSVG", "rank_nnSVG", "pval_nnSVG", "padj_nnSVG")

res_mOB_nnSVG_clusters <- 
  rowData(spe_mOB_nnSVG_clusters)[, c("LR_stat", "rank", "pval", "padj")]
colnames(res_mOB_nnSVG_clusters) <- 
  c("stat_nnSVG_clusters", "rank_nnSVG_clusters", "pval_nnSVG_clusters", "padj_nnSVG_clusters")

res_mOB_HVGs <- 
  rowData(spe_mOB_HVGs)[, c("bio", "rank", "p.value", "FDR")]
colnames(res_mOB_HVGs) <- 
  c("stat_HVGs", "rank_HVGs", "pval_HVGs", "padj_HVGs")

res_mOB_nnSVG <- as.data.frame(res_mOB_nnSVG)
res_mOB_nnSVG_clusters <- as.data.frame(res_mOB_nnSVG_clusters)
res_mOB_HVGs <- as.data.frame(res_mOB_HVGs)

res_mOB <- cbind(res_mOB_nnSVG, res_mOB_nnSVG_clusters, res_mOB_HVGs)


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of top n genes from method 1 (e.g. HVGs) 
# that are also in the set of top n genes from method 2 (e.g. nnSVG)

overlaps <- c(10, 20, 50, 100, 200, 500, 1000)
top_nnSVG <- rep(NA, length(top_HVGs))
top_nnSVG_clusters <- rep(NA, length(top_HVGs))

for (k in seq_along(overlaps)) {
  # select top gene names
  genes_k <- rownames(filter(res_mOB_HVGs, rank_HVGs <= overlaps[k]))
  
  # calculate overlaps
  top_nnSVG[k] <- nrow(filter(res_mOB_nnSVG[genes_k, ], 
                              rank_nnSVG <= overlaps[k]))
  top_nnSVG_clusters[k] <- nrow(filter(res_mOB_nnSVG_clusters[genes_k, ], 
                                       rank_nnSVG_clusters <= overlaps[k]))
}

# calculate proportions
df_overlaps <- data.frame(
  top_HVGs = overlaps, 
  nnSVG = top_nnSVG / overlaps, 
  nnSVG_clusters = top_nnSVG_clusters / overlaps
)

df_overlaps <- pivot_longer(
  df_overlaps, 
  cols = c("nnSVG", "nnSVG_clusters"), 
  names_to = "method", 
  values_to = "prop_overlap"
)

# plot
ggplot(df_overlaps, aes(x = top_HVGs, y = prop_overlap, 
                        group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_ochre(palette = "parliament") + 
  scale_x_continuous(breaks = top_HVGs, trans = "log10") + 
  labs(x = "top HVGs", 
       y = "proportion overlapping") + 
  ggtitle("mOB: Overlap with top HVGs") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank())

fn <- here("plots", "overlaps", "overlaps_mOB.pdf")
ggsave(fn, width = 6, height = 4)


# ----------------------
# calculate correlations
# ----------------------

ggplot(res_mOB, aes(x = rank_HVGs, y = rank_nnSVG)) + 
  geom_point() + 
  xlim(c(0, 200)) + 
  ylim(c(0, 200)) + 
  coord_fixed() + 
  theme_bw()

ggplot(res_mOB, aes(x = rank_nnSVG, y = rank_nnSVG_clusters)) + 
  geom_point() + 
  xlim(c(0, 200)) + 
  ylim(c(0, 200)) + 
  coord_fixed() + 
  theme_bw()

