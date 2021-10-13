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
library(ggrepel)
library(ggsci)


# ------------
# load results
# ------------

# SlideSeqHippo dataset

list_SlideSeqHippo <- list(
  #spe_SlideSeqHippo_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds"))), 
  spe_SlideSeqHippo_nnSVG_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  
  #spe_SlideSeqHippo_nnSVG_logcounts = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_logcounts.rds"))), 
  #spe_SlideSeqHippo_nnSVG_logcounts_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_logcounts_clusters.rds"))), 
  
  #spe_SlideSeqHippo_SPARKX = rowData(readRDS(here("outputs", "results", "SPARK-X", "spe_SlideSeqHippo_SPARK-X.rds"))), 
  spe_SlideSeqHippo_SPARKX_clusters = rowData(readRDS(here("outputs", "results", "SPARK-X", "spe_SlideSeqHippo_SPARK-X_clusters.rds"))), 
  
  #spe_SlideSeqHippo_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds"))), 
  #spe_SlideSeqHippo_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance.rds")))#, 
  spe_SlideSeqHippo_deviance_clusters = rowData(readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance_clusters.rds")))
)

for (i in seq_along(list_SlideSeqHippo)){
  stopifnot(all(list_SlideSeqHippo[[i]]$gene_name == list_SlideSeqHippo[[1]]$gene_name))
}

# method and dataset names

names_SlideSeqHippo <- gsub("^spe_", "", names(list_SlideSeqHippo))
names_SlideSeqHippo

names(list_SlideSeqHippo) <- names_SlideSeqHippo


# -------------
# extract ranks
# -------------

# add method names to columns
for (i in seq_along(list_SlideSeqHippo)) {
  colnames(list_SlideSeqHippo[[i]])[-1] <- paste0(colnames(list_SlideSeqHippo[[i]])[-1], "_", names_SlideSeqHippo[i])
}

res_SlideSeqHippo <- as.data.frame(cbind(
  list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]], 
  list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]][, -1], 
  list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]][, -1]
))


# ----------
# plot ranks
# ----------

# plot ranks of top n genes from each pair of methods

ranks <- c(20, 100, 1000)


for (k in seq_along(ranks)) {
  p <- ggplot(res_SlideSeqHippo, aes(x = rank_SlideSeqHippo_deviance_clusters, 
                                     y = rank_SlideSeqHippo_nnSVG_clusters, 
                                     label = gene_name)) + 
    geom_point() + 
    xlim(c(0, ranks[k])) + 
    ylim(c(0, ranks[k])) + 
    coord_fixed() + 
    ggtitle(paste0("Ranks of top ", ranks[k], " genes per method")) + 
    theme_bw()
  
  if (ranks[k] == 20) {
    p <- p + geom_text_repel(position = "dodge")
  }
  
  print(p)
  
  fn <- here("plots", "ranks", paste0("ranks_SlideSeqHippo_deviance_clusters_nnSVG_clusters_", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}


for (k in seq_along(ranks)) {
  p <- ggplot(res_SlideSeqHippo, aes(x = rank_SlideSeqHippo_deviance_clusters, 
                                     y = rank_SlideSeqHippo_SPARKX_clusters, 
                                     label = gene_name)) + 
    geom_point() + 
    xlim(c(0, ranks[k])) + 
    ylim(c(0, ranks[k])) + 
    coord_fixed() + 
    ggtitle(paste0("Ranks of top ", ranks[k], " genes per method")) + 
    theme_bw()
  
  if (ranks[k] == 20) {
    p <- p + geom_text_repel(position = "dodge")
  }
  
  print(p)
  
  fn <- here("plots", "ranks", paste0("ranks_SlideSeqHippo_deviance_clusters_SPARKX_clusters", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}

