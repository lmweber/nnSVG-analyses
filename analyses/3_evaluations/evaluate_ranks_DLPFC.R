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

# DLPFC dataset

list_DLPFC <- list(
  spe_DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  #spe_DLPFC_nnSVG_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG_clusters.rds"))), 
  
  spe_DLPFC_nnSVG_logcounts = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG_logcounts.rds"))), 
  #spe_DLPFC_nnSVG_logcounts_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG_logcounts_clusters.rds"))), 
  
  spe_DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARK-X", "spe_DLPFC_SPARK-X.rds"))), 
  
  spe_DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds"))), 
  spe_DLPFC_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_DLPFC_deviance.rds")))#, 
  #spe_DLPFC_deviance_clusters = rowData(readRDS(here("outputs", "results", "deviance", "spe_DLPFC_deviance_clusters.rds")))
)

for (i in seq_along(list_DLPFC)){
  stopifnot(all(rowData(list_DLPFC[[i]])$gene_name == rowData(list_DLPFC[[1]])$gene_name))
}

# method and dataset names

names_DLPFC <- gsub("^spe_", "", names(list_DLPFC))
names_DLPFC

names(list_DLPFC) <- names_DLPFC


# -------------
# extract ranks
# -------------

# add method names to columns
for (i in seq_along(list_DLPFC)) {
  colnames(list_DLPFC[[i]])[-1] <- paste0(colnames(list_DLPFC[[i]])[-1], "_", names_DLPFC[i])
}

res_DLPFC <- as.data.frame(cbind(
  list_DLPFC[["DLPFC_nnSVG"]], 
  list_DLPFC[["DLPFC_nnSVG_logcounts"]][, -1], 
  list_DLPFC[["DLPFC_SPARKX"]][, -1], 
  list_DLPFC[["DLPFC_HVGs"]][, -1], 
  list_DLPFC[["DLPFC_deviance"]][, -1]
))


# ----------
# plot ranks
# ----------

# plot ranks of top n genes from each pair of methods

ranks <- c(20, 100, 1000)


for (k in seq_along(ranks)) {
  p <- ggplot(res_DLPFC, aes(x = rank_DLPFC_HVGs, y = rank_DLPFC_nnSVG, label = gene_name_DLPFC_HVGs)) + 
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
  
  fn <- here("plots", "ranks", paste0("ranks_DLPFC_HVGs_nnSVG_", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}


for (k in seq_along(ranks)) {
  p <- ggplot(res_DLPFC, aes(x = rank_DLPFC_HVGs, y = rank_DLPFC_SPARKX, label = gene_name_DLPFC_HVGs)) + 
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
  
  fn <- here("plots", "ranks", paste0("ranks_DLPFC_HVGs_SPARKX_", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}


for (k in seq_along(ranks)) {
  p <- ggplot(res_DLPFC, aes(x = rank_DLPFC_nnSVG, y = rank_DLPFC_SPARKX, label = gene_name_DLPFC_HVGs)) + 
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
  
  fn <- here("plots", "ranks", paste0("ranks_DLPFC_nnSVG_SPARKX_", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}



