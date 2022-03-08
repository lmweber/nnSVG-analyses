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
  stopifnot(all(list_DLPFC[[i]]$gene_name == list_DLPFC[[1]])$gene_name)
}

# method and dataset names

names_DLPFC <- gsub("^spe_", "", names(list_DLPFC))
names_DLPFC

names(list_DLPFC) <- names_DLPFC


# ----------------------------------
# extract known SVGs in this dataset
# ----------------------------------

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

# known SVGs in this dataset
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
ix_known <- match(known_genes, res_DLPFC[, "gene_name_DLPFC_nnSVG"])
ix_known

res_known <- res_DLPFC[ix_known, ]


# ---------------
# plot statistics
# ---------------

# plot comparison of ranks for pairs of methods

ggplot(res_known, aes(x = rank_DLPFC_nnSVG, 
                      y = rank_DLPFC_SPARKX, 
                      label = gene_name_DLPFC_HVGs)) + 
  geom_point(color = "purple3") + 
  geom_text_repel(position = "dodge", color = "purple3") + 
  ggtitle(paste0("Known genes: comparison of ranks")) + 
  theme_bw()

fn <- here("plots", "known", paste0("known_genes_ranks_DLPFC_nnSVG_SPARKX"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)


# plot comparison of statistics / significance for pairs of methods

ggplot(res_known, aes(x = LR_stat_DLPFC_nnSVG, 
                      y = -log10(combinedPval_DLPFC_SPARKX), 
                      label = gene_name_DLPFC_HVGs)) + 
  geom_point(color = "purple3") + 
  geom_text_repel(position = "dodge", color = "purple3") + 
  ggtitle(paste0("Known genes: comparison of statistics / significance")) + 
  theme_bw()

fn <- here("plots", "known", paste0("known_genes_stats_DLPFC_nnSVG_SPARKX"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)


# plot rank vs effect size for nnSVG

ggplot(res_known, aes(x = rank_DLPFC_nnSVG, 
                      y = prop_sv_DLPFC_nnSVG, 
                      label = gene_name_DLPFC_HVGs)) + 
  geom_point(color = "blue") + 
  geom_text_repel(position = "dodge", color = "blue") + 
  ggtitle(paste0("Known genes: rank vs. effect size")) + 
  theme_bw()

fn <- here("plots", "known", paste0("rank_vs_effect_size_DLPFC_nnSVG"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)

