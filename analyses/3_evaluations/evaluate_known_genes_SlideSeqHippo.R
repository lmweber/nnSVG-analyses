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


# ----------------------------------
# extract known SVGs in this dataset
# ----------------------------------

# add method names to columns
for (i in seq_along(list_SlideSeqHippo)) {
  colnames(list_SlideSeqHippo[[i]])[-1] <- paste0(colnames(list_SlideSeqHippo[[i]])[-1], "_", names_SlideSeqHippo[i])
}

res_SlideSeqHippo <- as.data.frame(cbind(
  list_SlideSeqHippo[["SlideSeqHippo_nnSVG_clusters"]], 
  list_SlideSeqHippo[["SlideSeqHippo_SPARKX_clusters"]][, -1], 
  list_SlideSeqHippo[["SlideSeqHippo_deviance_clusters"]][, -1]
))

# known SVGs in this dataset
known_genes <- c("Rgs14", "Cpne9")
ix_known <- match(known_genes, res_SlideSeqHippo[, "gene_name"])
ix_known

res_known <- res_SlideSeqHippo[ix_known, ]


# ---------------
# plot statistics
# ---------------

# plot comparison of ranks for pairs of methods

ggplot(res_known, aes(x = rank_SlideSeqHippo_nnSVG_clusters, 
                      y = rank_SlideSeqHippo_SPARKX_clusters, 
                      label = gene_name)) + 
  geom_point(color = "purple3") + 
  geom_text_repel(position = "dodge", color = "purple3") + 
  ggtitle(paste0("Known genes: comparison of ranks")) + 
  theme_bw()

fn <- here("plots", "known", paste0("known_genes_ranks_SlideSeqHippo_nnSVG_SPARKX"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)


# plot comparison of statistics / significance for pairs of methods

ggplot(res_known, aes(x = LR_stat_SlideSeqHippo_nnSVG_clusters, 
                      y = -log10(combinedPval_SlideSeqHippo_SPARKX_clusters), 
                      label = gene_name)) + 
  geom_point(color = "purple3") + 
  geom_text_repel(position = "dodge", color = "purple3") + 
  ggtitle(paste0("Known genes: comparison of statistics / significance")) + 
  theme_bw()

fn <- here("plots", "known", paste0("known_genes_stats_SlideSeqHippo_nnSVG_SPARKX"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)


# plot rank vs effect size for nnSVG

ggplot(res_known, aes(x = rank_SlideSeqHippo_nnSVG_clusters, 
                      y = prop_sv_SlideSeqHippo_nnSVG_clusters, 
                      label = gene_name)) + 
  geom_point(color = "blue") + 
  geom_text_repel(position = "dodge", color = "blue") + 
  ggtitle(paste0("Known genes: rank vs. effect size")) + 
  theme_bw()

fn <- here("plots", "known", paste0("rank_vs_effect_size_SlideSeqHippo_nnSVG"))
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)

