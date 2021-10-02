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
library(ggplot2)


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

res_mOB <- as.data.frame(cbind(res_mOB_nnSVG, res_mOB_nnSVG_clusters, res_mOB_HVGs))


# ------------
# plot results
# ------------

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

