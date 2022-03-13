#################################
# Script to calculate evaluations
# Lukas Weber, Mar 2022
#################################

# data set: human DLPFC
# filtering: no filtering for nnSVG


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "evaluations", "humanDLPFC", "no_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

# note choice of filtering per method
res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG_noFilt.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX_noFilt.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method: no filtering for both nnSVG and SPARK-X

table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_nnSVG$gene_id == res_list$humanDLPFC_SPARKX$gene_id)

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_HVGs$gene_id == res_list$humanDLPFC_nnSVG$gene_id)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$humanDLPFC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

