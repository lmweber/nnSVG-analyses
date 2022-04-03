###################################################
# Script to calculate evaluations: null simulations
# Lukas Weber, Apr 2022
###################################################

# data set: mouse OB
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "null_sims"))


# ------------
# load results
# ------------

# note choice of filtering per method
res_list <- list(
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "null_sims", "spe_mouseOB_nnSVG_nullSim.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseOB_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-(1:2)]


# note filtering
dim(res_list$mouseOB_nnSVG)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseOB_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  labs(x = "p-values", 
       y = "frequency", 
       title = "nnSVG p-values: mouse OB", 
       subtitle = "null simulation") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseOB_nullSim")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4, height = 3.5)

