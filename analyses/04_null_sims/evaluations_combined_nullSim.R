###################################################
# Script to calculate evaluations: null simulations
# Lukas Weber, May 2022
###################################################

# data set: combined
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggnewscale)


# directory to save plots
dir_plots <- here(file.path("plots", "null_sims"))


# ------------
# load results
# ------------

res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "null_sims", "spe_humanDLPFC_nnSVG_nullSim.rds"))), 
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "null_sims", "spe_mouseOB_nnSVG_nullSim.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["mouseOB_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-1]


# note filtering
dim(res_list$humanDLPFC_nnSVG)
dim(res_list$mouseOB_nnSVG)


# -------------
# error control
# -------------

# calculate false positives

fp_prop_true <- c(0.05, 0.05)
names(fp_prop_true) <- c("humanDLPFC", "mouseOB")

fp_prop <- rep(NA, length(fp_prop_true))
names(fp_prop) <- names(fp_prop_true)

for (i in seq_along(fp_prop)) {
  fp_prop[i] <- mean(res_list[[i]]$pval_nnSVG <= fp_prop_true[i])
}

df_fpr <- 
  data.frame(
    dataset = names(fp_prop_true), 
    cutoff = unname(fp_prop_true), 
    proportion = unname(fp_prop)) %>% 
  mutate(dataset = factor(dataset))


ggplot(df_fpr, aes(x = dataset, y = proportion)) + 
  geom_point(pch = 4, size = 2, stroke = 1.5, color = "red") + 
  geom_hline(yintercept = 0.05, linetype = "dashed") + 
  ylim(c(0, 0.05)) + 
  labs(y = "proportion false SVGs") + 
  ggtitle("Error control") + 
  theme_bw() + 
  theme(legend.title = element_blank(), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "error_control_nnSVG_combined_nullSim")
ggsave(paste0(fn, ".pdf"), width = 3.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 3.5, height = 3.5)

