###################################################
# Script to calculate evaluations: null simulations
# Lukas Weber, May 2022
###################################################

# data set: mouse OB
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

# note choice of filtering per method
res_list <- list(
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "null_sims", "spe_mouseOB_nnSVG_nullSim.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseOB_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-1]


# note filtering
dim(res_list$mouseOB_nnSVG)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseOB_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  ylim(c(0, 2500)) + 
  labs(x = "p-values", 
       y = "frequency", 
       title = "nnSVG p-values: mouse OB", 
       subtitle = "null simulation") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseOB_nullSim")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4, height = 3.5)


# -------------
# error control
# -------------

# calculate false positives

fp_prop_true <- c(0.01, 0.02, 0.05, 0.1, 0.2)

fp_prop <- rep(NA, length(fp_prop_true))
names(fp_prop) <- fp_prop_true

for (i in seq_along(fp_prop)) {
  fp_prop[i] <- mean(df_pvals$pval_nnSVG <= fp_prop_true[i])
}

df_fpr <- data.frame(
  cutoff = as.numeric(names(fp_prop)), 
  proportion = unname(fp_prop), 
  type = "observed"
)
df_true <- data.frame(
  cutoff = as.numeric(names(fp_prop)), 
  proportion = as.numeric(names(fp_prop)), 
  type = "expected"
)

pal <- c("black", "red")

ggplot() + 
  geom_point(data = df_fpr, aes(x = cutoff, y = proportion), 
             color = "red", size = 2.5) + 
  geom_line(data = df_fpr, aes(x = cutoff, y = proportion, color = type), 
            linetype = "solid") + 
  geom_line(data = df_true, aes(x = cutoff, y = proportion, color = type), 
            linetype = "dashed") + 
  scale_color_manual(values = pal) + 
  coord_fixed() + 
  xlim(c(0, 0.2)) + 
  ylim(c(0, 0.2)) + 
  labs(x = "p-value cutoff", 
       y = "proportion false SVGs") + 
  ggtitle("Error control: mouse OB", 
          subtitle = "null simulation") + 
  theme_bw() + 
  theme(legend.title = element_blank(), 
        panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "error_control_nnSVG_mouseOB_nullSim")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4.5, height = 3.5)

