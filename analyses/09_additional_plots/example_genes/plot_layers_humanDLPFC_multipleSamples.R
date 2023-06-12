######################################################
# Script to plot ground truth layers: multiple samples
# Lukas Weber, updated Jun 2023
######################################################

library(SpatialExperiment)
library(here)
library(scran)
library(scater)
library(ggplot2)
library(dplyr)
library(tidyr)


# directory to save plots
dir_plots <- here(file.path("plots", "example_genes"))


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "multiple_samples", "spe_humanDLPFC_multipleSamples_preprocessed.rds")
spe <- readRDS(fn)

dim(spe)
table(colData(spe)$sample_id)
table(colData(spe)$subject)

# 12 samples
n_samples <- length(table(colData(spe)$sample_id))
n_samples

# 3 donors
n_donors <- length(table(colData(spe)$subject))
n_donors


# --------------
# generate plots
# --------------

df <- cbind(colData(spe), spatialCoords(spe)) %>% 
  as.data.frame() %>% 
  filter(in_tissue == 1)

# update factors
levs <- c(sort(as.character(unique(df$layer_guess_reordered))), "none")
df$ground_truth <- as.character(df$layer_guess_reordered)
df$ground_truth[is.na(df$ground_truth)] <- "none"
df$ground_truth <- factor(df$ground_truth, levels = levs)

df$donor_id <- as.factor(df$subject)

df$sample_id <- factor(df$sample_id, 
                       levels = unique(df$sample_id), 
                       labels = paste0("sample ", unique(df$sample_id)))


# layer color palette from paper
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
         "#FFD700", "#FF7F00", "#1A1A1A", "#666666")


# plot all layers
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = ground_truth)) + 
  facet_wrap(~ sample_id, nrow = 3) + 
  geom_point(size = 0.1) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = pal, name = "label") + 
  ggtitle("Layer labels: human DLPFC") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "ground_truth_humanDLPFC_multiple_samples")
ggsave(paste0(fn, ".png"), width = 8.25, height = 6.25)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S23B.RData"

save(df, file = here(dir_sd, fn_sd))

