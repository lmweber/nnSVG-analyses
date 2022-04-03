####################################
# Script to plot ground truth layers
# Lukas Weber, Apr 2022
####################################

library(SpatialExperiment)
library(STexampleData)
library(ggplot2)
library(dplyr)
library(tidyr)
library(here)


# directory to save plots
dir_plots <- here(file.path("plots", "example_genes"))


# ------------------
# Visium human DLPFC
# ------------------

# Maynard and Collado-Torres et al. (2021)

spe <- Visium_humanDLPFC()

# ground truth layer labels are stored in colData
colData(spe)

df <- cbind(colData(spe), spatialCoords(spe)) %>% 
  as.data.frame() %>% 
  filter(in_tissue == 1)

# add additional column of white matter (WM) vs. layers
df$ground_truth_collapsed <- as.numeric(df$ground_truth == "WM")
df$ground_truth_collapsed[is.na(df$ground_truth_collapsed)] <- 2

# update factors
levs <- c(sort(unique(df$ground_truth)), "none")
df$ground_truth[is.na(df$ground_truth)] <- "none"
df$ground_truth <- factor(df$ground_truth, levels = levs)

levs_collapsed <- c("Layers", "WM", "none")
df$ground_truth_collapsed <- factor(df$ground_truth_collapsed, labels = levs_collapsed)


# layer color palette from paper
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", 
         "#FFD700", "#FF7F00", "#1A1A1A", "#666666")


# plot all layers
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = ground_truth)) + 
  geom_point(size = 0.3) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = pal, name = "label") + 
  ggtitle("Ground truth labels: human DLPFC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "ground_truth_humanDLPFC")
ggsave(paste0(fn, ".png"), width = 4, height = 3.75)


# plot WM vs. layers
ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
               color = ground_truth_collapsed)) + 
  geom_point(size = 0.3) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_manual(values = c("firebrick1", "gray20", "gray50"), name = "label") + 
  ggtitle("Ground truth labels: human DLPFC") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "ground_truth_collapsed_humanDLPFC")
ggsave(paste0(fn, ".png"), width = 4, height = 3.75)

