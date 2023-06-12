###############################
# Script to plot cell locations
# Lukas Weber, updated Jun 2023
###############################

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

# Lohoff and Ghazanfar et al. (2021)

spe <- seqFISH_mouseEmbryo()

# cell locations are stored in spatialCoords
colData(spe)
head(spatialCoords(spe))

df <- as.data.frame(spatialCoords(spe))


# plot cell locations
ggplot(df, aes(x = x, y = y)) + 
  geom_point(size = 0.01, color = "gray20") + 
  coord_fixed() + 
  ggtitle("Cells: mouse embryo") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "cell_locations_mouseEmbryo")
ggsave(paste0(fn, ".png"), width = 3, height = 4.25)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S9A.RData"

save(df, file = here(dir_sd, fn_sd))

