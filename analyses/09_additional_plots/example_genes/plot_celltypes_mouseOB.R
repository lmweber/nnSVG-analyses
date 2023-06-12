#################################
# Script to plot cell type layers
# Lukas Weber, updated Jun 2023
#################################

library(SpatialExperiment)
library(STexampleData)
library(dplyr)
library(tidyr)
library(ggplot2)
library(here)


# directory to save plots
dir_plots <- here(file.path("plots", "example_genes"))


# -----------
# ST mouse OB
# -----------

# Stahl et al. (2016)

spe <- ST_mouseOB()

# cell type layer labels are stored in colData
colData(spe)

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe)))
levs <- c(sort(unique(df$layer)), "none")
df$layer[is.na(df$layer)] <- "none"
df$layer <- factor(df$layer, levels = levs)

pal <- unname(palette.colors(8, "Okabe-Ito"))
pal <- c(pal[c(2, 3, 7, 5, 6)], "gray80")


ggplot(df, aes(x = x, y = y, color = layer)) + 
  geom_point(size = 3.5) + 
  scale_color_manual(values = pal) + 
  coord_fixed() + 
  ggtitle("Mouse OB") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "mouseOB_layers")
ggsave(paste0(fn, ".png"), width = 5.25, height = 3.25)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S7A.RData"

save(df, file = here(dir_sd, fn_sd))

