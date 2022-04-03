############################
# Additional plots: top SVGs
# Lukas Weber, Apr 2022
############################

# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(scales)


# directory to save plots
dir_plots <- here(file.path("plots", "top_SVGs"))


# ------------
# load results
# ------------

# load as SPE object containing counts
spe <- readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))

# note filtering
dim(spe)


# -------------------------
# top SVGs from this method
# -------------------------

# select top n SVGs

n <- 20
ix_top <- order(rowData(spe)$rank)[seq_len(n)]

# check gene names
top_names <- rowData(spe)$gene_name[ix_top]
top_names

# counts for top n SVGs
counts_top <- counts(spe)[ix_top, ] %>% as.matrix() %>% t()
colnames(counts_top) <- top_names

dim(counts_top)
head(counts_top)

stopifnot(all(rownames(counts_top) == rownames(spatialCoords(spe))))


# plot top n SVGs

df <- as.data.frame(cbind(colData(spe), spatialCoords(spe), counts_top)) %>% 
  select(c("barcode_id", "in_tissue", top_names, 
           "pxl_col_in_fullres", "pxl_row_in_fullres")) %>% 
  filter(in_tissue == 1) %>% 
  pivot_longer(., cols = top_names, 
               names_to = "gene", values_to = "counts") %>% 
  mutate(gene = factor(gene, levels = top_names))

max_counts <- max(counts_top)
max_counts


ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = counts)) + 
  facet_wrap(~ gene, nrow = 4) + 
  geom_point(size = 0.05) + 
  coord_fixed() + 
  scale_y_reverse() + 
  scale_color_gradientn(trans = "sqrt", 
                        colors = c("gray90", "blue", "black"), 
                        values = rescale(c(0, 120, 195)), 
                        labels = c("0", "", "", "120", "")) + 
  ggtitle("Top SVGs: human DLPFC, nnSVG") + 
  theme_bw() + 
  guides(color = guide_colorbar(ticks = FALSE)) + 
  theme(panel.grid = element_blank(), 
        axis.title = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank())

fn <- file.path(dir_plots, "top_SVGs_humanDLPFC_nnSVG")
ggsave(paste0(fn, ".png"), width = 8.25, height = 8)

