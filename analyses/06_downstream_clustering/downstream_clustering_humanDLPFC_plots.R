#####################################################
# Script for downstream clustering comparison - plots
# Lukas Weber, updated Jun 2023
#####################################################

# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(mclust)
library(ggplot2)


dir_clustering <- here("outputs", "downstream_clustering")
dir_plots <- here("plots", "downstream_clustering")


# ------------
# load results
# ------------

fn <- here(dir_clustering, "res_downstream_clustering.rds")
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S15ABCDE.RData"

save(res_out, file = here(dir_sd, fn_sd))


# --------------
# match clusters
# --------------

# match clusters to ground truth layers for each method

match_nnSVG <- c(5, 8, 1, 7, 2, 6, 3, 4)
coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, levels = match_nnSVG)

match_SPARKX <- c(1, 6, 2, 8, 3, 4, 7, 5)
coldata_out$SPARKX$label <- factor(
  coldata_out$SPARKX$label, levels = match_SPARKX)

match_HVGs <- c(1, 8, 2, 4, 3, 6, 7, 5)
coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, levels = match_HVGs)

match_MoransI <- c(5, 2, 8, 3, 4, 7, 1, 6)
coldata_out$MoransI$label <- factor(
  coldata_out$MoransI$label, levels = match_MoransI)


# ----------
# spot plots
# ----------

# plot clustering

# LIBD layer colors palette
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")

for (i in seq_along(coldata_out)) {
  
  df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]]))
  
  ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                 color = label)) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human DLPFC: ", 
                              gsub("MoransI", "Moran's I", 
                                   gsub("SPARKX", "SPARK-X", names(coldata_out[i]))))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(dir_plots, 
                  paste0("clustering_humanDLPFC_", names(coldata_out[i])))
  ggsave(paste0(fn, ".pdf"), width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}


# -------------------
# adjusted Rand index
# -------------------

# calculate adjusted Rand index to evaluate clustering performance


# re-label factor levels combine clusters 7 and 8 into a single cluster
# representing white matter

coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

coldata_out$SPARKX$label <- factor(
  coldata_out$SPARKX$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))

coldata_out$MoransI$label <- factor(
  coldata_out$MoransI$label, labels = c(paste0("Layer", 1:6), rep("WM", 2)))


# calculate adjusted Rand index

ari_nnSVG <- adjustedRandIndex(coldata_out$nnSVG$label, 
                               coldata_out$nnSVG$ground_truth)

ari_SPARKX <- adjustedRandIndex(coldata_out$SPARKX$label, 
                                coldata_out$SPARKX$ground_truth)

ari_HVGs <- adjustedRandIndex(coldata_out$HVGs$label, 
                              coldata_out$HVGs$ground_truth)

ari_MoransI <- adjustedRandIndex(coldata_out$MoransI$label, 
                                 coldata_out$MoransI$ground_truth)


# plot adjusted Rand index

df <- data.frame(
  method = c("nnSVG", "SPARKX", "HVGs", "MoransI"), 
  ARI = c(ari_nnSVG, ari_SPARKX, ari_HVGs, ari_MoransI)
)
df$method <- factor(df$method, 
                    levels = df$method, 
                    labels = gsub("MoransI", "Moran's I", 
                                  gsub("SPARKX", "SPARK-X", df$method)))


pal_methods <- c("blue3", "deepskyblue2", "darkorange", "firebrick3")


ggplot(df, aes(x = method, y = ARI, shape = method, color = method)) + 
  geom_point(stroke = 1.5, size = 2) + 
  scale_shape_manual(values = c(4, 3, 1, 2)) + 
  scale_color_manual(values = pal_methods) + 
  ylim(c(0, 1)) + 
  ggtitle("Downstream clustering performance") + 
  labs(y = "Adjusted Rand Index") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())

fn <- file.path(dir_plots, "summary_clustering_performance_ARI")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 3.1)
ggsave(paste0(fn, ".png"), width = 4.75, height = 3.1)

