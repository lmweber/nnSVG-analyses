########################################
# Script to calculate evaluation metrics
# Lukas Weber, Oct 2021
########################################

# to run in interactive session on cluster:
# module load conda_R/4.1.x
# Rscript filename.R


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsci)


# ------------
# load results
# ------------

# SlideSeqHippo dataset

list_SlideSeqHippo <- list(
  #spe_SlideSeqHippo_nnSVG = readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds")), 
  spe_SlideSeqHippo_nnSVG_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds")), 
  
  #spe_SlideSeqHippo_nnSVG_logcounts = readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_logcounts.rds")), 
  #spe_SlideSeqHippo_nnSVG_logcounts_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_logcounts_clusters.rds")), 
  
  #spe_SlideSeqHippo_SPARKX = readRDS(here("outputs", "results", "SPARK-X", "spe_SlideSeqHippo_SPARK-X.rds")), 
  spe_SlideSeqHippo_SPARKX_clusters = readRDS(here("outputs", "results", "SPARK-X", "spe_SlideSeqHippo_SPARK-X_clusters.rds")), 
  
  #spe_SlideSeqHippo_HVGs = readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds"))#, 
  #spe_SlideSeqHippo_deviance = readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance.rds"))#, 
  spe_SlideSeqHippo_deviance_clusters = readRDS(here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance_clusters.rds"))
)

for (i in seq_along(list_SlideSeqHippo)){
  stopifnot(all(rowData(list_SlideSeqHippo[[i]])$gene_name == rowData(list_SlideSeqHippo[[1]])$gene_name))
}

# method and dataset names

names_SlideSeqHippo <- gsub("^spe_", "", names(list_SlideSeqHippo))
names_SlideSeqHippo

names(list_SlideSeqHippo) <- names_SlideSeqHippo


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of top n genes from method 1 (e.g. HVGs)
# that are also in the set of top n genes from method 2 (e.g. nnSVG)

# overlap sizes
overlaps <- c(10, 20, 50, 100, 200, 500, 1000)

# function to calculate overlaps for a pair of methods
calc_overlaps <- function(method1, method2) {
  
  res_method1 <- rowData(list_SlideSeqHippo[[method1]])
  res_method2 <- rowData(list_SlideSeqHippo[[method2]])
  
  top_method2 <- rep(NA, length(overlaps))
  
  for (k in seq_along(overlaps)) {
    # select top gene names from method 1
    genes_k <- rownames(filter(as.data.frame(res_method1), 
                               rank <= overlaps[k]))
    # calculate overlaps
    top_method2[k] <- nrow(filter(as.data.frame(res_method2[genes_k, ]), 
                                  rank <= overlaps[k]))
  }
  
  # calculate proportions
  top_method2 / overlaps
}


# use function to calculate overlaps for sets of methods

# df_overlaps_SlideSeqHippo_HVGs <- data.frame(
#   top_n = overlaps, 
#   #nnSVG = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_nnSVG"), 
#   nnSVG_clusters = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_nnSVG_clusters"), 
#   #nnSVG_logcounts = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_nnSVG_logcounts"), 
#   #nnSVG_logcounts_clusters = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_nnSVG_logcounts_clusters"), 
#   #SPARKX = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_SPARKX"), 
#   SPARKX_clusters = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_SPARKX_clusters")
# )

# df_overlaps_SlideSeqHippo_deviance <- data.frame(
#   top_n = overlaps, 
#   nnSVG = calc_overlaps("SlideSeqHippo_deviance", "SlideSeqHippo_nnSVG"), 
#   #nnSVG_logcounts = calc_overlaps("SlideSeqHippo_deviance", "SlideSeqHippo_nnSVG_logcounts"), 
#   SPARKX = calc_overlaps("SlideSeqHippo_deviance", "SlideSeqHippo_SPARKX")
# )

# df_overlaps_SlideSeqHippo_HVGs_vs_deviance <- data.frame(
#   top_n = overlaps, 
#   deviance = calc_overlaps("SlideSeqHippo_HVGs", "SlideSeqHippo_deviance")
# )

df_overlaps_SlideSeqHippo_deviance_clusters <- data.frame(
  top_n = overlaps,
  nnSVG_clusters = calc_overlaps("SlideSeqHippo_deviance_clusters", "SlideSeqHippo_nnSVG_clusters"),
  #nnSVG_logcounts_clusters = calc_overlaps("SlideSeqHippo_deviance_clusters", "SlideSeqHippo_nnSVG_logcounts_clusters"), 
  SPARKX_clusters = calc_overlaps("SlideSeqHippo_deviance_clusters", "SlideSeqHippo_SPARKX_clusters")
)


# create data frames for plotting

# df_overlaps_SlideSeqHippo_HVGs <- pivot_longer(
#   df_overlaps_SlideSeqHippo_HVGs, 
#   cols = colnames(df_overlaps_SlideSeqHippo_HVGs)[-1], 
#   names_to = "method", 
#   values_to = "prop_overlap"
# )

# df_overlaps_SlideSeqHippo_deviance <- pivot_longer(
#   df_overlaps_SlideSeqHippo_deviance, 
#   cols = colnames(df_overlaps_SlideSeqHippo_deviance)[-1], 
#   names_to = "method", 
#   values_to = "prop_overlap"
# )

# df_overlaps_SlideSeqHippo_HVGs_vs_deviance <- pivot_longer(
#   df_overlaps_SlideSeqHippo_HVGs_vs_deviance, 
#   cols = colnames(df_overlaps_SlideSeqHippo_HVGs_vs_deviance)[-1], 
#   names_to = "method", 
#   values_to = "prop_overlap"
# )

df_overlaps_SlideSeqHippo_deviance_clusters <- pivot_longer(
  df_overlaps_SlideSeqHippo_deviance_clusters,
  cols = colnames(df_overlaps_SlideSeqHippo_deviance_clusters)[-1],
  names_to = "method",
  values_to = "prop_overlap"
)


# function to generate plots

plot_overlaps <- function(df, name_x) {
  ggplot(df, aes(x = top_n, y = prop_overlap, 
                 group = method, color = method)) + 
    geom_line() + 
    geom_point() + 
    scale_color_startrek() + 
    scale_x_continuous(breaks = df$top_n, trans = "log10") + 
    ylim(c(0, 1)) + 
    labs(x = paste0("top n genes: ", name_x), 
         y = "proportion overlapping") + 
    ggtitle(paste0("SlideSeqHippo: Overlap with top ", name_x)) + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank())
}


# generate and save plots

# plot_overlaps(df_overlaps_SlideSeqHippo_HVGs, "HVGs")
# fn <- here("plots", "overlaps", "overlaps_SlideSeqHippo_HVGs")
# ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
# ggsave(paste0(fn, ".png"), width = 6, height = 4)

# plot_overlaps(df_overlaps_SlideSeqHippo_deviance, "deviance")
# fn <- here("plots", "overlaps", "overlaps_SlideSeqHippo_deviance")
# ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
# ggsave(paste0(fn, ".png"), width = 6, height = 4)

# plot_overlaps(df_overlaps_SlideSeqHippo_HVGs_vs_deviance, "HVGs")
# fn <- here("plots", "overlaps", "overlaps_SlideSeqHippo_HVGs_vs_deviance")
# ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
# ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_SlideSeqHippo_deviance_clusters, "deviance_clusters")
fn <- here("plots", "overlaps", "overlaps_SlideSeqHippo_deviance_clusters")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

