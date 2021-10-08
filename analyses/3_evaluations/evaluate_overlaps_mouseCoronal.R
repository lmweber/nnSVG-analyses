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

# mouseCoronal dataset

list_mouseCoronal <- list(
  spe_mouseCoronal_nnSVG = readRDS(here("outputs", "results", "nnSVG", "spe_mouseCoronal_nnSVG.rds")), 
  spe_mouseCoronal_nnSVG_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_mouseCoronal_nnSVG_clusters.rds")), 
  
  spe_mouseCoronal_nnSVG_logcounts = readRDS(here("outputs", "results", "nnSVG", "spe_mouseCoronal_nnSVG_logcounts.rds")), 
  spe_mouseCoronal_nnSVG_logcounts_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_mouseCoronal_nnSVG_logcounts_clusters.rds")), 
  
  spe_mouseCoronal_HVGs = readRDS(here("outputs", "results", "HVGs", "spe_mouseCoronal_HVGs.rds")), 
  
  spe_mouseCoronal_deviance = readRDS(here("outputs", "results", "deviance", "spe_mouseCoronal_deviance.rds")), 
  spe_mouseCoronal_deviance_clusters = readRDS(here("outputs", "results", "deviance", "spe_mouseCoronal_deviance_clusters.rds"))
)

for (i in seq_along(list_mouseCoronal)){
  stopifnot(all(rowData(list_mouseCoronal[[i]])$gene_name == rowData(list_mouseCoronal[[1]])$gene_name))
}

# method and dataset names

names_mouseCoronal <- gsub("^spe_", "", names(list_mouseCoronal))
names_mouseCoronal

names(list_mouseCoronal) <- names_mouseCoronal


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of top n genes from method 1 (e.g. HVGs)
# that are also in the set of top n genes from method 2 (e.g. nnSVG)

# overlap sizes
overlaps <- c(10, 20, 50, 100, 200, 500, 1000)

# function to calculate overlaps for a pair of methods
calc_overlaps <- function(method1, method2) {
  
  res_method1 <- rowData(list_mouseCoronal[[method1]])
  res_method2 <- rowData(list_mouseCoronal[[method2]])
  
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

df_overlaps_mouseCoronal_HVGs <- data.frame(
  top_n = overlaps, 
  nnSVG = calc_overlaps("mouseCoronal_HVGs", "mouseCoronal_nnSVG"), 
  nnSVG_clusters = calc_overlaps("mouseCoronal_HVGs", "mouseCoronal_nnSVG_clusters"), 
  nnSVG_logcounts = calc_overlaps("mouseCoronal_HVGs", "mouseCoronal_nnSVG_logcounts"), 
  nnSVG_logcounts_clusters = calc_overlaps("mouseCoronal_HVGs", "mouseCoronal_nnSVG_logcounts_clusters")
)

df_overlaps_mouseCoronal_deviance <- data.frame(
  top_n = overlaps, 
  nnSVG = calc_overlaps("mouseCoronal_deviance", "mouseCoronal_nnSVG"), 
  nnSVG_logcounts = calc_overlaps("mouseCoronal_deviance", "mouseCoronal_nnSVG_logcounts")
)

df_overlaps_mouseCoronal_HVGs_vs_deviance <- data.frame(
  top_n = overlaps, 
  deviance = calc_overlaps("mouseCoronal_HVGs", "mouseCoronal_deviance")
)

df_overlaps_mouseCoronal_deviance_clusters <- data.frame(
  top_n = overlaps, 
  nnSVG_clusters = calc_overlaps("mouseCoronal_deviance_clusters", "mouseCoronal_nnSVG_clusters"), 
  nnSVG_logcounts_clusters = calc_overlaps("mouseCoronal_deviance_clusters", "mouseCoronal_nnSVG_logcounts_clusters")
)


# create data frames for plotting

df_overlaps_mouseCoronal_HVGs <- pivot_longer(
  df_overlaps_mouseCoronal_HVGs, 
  cols = colnames(df_overlaps_mouseCoronal_HVGs)[-1], 
  names_to = "method", 
  values_to = "prop_overlap"
)

df_overlaps_mouseCoronal_deviance <- pivot_longer(
  df_overlaps_mouseCoronal_deviance, 
  cols = colnames(df_overlaps_mouseCoronal_deviance)[-1], 
  names_to = "method", 
  values_to = "prop_overlap"
)

df_overlaps_mouseCoronal_HVGs_vs_deviance <- pivot_longer(
  df_overlaps_mouseCoronal_HVGs_vs_deviance, 
  cols = colnames(df_overlaps_mouseCoronal_HVGs_vs_deviance)[-1], 
  names_to = "method", 
  values_to = "prop_overlap"
)

df_overlaps_mouseCoronal_deviance_clusters <- pivot_longer(
  df_overlaps_mouseCoronal_deviance_clusters, 
  cols = colnames(df_overlaps_mouseCoronal_deviance_clusters)[-1], 
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
    labs(x = paste0("top n genes: ", name_x), 
         y = "proportion overlapping") + 
    ggtitle(paste0("mouseCoronal: Overlap with top ", name_x)) + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank())
}


# generate and save plots

plot_overlaps(df_overlaps_mouseCoronal_HVGs, "HVGs")
fn <- here("plots", "overlaps", "overlaps_mouseCoronal_HVGs")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_mouseCoronal_deviance, "deviance")
fn <- here("plots", "overlaps", "overlaps_mouseCoronal_deviance")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_mouseCoronal_HVGs_vs_deviance, "HVGs")
fn <- here("plots", "overlaps", "overlaps_mouseCoronal_HVGs_vs_deviance")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_mouseCoronal_deviance_clusters, "deviance_clusters")
fn <- here("plots", "overlaps", "overlaps_mouseCoronal_deviance_clusters")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

