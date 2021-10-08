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

# mOB dataset

list_mOB <- list(
  spe_mOB_nnSVG = readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG.rds")), 
  spe_mOB_nnSVG_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG_clusters.rds")), 
  
  spe_mOB_nnSVG_logcounts = readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG_logcounts.rds")), 
  spe_mOB_nnSVG_logcounts_clusters = readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG_logcounts_clusters.rds")), 
  
  spe_mOB_HVGs = readRDS(here("outputs", "results", "HVGs", "spe_mOB_HVGs.rds")), 
  
  spe_mOB_deviance = readRDS(here("outputs", "results", "deviance", "spe_mOB_deviance.rds")), 
  spe_mOB_deviance_clusters = readRDS(here("outputs", "results", "deviance", "spe_mOB_deviance_clusters.rds"))
)

for (i in seq_along(list_mOB)){
  stopifnot(all(rowData(list_mOB[[i]])$gene_name == rowData(list_mOB[[1]])$gene_name))
}

# method and dataset names

names_mOB <- gsub("^spe_", "", names(list_mOB))
names_mOB

names(list_mOB) <- names_mOB


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of top n genes from method 1 (e.g. HVGs)
# that are also in the set of top n genes from method 2 (e.g. nnSVG)

# overlap sizes
overlaps <- c(10, 20, 50, 100, 200, 500, 1000)

# function to calculate overlaps for a pair of methods
calc_overlaps <- function(method1, method2) {
  
  res_method1 <- rowData(list_mOB[[method1]])
  res_method2 <- rowData(list_mOB[[method2]])
  
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

df_overlaps_mOB_HVGs <- data.frame(
  top_n = overlaps, 
  nnSVG = calc_overlaps("mOB_HVGs", "mOB_nnSVG"), 
  nnSVG_clusters = calc_overlaps("mOB_HVGs", "mOB_nnSVG_clusters"), 
  nnSVG_logcounts = calc_overlaps("mOB_HVGs", "mOB_nnSVG_logcounts"), 
  nnSVG_logcounts_clusters = calc_overlaps("mOB_HVGs", "mOB_nnSVG_logcounts_clusters")
)

df_overlaps_mOB_deviance <- data.frame(
  top_n = overlaps, 
  nnSVG = calc_overlaps("mOB_deviance", "mOB_nnSVG"), 
  nnSVG_clusters = calc_overlaps("mOB_deviance", "mOB_nnSVG_clusters"), 
  nnSVG_logcounts = calc_overlaps("mOB_deviance", "mOB_nnSVG_logcounts"), 
  nnSVG_logcounts_clusters = calc_overlaps("mOB_deviance", "mOB_nnSVG_logcounts_clusters")
)

df_overlaps_mOB_HVGs_vs_deviance <- data.frame(
  top_n = overlaps, 
  deviance = calc_overlaps("mOB_HVGs", "mOB_deviance"), 
  deviance_clusters = calc_overlaps("mOB_HVGs", "mOB_deviance_clusters")
)


# create data frames for plotting

df_overlaps_mOB_HVGs <- pivot_longer(
  df_overlaps_mOB_HVGs, 
  cols = colnames(df_overlaps_mOB_HVGs)[-1], 
  names_to = "method", 
  values_to = "prop_overlap"
)

df_overlaps_mOB_deviance <- pivot_longer(
  df_overlaps_mOB_deviance, 
  cols = colnames(df_overlaps_mOB_deviance)[-1], 
  names_to = "method", 
  values_to = "prop_overlap"
)

df_overlaps_mOB_HVGs_vs_deviance <- pivot_longer(
  df_overlaps_mOB_HVGs_vs_deviance, 
  cols = colnames(df_overlaps_mOB_HVGs_vs_deviance)[-1], 
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
    labs(x = paste0("top n ", name_x), 
         y = "proportion overlapping") + 
    ggtitle(paste0("mOB: Overlap with top ", name_x)) + 
    theme_bw() + 
    theme(panel.grid.minor.x = element_blank())
}


# generate and save plots

plot_overlaps(df_overlaps_mOB_HVGs, "HVGs")
fn <- here("plots", "overlaps", "overlaps_mOB_HVGs")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_mOB_deviance, "deviance")
fn <- here("plots", "overlaps", "overlaps_mOB_deviance")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

plot_overlaps(df_overlaps_mOB_HVGs_vs_deviance, "HVGs_vs_deviance")
fn <- here("plots", "overlaps", "overlaps_mOB_HVGs_vs_deviance")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)



# ----------------------
# calculate correlations
# ----------------------

# correlations are calculated as Spearman rank correlations between ranks of 
# genes within the top n ranks from both method 1 (e.g. HVGs) and method 2 
# (e.g. nnSVG); note this is usually less than n genes in total

correlations <- c(100, 500, 1000)
cor_HVGs_nnSVG <- rep(NA, length(correlations))
cor_HVGs_nnSVG_clusters <- rep(NA, length(correlations))
cor_nnSVG_nnSVG_clusters <- rep(NA, length(correlations))

# subset and calculate correlations for each combination of methods
for (k in seq_along(correlations)) {
  res_sub <- 
    filter(res_mOB, rank_HVGs <= correlations[k] & rank_nnSVG <= correlations[k])
  cor_HVGs_nnSVG[k] <- 
    cor(res_sub$rank_HVGs, res_sub$rank_nnSVG, method = "spearman")
  
  res_sub <- 
    filter(res_mOB, rank_HVGs <= correlations[k] & rank_nnSVG_clusters <= correlations[k])
  cor_HVGs_nnSVG_clusters[k] <- 
    cor(res_sub$rank_HVGs, res_sub$rank_nnSVG_clusters, method = "spearman")
  
  res_sub <- 
    filter(res_mOB, rank_nnSVG <= correlations[k] & rank_nnSVG_clusters <= correlations[k])
  cor_nnSVG_nnSVG_clusters[k] <- 
    cor(res_sub$rank_nnSVG, res_sub$rank_nnSVG_clusters, method = "spearman")
}

# summary data frame
df_correlations <- data.frame(
  top_n = correlations, 
  HVGs_nnSVG = cor_HVGs_nnSVG, 
  HVGs_nnSVGclusters = cor_HVGs_nnSVG_clusters, 
  nnSVG_nnSVGclusters = cor_nnSVG_nnSVG_clusters
)

df_correlations <- pivot_longer(
  df_correlations, 
  cols = colnames(df_correlations)[-1], 
  names_to = "methods", 
  values_to = "correlation"
)

# plot
ggplot(df_correlations, aes(x = top_n, y = correlation, 
                            group = methods, color = methods)) + 
  geom_line() + 
  geom_point() + 
  scale_color_ochre(palette = "nolan_ned") + 
  scale_x_continuous(breaks = df_correlations$top_n) + 
  labs(x = "genes within top n genes from both methods", 
       y = "correlation between ranks") + 
  ggtitle("mOB: Correlations between top genes") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank())

fn <- here("plots", "correlations", "correlations_mOB")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)



# ----------
# plot ranks
# ----------

# plot ranks of top n genes from each pair of methods

ranks <- c(20, 100, 1000)

for (k in seq_along(ranks)) {
  p <- ggplot(res_mOB, aes(x = rank_HVGs, y = rank_nnSVG, label = gene_name)) + 
    geom_point() + 
    xlim(c(0, ranks[k])) + 
    ylim(c(0, ranks[k])) + 
    coord_fixed() + 
    ggtitle(paste0("Ranks of top ", ranks[k], " genes per method")) + 
    theme_bw()
  
  if (ranks[k] == 20) {
    p <- p + geom_text_repel(position = "dodge")
  }
  
  print(p)
  
  fn <- here("plots", "ranks", paste0("ranks_mOB_HVGs_nnSVG_", ranks[k]))
  ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4.5)
  ggsave(paste0(fn, ".png"), width = 4.5, height = 4.5)
}

