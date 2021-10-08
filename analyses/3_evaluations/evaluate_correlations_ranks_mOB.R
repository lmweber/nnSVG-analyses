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

