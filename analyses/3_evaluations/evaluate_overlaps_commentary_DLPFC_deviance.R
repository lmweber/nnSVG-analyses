########################################
# Script to calculate evaluation metrics
# Lukas Weber, Nov 2021
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

# DLPFC dataset
res_DLPFC_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_DLPFC_deviance.rds")))
res_DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds")))
res_DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds")))

# check row order of gene names is consistent
stopifnot(all(res_DLPFC_deviance$gene_id == res_DLPFC_nnSVG$gene_id))
stopifnot(all(res_DLPFC_deviance$gene_id == res_DLPFC_SPARKX$gene_id))


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of significant SVGs from a given method 
# that are also in the set of top n deviance (and vice versa)

# identify top n deviance
top10deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 10]
top20deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 20]
top100deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 100]
top200deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 200]
top1000deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 1000]
top2000deviance <- res_DLPFC_deviance$gene_id[res_DLPFC_deviance$rank <= 2000]

# number of significant SVGs per method
n_signnSVG <- sum(res_DLPFC_nnSVG$padj <= 0.05)
n_sigSPARKX <- sum(res_DLPFC_SPARKX$adjustedPval <= 0.05)

# identify significant SVGs
signnSVG <- res_DLPFC_nnSVG$gene_id[res_DLPFC_nnSVG$padj <= 0.05]
sigSPARKX <- res_DLPFC_SPARKX$gene_id[res_DLPFC_SPARKX$adjustedPval <= 0.05]

# identify top n SVGs
top10nnSVG <- res_DLPFC_nnSVG$gene_id[res_DLPFC_nnSVG$rank <= 10]
top100nnSVG <- res_DLPFC_nnSVG$gene_id[res_DLPFC_nnSVG$rank <= 100]
top1000nnSVG <- res_DLPFC_nnSVG$gene_id[res_DLPFC_nnSVG$rank <= 1000]

top10SPARKX <- res_DLPFC_SPARKX$gene_id[res_DLPFC_SPARKX$rank <= 10]
top100SPARKX <- res_DLPFC_SPARKX$gene_id[res_DLPFC_SPARKX$rank <= 100]
top1000SPARKX <- res_DLPFC_SPARKX$gene_id[res_DLPFC_SPARKX$rank <= 1000]


# calculate overlaps for each method


# SIGNIFICANT SVGS

# deviance within set of SVGs (since n_SVGs >> n_deviance)
prop_overlap_top10deviance_signnSVG <- mean(top10deviance %in% signnSVG)
prop_overlap_top100deviance_signnSVG <- mean(top100deviance %in% signnSVG)
prop_overlap_top1000deviance_signnSVG <- mean(top1000deviance %in% signnSVG)

prop_overlap_top10deviance_sigSPARKX <- mean(top10deviance %in% sigSPARKX)
prop_overlap_top100deviance_sigSPARKX <- mean(top100deviance %in% sigSPARKX)
prop_overlap_top1000deviance_sigSPARKX <- mean(top1000deviance %in% sigSPARKX)


# TOP-RANKED SVGS

# SVGs within set of deviance (since n_deviance >> n_SVGs for most of these)
prop_overlap_top10nnSVG_top10deviance <- mean(top10nnSVG %in% top10deviance)
prop_overlap_top100nnSVG_top100deviance <- mean(top100nnSVG %in% top100deviance)
prop_overlap_top1000nnSVG_top1000deviance <- mean(top1000nnSVG %in% top1000deviance)

prop_overlap_top10nnSVG_top20deviance <- mean(top10nnSVG %in% top20deviance)
prop_overlap_top100nnSVG_top200deviance <- mean(top100nnSVG %in% top200deviance)
prop_overlap_top1000nnSVG_top2000deviance <- mean(top1000nnSVG %in% top2000deviance)


prop_overlap_top10SPARKX_top10deviance <- mean(top10SPARKX %in% top10deviance)
prop_overlap_top100SPARKX_top100deviance <- mean(top100SPARKX %in% top100deviance)
prop_overlap_top1000SPARKX_top1000deviance <- mean(top1000SPARKX %in% top1000deviance)

prop_overlap_top10SPARKX_top20deviance <- mean(top10SPARKX %in% top20deviance)
prop_overlap_top100SPARKX_top200deviance <- mean(top100SPARKX %in% top200deviance)
prop_overlap_top1000SPARKX_top2000deviance <- mean(top1000SPARKX %in% top2000deviance)


# ---------------
# plots: overlaps
# ---------------

# top deviance within significant SVGs

df_topdeviance_sigSVGs <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10deviance_signnSVG, prop_overlap_top100deviance_signnSVG, prop_overlap_top1000deviance_signnSVG), 
  SPARKX = c(prop_overlap_top10deviance_sigSPARKX, prop_overlap_top100deviance_sigSPARKX, prop_overlap_top1000deviance_sigSPARKX)
))
colnames(df_topdeviance_sigSVGs) <- c(10, 100, 1000)
df_topdeviance_sigSVGs$method <- as.factor(rownames(df_topdeviance_sigSVGs))
df_topdeviance_sigSVGs <- pivot_longer(
  df_topdeviance_sigSVGs, cols = c("10", "100", "1000"), names_to = "n_deviance", values_to = "proportion"
)
df_topdeviance_sigSVGs$n_deviance <- as.numeric(df_topdeviance_sigSVGs$n_deviance)

ggplot(df_topdeviance_sigSVGs, aes(x = n_deviance, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(10, 100, 1000), labels = c(10, 100, 1000), limits = c(10, 1000), trans = "sqrt") + 
  xlab("number of top deviance") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top deviance within significant SVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topdeviance_sigSVGs")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)



# top SVGs within top deviance

df_topSVGs_topdeviance <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10nnSVG_top10deviance, prop_overlap_top100nnSVG_top100deviance, prop_overlap_top1000nnSVG_top1000deviance), 
  SPARKX = c(prop_overlap_top10SPARKX_top10deviance, prop_overlap_top100SPARKX_top100deviance, prop_overlap_top1000SPARKX_top1000deviance)
))
colnames(df_topSVGs_topdeviance) <- c(10, 100, 1000)
df_topSVGs_topdeviance$method <- as.factor(rownames(df_topSVGs_topdeviance))
df_topSVGs_topdeviance <- pivot_longer(
  df_topSVGs_topdeviance, cols = c("10", "100", "1000"), names_to = "n_deviance", values_to = "proportion"
)
df_topSVGs_topdeviance$n_deviance <- as.numeric(df_topSVGs_topdeviance$n_deviance)

ggplot(df_topSVGs_topdeviance, aes(x = n_deviance, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(10, 100, 1000), labels = c(10, 100, 1000), limits = c(10, 1000), trans = "sqrt") + 
  xlab("number of top deviance") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top SVGs within top deviance") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topSVGs_topdeviance")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# top SVGs within top deviance (x2)

df_topSVGs_topdeviancex2 <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10nnSVG_top20deviance, prop_overlap_top100nnSVG_top200deviance, prop_overlap_top1000nnSVG_top2000deviance), 
  SPARKX = c(prop_overlap_top10SPARKX_top20deviance, prop_overlap_top100SPARKX_top200deviance, prop_overlap_top1000SPARKX_top2000deviance)
))
colnames(df_topSVGs_topdeviancex2) <- c(20, 200, 2000)
df_topSVGs_topdeviancex2$method <- as.factor(rownames(df_topSVGs_topdeviancex2))
df_topSVGs_topdeviancex2 <- pivot_longer(
  df_topSVGs_topdeviancex2, cols = c("20", "200", "2000"), names_to = "n_deviance", values_to = "proportion"
)
df_topSVGs_topdeviancex2$n_deviance <- as.numeric(df_topSVGs_topdeviancex2$n_deviance)

ggplot(df_topSVGs_topdeviancex2, aes(x = n_deviance, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(20, 200, 2000), labels = c(20, 200, 2000), limits = c(20, 2000), trans = "sqrt") + 
  xlab("number of top deviance") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top SVGs within top deviance (x2)") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topSVGs_topdeviancex2")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# ------------------------
# plots: rank correlations
# ------------------------

df_ranks <- data.frame(
  deviance = res_DLPFC_deviance$rank, 
  nnSVG = res_DLPFC_nnSVG$rank, 
  SPARKX = res_DLPFC_SPARKX$rank
)
df_ranks <- pivot_longer(
  df_ranks, cols = c("nnSVG", "SPARKX"), names_to = "method", values_to = "SVGs"
)
df_ranks$method <- as.factor(df_ranks$method)


# top 10
ggplot(df_ranks, aes(x = deviance, y = SVGs, color = method, shape = method)) + 
  geom_point(size = 1.5, stroke = 1) + 
  scale_color_startrek() + 
  coord_fixed() + 
  scale_x_continuous(breaks = seq(0, 10, length.out = 6), limits = c(0, 10)) + 
  scale_y_continuous(breaks = seq(0, 10, length.out = 6), limits = c(0, 10)) + 
  scale_shape_manual(values = c(1, 4)) + 
  xlab("rank deviance") + 
  ylab("rank SVGs") + 
  ggtitle("Ranks of top 10 deviance and SVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "ranks_scatter_top10_deviance")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4.25)
ggsave(paste0(fn, ".png"), width = 5, height = 4.25)


# top 100
ggplot(df_ranks, aes(x = deviance, y = SVGs, color = method, shape = method)) + 
  geom_point(size = 1.5, stroke = 1) + 
  scale_color_startrek() + 
  coord_fixed() + 
  scale_x_continuous(breaks = seq(0, 100, length.out = 6), limits = c(0, 100)) + 
  scale_y_continuous(breaks = seq(0, 100, length.out = 6), limits = c(0, 100)) + 
  scale_shape_manual(values = c(1, 4)) + 
  xlab("rank deviance") + 
  ylab("rank SVGs") + 
  ggtitle("Ranks of top 100 deviance and SVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "ranks_scatter_top100_deviance")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4.25)
ggsave(paste0(fn, ".png"), width = 5, height = 4.25)


# top 1000
ggplot(df_ranks, aes(x = deviance, y = SVGs, color = method, shape = method)) + 
  geom_point(size = 1.5, stroke = 1) + 
  scale_color_startrek() + 
  coord_fixed() + 
  scale_x_continuous(breaks = seq(0, 1000, length.out = 6), limits = c(0, 1000)) + 
  scale_y_continuous(breaks = seq(0, 1000, length.out = 6), limits = c(0, 1000)) + 
  scale_shape_manual(values = c(1, 4)) + 
  xlab("rank deviance") + 
  ylab("rank SVGs") + 
  ggtitle("Ranks of top 1000 deviance and SVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "ranks_scatter_top1000_deviance")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4.25)
ggsave(paste0(fn, ".png"), width = 5, height = 4.25)

