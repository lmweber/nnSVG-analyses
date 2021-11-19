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
res_DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds")))
res_DLPFC_deviance = rowData(readRDS(here("outputs", "results", "deviance", "spe_DLPFC_deviance.rds")))
res_DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds")))
res_DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds")))

# check row order of gene names is consistent
stopifnot(all(res_DLPFC_HVGs$gene_id == res_DLPFC_deviance$gene_id))
stopifnot(all(res_DLPFC_HVGs$gene_id == res_DLPFC_nnSVG$gene_id))
stopifnot(all(res_DLPFC_HVGs$gene_id == res_DLPFC_SPARKX$gene_id))


# ------------------
# calculate overlaps
# ------------------

# overlaps calculated as proportion of significant SVGs from a given method 
# that are also in the set of top n HVGs (and vice versa)

# identify top n HVGs
top10HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 10]
top20HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 20]
top100HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 100]
top200HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 200]
top1000HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 1000]
top2000HVGs <- res_DLPFC_HVGs$gene_id[res_DLPFC_HVGs$rank <= 2000]

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

# HVGs within set of SVGs (since n_SVGs >> n_HVGs)
prop_overlap_top10HVGs_signnSVG <- mean(top10HVGs %in% signnSVG)
prop_overlap_top100HVGs_signnSVG <- mean(top100HVGs %in% signnSVG)
prop_overlap_top1000HVGs_signnSVG <- mean(top1000HVGs %in% signnSVG)

prop_overlap_top10HVGs_sigSPARKX <- mean(top10HVGs %in% sigSPARKX)
prop_overlap_top100HVGs_sigSPARKX <- mean(top100HVGs %in% sigSPARKX)
prop_overlap_top1000HVGs_sigSPARKX <- mean(top1000HVGs %in% sigSPARKX)


# TOP-RANKED SVGS

# SVGs within set of HVGs (since n_HVGs >> n_SVGs for most of these)
prop_overlap_top10nnSVG_top10HVGs <- mean(top10nnSVG %in% top10HVGs)
prop_overlap_top100nnSVG_top100HVGs <- mean(top100nnSVG %in% top100HVGs)
prop_overlap_top1000nnSVG_top1000HVGs <- mean(top1000nnSVG %in% top1000HVGs)

prop_overlap_top10nnSVG_top20HVGs <- mean(top10nnSVG %in% top20HVGs)
prop_overlap_top100nnSVG_top200HVGs <- mean(top100nnSVG %in% top200HVGs)
prop_overlap_top1000nnSVG_top2000HVGs <- mean(top1000nnSVG %in% top2000HVGs)


prop_overlap_top10SPARKX_top10HVGs <- mean(top10SPARKX %in% top10HVGs)
prop_overlap_top100SPARKX_top100HVGs <- mean(top100SPARKX %in% top100HVGs)
prop_overlap_top1000SPARKX_top1000HVGs <- mean(top1000SPARKX %in% top1000HVGs)

prop_overlap_top10SPARKX_top20HVGs <- mean(top10SPARKX %in% top20HVGs)
prop_overlap_top100SPARKX_top200HVGs <- mean(top100SPARKX %in% top200HVGs)
prop_overlap_top1000SPARKX_top2000HVGs <- mean(top1000SPARKX %in% top2000HVGs)


# -------------
# summary plots
# -------------

# top HVGs within significant SVGs

df_topHVGs_sigSVGs <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10HVGs_signnSVG, prop_overlap_top100HVGs_signnSVG, prop_overlap_top1000HVGs_signnSVG), 
  SPARKX = c(prop_overlap_top10HVGs_sigSPARKX, prop_overlap_top100HVGs_sigSPARKX, prop_overlap_top1000HVGs_sigSPARKX)
))
colnames(df_topHVGs_sigSVGs) <- c(10, 100, 1000)
df_topHVGs_sigSVGs$method <- as.factor(rownames(df_topHVGs_sigSVGs))
df_topHVGs_sigSVGs <- pivot_longer(
  df_topHVGs_sigSVGs, cols = c("10", "100", "1000"), names_to = "n_HVGs", values_to = "proportion"
)
df_topHVGs_sigSVGs$n_HVGs <- as.numeric(df_topHVGs_sigSVGs$n_HVGs)

ggplot(df_topHVGs_sigSVGs, aes(x = n_HVGs, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(10, 100, 1000), labels = c(10, 100, 1000), limits = c(10, 1000), trans = "sqrt") + 
  xlab("number of top HVGs") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top HVGs within significant SVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topHVGs_sigSVGs")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)



# top SVGs within top HVGs

df_topSVGs_topHVGs <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10nnSVG_top10HVGs, prop_overlap_top100nnSVG_top100HVGs, prop_overlap_top1000nnSVG_top1000HVGs), 
  SPARKX = c(prop_overlap_top10SPARKX_top10HVGs, prop_overlap_top100SPARKX_top100HVGs, prop_overlap_top1000SPARKX_top1000HVGs)
))
colnames(df_topSVGs_topHVGs) <- c(10, 100, 1000)
df_topSVGs_topHVGs$method <- as.factor(rownames(df_topSVGs_topHVGs))
df_topSVGs_topHVGs <- pivot_longer(
  df_topSVGs_topHVGs, cols = c("10", "100", "1000"), names_to = "n_HVGs", values_to = "proportion"
)
df_topSVGs_topHVGs$n_HVGs <- as.numeric(df_topSVGs_topHVGs$n_HVGs)

ggplot(df_topSVGs_topHVGs, aes(x = n_HVGs, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(10, 100, 1000), labels = c(10, 100, 1000), limits = c(10, 1000), trans = "sqrt") + 
  xlab("number of top HVGs") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top SVGs within top HVGs") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topSVGs_topHVGs")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# top SVGs within top HVGs (x2)

df_topSVGs_topHVGsx2 <- as.data.frame(rbind(
  nnSVG = c(prop_overlap_top10nnSVG_top20HVGs, prop_overlap_top100nnSVG_top200HVGs, prop_overlap_top1000nnSVG_top2000HVGs), 
  SPARKX = c(prop_overlap_top10SPARKX_top20HVGs, prop_overlap_top100SPARKX_top200HVGs, prop_overlap_top1000SPARKX_top2000HVGs)
))
colnames(df_topSVGs_topHVGsx2) <- c(20, 200, 2000)
df_topSVGs_topHVGsx2$method <- as.factor(rownames(df_topSVGs_topHVGsx2))
df_topSVGs_topHVGsx2 <- pivot_longer(
  df_topSVGs_topHVGsx2, cols = c("20", "200", "2000"), names_to = "n_HVGs", values_to = "proportion"
)
df_topSVGs_topHVGsx2$n_HVGs <- as.numeric(df_topSVGs_topHVGsx2$n_HVGs)

ggplot(df_topSVGs_topHVGsx2, aes(x = n_HVGs, y = proportion, group = method, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_startrek() + 
  ylim(c(0, 1)) + 
  scale_x_continuous(breaks = c(20, 200, 2000), labels = c(20, 200, 2000), limits = c(20, 2000), trans = "sqrt") + 
  xlab("number of top HVGs") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap top SVGs within top HVGs (x2)") + 
  theme_bw()

fn <- here("plots", "commentary", "overlap_topSVGs_topHVGsx2")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)

