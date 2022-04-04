#################################
# Script to calculate evaluations
# Lukas Weber, Apr 2022
#################################

# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "effect_sizes"))


# ------------
# load results
# ------------

# note choice of filtering per method
res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]


table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)


# ------------
# effect sizes
# ------------

# known SVGs in this dataset
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")

# layer-specific marker genes from manually guided analyses in spatialLIBD
# see original code at:
# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/SpatialDE_clustering.Rmd
sig_genes <- read_csv(here("inputs", "spatialLIBD", "sig_genes.csv"))
markers <- sig_genes[, c("ensembl", "gene")]
colnames(markers) <- c("gene_id", "gene_name")
dim(markers)
# remove duplicates (genes from multiple layers)
markers <- distinct(markers)
markers_names <- sort(markers$gene_name)
length(markers_names)
head(markers_names)
# combined set of 6 known and 198 marker (201 total)
markers_names_all <- unique(c(known_genes, markers_names))
length(markers_names_all)

# note: numbers of genes
# 198 manual marker genes
length(markers_names)
# 195 additional manual marker genes excluding overlap with 6 known genes
length(markers_names) - sum(known_genes %in% markers_names)
# 201 total manual marker genes and known genes (i.e. 3 overlapping)
length(markers_names_all)


# set up data frames

df_nnSVG <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% markers_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)

df_effect <- 
  as.data.frame(df_nnSVG) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)

# calculate adjusted effect size
df_adj_effect <- 
  full_join(as.data.frame(res_list$humanDLPFC_nnSVG), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(adj_effect_nnSVG = (sigma.sq_nnSVG - tech_HVGs) / (sigma.sq_nnSVG + tau.sq_nnSVG - tech_HVGs)) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% markers_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)


# variance vs. mean
ggplot(df_effect, 
       aes(x = mean_nnSVG, y = var_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_effect %>% filter(is_known), 
             aes(shape = is_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "known") + 
  geom_text_repel(
    data = df_effect %>% filter(is_known), 
    aes(label = gene_name), color = "red", size = 3, nudge_x = -0.1, nudge_y = 0.3) + 
  labs(x = "mean logcounts", 
       y = "variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "effectSize_var_vs_mean_humanDLPFC")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# spatial variance vs. mean
ggplot(df_effect, 
       aes(x = mean_nnSVG, y = sigma.sq_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_effect %>% filter(is_known), 
             aes(shape = is_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "known") + 
  geom_text_repel(
    data = df_effect %>% filter(is_known), 
    aes(label = gene_name), color = "red", size = 3, nudge_x = -0.1, nudge_y = 0.3) + 
  labs(x = "mean logcounts", 
       y = "spatial variance (sigma^2)", 
       color = "LR statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "effectSize_spatialVar_vs_mean_humanDLPFC")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# adjusted effect size vs. mean
ggplot(df_adj_effect, 
       aes(x = mean_nnSVG, y = adj_effect_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_adj_effect %>% filter(is_known), 
             aes(shape = is_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "known") + 
  geom_text_repel(
    data = df_adj_effect %>% filter(is_known), 
    aes(label = gene_name), color = "red", size = 3, nudge_x = 0, nudge_y = 1000) + 
  labs(x = "mean logcounts", 
       y = "adjusted effect size", 
       color = "LR statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "adjEffectSize_vs_mean_humanDLPFC")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)

