#################################
# Script to calculate evaluations
# Lukas Weber, Mar 2022
#################################

# data set: human DLPFC
# gene filtering: same for both nnSVG and SPARK-X


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "evaluations", "humanDLPFC", "same_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

# note choice of filtering per method
res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_nofilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# note filtering per method: same for both nnSVG and SPARK-X

table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_nnSVG$gene_id == res_list$humanDLPFC_SPARKX$gene_id)

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: MOBP, PCP4, SNAP25, HBB, IGKC, NPY

known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")

df_known <- 
  full_join(as.data.frame(res_list$humanDLPFC_nnSVG), 
            as.data.frame(res_list$humanDLPFC_SPARKX), 
            by = c("gene_id", "gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_SPARKX", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("nnSVG", "SPARKX", "HVGs"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))


# plot ranks
ggplot(as.data.frame(df_known), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("blue3", "maroon", "darkorange")) + 
  scale_y_log10(limits = c(3, 6000)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.3, size = 1.75, segment.color = NA, show.legend = FALSE) + 
  annotate("text", label = "large length scale", x = 2, y = 6000, size = 4) + 
  annotate("text", label = "small length scale", x = 5, y = 6000, size = 4) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Example SVGs: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "example_SVGs_ranks_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# --------------------------------------------
# nnSVG: likelihood ratio statistics vs. ranks
# --------------------------------------------

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


# plot likelihood ratio statistics vs. ranks

df_nnSVG <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% markers_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_nnSVG <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  filter(padj_nnSVG <= 0.05) %>% 
  summarize(max = max(rank_nnSVG)) %>% 
  unlist

padj_cutoff_nnSVG

# number of significant SVGs
table(df_nnSVG$padj_nnSVG <= 0.05)

# number of markers and known in nnSVG gene set (i.e. pass filtering)
table(df_nnSVG$gene_name %in% markers_names_all)
# number of marker genes identified as significant
table(filter(df_nnSVG, is_marker)$padj_nnSVG <= 0.05)
# number of combined set identified as significant
table(filter(df_nnSVG, is_marker_or_known)$padj_nnSVG <= 0.05)


# highlighting known genes
ggplot(as.data.frame(df_nnSVG), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_nnSVG, gene_name %in% known_genes), 
                  nudge_x = 80, nudge_y = 350, size = 3, color = "firebrick3") + 
  xlim(c(0, 1000)) + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_top1000_nnSVG_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# highlighting known and marker genes
# set seed for overlapping geom_text_repel
set.seed(1)
ggplot(as.data.frame(df_nnSVG), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  ylim(c(-1, 6200)) + 
  geom_point(data = filter(df_nnSVG, gene_name %in% markers_names_all), 
             pch = 1, size = 2, color = "firebrick3") + 
  geom_point(data = filter(df_nnSVG, gene_name %in% known_genes), 
             pch = 1, size = 2, stroke = 0.75, color = "black") + 
  geom_text_repel(data = filter(df_nnSVG, gene_name %in% known_genes), 
                  nudge_x = 500, nudge_y = 800, size = 3, 
                  segment.color = "black", color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_nnSVG, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_nnSVG, ")"), 
           x = 2850, y = 5000, size = 3, color = "darkorange2") + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_markers_nnSVG_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

# plot combined p-values vs. ranks

df_SPARKX <- 
  as.data.frame(res_list$humanDLPFC_SPARKX) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% markers_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$humanDLPFC_SPARKX) %>% 
  filter(adjustedPval_SPARKX <= 0.05) %>% 
  summarize(max = max(rank_SPARKX)) %>% 
  unlist

padj_cutoff_SPARKX

# number of significant SVGs
table(df_SPARKX$adjustedPval_SPARKX <= 0.05)

# number of manual and known in SPARK-X gene set
table(df_SPARKX$gene_name %in% markers_names_all)
# number of marker genes identified as significant
table(filter(df_SPARKX, is_marker)$adjustedPval_SPARKX <= 0.05)
# number of combined set identified as significant
table(filter(df_SPARKX, is_marker_or_known)$adjustedPval_SPARKX <= 0.05)


# highlighting combined set of known and layer-specific markers
ggplot(as.data.frame(df_SPARKX), 
       aes(x = rank_SPARKX, y = -log10(combinedPval_SPARKX), label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_SPARKX, gene_name %in% markers_names_all), 
             pch = 1, size = 2, color = "firebrick3") + 
  geom_point(data = filter(df_SPARKX, gene_name %in% known_genes), 
             pch = 1, size = 2, stroke = 0.75, color = "black") + 
  geom_text_repel(data = filter(df_SPARKX, gene_name %in% known_genes), 
                  nudge_x = 200, nudge_y = 40, size = 3, 
                  segment.color = "black", color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 2800, y = 250, size = 3, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------
# effect sizes
# ------------

# LR statistic vs. effect size

df_effect <- 
  as.data.frame(df_nnSVG) %>% 
  mutate(is_marker_or_known = is_marker | is_known) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)


# effect size vs. mean
ggplot(df_effect, 
       aes(x = mean_nnSVG, y = prop_sv_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_effect %>% filter(is_known), 
             aes(shape = is_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "known") + 
  geom_text_repel(
    data = df_effect %>% filter(is_known), 
    aes(label = gene_name), color = "red", size = 3, nudge_x = 0.4, nudge_y = 0.05) + 
  ylim(c(0, 1)) + 
  labs(x = "mean logcounts", 
       y = "proportion spatial variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "effect_size_nnSVG_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# --------------------------------------------------------------------
# proportion overlap between top n SVGs for each method and top n HVGs
# --------------------------------------------------------------------

# overlap sizes
overlaps <- c(10, 20, 50, 100, 200)

# function to calculate overlaps for a pair of methods, i.e. proportion of top n
# genes from method 1 (e.g. HVGs) that are also in the set of top n genes from
# method 2 (e.g. nnSVG)
calc_overlaps <- function(method1, method2) {
  
  res_method1 <- res_list[[method1]]
  res_method2 <- res_list[[method2]]
  
  # remove method names from column names
  colnames(res_method1)[-(1:2)] <- gsub("_.*$", "", colnames(res_method1)[-(1:2)])
  colnames(res_method2)[-(1:2)] <- gsub("_.*$", "", colnames(res_method2)[-(1:2)])
  
  top_method2 <- rep(NA, length(overlaps))
  
  for (k in seq_along(overlaps)) {
    # select top gene IDs from method 1
    genes_k <- rownames(filter(as.data.frame(res_method1), 
                               rank <= overlaps[k]))
    # calculate overlaps
    top_method2[k] <- nrow(filter(as.data.frame(res_method2[genes_k, ]), 
                                  rank <= overlaps[k]))
  }
  
  # calculate proportions
  return(top_method2 / overlaps)
}

# use function to calculate overlaps

df_overlaps <- data.frame(
  top_n = overlaps, 
  dataset = "humanDLPFC", 
  nnSVG = calc_overlaps("humanDLPFC_HVGs", "humanDLPFC_nnSVG"), 
  SPARKX = calc_overlaps("humanDLPFC_HVGs", "humanDLPFC_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps, cols = c("nnSVG", "SPARKX"), 
               names_to = "method", values_to = "proportion")


# plot overlaps
ggplot(as.data.frame(df_overlaps), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  geom_line(lwd = 0.75) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("blue3", "maroon")) + 
  scale_x_continuous(breaks = overlaps, trans = "log10") + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap SVGs and HVGs: human DLPFC") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "overlaps_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

df_ranks_nnSVG_HVGs <- 
  full_join(as.data.frame(res_list$humanDLPFC_nnSVG), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_nnSVG) %>% 
  mutate(method = "nnSVG") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$humanDLPFC_SPARKX), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_SPARKX) %>% 
  mutate(method = "SPARKX") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks <- full_join(df_ranks_nnSVG_HVGs, df_ranks_SPARKX_HVGs)


# calculate Spearman correlations
cor_nnSVG <- cor(df_ranks_nnSVG_HVGs$rank_method, 
                 df_ranks_nnSVG_HVGs$rank_HVGs, method = "spearman")
cor_SPARKX <- cor(df_ranks_SPARKX_HVGs$rank_method, 
                  df_ranks_SPARKX_HVGs$rank_HVGs, method = "spearman")

ann_text <- data.frame(
  x = 820, 
  y = 50, 
  label = paste0("cor = ", c(round(cor_nnSVG, 2), round(cor_SPARKX, 2))), 
  method = factor(c("nnSVG", "SPARKX"), levels = c("nnSVG", "SPARKX"))
)


# plot comparisons of ranks
ggplot(as.data.frame(df_ranks), 
       aes(x = rank_HVGs, y = rank_method, color = method)) + 
  facet_wrap(~ method) + 
  geom_point() + 
  geom_text(data = ann_text, aes(x = x, y = y, label = label), 
            size = 5, color = "black") + 
  scale_color_manual(values = c("blue3", "maroon")) + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  xlab("rank HVGs") + 
  ylab("rank SVGs") + 
  ggtitle("Ranks SVGs and HVGs: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$humanDLPFC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# --------------------------------------
# length scale (bandwidth) distributions
# --------------------------------------

df_bandwidth <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)

phis <- df_bandwidth$phi_nnSVG
names(phis) <- df_bandwidth$gene_name
phis_known <- phis[known_genes]

ls_known <- 1 / phis_known
ls_known

ann_text <- data.frame(
  x = unname(ls_known), 
  y =  c(0.25, 2.25, 1.25, 6.2, 6.15, 5), 
  label = paste0(names(ls_known), " = ", format(round(ls_known, 3), nsmall = 3))
)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth), aes(x = l_nnSVG)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  geom_point(data = ann_text, aes(x = x, y = y), color = "red", size = 2) + 
  geom_text_repel(data = ann_text, aes(x = x, y = y, label = label), 
                  nudge_x = 0.15, nudge_y = 0.1, color = "red", size = 3) + 
  xlab("estimated length scale") + 
  ylab("density") + 
  ggtitle("nnSVG length scales: human DLPFC") + 
  theme_bw()

fn <- file.path(dir_plots, "lengthscales_nnSVG_humanDLPFC_sameFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

