#################################
# Script to calculate evaluations
# Lukas Weber, updated Jun 2023
#################################

# data set: human DLPFC
# filtering: no filtering of low-expressed genes (i.e. default for SPARK-X)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "evaluations", "humanDLPFC", "no_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG_noFilt.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX_noFilt.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds"))), 
  humanDLPFC_MoransI = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["humanDLPFC_MoransI"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_MoransI"]]), "_MoransI")[-(1:2)]


# note filtering per method: no filtering for either nnSVG or SPARK-X
table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_SPARKX$gene_id == res_list$humanDLPFC_nnSVG$gene_id)

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)

table(res_list$humanDLPFC_MoransI$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$humanDLPFC_MoransI$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S2ACF_S4A_S19D.RData"

save(res_list, file = here(dir_sd, fn_sd))


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
  full_join(., 
            as.data.frame(res_list$humanDLPFC_MoransI), 
            by = c("gene_id", "gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_SPARKX", "rank_HVGs", "rank_MoransI"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("MoransI", "Moran's I", 
                              gsub("SPARKX", "SPARK-X", 
                                   gsub("^rank_", "", method))), 
                         levels = c("nnSVG", "SPARK-X", "HVGs", "Moran's I"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))


# plot ranks
ggplot(as.data.frame(df_known), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 3, 1, 2)) + 
  scale_color_manual(values = c("blue3", "deepskyblue2", "darkorange", "firebrick3")) + 
  scale_y_log10(limits = c(3, 35000)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  annotate("text", label = "large length scale", x = 2, y = 35000, size = 4) + 
  annotate("text", label = "small length scale", x = 5, y = 35000, size = 4) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: human DLPFC (no filtering)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

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
                  nudge_x = 2000, nudge_y = 20, size = 3, 
                  segment.color = "black", color = "firebrick3", 
                  fontface = "italic") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 15000, y = 250, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: human DLPFC (no filtering)") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


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
    top_method2[k] <- nrow(filter(as.data.frame(na.omit(res_method2[genes_k, ])), 
                                  rank <= overlaps[k]))
  }
  
  # calculate proportions
  return(top_method2 / overlaps)
}

# use function to calculate overlaps

df_overlaps <- data.frame(
  top_n = overlaps, 
  dataset = "humanDLPFC", 
  SPARKX = calc_overlaps("humanDLPFC_HVGs", "humanDLPFC_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps, cols = "SPARKX", 
               names_to = "method", values_to = "proportion") %>% 
  mutate(method = gsub("SPARKX", "SPARK-X", method))


# plot overlaps
ggplot(as.data.frame(df_overlaps), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  geom_line(lwd = 0.75) + 
  geom_point(size = 2) + 
  scale_color_manual(values = "deepskyblue2") + 
  scale_x_continuous(breaks = overlaps, trans = "log10") + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap SVGs and HVGs: human DLPFC (no filtering)") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "overlaps_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

# SPARK-X

df_ranks_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$humanDLPFC_SPARKX), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_baseline = rank_HVGs) %>% 
  mutate(baseline = "HVGs") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX_MoransI <- 
  full_join(as.data.frame(res_list$humanDLPFC_SPARKX), 
            as.data.frame(res_list$humanDLPFC_MoransI), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_baseline = rank_MoransI) %>% 
  mutate(baseline = "Moran's I") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX <- full_join(df_ranks_SPARKX_HVGs, df_ranks_SPARKX_MoransI)


# calculate Spearman correlations

cor_SPARKX_HVGs <- cor(df_ranks_SPARKX_HVGs$rank_SPARKX, 
                       df_ranks_SPARKX_HVGs$rank_baseline, method = "spearman")
cor_SPARKX_MoransI <- cor(df_ranks_SPARKX_MoransI$rank_SPARKX, 
                          df_ranks_SPARKX_MoransI$rank_baseline, method = "spearman")


ann_text_SPARKX <- data.frame(
  x = 220, 
  y = 950, 
  label = paste0("cor = ", c(round(cor_SPARKX_HVGs, 2), round(cor_SPARKX_MoransI, 2))), 
  baseline = factor(c("HVGs", "Moran's I"), levels = c("HVGs", "Moran's I"))
)


# plot comparisons of ranks

# SPARK-X
# seed for placement of text labels
set.seed(1)
ggplot(as.data.frame(df_ranks_SPARKX), 
       aes(x = rank_baseline, y = rank_SPARKX, color = baseline)) + 
  facet_wrap(~ baseline) + 
  geom_point(size = 0.75) + 
  geom_text(data = ann_text_SPARKX, aes(x = x, y = y, label = label), 
            size = 3.75, color = "black") + 
  scale_color_manual(values = c("darkorange", "firebrick3")) + 
  geom_point(
    data = df_ranks_SPARKX %>% filter(gene_name %in% known_genes), 
    color = "black", size = 1) + 
  geom_text_repel(
    data = df_ranks_SPARKX %>% filter(gene_name %in% known_genes), 
    aes(label = gene_name), color = "black", size = 3.25, fontface = "italic", 
    nudge_x = 100, nudge_y = 100, box.padding = 0.5) + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  xlab("rank baseline") + 
  ylab("rank SPARK-X") + 
  ggtitle("Ranks SPARK-X vs. baselines: human DLPFC (no filt.)") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_SPARKX_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 2.75)
ggsave(paste0(fn, ".png"), width = 5.25, height = 2.75)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$humanDLPFC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  ylim(c(0, 10500)) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: human DLPFC", 
          subtitle = "no filtering") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4, height = 3.5)


# --------------------------------------
# length scale (bandwidth) distributions
# --------------------------------------

df_bandwidth <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)


# boxplot of ranks for genes with very low bandwidth l

threshold <- 0.01
n_small <- sum(df_bandwidth$l_nnSVG < threshold)
n_large <- sum(df_bandwidth$l_nnSVG >= threshold)
levs <- paste0(c("small (n = ", "large (n = "), c(n_small, n_large), c(")", ")"))

df <- df_bandwidth %>% 
  select(c("gene_id", "gene_name", "l_nnSVG", "rank_nnSVG")) %>% 
  mutate(small_lengthscale = 
           factor(ifelse(l_nnSVG < threshold, levs[1], levs[2]), levels = levs))

# number of genes and range of ranks
table(df$small_lengthscale)
summary(df$rank_nnSVG[df$small_lengthscale == levs[1]])

ggplot(as.data.frame(df), aes(x = small_lengthscale, y = rank_nnSVG, 
                              fill = small_lengthscale)) + 
  geom_violin(width = 0.5) + 
  geom_boxplot(width = 0.1, color = "black", alpha = 0.2, show.legend = FALSE) + 
  scale_fill_manual(values = c("purple3", "darkgoldenrod1")) + 
  ylim(c(0, 3396)) + 
  labs(x = "length scale", 
       y = "rank", 
       fill = "length scale") + 
  ggtitle("nnSVG length scales: human DLPFC", 
          subtitle = "Small (< 0.01) vs. large (>= 0.01) length scales") + 
  theme_bw()

fn <- file.path(dir_plots, "lengthscales_nnSVG_smallVsLarge_humanDLPFC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)

