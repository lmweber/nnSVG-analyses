#################################
# Script to calculate evaluations
# Lukas Weber, Mar 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_humanDLPFC_SPARKX_nofilt.rds"))), 
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_humanDLPFC_HVGs_nofilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# nnSVG includes additional gene filtering

table(res_list$humanDLPFC_SPARKX$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_nnSVG$gene_id %in% res_list$humanDLPFC_SPARKX$gene_id)

table(res_list$humanDLPFC_HVGs$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
all(res_list$humanDLPFC_nnSVG$gene_id %in% res_list$humanDLPFC_HVGs$gene_id)


# ---------------------------------
# known SVGs in human DLPFC dataset
# ---------------------------------

# known SVGs: MOBP, PCP4, SNAP25, HBB, IGKC, NPY

known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")

df_known_humanDLPFC <- 
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
ggplot(as.data.frame(df_known_humanDLPFC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("blue3", "maroon", "darkorange")) + 
  scale_y_log10(limits = c(3, 35000)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.3, size = 1.75, segment.color = NA, show.legend = FALSE) + 
  annotate("text", label = "large bandwidth", x = 2, y = 35000, size = 4) + 
  annotate("text", label = "small bandwidth", x = 5, y = 35000, size = 4) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Example SVGs: human DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_ranks_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# --------------------------------------------
# nnSVG: likelihood ratio statistics vs. ranks
# --------------------------------------------

# layer-specific marker genes from manually guided analyses in spatialLIBD
# see original code at:
# https://github.com/LieberInstitute/HumanPilot/blob/master/Analysis/SpatialDE_clustering.Rmd
sig_genes <- read_csv(here("inputs", "spatialLIBD", "sig_genes.csv"))
manual_genes <- sig_genes[, c("ensembl", "gene")]
colnames(manual_genes) <- c("gene_id", "gene_name")
dim(manual_genes)
# remove duplicates (genes from multiple layers)
manual_genes <- distinct(manual_genes)
manual_gene_names <- sort(manual_genes$gene_name)
length(manual_gene_names)
head(manual_gene_names)

# combined set of 6 known and 198 marker (201 total)
manual_genes_all <- unique(c(known_genes, manual_gene_names))
length(manual_genes_all)


# plot likelihood ratio statistics vs. ranks

df_nnSVG_humanDLPFC <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% manual_gene_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)

# rank at p-value = 0.05 cutoff
padj_cutoff_nnSVG <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  filter(padj_nnSVG <= 0.05) %>% 
  summarize(max = max(rank_nnSVG)) %>% 
  unlist

padj_cutoff_nnSVG

# number of 190 marker genes identified as significant
table(filter(df_nnSVG_humanDLPFC, is_marker)$padj_nnSVG <= 0.05)
# number of combined set identified as significant
table(filter(df_nnSVG_humanDLPFC, is_marker_or_known)$padj_nnSVG <= 0.05)


# highlighting 6 known genes
ggplot(as.data.frame(df_nnSVG_humanDLPFC), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG_humanDLPFC, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_nnSVG_humanDLPFC, gene_name %in% known_genes), 
                  nudge_x = 80, nudge_y = 350, size = 3, color = "firebrick3") + 
  xlim(c(0, 1000)) + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: human DLPFC, example SVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "LR_stat_ranks_6known_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# highlighting 6 known genes and 193 combined set of known and layer-specific markers
# set seed for overlapping geom_text_repel
set.seed(1)
ggplot(as.data.frame(df_nnSVG_humanDLPFC), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  ylim(c(-1, 6200)) + 
  geom_point(data = filter(df_nnSVG_humanDLPFC, gene_name %in% manual_genes_all), 
             pch = 1, size = 2, color = "firebrick3") + 
  geom_point(data = filter(df_nnSVG_humanDLPFC, gene_name %in% known_genes), 
             pch = 1, size = 2, stroke = 0.75, color = "black") + 
  geom_text_repel(data = filter(df_nnSVG_humanDLPFC, gene_name %in% known_genes), 
                  nudge_x = 500, nudge_y = 800, size = 3, 
                  segment.color = "black", color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_nnSVG, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_nnSVG, ")"), 
           x = 2850, y = 5000, size = 3, color = "darkorange2") + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: human DLPFC, example SVGs and markers") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "LR_stat_ranks_193knownAndMarkers_nnSVG_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------------------------------
# SPARK-X: adjusted p-values vs. ranks
# ------------------------------------

# plot adjusted p-values vs. ranks

df_SPARKX_humanDLPFC <- 
  as.data.frame(res_list$humanDLPFC_SPARKX) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% manual_gene_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)

# rank at p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$humanDLPFC_SPARKX) %>% 
  filter(adjustedPval_SPARKX <= 0.05) %>% 
  summarize(max = max(rank_SPARKX)) %>% 
  unlist

padj_cutoff_SPARKX

# number of 190 marker genes identified as significant
table(filter(df_SPARKX_humanDLPFC, is_marker)$adjustedPval_SPARKX <= 0.05)
# number of combined set identified as significant
table(filter(df_SPARKX_humanDLPFC, is_marker_or_known)$adjustedPval_SPARKX <= 0.05)


# highlighting 193 (combined set of known and layer-specific marker)
ggplot(as.data.frame(df_SPARKX_humanDLPFC), 
       aes(x = rank_SPARKX, y = -log10(adjustedPval_SPARKX), label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_SPARKX_humanDLPFC, gene_name %in% manual_genes_all), 
             pch = 1, size = 2, color = "firebrick3") + 
  geom_point(data = filter(df_SPARKX_humanDLPFC, gene_name %in% known_genes), 
             pch = 1, size = 2, stroke = 0.75, color = "black") + 
  geom_text_repel(data = filter(df_SPARKX_humanDLPFC, gene_name %in% known_genes), 
                  nudge_x = 2000, nudge_y = 20, size = 3, 
                  segment.color = "black", color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value\n = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 14000, y = 230, size = 4, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(adjusted p-value)") + 
  ggtitle("SPARK-X: human DLPFC, example SVGs and markers") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "adjPvals_ranks_193knownAndMarkers_SPARKX_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------
# effect sizes
# ------------

# LR statistic vs. effect size

df_effect <- 
  as.data.frame(df_nnSVG_humanDLPFC) %>% 
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

fn <- here(file.path("plots", "evaluations", "effect_size_mean_humanDLPFC"))
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

df_overlaps_humanDLPFC <- data.frame(
  top_n = overlaps, 
  dataset = "humanDLPFC", 
  nnSVG = calc_overlaps("humanDLPFC_HVGs", "humanDLPFC_nnSVG"), 
  SPARKX = calc_overlaps("humanDLPFC_HVGs", "humanDLPFC_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps_humanDLPFC, cols = c("nnSVG", "SPARKX"), 
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

fn <- here(file.path("plots", "evaluations", "prop_overlap_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

df_ranks_humanDLPFC_nnSVG_HVGs <- 
  full_join(as.data.frame(res_list$humanDLPFC_nnSVG), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_nnSVG) %>% 
  mutate(method = "nnSVG") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks_humanDLPFC_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$humanDLPFC_SPARKX), 
            as.data.frame(res_list$humanDLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_SPARKX) %>% 
  mutate(method = "SPARKX") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks_humanDLPFC <- full_join(df_ranks_humanDLPFC_nnSVG_HVGs, df_ranks_humanDLPFC_SPARKX_HVGs)


# calculate Spearman correlations
cor_nnSVG <- cor(df_ranks_humanDLPFC_nnSVG_HVGs$rank_method, 
                 df_ranks_humanDLPFC_nnSVG_HVGs$rank_HVGs, method = "spearman")
cor_SPARKX <- cor(df_ranks_humanDLPFC_SPARKX_HVGs$rank_method, 
                  df_ranks_humanDLPFC_SPARKX_HVGs$rank_HVGs, method = "spearman")

ann_text <- data.frame(
  x = 820, 
  y = 50, 
  label = paste0("cor = ", c(round(cor_nnSVG, 2), round(cor_SPARKX, 2))), 
  method = factor(c("nnSVG", "SPARKX"), levels = c("nnSVG", "SPARKX"))
)


# plot comparisons of ranks
ggplot(as.data.frame(df_ranks_humanDLPFC), 
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
  ggtitle("Ranks SVGs vs. HVGs: human DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "ranks_scatter_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals_humanDLPFC <- as.data.frame(res_list$humanDLPFC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals_humanDLPFC), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("P-values: nnSVG, human DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "pvals_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# -----------------------
# bandwidth distributions
# -----------------------

df_bandwidth_humanDLPFC <- 
  as.data.frame(res_list$humanDLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)

phis <- df_bandwidth_humanDLPFC$phi_nnSVG
names(phis) <- df_bandwidth_humanDLPFC$gene_name
phis_known <- phis[known_genes]

ls_known <- 1 / phis_known
ls_known

ann_text <- data.frame(
  x = unname(ls_known), 
  y =  c(0.25, 2.25, 1.25, 6.2, 6.15, 5), 
  label = paste0(names(ls_known), " = ", format(round(ls_known, 3), nsmall = 3))
)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth_humanDLPFC), aes(x = l_nnSVG)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  geom_point(data = ann_text, aes(x = x, y = y), color = "red", size = 2) + 
  geom_text_repel(data = ann_text, aes(x = x, y = y, label = label), 
                  nudge_x = 0.15, nudge_y = 0.1, color = "red", size = 3) + 
  xlab("bandwidth (1 / phi)") + 
  ylab("density") + 
  ggtitle("Bandwidths: nnSVG, human DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "bandwidths_humanDLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

