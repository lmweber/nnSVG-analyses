#################################
# Script to calculate evaluations
# Lukas Weber, Feb 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

res_list <- list(
  DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds"))), 
  DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["DLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["DLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["DLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# ---------------------------
# known SVGs in DLPFC dataset
# ---------------------------

# known SVGs: SNAP25, MOBP, PCP4, HBB, IGKC, NPY

known_genes <- c("SNAP25", "MOBP", "PCP4", "HBB", "IGKC", "NPY")

all(res_list$DLPFC_nnSVG$gene_id == res_list$DLPFC_SPARKX$gene_id)
all(res_list$DLPFC_nnSVG$gene_id == res_list$DLPFC_HVGs$gene_id)

df_known_DLPFC <- 
  full_join(as.data.frame(res_list$DLPFC_nnSVG), 
            as.data.frame(res_list$DLPFC_SPARKX), 
            by = c("gene_id", "gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$DLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_SPARKX", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("nnSVG", "SPARKX", "HVGs"))) %>% 
  mutate(gene_name = factor(gene_name, levels = c("MOBP", "PCP4", "SNAP25", 
                                                  "HBB", "IGKC", "NPY")))


# plot ranks
ggplot(as.data.frame(df_known_DLPFC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("blue3", "maroon", "darkorange")) + 
  scale_y_log10() + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.25, size = 2, show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Known SVGs: DLPFC dataset") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_ranks_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# -------------------------------------
# likelihood ratio statistics vs. ranks
# -------------------------------------

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


# plot likelihood ratio statistics vs. ranks

df_nnSVG_DLPFC <- 
  as.data.frame(res_list$DLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% manual_gene_names)

# rank at p-value = 0.05 cutoff
padj_cutoff <- 
  as.data.frame(res_list$DLPFC_nnSVG) %>% 
  filter(padj_nnSVG <= 0.05) %>% 
  summarize(max(rank_nnSVG)) %>% 
  unlist

# number of 190 marker genes identified as significant
table(filter(df_nnSVG_DLPFC, is_marker)$padj_nnSVG <= 0.05)


# highlighting known 6 genes
ggplot(as.data.frame(df_nnSVG_DLPFC), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG_DLPFC, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_nnSVG_DLPFC, gene_name %in% known_genes), 
                  nudge_x = 80, nudge_y = 350, size = 3, color = "firebrick3") + 
  xlim(c(0, 1000)) + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: DLPFC, known SVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "LR_stat_ranks_6known_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# highlighting 190 layer-specific marker genes
ggplot(as.data.frame(df_nnSVG_DLPFC), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG_DLPFC, gene_name %in% manual_gene_names), 
             pch = 1, size = 2, color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff, linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = "adjusted p-value  = 0.05", 
           x = 6750, y = 7000, size = 4, color = "darkorange2") + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: DLPFC, layer-specific marker genes") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "LR_stat_ranks_190markers_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------
# effect sizes
# ------------

df_effect <- 
  as.data.frame(df_nnSVG_DLPFC) %>% 
  mutate(is_marker_or_known = is_marker | is_known) %>% 
  filter(rank_nnSVG <= 1000)

ggplot(df_effect, 
       aes(x = prop_sv_nnSVG, y = LR_stat_nnSVG, 
           color = is_marker_or_known, shape = is_marker_or_known)) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red"), name = "marker or\nknown") + 
  scale_shape_manual(values = c(1, 19), name = "marker or\nknown") + 
  geom_text_repel(data = df_effect %>% filter(is_known), 
                  aes(label = gene_name), nudge_y = 2000, show.legend = FALSE) + 
  labs(x = "proportion of spatial variance", 
       y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_sizes_DLPFC"))
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

df_overlaps_DLPFC <- data.frame(
  top_n = overlaps, 
  dataset = "DLPFC", 
  nnSVG = calc_overlaps("DLPFC_HVGs", "DLPFC_nnSVG"), 
  SPARKX = calc_overlaps("DLPFC_HVGs", "DLPFC_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps_DLPFC, cols = c("nnSVG", "SPARKX"), 
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
  ggtitle("Overlap SVGs and HVGs: DLPFC") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- here(file.path("plots", "evaluations", "prop_overlap_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

df_ranks_DLPFC_nnSVG_HVGs <- 
  full_join(as.data.frame(res_list$DLPFC_nnSVG), 
            as.data.frame(res_list$DLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_nnSVG) %>% 
  mutate(method = "nnSVG") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks_DLPFC_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$DLPFC_SPARKX), 
            as.data.frame(res_list$DLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(rank_method = rank_SPARKX) %>% 
  mutate(method = "SPARKX") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_id", "gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks_DLPFC <- full_join(df_ranks_DLPFC_nnSVG_HVGs, df_ranks_DLPFC_SPARKX_HVGs)


# calculate Spearman correlations
cor_nnSVG <- cor(df_ranks_DLPFC_nnSVG_HVGs$rank_method, 
                 df_ranks_DLPFC_nnSVG_HVGs$rank_HVGs, method = "spearman")
cor_SPARKX <- cor(df_ranks_DLPFC_SPARKX_HVGs$rank_method, 
                  df_ranks_DLPFC_SPARKX_HVGs$rank_HVGs, method = "spearman")

ann_text <- data.frame(
  x = 800, 
  y = 50, 
  label = paste0("cor = ", c(round(cor_nnSVG, 2), round(cor_SPARKX, 2))), 
  method = factor(c("nnSVG", "SPARKX"), levels = c("nnSVG", "SPARKX"))
)


# plot comparisons of ranks
ggplot(as.data.frame(df_ranks_DLPFC), 
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
  ggtitle("Ranks SVGs vs. HVGs: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "ranks_scatter_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals_DLPFC <- as.data.frame(res_list$DLPFC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals_DLPFC), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  xlab("p-values") + 
  ggtitle("P-value distribution: nnSVG, DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "pvals_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# -----------------------
# bandwidth distributions
# -----------------------

df_bandwidth_DLPFC <- 
  as.data.frame(res_list$DLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)

phis <- df_bandwidth_DLPFC$phi_nnSVG
names(phis) <- df_bandwidth_DLPFC$gene_name
phis_known <- phis[known_genes]

ls_known <- 1 / phis_known
ls_known

ann_text <- data.frame(
  x = unname(ls_known), 
  y =  c(0.5, 0.1, 1.5, 8, 9, 13), 
  label = paste0(names(ls_known), " = ", format(round(ls_known, 3), nsmall = 3))
)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth_DLPFC), aes(x = l_nnSVG)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  geom_point(data = ann_text, aes(x = x, y = y), color = "red", size = 2) + 
  geom_text_repel(data = ann_text, aes(x = x, y = y, label = label), 
                  nudge_x = 0.13, nudge_y = 1, color = "red", size = 3) + 
  xlab("bandwidth") + 
  ylab("density") + 
  ggtitle("Bandwidth: nnSVG, DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "bandwidths_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

