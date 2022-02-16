#################################
# Script to calculate evaluations
# Lukas Weber, Feb 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
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
  geom_point(stroke = 2, size = 2) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("royalblue3", "maroon", "darkorange")) + 
  scale_y_log10(limits = range(df_known_DLPFC$rank)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.3, size = 3.5, segment.color = NA, show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("DLPFC dataset: known SVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_DLPFC"))
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
  geom_line() + 
  geom_point() + 
  scale_x_continuous(breaks = overlaps, trans = "log10") + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlapping") + 
  ggtitle("Proportion overlap between top n SVGs and HVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "prop_overlap_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

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


# plot comparisons of ranks
ggplot(as.data.frame(df_ranks_DLPFC), 
       aes(x = rank_HVGs, y = rank_method, color = method)) + 
  facet_wrap(~ method) + 
  geom_point() + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  xlab("rank HVGs") + 
  ylab("rank method for SVGs") + 
  ggtitle("Comparison of ranks SVGs vs. HVGs: DLPFC dataset") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "ranks_scatter_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)


# ---------------------
# correlations of ranks
# ---------------------

# Spearman correlations comparing ranks against HVGs

df_correlations_DLPFC <- data.frame(
  nnSVG = cor(df_ranks_DLPFC_nnSVG_HVGs$rank_method, df_ranks_DLPFC_nnSVG_HVGs$rank_HVGs, method = "spearman"), 
  SPARKX = cor(df_ranks_DLPFC_SPARKX_HVGs$rank_method, df_ranks_DLPFC_SPARKX_HVGs$rank_HVGs, method = "spearman"), 
  dataset = "DLPFC"
)

df_correlations <- 
  pivot_longer(df_correlations_DLPFC, cols = c("nnSVG", "SPARKX"), 
               names_to = "method", values_to = "correlation")


# plot correlations
ggplot(as.data.frame(df_correlations), 
       aes(x = correlation, y = dataset, color = method)) + 
  geom_point(pch = 3, stroke = 2) + 
  xlim(c(0, 1)) + 
  xlab("Spearman correlation") + 
  ggtitle("Correlation between ranks top 1000 SVGs vs. HVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "correlations_ranks_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals_DLPFC <- 
  as.data.frame(res_list$DLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes)

# plot p-values
ggplot(df_pvals_DLPFC, aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "mediumturquoise") + 
  xlab("adjusted p-values") + 
  ggtitle("P-values: nnSVG, DLPFC dataset") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "pval_DLPFC"))
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

# plot phi
ggplot(df_bandwidth_DLPFC, aes(x = phi_nnSVG)) + 
  geom_histogram(color = "royalblue", fill = "royalblue") + 
  sapply(phis_known, function(xint) geom_vline(xintercept = xint, color = "red")) + 
  geom_text(data = df_bandwidth_DLPFC %>% filter(gene_name %in% known_genes), 
            aes(y = phi_nnSVG * 50, label = gene_name), size = 3) + 
  xlab("phi") + 
  ggtitle("Estimated bandwidth: nnSVG, DLPFC dataset") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "phi_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# plot l (inverse bandwidth, on same scale as distance scaled to 0 to 1)
ggplot(df_bandwidth_DLPFC, aes(x = l_nnSVG)) + 
  geom_histogram(bins = 50, color = "royalblue", fill = "royalblue") + 
  sapply(ls_known, function(xint) geom_vline(xintercept = xint, color = "red")) + 
  geom_text(data = df_bandwidth_DLPFC %>% filter(gene_name %in% known_genes), 
            aes(y = l_nnSVG * 3000, label = gene_name), size = 3) + 
  xlim(c(0.01, 1)) + 
  xlab("l") + 
  ggtitle("Estimated bandwidth: nnSVG, DLPFC dataset") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "bandwidth_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

