#################################
# Script to calculate evaluations
# Lukas Weber, updated Jun 2023
#################################

# data set: mouse OB
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
dir_plots <- here(file.path("plots", "evaluations", "mouseOB", "no_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
res_list <- list(
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG_noFilt.rds"))), 
  mouseOB_SPARKX = rowData(readRDS(here("outputs", "results", "spe_mouseOB_SPARKX_noFilt.rds"))), 
  mouseOB_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseOB_HVGs_noFilt.rds"))), 
  mouseOB_MoransI = rowData(readRDS(here("outputs", "results", "spe_mouseOB_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseOB_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseOB_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseOB_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseOB_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseOB_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseOB_MoransI"]])[-1] <- paste0(colnames(res_list[["mouseOB_MoransI"]]), "_MoransI")[-1]


# note filtering per method: no filtering for either nnSVG or SPARK-X
table(res_list$mouseOB_SPARKX$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
all(res_list$mouseOB_SPARKX$gene_name == res_list$mouseOB_nnSVG$gene_name)

table(res_list$mouseOB_HVGs$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
table(res_list$mouseOB_HVGs$gene_name %in% res_list$mouseOB_SPARKX$gene_name)

table(res_list$mouseOB_MoransI$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
table(res_list$mouseOB_MoransI$gene_name %in% res_list$mouseOB_SPARKX$gene_name)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S19E.RData"

save(res_list, file = here(dir_sd, fn_sd))


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: Penk, Doc2g, Kctd12, Slc17a7, Cdhr1, Sv2b, Shisa3

known_genes <- c("Penk", "Doc2g", "Kctd12", "Slc17a7", "Cdhr1", "Sv2b", "Shisa3")

df_known <- 
  full_join(as.data.frame(res_list$mouseOB_nnSVG), 
            as.data.frame(res_list$mouseOB_SPARKX), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseOB_HVGs), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseOB_MoransI), 
            by = c("gene_name")) %>% 
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
  scale_y_log10(limits = c(6, 17000)) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouseOB (no filtering)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

# plot combined p-values vs. ranks

df_SPARKX <- 
  as.data.frame(res_list$mouseOB_SPARKX)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$mouseOB_SPARKX) %>% 
  filter(adjustedPval_SPARKX <= 0.05) %>% 
  summarize(max = max(rank_SPARKX)) %>% 
  unlist

padj_cutoff_SPARKX

# number of significant SVGs
table(df_SPARKX$adjustedPval_SPARKX <= 0.05)


# all genes
ggplot(as.data.frame(df_SPARKX), 
       aes(x = rank_SPARKX, y = -log10(combinedPval_SPARKX))) + 
  geom_line(color = "navy") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 5500, y = 11, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: mouse OB (no filtering)") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_mouseOB_noFilt")
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
  colnames(res_method1)[-1] <- gsub("_.*$", "", colnames(res_method1)[-1])
  colnames(res_method2)[-1] <- gsub("_.*$", "", colnames(res_method2)[-1])
  
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
  dataset = "mouseOB", 
  SPARKX = calc_overlaps("mouseOB_HVGs", "mouseOB_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps, "SPARKX", 
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
  ggtitle("Overlap SVGs and HVGs: mouse OB (no filtering)") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "overlaps_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

# SPARK-X

df_ranks_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$mouseOB_SPARKX), 
            as.data.frame(res_list$mouseOB_HVGs), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_HVGs) %>% 
  mutate(baseline = "HVGs") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX_MoransI <- 
  full_join(as.data.frame(res_list$mouseOB_SPARKX), 
            as.data.frame(res_list$mouseOB_MoransI), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_MoransI) %>% 
  mutate(baseline = "Moran's I") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX <- full_join(df_ranks_SPARKX_HVGs, df_ranks_SPARKX_MoransI)


# calculate Spearman correlations

cor_SPARKX_HVGs <- cor(df_ranks_SPARKX_HVGs$rank_SPARKX, 
                       df_ranks_SPARKX_HVGs$rank_baseline, method = "spearman")
cor_SPARKX_MoransI <- cor(df_ranks_SPARKX_MoransI$rank_SPARKX, 
                          df_ranks_SPARKX_MoransI$rank_baseline, method = "spearman")


ann_text_SPARKX <- data.frame(
  x = 800, 
  y = 50, 
  label = paste0("cor = ", c(round(cor_SPARKX_HVGs, 2), round(cor_SPARKX_MoransI, 2))), 
  baseline = factor(c("HVGs", "Moran's I"), levels = c("HVGs", "Moran's I"))
)


# plot comparisons of ranks

# SPARK-X
ggplot(as.data.frame(df_ranks_SPARKX), 
       aes(x = rank_baseline, y = rank_SPARKX, color = baseline)) + 
  facet_wrap(~ baseline) + 
  geom_point(size = 0.75) + 
  geom_text(data = ann_text_SPARKX, aes(x = x, y = y, label = label), 
            size = 3.75, color = "black") + 
  scale_color_manual(values = c("darkorange", "firebrick3")) + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  xlab("rank baseline") + 
  ylab("rank SPARK-X") + 
  ggtitle("Ranks SPARK-X vs. baselines: mouse OB (no filtering)") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_SPARKX_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 2.75)
ggsave(paste0(fn, ".png"), width = 5.25, height = 2.75)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseOB_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  ylim(c(0, 8500)) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: mouse OB", 
          subtitle = "no filtering") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4, height = 3.5)

