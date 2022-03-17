#################################
# Script to calculate evaluations
# Lukas Weber, Mar 2022
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

# scalable methods: nnSVG, SPARK-X, HVGs

# note choice of filtering per method
res_list <- list(
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG_noFilt.rds"))), 
  mouseOB_SPARKX = rowData(readRDS(here("outputs", "results", "spe_mouseOB_SPARKX_noFilt.rds"))), 
  mouseOB_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseOB_HVGs_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseOB_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseOB_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseOB_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseOB_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseOB_HVGs"]]), "_HVGs")[-1]


# note filtering per method: no filtering for either nnSVG or SPARK-X
table(res_list$mouseOB_SPARKX$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
all(res_list$mouseOB_SPARKX$gene_name == res_list$mouseOB_nnSVG$gene_name)

table(res_list$mouseOB_HVGs$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
table(res_list$mouseOB_HVGs$gene_name %in% res_list$mouseOB_SPARKX$gene_name)


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
           x = 5000, y = 11, size = 3, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: mouse OB") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_mouseOB_noFilt")
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
               names_to = "method", values_to = "proportion")


# plot overlaps
ggplot(as.data.frame(df_overlaps), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  geom_line(lwd = 0.75) + 
  geom_point(size = 2) + 
  scale_color_manual(values = "maroon") + 
  scale_x_continuous(breaks = overlaps, trans = "log10") + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap SVGs and HVGs: mouse OB") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "overlaps_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

df_ranks_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$mouseOB_SPARKX), 
            as.data.frame(res_list$mouseOB_HVGs), 
            by = "gene_name") %>% 
  mutate(rank_method = rank_SPARKX) %>% 
  mutate(method = "SPARKX") %>% 
  filter(rank_method <= 1000) %>% 
  filter(rank_HVGs <= 1000) %>% 
  select(c("gene_name", "rank_HVGs", "rank_method", "method"))

df_ranks <- df_ranks_SPARKX_HVGs


# calculate Spearman correlations
cor_SPARKX <- cor(df_ranks_SPARKX_HVGs$rank_method, 
                  df_ranks_SPARKX_HVGs$rank_HVGs, method = "spearman")

ann_text <- data.frame(
  x = 820, 
  y = 50, 
  label = paste0("cor = ", round(cor_SPARKX, 2)), 
  method = factor("SPARKX")
)


# plot comparisons of ranks
ggplot(as.data.frame(df_ranks), 
       aes(x = rank_HVGs, y = rank_method, color = method)) + 
  facet_wrap(~ method) + 
  geom_point() + 
  geom_text(data = ann_text, aes(x = x, y = y, label = label), 
            size = 5, color = "black") + 
  scale_color_manual(values = "maroon") + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  xlab("rank HVGs") + 
  ylab("rank SVGs") + 
  ggtitle("Ranks SVGs and HVGs: mouse OB") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 4)
ggsave(paste0(fn, ".png"), width = 4.75, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseOB_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: mouse OB") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseOB_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

