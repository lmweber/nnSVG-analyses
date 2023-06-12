#################################
# Script to calculate evaluations
# Lukas Weber, updated Jun 2023
#################################

# data set: mouse embryo
# filtering: no filtering of low-expressed genes (note: this dataset contains a targeted gene set)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "evaluations", "mouseEmbryo", "no_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
res_list <- list(
  mouseEmbryo_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_nnSVG_noFilt.rds"))), 
  mouseEmbryo_SPARKX = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_SPARKX_noFilt.rds"))), 
  mouseEmbryo_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_HVGs_noFilt.rds"))), 
  mouseEmbryo_MoransI = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseEmbryo_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseEmbryo_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseEmbryo_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseEmbryo_MoransI"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_MoransI"]]), "_MoransI")[-1]


# note filtering per method: no filtering

table(res_list$mouseEmbryo_SPARKX$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)
all(res_list$mouseEmbryo_nnSVG$gene_name == res_list$mouseEmbryo_SPARKX$gene_name)

table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)
table(res_list$mouseEmbryo_HVGs$gene_name %in% res_list$mouseEmbryo_SPARKX$gene_name)

table(res_list$mouseEmbryo_MoransI$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)
table(res_list$mouseEmbryo_MoransI$gene_name %in% res_list$mouseEmbryo_SPARKX$gene_name)


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: Ttn, Popdc2, Hand1, Gata5, Six3, Lhx2, Otx2, Pou3f1, Sox2, Foxf1, Foxa1, Cldn4

known_genes <- c("Ttn", "Popdc2", "Hand1", "Gata5", "Six3", "Lhx2", 
                 "Otx2", "Pou3f1", "Sox2", "Foxf1", "Foxa1", "Cldn4")

df_known <- 
  full_join(as.data.frame(res_list$mouseEmbryo_nnSVG), 
            as.data.frame(res_list$mouseEmbryo_SPARKX), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseEmbryo_HVGs), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseEmbryo_MoransI), 
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
  scale_y_log10(limits = c(3, 300)) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse embryo") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4)


# --------------------------------------------
# nnSVG: likelihood ratio statistics vs. ranks
# --------------------------------------------

# plot likelihood ratio statistics vs. ranks

df_nnSVG <- 
  as.data.frame(res_list$mouseEmbryo_nnSVG)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_nnSVG <- 
  as.data.frame(res_list$mouseEmbryo_nnSVG) %>% 
  filter(padj_nnSVG <= 0.05) %>% 
  summarize(max = max(rank_nnSVG)) %>% 
  unlist

padj_cutoff_nnSVG

# number of significant SVGs
table(df_nnSVG$padj_nnSVG <= 0.05)


# all genes
ggplot(as.data.frame(df_nnSVG), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG)) + 
  geom_line(color = "navy") + 
  geom_vline(xintercept = padj_cutoff_nnSVG, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_nnSVG, ")"), 
           x = 240, y = 7500, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_all_nnSVG_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

# plot combined p-values vs. ranks

df_SPARKX <- 
  as.data.frame(res_list$mouseEmbryo_SPARKX)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$mouseEmbryo_SPARKX) %>% 
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
           x = 240, y = 200, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ------------
# effect sizes
# ------------

# LR statistic vs. effect size

df_effect <- 
  as.data.frame(df_nnSVG) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)


# effect size vs. mean
ggplot(df_effect, 
       aes(x = mean_nnSVG, y = prop_sv_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  ylim(c(0, 1)) + 
  labs(x = "mean logcounts", 
       y = "proportion spatial variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "effect_size_nnSVG_mouseEmbryo_noFilt")
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
  dataset = "mouseEmbryo", 
  nnSVG = calc_overlaps("mouseEmbryo_HVGs", "mouseEmbryo_nnSVG"), 
  SPARKX = calc_overlaps("mouseEmbryo_HVGs", "mouseEmbryo_SPARKX")
)

df_overlaps <- 
  pivot_longer(df_overlaps, cols = c("nnSVG", "SPARKX"), 
               names_to = "method", values_to = "proportion") %>% 
  mutate(method = gsub("SPARKX", "SPARK-X", method))


# plot overlaps
ggplot(as.data.frame(df_overlaps), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  geom_line(lwd = 0.75) + 
  geom_point(size = 2) + 
  scale_color_manual(values = c("blue3", "deepskyblue2")) + 
  scale_x_continuous(breaks = overlaps, trans = "log10") + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlap") + 
  ggtitle("Overlap SVGs and HVGs: mouse embryo") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "overlaps_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------------
# scatterplot comparing ranks
# ---------------------------

# top 1000 ranked genes from each method

# nnSVG

df_ranks_nnSVG_HVGs <- 
  full_join(as.data.frame(res_list$mouseEmbryo_nnSVG), 
            as.data.frame(res_list$mouseEmbryo_HVGs), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_HVGs) %>% 
  mutate(baseline = "HVGs") %>% 
  filter(rank_nnSVG <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_nnSVG", "rank_baseline", "baseline"))

df_ranks_nnSVG_MoransI <- 
  full_join(as.data.frame(res_list$mouseEmbryo_nnSVG), 
            as.data.frame(res_list$mouseEmbryo_MoransI), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_MoransI) %>% 
  mutate(baseline = "Moran's I") %>% 
  filter(rank_nnSVG <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_nnSVG", "rank_baseline", "baseline"))

df_ranks_nnSVG <- full_join(df_ranks_nnSVG_HVGs, df_ranks_nnSVG_MoransI)


# SPARK-X

df_ranks_SPARKX_HVGs <- 
  full_join(as.data.frame(res_list$mouseEmbryo_SPARKX), 
            as.data.frame(res_list$mouseEmbryo_HVGs), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_HVGs) %>% 
  mutate(baseline = "HVGs") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX_MoransI <- 
  full_join(as.data.frame(res_list$mouseEmbryo_SPARKX), 
            as.data.frame(res_list$mouseEmbryo_MoransI), 
            by = c("gene_name")) %>% 
  mutate(rank_baseline = rank_MoransI) %>% 
  mutate(baseline = "Moran's I") %>% 
  filter(rank_SPARKX <= 1000) %>% 
  filter(rank_baseline <= 1000) %>% 
  select(c("gene_name", "rank_SPARKX", "rank_baseline", "baseline"))

df_ranks_SPARKX <- full_join(df_ranks_SPARKX_HVGs, df_ranks_SPARKX_MoransI)


# calculate Spearman correlations

cor_nnSVG_HVGs <- cor(df_ranks_nnSVG_HVGs$rank_nnSVG, 
                      df_ranks_nnSVG_HVGs$rank_baseline, method = "spearman")
cor_nnSVG_MoransI <- cor(df_ranks_nnSVG_MoransI$rank_nnSVG, 
                         df_ranks_nnSVG_MoransI$rank_baseline, method = "spearman")

cor_SPARKX_HVGs <- cor(df_ranks_SPARKX_HVGs$rank_SPARKX, 
                       df_ranks_SPARKX_HVGs$rank_baseline, method = "spearman")
cor_SPARKX_MoransI <- cor(df_ranks_SPARKX_MoransI$rank_SPARKX, 
                          df_ranks_SPARKX_MoransI$rank_baseline, method = "spearman")


ann_text_nnSVG <- data.frame(
  x = 275, 
  y = 25, 
  label = paste0("cor = ", c(round(cor_nnSVG_HVGs, 2), round(cor_nnSVG_MoransI, 2))), 
  baseline = factor(c("HVGs", "Moran's I"), levels = c("HVGs", "Moran's I"))
)

ann_text_SPARKX <- data.frame(
  x = 275, 
  y = 25, 
  label = paste0("cor = ", c(round(cor_SPARKX_HVGs, 2), round(cor_SPARKX_MoransI, 2))), 
  baseline = factor(c("HVGs", "Moran's I"), levels = c("HVGs", "Moran's I"))
)


# plot comparisons of ranks

# nnSVG
ggplot(as.data.frame(df_ranks_nnSVG), 
       aes(x = rank_baseline, y = rank_nnSVG, color = baseline)) + 
  facet_wrap(~ baseline) + 
  geom_point(size = 0.75) + 
  geom_text(data = ann_text_nnSVG, aes(x = x, y = y, label = label), 
            size = 3.75, color = "black") + 
  scale_color_manual(values = c("darkorange", "firebrick3")) + 
  coord_fixed() + 
  xlim(c(0, 351)) + 
  ylim(c(0, 351)) + 
  xlab("rank baseline") + 
  ylab("rank nnSVG") + 
  ggtitle("Ranks nnSVG vs. baselines: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_nnSVG_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 2.75)
ggsave(paste0(fn, ".png"), width = 5.25, height = 2.75)


# SPARK-X
ggplot(as.data.frame(df_ranks_SPARKX), 
       aes(x = rank_baseline, y = rank_SPARKX, color = baseline)) + 
  facet_wrap(~ baseline) + 
  geom_point(size = 0.75) + 
  geom_text(data = ann_text_SPARKX, aes(x = x, y = y, label = label), 
            size = 3.75, color = "black") + 
  scale_color_manual(values = c("darkorange", "firebrick3")) + 
  coord_fixed() + 
  xlim(c(0, 351)) + 
  ylim(c(0, 351)) + 
  xlab("rank baseline") + 
  ylab("rank SPARK-X") + 
  ggtitle("Ranks SPARK-X vs. baselines: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_SPARKX_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 2.75)
ggsave(paste0(fn, ".png"), width = 5.25, height = 2.75)


# --------------------------------------
# length scale (bandwidth) distributions
# --------------------------------------

df_bandwidth <- 
  as.data.frame(res_list$mouseEmbryo_nnSVG) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth), aes(x = l_nnSVG)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  xlab("estimated length scale") + 
  ylab("density") + 
  ggtitle("nnSVG length scales: mouse embryo") + 
  theme_bw()

fn <- file.path(dir_plots, "lengthscales_nnSVG_mouseEmbryo_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S7BCDEFGHI.rds"

saveRDS(res_list, here(dir_sd, fn_sd))

