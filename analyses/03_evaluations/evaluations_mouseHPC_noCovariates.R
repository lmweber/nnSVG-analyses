#################################
# Script to calculate evaluations
# Lukas Weber, May 2022
#################################

# data set: mouse HPC
# methods: nnSVG and SPARK-X without covariates
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
dir_plots <- here(file.path("plots", "evaluations", "mouseHPC", "no_covariates"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
res_list <- list(
  mouseHPC_nnSVG_noCovariates = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG_noCovariates.rds"))), 
  mouseHPC_SPARKX_noCovariates = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_SPARKX_noCovariates.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_HVGs_noFilt.rds"))), 
  mouseHPC_MoransI = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseHPC_nnSVG_noCovariates"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG_noCovariates"]]), "_nnSVG_noCovariates")[-1]
colnames(res_list[["mouseHPC_SPARKX_noCovariates"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX_noCovariates"]]), "_SPARKX_noCovariates")[-1]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseHPC_MoransI"]])[-1] <- paste0(colnames(res_list[["mouseHPC_MoransI"]]), "_MoransI")[-1]


# note filtering per method: same for both nnSVG and SPARK-X
table(res_list$mouseHPC_SPARKX_noCovariates$gene_name %in% res_list$mouseHPC_nnSVG_noCovariates$gene_name)
all(res_list$mouseHPC_nnSVG_noCovariates$gene_name == res_list$mouseHPC_SPARKX_noCovariates$gene_name)

table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_nnSVG_noCovariates$gene_name)
table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_SPARKX_noCovariates$gene_name)

table(res_list$mouseHPC_MoransI$gene_name %in% res_list$mouseHPC_nnSVG_noCovariates$gene_name)
table(res_list$mouseHPC_MoransI$gene_name %in% res_list$mouseHPC_SPARKX_noCovariates$gene_name)


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: Cpne9, Rgs14

known_genes <- c("Cpne9", "Rgs14")

df_known <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG_noCovariates), 
            as.data.frame(res_list$mouseHPC_SPARKX_noCovariates), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_HVGs), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_MoransI), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG_noCovariates", "rank_SPARKX_noCovariates", "rank_HVGs", "rank_MoransI"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("_noCovariates", " (no covariates)", 
                              gsub("MoransI", "Moran's I", 
                                   gsub("SPARKX", "SPARK-X", 
                                        gsub("^rank_", "", method)))), 
                         levels = c("nnSVG (no covariates)", "SPARK-X (no covariates)", "HVGs", "Moran's I"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))


# plot ranks
ggplot(as.data.frame(df_known), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 3, 1, 2)) + 
  scale_color_manual(values = c("blue3", "deepskyblue2", "darkorange", "firebrick3")) + 
  scale_y_log10(limits = c(10, 25000)) + 
  geom_text_repel(nudge_x = 0.2, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouseHPC (no covariates)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# --------------------------------------------
# nnSVG: likelihood ratio statistics vs. ranks
# --------------------------------------------

# plot likelihood ratio statistics vs. ranks

df_nnSVG <- 
  as.data.frame(res_list$mouseHPC_nnSVG_noCovariates) %>% 
  mutate(is_known = gene_name %in% known_genes)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_nnSVG <- 
  as.data.frame(res_list$mouseHPC_nnSVG_noCovariates) %>% 
  filter(padj_nnSVG_noCovariates <= 0.05) %>% 
  summarize(max = max(rank_nnSVG_noCovariates)) %>% 
  unlist

padj_cutoff_nnSVG


# highlighting known genes
ggplot(as.data.frame(df_nnSVG), 
       aes(x = rank_nnSVG_noCovariates, y = LR_stat_nnSVG_noCovariates, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_nnSVG, gene_name %in% known_genes), 
                  nudge_x = 500, nudge_y = 1200, size = 3, 
                  color = "firebrick3", fontface = "italic") + 
  geom_vline(xintercept = padj_cutoff_nnSVG, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_nnSVG, ")"), 
           x = 5500, y = 7000, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: mouseHPC (no covariates)") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_nnSVG_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

# plot combined p-values vs. ranks

df_SPARKX <- 
  as.data.frame(res_list$mouseHPC_SPARKX_noCovariates) %>% 
  mutate(is_known = gene_name %in% known_genes)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$mouseHPC_SPARKX_noCovariates) %>% 
  filter(adjustedPval_SPARKX_noCovariates <= 0.05) %>% 
  summarize(max = max(rank_SPARKX_noCovariates)) %>% 
  unlist

padj_cutoff_SPARKX


# highlighting known genes
ggplot(as.data.frame(df_SPARKX), 
       aes(x = rank_SPARKX_noCovariates, y = -log10(combinedPval_SPARKX_noCovariates), label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_SPARKX, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_SPARKX, gene_name %in% known_genes), 
                  nudge_x = 700, nudge_y = 20, size = 3, 
                  color = "firebrick3", fontface = "italic") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 7000, y = 170, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: mouseHPC (no covariates)") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ------------
# effect sizes
# ------------

# LR statistic vs. effect size

df_effect <- 
  as.data.frame(df_nnSVG) %>% 
  filter(rank_nnSVG_noCovariates <= 1000)


# effect size vs. mean
ggplot(df_effect, 
       aes(x = mean_nnSVG_noCovariates, y = prop_sv_nnSVG_noCovariates, color = LR_stat_nnSVG_noCovariates)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_effect %>% filter(is_known), 
             aes(shape = is_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "known") + 
  geom_text_repel(
    data = df_effect %>% filter(is_known), 
    aes(label = gene_name), color = "red", size = 3, fontface = "italic", 
    nudge_x = 0.4, nudge_y = 0.25) + 
  ylim(c(0, 1)) + 
  labs(x = "mean logcounts", 
       y = "proportion spatial variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: mouse HPC (no covariates)") + 
  theme_bw()

fn <- file.path(dir_plots, "effect_size_nnSVG_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseHPC_nnSVG_noCovariates)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG_noCovariates)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  ylim(c(0, 13750)) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: mouse HPC", 
          subtitle = "no covariates") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.25)
ggsave(paste0(fn, ".png"), width = 4, height = 3.25)


# --------------------------------------
# length scale (bandwidth) distributions
# --------------------------------------

df_bandwidth <- 
  as.data.frame(res_list$mouseHPC_nnSVG_noCovariates) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG_noCovariates = 1 / phi_nnSVG_noCovariates)

phis <- df_bandwidth$phi_nnSVG_noCovariates
names(phis) <- df_bandwidth$gene_name
phis_known <- phis[known_genes]

ls_known <- 1 / phis_known
ls_known

ann_text <- data.frame(
  x = unname(ls_known), 
  y =  c(1.6, 1.0), 
  label = paste0(names(ls_known), " = ", round(ls_known, 3))
)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth), aes(x = l_nnSVG_noCovariates)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  geom_point(data = ann_text, aes(x = x, y = y), color = "red", size = 2) + 
  geom_text_repel(data = ann_text, aes(x = x, y = y, label = label), 
                  color = "red", size = 3, fontface = "italic", 
                  nudge_x = 0.13, nudge_y = 3) + 
  xlab("estimated length scale") + 
  ylab("density") + 
  ggtitle("nnSVG length scales: mouseHPC (no covariates)") + 
  theme_bw()

fn <- file.path(dir_plots, "lengthscales_nnSVG_mouseHPC_noCovariates")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)

