#################################
# Script to calculate evaluations
# Lukas Weber, Mar 2022
#################################

# data set: mouse HPC
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
dir_plots <- here(file.path("plots", "evaluations", "mouseHPC", "with_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

# note choice of filtering per method
res_list <- list(
  mouseHPC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG.rds"))), 
  mouseHPC_nnSVG_noCovariates = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG_noCovariates.rds"))), 
  mouseHPC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_SPARKX.rds"))), 
  mouseHPC_SPARKX_noCovariates = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_SPARKX_noCovariates.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_HVGs_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseHPC_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseHPC_nnSVG_noCovariates"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG_noCovariates"]]), "_nnSVG_noCovariates")[-1]
colnames(res_list[["mouseHPC_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseHPC_SPARKX_noCovariates"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX_noCovariates"]]), "_SPARKX_noCovariates")[-1]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]


# note filtering per method: same for both nnSVG and SPARK-X

table(res_list$mouseHPC_SPARKX$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
all(res_list$mouseHPC_nnSVG$gene_name == res_list$mouseHPC_SPARKX$gene_name)

table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_SPARKX$gene_name)


table(res_list$mouseHPC_SPARKX_noCovariates$gene_name %in% res_list$mouseHPC_nnSVG_noCovariates$gene_name)
all(res_list$mouseHPC_nnSVG_noCovariates$gene_name == res_list$mouseHPC_SPARKX_noCovariates$gene_name)

table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_nnSVG_noCovariates$gene_name)
table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_SPARKX_noCovariates$gene_name)


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: Cpne9, Rgs14

known_genes <- c("Cpne9", "Rgs14")

# with covariates
df_known_mouseHPC <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG), 
            as.data.frame(res_list$mouseHPC_SPARKX), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_HVGs), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_SPARKX", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("nnSVG", "SPARKX", "HVGs"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))

# plot ranks
ggplot(as.data.frame(df_known_mouseHPC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("blue3", "maroon", "darkorange")) + 
  scale_y_log10(limits = c(19, 8000)) + 
  geom_text_repel(nudge_x = 0.2, size = 2, segment.color = NA, show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Example SVGs: mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_ranks_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


# without covariates
df_known_mouseHPC_nocovariate <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG_nocovariate), 
            as.data.frame(res_list$mouseHPC_SPARKX_nocovariate), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_HVGs), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG_nocovariate", "rank_SPARKX_nocovariate", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("nnSVG_nocovariate", "SPARKX_nocovariate", "HVGs"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))

# plot ranks
ggplot(as.data.frame(df_known_mouseHPC_nocovariate), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 3, 1)) + 
  scale_color_manual(values = c("deepskyblue", "maroon2", "darkorange")) + 
  scale_y_log10(limits = c(19, 8000)) + 
  geom_text_repel(nudge_x = 0.2, size = 2, segment.color = NA, show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Example SVGs: mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_ranks_mouseHPC_nocovariate"))
ggsave(paste0(fn, ".pdf"), width = 5.75, height = 4)
ggsave(paste0(fn, ".png"), width = 5.75, height = 4)


# --------------------------------------------
# nnSVG: likelihood ratio statistics vs. ranks
# --------------------------------------------

# plot likelihood ratio statistics vs. ranks

df_nnSVG_mouseHPC <- 
  as.data.frame(res_list$mouseHPC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes)

# rank at p-value = 0.05 cutoff
padj_cutoff_nnSVG <- 
  as.data.frame(res_list$mouseHPC_nnSVG) %>% 
  filter(padj_nnSVG <= 0.05) %>% 
  summarize(max = max(rank_nnSVG)) %>% 
  unlist

padj_cutoff_nnSVG


# highlighting 2 known genes
ggplot(as.data.frame(df_nnSVG_mouseHPC), 
       aes(x = rank_nnSVG, y = LR_stat_nnSVG, label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_nnSVG_mouseHPC, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_nnSVG_mouseHPC, gene_name %in% known_genes), 
                  nudge_x = 250, nudge_y = 500, size = 3, color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_nnSVG, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value\n = 0.05\n(rank ", padj_cutoff_nnSVG, ")"), 
           x = 2650, y = 4000, size = 4, color = "darkorange2") + 
  xlim(c(0, 3000)) + 
  labs(x = "rank", y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: mouseHPC, example SVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "LR_stat_ranks_2known_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------------------------------
# SPARK-X: adjusted p-values vs. ranks
# ------------------------------------

# plot adjusted p-values vs. ranks

df_SPARKX_mouseHPC <- 
  as.data.frame(res_list$mouseHPC_SPARKX) %>% 
  mutate(is_known = gene_name %in% known_genes)

# rank at p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$mouseHPC_SPARKX) %>% 
  filter(adjustedPval_SPARKX <= 0.05) %>% 
  summarize(max = max(rank_SPARKX)) %>% 
  unlist

padj_cutoff_SPARKX


# highlighting 193 (combined set of known and layer-specific marker)
ggplot(as.data.frame(df_SPARKX_mouseHPC), 
       aes(x = rank_SPARKX, y = -log10(adjustedPval_SPARKX), label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_SPARKX_mouseHPC, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_SPARKX_mouseHPC, gene_name %in% known_genes), 
                  nudge_x = 500, nudge_y = 15, size = 3, color = "firebrick3") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value\n = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 3400, y = 105, size = 4, color = "darkorange2") + 
  #xlim(c(0, 3000)) + 
  labs(x = "rank", y = "-log10(adjusted p-value)") + 
  ggtitle("SPARK-X: mouseHPC, example SVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "adjPvals_ranks_2known_SPARKX_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ------------
# effect sizes
# ------------

df_effect <- 
  as.data.frame(df_nnSVG_mouseHPC) %>% 
  filter(rank_nnSVG <= 1000)

ggplot(df_effect, 
       aes(x = prop_sv_nnSVG, y = LR_stat_nnSVG, 
           color = is_known, shape = is_known)) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red"), name = "example SVG") + 
  scale_shape_manual(values = c(1, 19), name = "example SVG") + 
  geom_text_repel(data = df_effect %>% filter(is_known), 
                  aes(label = gene_name), nudge_y = 1000, show.legend = FALSE) + 
  labs(x = "proportion of spatial variance", 
       y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_sizes_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# ---------------------
# p-value distributions
# ---------------------

df_pvals_mouseHPC <- as.data.frame(res_list$mouseHPC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals_mouseHPC), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  xlab("p-values") + 
  ggtitle("P-values: nnSVG, mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "pvals_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# -----------------------
# bandwidth distributions
# -----------------------

df_bandwidth_mouseHPC <- 
  as.data.frame(res_list$mouseHPC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG)

phis <- df_bandwidth_mouseHPC$phi_nnSVG
names(phis) <- df_bandwidth_mouseHPC$gene_name
phis_known <- phis[known_genes]

ls_known <- 1 / phis_known
ls_known

ann_text <- data.frame(
  x = unname(ls_known), 
  y =  c(5, 5), 
  label = paste0(names(ls_known), " = ", round(ls_known, 2))
)

# plot bandwidth l (inverse of phi, scaled to distances 0 to 1)
ggplot(as.data.frame(df_bandwidth_mouseHPC), aes(x = l_nnSVG)) + 
  geom_density(color = "black", fill = "skyblue") + 
  xlim(c(0, 1)) + 
  #geom_point(data = ann_text, aes(x = x, y = y), color = "red", size = 2) + 
  #geom_text_repel(data = ann_text, aes(x = x, y = y, label = label), 
  #                nudge_x = 0.13, nudge_y = 1, color = "red", size = 3) + 
  xlab("bandwidth (1 / phi)") + 
  ylab("density") + 
  ggtitle("Bandwidths: nnSVG, mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "bandwidths_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

