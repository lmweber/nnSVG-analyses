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
  DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds"))), 
  DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["DLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["DLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["DLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_HVGs"]]), "_HVGs")[-(1:2)]


# ------------
# effect sizes
# ------------

# known genes
known_genes <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")

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


df_nnSVG_DLPFC <- 
  as.data.frame(res_list$DLPFC_nnSVG) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% manual_gene_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known)


df_effect <- 
  as.data.frame(df_nnSVG_DLPFC) %>% 
  mutate(is_marker_or_known = is_marker | is_known) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)


# LR statistic vs. effect size

ggplot(df_effect, 
       aes(x = prop_sv_nnSVG, y = LR_stat_nnSVG, 
           color = is_marker_or_known, shape = is_marker_or_known)) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red"), name = "example SVG\nor marker") + 
  scale_shape_manual(values = c(1, 19), name = "example SVG\nor marker") + 
  geom_text_repel(data = df_effect %>% filter(is_known), 
                  aes(label = gene_name), nudge_y = 2000, show.legend = FALSE) + 
  labs(x = "proportion of spatial variance", 
       y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_sizes_LRstat_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# rank vs. effect size

ggplot(df_effect, 
       aes(x = prop_sv_nnSVG, y = rank_nnSVG, 
           color = is_marker_or_known, shape = is_marker_or_known)) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red"), name = "example SVG\nor marker") + 
  scale_shape_manual(values = c(1, 19), name = "example SVG\nor marker") + 
  geom_text_repel(data = df_effect %>% filter(is_known), 
                  aes(label = gene_name), nudge_y = 100, show.legend = FALSE) + 
  labs(x = "proportion of spatial variance", 
       y = "rank") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_sizes_rank_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# bandwidth vs. effect size

ggplot(df_effect, 
       aes(x = prop_sv_nnSVG, y = l_nnSVG, 
           color = is_marker_or_known, shape = is_marker_or_known)) + 
  geom_point() + 
  ylim(c(0, 1)) + 
  scale_color_manual(values = c("black", "red"), name = "example SVG\nor marker") + 
  scale_shape_manual(values = c(1, 19), name = "example SVG\nor marker") + 
  geom_text_repel(data = df_effect %>% filter(is_known), 
                  aes(label = gene_name), show.legend = FALSE) + 
  labs(x = "proportion of spatial variance", 
       y = "bandwidth (l = 1/phi)") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_sizes_bandwidth_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# LR stat vs. adjusted effect size (subtracting HVGs trend)

# calculate adjusted effect size
df_adj_effect <- 
  full_join(as.data.frame(res_list$DLPFC_nnSVG), 
            as.data.frame(res_list$DLPFC_HVGs), 
            by = c("gene_id", "gene_name")) %>% 
  mutate(adj_effect_nnSVG = (sigma.sq_nnSVG - tech_HVGs) / (sigma.sq_nnSVG + tau.sq_nnSVG - tech_HVGs)) %>% 
  mutate(is_known = gene_name %in% known_genes) %>% 
  mutate(is_marker = gene_name %in% manual_gene_names) %>% 
  mutate(is_marker_or_known = is_marker | is_known) %>% 
  mutate(l_nnSVG = 1 / phi_nnSVG) %>% 
  filter(rank_nnSVG <= 1000)

ggplot(df_adj_effect, 
       aes(x = adj_effect_nnSVG, y = LR_stat_nnSVG, 
           color = is_marker_or_known, shape = is_marker_or_known)) + 
  geom_point() + 
  scale_color_manual(values = c("black", "red"), name = "example SVG\nor marker") + 
  scale_shape_manual(values = c(1, 19), name = "example SVG\nor marker") + 
  geom_text_repel(data = df_adj_effect %>% filter(is_known), 
                  aes(label = gene_name), show.legend = FALSE) + 
  labs(x = "adjusted effect size", 
       y = "likelihood ratio statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "adj_effect_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# adjusted effect size (subtracting HVGs trend) vs. mean

ggplot(df_adj_effect, 
       aes(x = mean_nnSVG, y = adj_effect_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.8) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_adj_effect %>% filter(is_marker_or_known), 
             color = "red", size = 0.8) + 
  geom_text_repel(data = df_adj_effect %>% filter(is_known), 
                  aes(label = gene_name), color = "red", show.legend = FALSE) + 
  labs(x = "mean", 
       y = "adjusted effect size", 
       color = "LR statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "adj_effect_colors_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# effect size vs. mean

ggplot(df_adj_effect, 
       aes(x = mean_nnSVG, y = prop_sv_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.8) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_adj_effect %>% filter(is_marker_or_known), 
             color = "red", size = 0.8) + 
  geom_text_repel(data = df_adj_effect %>% filter(is_known), 
                  aes(label = gene_name), color = "red", show.legend = FALSE) + 
  ylim(c(0, 1)) + 
  labs(x = "mean", 
       y = "proportion spatial variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "prop_SV_colors_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)


# effect size vs. mean (formatted)

ggplot(df_adj_effect, 
       aes(x = mean_nnSVG, y = prop_sv_nnSVG, color = LR_stat_nnSVG)) + 
  geom_point(size = 0.75) + 
  scale_color_viridis(trans = "log10") + 
  geom_point(data = df_adj_effect %>% filter(is_marker_or_known), 
             aes(shape = is_marker_or_known), color = "red", size = 0.8) + 
  scale_shape_manual(values = 1, name = "markers") + 
  geom_text_repel(
    data = df_adj_effect %>% filter(is_known | (is_marker_or_known & (mean_nnSVG > 0.6 & prop_sv_nnSVG > 0.48))), 
    aes(label = gene_name), color = "red", size = 3, nudge_x = 0.5, nudge_y = 0.07) + 
  ylim(c(0, 1)) + 
  labs(x = "mean logcounts", 
       y = "proportion spatial variance", 
       color = "LR statistic") + 
  ggtitle("nnSVG: DLPFC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "effect_size_mean_DLPFC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

