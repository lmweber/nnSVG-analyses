#################################
# Script to calculate evaluations
# Lukas Weber, Feb 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------
# load results
# ------------

# DLPFC dataset

res_list <- list(
  hippo_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds"))), 
  hippo_nnSVG_clusters = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  hippo_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX.rds"))), 
  hippo_SPARKX_clusters = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX_clusters.rds"))), 
  hippo_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["hippo_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["hippo_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["hippo_nnSVG_clusters"]])[-(1:2)] <- paste0(colnames(res_list[["hippo_nnSVG_clusters"]]), "_nnSVG_clusters")[-(1:2)]
colnames(res_list[["hippo_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["hippo_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["hippo_SPARKX_clusters"]])[-(1:2)] <- paste0(colnames(res_list[["hippo_SPARKX_clusters"]]), "_SPARKX_clusters")[-(1:2)]
colnames(res_list[["hippo_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["hippo_HVGs"]]), "_HVGs")[-(1:2)]


# -----------------------------------
# known SVGs in SlideSeqHippo dataset
# -----------------------------------

# known SVGs within spatial domains: Rgs14, Cpne9

known_genes <- c("Rgs14", "Cpne9")

all(res_list$hippo_nnSVG$gene_name == res_list$hippo_nnSVG_clusters$gene_name)
all(res_list$hippo_nnSVG$gene_name == res_list$hippo_SPARKX$gene_name)
all(res_list$hippo_nnSVG$gene_name == res_list$hippo_SPARKX_clusters$gene_name)
all(res_list$hippo_nnSVG$gene_name == res_list$hippo_HVGs$gene_name)

df_known_hippo <- 
  full_join(as.data.frame(res_list$hippo_nnSVG), 
            as.data.frame(res_list$hippo_nnSVG_clusters), 
            by = "gene_name") %>% 
  full_join(., 
            as.data.frame(res_list$hippo_SPARKX), 
            by = "gene_name") %>% 
  full_join(., 
            as.data.frame(res_list$hippo_SPARKX_clusters), 
            by = "gene_name") %>% 
  full_join(., 
            as.data.frame(res_list$hippo_HVGs), 
            by = "gene_name") %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_clusters", "rank_SPARKX", "rank_SPARKX_clusters", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = gsub("^rank_", "", method))


# plot ranks
ggplot(as.data.frame(df_known_hippo), 
       aes(x = rank, y = gene_name, group = method, color = method)) + 
  geom_point(pch = 4, stroke = 2) + 
  scale_x_log10() + 
  scale_color_manual(values = c("maroon", "#F8766D", "#C77CFF", "#00BFC4", "#7CAE00")) + 
  ggtitle("SlideSeqHippo dataset, known SVGs") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

fn <- here(file.path("plots", "evaluations", "hippo_known_SVGs_ranks"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

