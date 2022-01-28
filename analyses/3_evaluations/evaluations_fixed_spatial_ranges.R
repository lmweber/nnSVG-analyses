#################################
# Script to calculate evaluations
# Lukas Weber, Jan 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsci)


# ------------
# load results
# ------------

# DLPFC and mOB datasets

res_list <- list(
  DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds"))), 
  DLPFC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_DLPFC_HVGs.rds"))), 
  mOB_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG.rds"))), 
  mOB_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_mOB_SPARKX.rds"))), 
  mOB_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_mOB_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["DLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["DLPFC_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["DLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["DLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["mOB_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["mOB_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["mOB_SPARKX"]])[-(1:2)] <- paste0(colnames(res_list[["mOB_SPARKX"]]), "_SPARKX")[-(1:2)]
colnames(res_list[["mOB_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["mOB_HVGs"]]), "_HVGs")[-(1:2)]


# ---------------------------
# known SVGs in DLPFC dataset
# ---------------------------

# known SVGs with fixed spatial ranges: SNAP25, MOBP, PCP4

known_genes <- c("SNAP25", "MOBP", "PCP4")

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
  mutate(method = gsub("^rank_", "", method))


# plot ranks
ggplot(as.data.frame(df_known_DLPFC), 
       aes(x = rank, y = gene_name, group = method, color = method)) + 
  geom_point(pch = 4, stroke = 2) + 
  scale_x_log10() + 
  ggtitle("DLPFC dataset, known SVGs") + 
  theme_bw() + 
  theme(axis.title.y = element_blank())

fn <- here(file.path("plots", "evaluations", "DLPFC_known_SVGs_ranks"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

