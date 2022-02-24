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
  mouseHPC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds"))), 
  mouseHPC_nnSVG_covariate = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  mouseHPC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX.rds"))), 
  mouseHPC_SPARKX_covariate = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX_clusters.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseHPC_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseHPC_nnSVG_covariate"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG_covariate"]]), "_nnSVG_covariate")[-1]
colnames(res_list[["mouseHPC_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseHPC_SPARKX_covariate"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX_covariate"]]), "_SPARKX_covariate")[-1]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]


# ------------------------------
# known SVGs in mouseHPC dataset
# ------------------------------

# known SVGs: Cpne9, Rgs14

known_genes <- c("Cpne9", "Rgs14")

all(res_list$mouseHPC_nnSVG$gene_name == res_list$mouseHPC_nnSVG_covariate$gene_name)
all(res_list$mouseHPC_nnSVG$gene_name == res_list$mouseHPC_SPARKX$gene_name)
all(res_list$mouseHPC_nnSVG$gene_name == res_list$mouseHPC_SPARKX_covariate$gene_name)
all(res_list$mouseHPC_nnSVG$gene_name == res_list$mouseHPC_HVGs$gene_name)

df_known_mouseHPC <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG), 
            as.data.frame(res_list$mouseHPC_nnSVG_covariate), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_SPARKX), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_SPARKX_covariate), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_HVGs), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_covariate", "rank_SPARKX", "rank_SPARKX_covariate", "rank_HVGs"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("nnSVG", "nnSVG_covariate", "SPARKX", "SPARKX_covariate", "HVGs"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes))


# plot ranks
# set seed for ggrepel label positioning
set.seed(1)
ggplot(as.data.frame(df_known_mouseHPC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.5, size = 1.75) + 
  scale_shape_manual(values = c(4, 4, 3, 3, 1)) + 
  scale_color_manual(values = c("blue3", "deepskyblue", "maroon", "maroon1", "darkorange")) + 
  scale_y_log10() + 
  geom_text_repel(nudge_x = 0.2, size = 2, segment.color = NA, show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Example SVGs: mouseHPC") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "known_genes_ranks_mouseHPC"))
ggsave(paste0(fn, ".pdf"), width = 5.25, height = 4)
ggsave(paste0(fn, ".png"), width = 5.25, height = 4)

