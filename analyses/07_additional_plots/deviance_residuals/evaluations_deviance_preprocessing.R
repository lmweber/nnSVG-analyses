#################################
# Script to calculate evaluations
# Lukas Weber, May 2022
#################################

# deviance residuals preprocessing for nnSVG


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "deviance_preprocessing"))


# ------------
# load results
# ------------

res_list <- list(
  humanDLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  humanDLPFC_nnSVG_devResid = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG_devResid.rds"))), 
  mouseHPC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG.rds"))), 
  mouseHPC_nnSVG_devResid = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG_devResid.rds"))), 
  mouseOB_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG.rds"))), 
  mouseOB_nnSVG_devResid = rowData(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG_devResid.rds"))), 
  mouseEmbryo_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_nnSVG_noFilt.rds"))), 
  mouseEmbryo_nnSVG_devResid = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_nnSVG_devResid.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_nnSVG"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG"]]), "_nnSVG")[-(1:2)]
colnames(res_list[["humanDLPFC_nnSVG_devResid"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_nnSVG_devResid"]]), "_nnSVG_devResid")[-(1:2)]
colnames(res_list[["mouseHPC_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseHPC_nnSVG_devResid"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG_devResid"]]), "_nnSVG_devResid")[-1]
colnames(res_list[["mouseOB_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseOB_nnSVG_devResid"]])[-1] <- paste0(colnames(res_list[["mouseOB_nnSVG_devResid"]]), "_nnSVG_devResid")[-1]
colnames(res_list[["mouseEmbryo_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseEmbryo_nnSVG_devResid"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_nnSVG_devResid"]]), "_nnSVG_devResid")[-1]


# note filtering per method
table(res_list$humanDLPFC_nnSVG_devResid$gene_id %in% res_list$humanDLPFC_nnSVG$gene_id)
table(res_list$mouseHPC_nnSVG_devResid$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
table(res_list$mouseOB_nnSVG_devResid$gene_name %in% res_list$mouseOB_nnSVG$gene_name)
table(res_list$mouseEmbryo_nnSVG_devResid$gene_name %in% res_list$mouseEmbryo_nnSVG$gene_name)


# --------------------------
# known SVGs in each dataset
# --------------------------

known_genes_humanDLPFC <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
known_genes_mouseHPC <- c("Cpne9", "Rgs14")
known_genes_mouseOB <- c("Penk", "Doc2g", "Kctd12", "Slc17a7", "Cdhr1", "Sv2b", "Shisa3")
known_genes_mouseEmbryo <- c("Ttn", "Popdc2", "Hand1", "Gata5", "Six3", "Lhx2", 
                             "Otx2", "Pou3f1", "Sox2", "Foxf1", "Foxa1", "Cldn4")


df_known_humanDLPFC <- 
  full_join(as.data.frame(res_list$humanDLPFC_nnSVG), 
            as.data.frame(res_list$humanDLPFC_nnSVG_devResid), 
            by = c("gene_id", "gene_name")) %>% 
  filter(gene_name %in% known_genes_humanDLPFC) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_devResid"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("nnSVG_devResid", "nnSVG (dev resid)", 
                              gsub("^rank_", "", method)), 
                         levels = c("nnSVG", "nnSVG (dev resid)"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_humanDLPFC))

df_known_mouseHPC <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG), 
            as.data.frame(res_list$mouseHPC_nnSVG_devResid), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseHPC) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_devResid"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("nnSVG_devResid", "nnSVG (dev resid)", 
                              gsub("^rank_", "", method)), 
                         levels = c("nnSVG", "nnSVG (dev resid)"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseHPC))

df_known_mouseOB <- 
  full_join(as.data.frame(res_list$mouseOB_nnSVG), 
            as.data.frame(res_list$mouseOB_nnSVG_devResid), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseOB) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_devResid"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("nnSVG_devResid", "nnSVG (dev resid)", 
                              gsub("^rank_", "", method)), 
                         levels = c("nnSVG", "nnSVG (dev resid)"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseOB))

df_known_mouseEmbryo <- 
  full_join(as.data.frame(res_list$mouseEmbryo_nnSVG), 
            as.data.frame(res_list$mouseEmbryo_nnSVG_devResid), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseEmbryo) %>% 
  pivot_longer(c("rank_nnSVG", "rank_nnSVG_devResid"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("nnSVG_devResid", "nnSVG (dev resid)", 
                              gsub("^rank_", "", method)), 
                         levels = c("nnSVG", "nnSVG (dev resid)"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseEmbryo))


# ----------------------
# plots for each dataset
# ----------------------

ggplot(as.data.frame(df_known_humanDLPFC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 4)) + 
  scale_color_manual(values = c("blue3", "green4")) + 
  scale_y_log10(limits = c(3, 650)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  annotate("text", label = "large length scale", x = 2, y = 650, size = 4) + 
  annotate("text", label = "small length scale", x = 5, y = 650, size = 4) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: human DLPFC") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_humanDLPFC_devResid")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4)


ggplot(as.data.frame(df_known_mouseHPC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 4)) + 
  scale_color_manual(values = c("blue3", "green4")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse HPC") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseHPC_devResid")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


ggplot(as.data.frame(df_known_mouseOB), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 4)) + 
  scale_color_manual(values = c("blue3", "green4")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse OB") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseOB_devResid")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)


ggplot(as.data.frame(df_known_mouseEmbryo), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(4, 4)) + 
  scale_color_manual(values = c("blue3", "green4")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse embryo") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseEmbryo_devResid")
ggsave(paste0(fn, ".pdf"), width = 8, height = 4)
ggsave(paste0(fn, ".png"), width = 8, height = 4)

