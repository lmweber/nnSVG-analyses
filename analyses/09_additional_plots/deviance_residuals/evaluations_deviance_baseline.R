#################################
# Script to calculate evaluations
# Lukas Weber, updated Jun 2023
#################################

# deviance residuals baseline method


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)
library(viridis)


# directory to save plots
dir_plots <- here(file.path("plots", "deviance_baseline"))


# ------------
# load results
# ------------

res_list <- list(
  humanDLPFC_HVGs = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_HVGs_noFilt.rds"))), 
  humanDLPFC_deviance = rowData(readRDS(here("outputs", "results", "spe_humanDLPFC_deviance_noFilt.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_HVGs_noFilt.rds"))), 
  mouseHPC_deviance = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_deviance_noFilt.rds"))), 
  mouseOB_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseOB_HVGs_noFilt.rds"))), 
  mouseOB_deviance = rowData(readRDS(here("outputs", "results", "spe_mouseOB_deviance_noFilt.rds"))), 
  mouseEmbryo_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_HVGs_noFilt.rds"))), 
  mouseEmbryo_deviance = rowData(readRDS(here("outputs", "results", "spe_mouseEmbryo_deviance_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["humanDLPFC_HVGs"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_HVGs"]]), "_HVGs")[-(1:2)]
colnames(res_list[["humanDLPFC_deviance"]])[-(1:2)] <- paste0(colnames(res_list[["humanDLPFC_deviance"]]), "_deviance")[-(1:2)]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseHPC_deviance"]])[-1] <- paste0(colnames(res_list[["mouseHPC_deviance"]]), "_deviance")[-1]
colnames(res_list[["mouseOB_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseOB_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseOB_deviance"]])[-1] <- paste0(colnames(res_list[["mouseOB_deviance"]]), "_deviance")[-1]
colnames(res_list[["mouseEmbryo_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseEmbryo_deviance"]])[-1] <- paste0(colnames(res_list[["mouseEmbryo_deviance"]]), "_deviance")[-1]


# note filtering per method
table(res_list$humanDLPFC_deviance$gene_id %in% res_list$humanDLPFC_HVGs$gene_id)
table(res_list$mouseHPC_deviance$gene_name %in% res_list$mouseHPC_HVGs$gene_name)
table(res_list$mouseOB_deviance$gene_name %in% res_list$mouseOB_HVGs$gene_name)
table(res_list$mouseEmbryo_deviance$gene_name %in% res_list$mouseEmbryo_HVGs$gene_name)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S21ABCD.RData"

save(res_list, file = here(dir_sd, fn_sd))


# --------------------------
# known SVGs in each dataset
# --------------------------

known_genes_humanDLPFC <- c("MOBP", "PCP4", "SNAP25", "HBB", "IGKC", "NPY")
known_genes_mouseHPC <- c("Cpne9", "Rgs14")
known_genes_mouseOB <- c("Penk", "Doc2g", "Kctd12", "Slc17a7", "Cdhr1", "Sv2b", "Shisa3")
known_genes_mouseEmbryo <- c("Ttn", "Popdc2", "Hand1", "Gata5", "Six3", "Lhx2", 
                             "Otx2", "Pou3f1", "Sox2", "Foxf1", "Foxa1", "Cldn4")


df_known_humanDLPFC <- 
  full_join(as.data.frame(res_list$humanDLPFC_HVGs), 
            as.data.frame(res_list$humanDLPFC_deviance), 
            by = c("gene_id", "gene_name")) %>% 
  filter(gene_name %in% known_genes_humanDLPFC) %>% 
  pivot_longer(c("rank_HVGs", "rank_deviance"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("HVGs", "deviance"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_humanDLPFC))

df_known_mouseHPC <- 
  full_join(as.data.frame(res_list$mouseHPC_HVGs), 
            as.data.frame(res_list$mouseHPC_deviance), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseHPC) %>% 
  pivot_longer(c("rank_HVGs", "rank_deviance"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("HVGs", "deviance"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseHPC))

df_known_mouseOB <- 
  full_join(as.data.frame(res_list$mouseOB_HVGs), 
            as.data.frame(res_list$mouseOB_deviance), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseOB) %>% 
  pivot_longer(c("rank_HVGs", "rank_deviance"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("HVGs", "deviance"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseOB))

df_known_mouseEmbryo <- 
  full_join(as.data.frame(res_list$mouseEmbryo_HVGs), 
            as.data.frame(res_list$mouseEmbryo_deviance), 
            by = c("gene_name")) %>% 
  filter(gene_name %in% known_genes_mouseEmbryo) %>% 
  pivot_longer(c("rank_HVGs", "rank_deviance"), 
               names_to = "method", 
               values_to = "rank") %>% 
  mutate(method = factor(gsub("^rank_", "", method), 
                         levels = c("HVGs", "deviance"))) %>% 
  mutate(gene_name = factor(gene_name, levels = known_genes_mouseEmbryo))


# ----------------------
# plots for each dataset
# ----------------------

ggplot(as.data.frame(df_known_humanDLPFC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(1, 1)) + 
  scale_color_manual(values = c("darkorange", "forestgreen")) + 
  scale_y_log10(limits = c(7, 120)) + 
  geom_vline(xintercept = 3.5, linetype = "dashed", color = "gray50") + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  annotate("text", label = "large length scale", x = 2, y = 120, size = 4) + 
  annotate("text", label = "small length scale", x = 5, y = 120, size = 4) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: human DLPFC") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_humanDLPFC_deviance")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4)
ggsave(paste0(fn, ".png"), width = 5, height = 4)


ggplot(as.data.frame(df_known_mouseHPC), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(1, 1)) + 
  scale_color_manual(values = c("darkorange", "forestgreen")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse HPC") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseHPC_deviance")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


ggplot(as.data.frame(df_known_mouseOB), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(1, 1)) + 
  scale_color_manual(values = c("darkorange", "forestgreen")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse OB") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseOB_deviance")
ggsave(paste0(fn, ".pdf"), width = 5.5, height = 4)
ggsave(paste0(fn, ".png"), width = 5.5, height = 4)


ggplot(as.data.frame(df_known_mouseEmbryo), 
       aes(x = gene_name, y = rank, group = method, color = method, 
           shape = method, label = rank)) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(1, 1)) + 
  scale_color_manual(values = c("darkorange", "forestgreen")) + 
  geom_text_repel(nudge_x = 0.35, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouse embryo") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseEmbryo_deviance")
ggsave(paste0(fn, ".pdf"), width = 7.5, height = 4)
ggsave(paste0(fn, ".png"), width = 7.5, height = 4)

