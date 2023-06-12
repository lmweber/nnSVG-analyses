#################################
# Script to calculate evaluations
# Lukas Weber, updated Jun 2023
#################################

# data set: mouse HPC
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
dir_plots <- here(file.path("plots", "evaluations", "mouseHPC", "no_filtering"))


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs, Moran's I

# note choice of filtering per method
res_list <- list(
  mouseHPC_nnSVG = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG_noFilt.rds"))), 
  mouseHPC_SPARKX = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_SPARKX_noFilt.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_HVGs_noFilt.rds"))), 
  mouseHPC_MoransI = rowData(readRDS(here("outputs", "results", "spe_mouseHPC_MoransI_noFilt.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseHPC_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseHPC_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]
colnames(res_list[["mouseHPC_MoransI"]])[-1] <- paste0(colnames(res_list[["mouseHPC_MoransI"]]), "_MoransI")[-1]


# note filtering per method: no filtering for either nnSVG or SPARK-X
table(res_list$mouseHPC_SPARKX$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
all(res_list$mouseHPC_SPARKX$gene_name == res_list$mouseHPC_nnSVG$gene_name)

table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_SPARKX$gene_name)

table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_nnSVG$gene_name)
table(res_list$mouseHPC_HVGs$gene_name %in% res_list$mouseHPC_SPARKX$gene_name)


# --------------------------
# known SVGs in this dataset
# --------------------------

# known SVGs: Cpne9, Rgs14

known_genes <- c("Cpne9", "Rgs14")

df_known <- 
  full_join(as.data.frame(res_list$mouseHPC_nnSVG), 
            as.data.frame(res_list$mouseHPC_SPARKX), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_HVGs), 
            by = c("gene_name")) %>% 
  full_join(., 
            as.data.frame(res_list$mouseHPC_MoransI), 
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
  scale_y_log10(limits = c(10, 25000)) + 
  geom_text_repel(nudge_x = 0.2, size = 2, segment.color = NA, box.padding = 0.1, 
                  show.legend = FALSE) + 
  labs(x = "gene", y = "rank") + 
  ggtitle("Selected SVGs: mouseHPC (no filtering)") + 
  theme_bw() + 
  theme(axis.text.x = element_text(face = "italic"))

fn <- file.path(dir_plots, "example_SVGs_ranks_mouseHPC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 4.5, height = 4)
ggsave(paste0(fn, ".png"), width = 4.5, height = 4)


# ------------------------------------
# SPARK-X: combined p-values vs. ranks
# ------------------------------------

# plot combined p-values vs. ranks

df_SPARKX <- 
  as.data.frame(res_list$mouseHPC_SPARKX) %>% 
  mutate(is_known = gene_name %in% known_genes)

# rank at adjusted p-value = 0.05 cutoff
padj_cutoff_SPARKX <- 
  as.data.frame(res_list$mouseHPC_SPARKX) %>% 
  filter(adjustedPval_SPARKX <= 0.05) %>% 
  summarize(max = max(rank_SPARKX)) %>% 
  unlist

padj_cutoff_SPARKX


# highlighting known genes
ggplot(as.data.frame(df_SPARKX), 
       aes(x = rank_SPARKX, y = -log10(combinedPval_SPARKX), label = gene_name)) + 
  geom_line(color = "navy") + 
  geom_point(data = filter(df_SPARKX, gene_name %in% known_genes), 
             size = 2, color = "firebrick3") + 
  geom_text_repel(data = filter(df_SPARKX, gene_name %in% known_genes), 
                  nudge_x = 1800, nudge_y = 12, size = 3, 
                  color = "firebrick3", fontface = "italic") + 
  geom_vline(xintercept = padj_cutoff_SPARKX, 
             linetype = "dashed", color = "darkorange2") + 
  annotate("text", label = paste0("adjusted p-value = 0.05\n(rank ", padj_cutoff_SPARKX, ")"), 
           x = 7000, y = 120, size = 3.5, color = "darkorange2") + 
  labs(x = "rank", y = "-log10(combined p-value)") + 
  ggtitle("SPARK-X: mouseHPC (no filtering)") + 
  theme_bw()

fn <- file.path(dir_plots, "stat_vs_rank_SPARKX_mouseHPC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 5, height = 3.75)
ggsave(paste0(fn, ".png"), width = 5, height = 3.75)


# ---------------------
# p-value distributions
# ---------------------

df_pvals <- as.data.frame(res_list$mouseHPC_nnSVG)

# plot p-values
ggplot(as.data.frame(df_pvals), aes(x = pval_nnSVG)) + 
  geom_histogram(color = "black", fill = "blue3", bins = 30) + 
  ylim(c(0, 13750)) + 
  labs(x = "p-values", 
       y = "frequency") + 
  ggtitle("nnSVG p-values: mouse HPC", 
          subtitle = "no filtering") + 
  theme_bw()

fn <- file.path(dir_plots, "pvals_nnSVG_mouseHPC_noFilt")
ggsave(paste0(fn, ".pdf"), width = 4, height = 3.5)
ggsave(paste0(fn, ".png"), width = 4, height = 3.5)


# ---------------------------------------------
# save source data file for publication figures
# ---------------------------------------------

dir_sd <- here("outputs", "source_data")
fn_sd <- "Source_Data_Figs_S11AC_S19D.rds"

saveRDS(res_list, here(dir_sd, fn_sd))

