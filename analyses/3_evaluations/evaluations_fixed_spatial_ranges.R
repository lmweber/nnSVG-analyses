#################################
# Script to calculate evaluations
# Lukas Weber, Jan 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


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


# ------------------------------------------------------
# overlaps: top n HVGs within top n SVGs for each method
# ------------------------------------------------------

# overlap sizes
overlaps <- c(10, 20, 50, 100, 200)

# function to calculate overlaps for a pair of methods, i.e. proportion of top n
# genes from method 1 (e.g. HVGs) that are also in the set of top n genes from
# method 2 (e.g. nnSVG)
calc_overlaps <- function(method1, method2) {
  
  res_method1 <- res_list[[method1]]
  res_method2 <- res_list[[method2]]
  
  # remove method names from column names
  colnames(res_method1)[-(1:2)] <- gsub("_.*$", "", colnames(res_method1)[-(1:2)])
  colnames(res_method2)[-(1:2)] <- gsub("_.*$", "", colnames(res_method2)[-(1:2)])
  
  top_method2 <- rep(NA, length(overlaps))
  
  for (k in seq_along(overlaps)) {
    # select top gene IDs from method 1
    genes_k <- rownames(filter(as.data.frame(res_method1), 
                               rank <= overlaps[k]))
    # calculate overlaps
    top_method2[k] <- nrow(filter(as.data.frame(res_method2[genes_k, ]), 
                                  rank <= overlaps[k]))
  }
  
  # calculate proportions
  return(top_method2 / overlaps)
}

# use function to calculate overlaps

df_overlaps_DLPFC <- data.frame(
  top_n = overlaps, 
  dataset = "DLPFC", 
  nnSVG = calc_overlaps("DLPFC_HVGs", "DLPFC_nnSVG"), 
  SPARKX = calc_overlaps("DLPFC_HVGs", "DLPFC_SPARKX")
)

df_overlaps_mOB <- data.frame(
  top_n = overlaps, 
  dataset = "mOB", 
  nnSVG = calc_overlaps("mOB_HVGs", "mOB_nnSVG"), 
  SPARKX = calc_overlaps("mOB_HVGs", "mOB_SPARKX")
)

df_overlaps <- 
  rbind(df_overlaps_DLPFC, df_overlaps_mOB) %>% 
  pivot_longer(., cols = c("nnSVG", "SPARKX"), 
               names_to = "method", values_to = "proportion")


# plot overlaps
ggplot(as.data.frame(df_overlaps), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  facet_wrap(~ dataset) + 
  geom_line() + 
  geom_point() + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlapping") + 
  scale_x_log10() + 
  ggtitle("Proportion of top n HVGs within top n SVGs per method") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "prop_overlap_HVGs_in_method"))
ggsave(paste0(fn, ".pdf"), width = 7, height = 4)
ggsave(paste0(fn, ".png"), width = 7, height = 4)


# -------------------------------------------------
# overlaps: top n SVGs from each method within HVGs
# -------------------------------------------------

# i.e. vice versa compared to above

# use function to calculate overlaps

df_overlaps_DLPFC_viceversa <- data.frame(
  top_n = overlaps, 
  dataset = "DLPFC", 
  nnSVG = calc_overlaps("DLPFC_nnSVG", "DLPFC_HVGs"), 
  SPARKX = calc_overlaps("DLPFC_SPARKX", "DLPFC_HVGs")
)

df_overlaps_mOB_viceversa <- data.frame(
  top_n = overlaps, 
  dataset = "mOB", 
  nnSVG = calc_overlaps("mOB_nnSVG", "mOB_HVGs"), 
  SPARKX = calc_overlaps("mOB_SPARKX", "mOB_HVGs")
)

df_overlaps_viceversa <- 
  rbind(df_overlaps_DLPFC_viceversa, df_overlaps_mOB_viceversa) %>% 
  pivot_longer(., cols = c("nnSVG", "SPARKX"), 
               names_to = "method", values_to = "proportion")


# plot overlaps
ggplot(as.data.frame(df_overlaps_viceversa), 
       aes(x = top_n, y = proportion, group = method, color = method)) + 
  facet_wrap(~ dataset) + 
  geom_line() + 
  geom_point() + 
  ylim(c(0, 1)) + 
  xlab("top n genes") + 
  ylab("proportion overlapping") + 
  scale_x_log10() + 
  ggtitle("Proportion of top n SVGs per method within top n HVGs") + 
  theme_bw()

fn <- here(file.path("plots", "evaluations", "prop_overlap_method_in_HVGs_viceversa"))
ggsave(paste0(fn, ".pdf"), width = 7, height = 4)
ggsave(paste0(fn, ".png"), width = 7, height = 4)

