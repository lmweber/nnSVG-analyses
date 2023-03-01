###############################
# Multiple samples: evaluations
# Lukas Weber, Mar 2023
###############################

# dataset: Visium human DLPFC (multiple samples)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(reshape2)
library(ggplot2)


# directory to save plots
dir_plots <- here(file.path("plots", "multiple_samples"))


# ------------
# load results
# ------------

res_list <- readRDS(here("outputs", "multiple_samples", "res_humanDLPFC_nnSVG_multipleSamples.rds"))

length(res_list)

sample_ids <- names(res_list)
sample_ids


# ----------------------------------------
# extract ranks and calculate correlations
# ----------------------------------------

# extract ranks per sample and calculate inter-sample correlations

df_ranks <- data.frame(
  gene_name = res_list[[1]][, "gene_name"], 
  sample = sample_ids[1], 
  rank = res_list[[1]][, "rank"]
)
for (s in 2:length(sample_ids)) {
  df_ranks_s <- data.frame(
    gene_name = res_list[[s]][, "gene_name"], 
    sample = sample_ids[s], 
    rank = res_list[[s]][, "rank"]
  )
  df_ranks <- rbind(df_ranks, df_ranks_s)
}
df_ranks$sample <- factor(df_ranks$sample, levels = sample_ids)


# calculate correlations

mat_cors <- matrix(NA, nrow = length(sample_ids), ncol = length(sample_ids))
rownames(mat_cors) <- colnames(mat_cors) <- sample_ids

for (i in 1:(length(sample_ids) - 1)) {
  for (j in (min(i + 1, length(sample_ids)):length(sample_ids))) {
    df_ranks_i <- df_ranks[df_ranks$sample == sample_ids[i], ]
    df_ranks_j <- df_ranks[df_ranks$sample == sample_ids[j], ]
    df_ranks_ij <- full_join(df_ranks_i, df_ranks_j, by = "gene_name")
    df_ranks_ij <- na.omit(df_ranks_ij)
    # calculate Spearman correlation
    cor_ij <- cor(df_ranks_ij$rank.x, df_ranks_ij$rank.y, method = "spearman")
    # store in matrix
    mat_cors[i, j] <- cor_ij
  }
}

mat_cors


# --------------------------
# plot scatterplots of ranks
# --------------------------

# compare ranks to sample 151673 (sample used in main results)

# sample with highest correlation: 151509
# sample with lowest correlation: 151669

df_ranks_151673 <- df_ranks[df_ranks$sample == "151673", c("gene_name", "rank")]
df_ranks_151509 <- df_ranks[df_ranks$sample == "151509", c("gene_name", "rank")]
df_ranks_151669 <- df_ranks[df_ranks$sample == "151669", c("gene_name", "rank")]

df_ranks_highest <- na.omit(full_join(df_ranks_151673, df_ranks_151509, 
                                      by = "gene_name", 
                                      suffix = c("_151673", "_151509")))
df_ranks_lowest <- na.omit(full_join(df_ranks_151673, df_ranks_151669, 
                                     by = "gene_name", 
                                     suffix = c("_151673", "_151669")))

cor_highest <- round(mat_cors["151509", "151673"], 3)
cor_lowest <- round(mat_cors["151669", "151673"], 3)


# plot scatterplots


# highest correlation
ggplot(as.data.frame(df_ranks_highest), 
       aes(x = rank_151673, y = rank_151509)) + 
  geom_point(size = 0.75, color = "navy") + 
  geom_text(label = paste0("cor = ", cor_highest), 
            x = 850, y = 25, color = "maroon") + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  ggtitle("Ranks vs. sample 151673", 
          subtitle = "Human DLPFC dataset") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_scatter_multipleSamples_highestCor")
ggsave(paste0(fn, ".pdf"), width = 4.25, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.25, height = 4.5)


# lowest correlation
ggplot(as.data.frame(df_ranks_lowest), 
       aes(x = rank_151673, y = rank_151669)) + 
  geom_point(size = 0.75, color = "navy") + 
  geom_text(label = paste0("cor = ", cor_lowest), 
            x = 850, y = 25, color = "maroon") + 
  coord_fixed() + 
  xlim(c(0, 1000)) + 
  ylim(c(0, 1000)) + 
  ggtitle("Ranks vs. sample 151673", 
          subtitle = "Human DLPFC dataset") + 
  theme_bw()

fn <- file.path(dir_plots, "ranks_scatter_multipleSamples_lowestCor")
ggsave(paste0(fn, ".pdf"), width = 4.25, height = 4.5)
ggsave(paste0(fn, ".png"), width = 4.25, height = 4.5)


# -------------------------------
# plot correlation matrix heatmap
# -------------------------------

# reshape correlation matrix for plotting
mat_cors_plot <- mat_cors
diag(mat_cors_plot) <- 1

melted_cormat <- melt(mat_cors_plot, value.name = "cor", na.rm = TRUE)
melted_cormat$Var1 <- factor(melted_cormat$Var1, levels = sample_ids)
melted_cormat$Var2 <- factor(melted_cormat$Var2, levels = sample_ids)

min_cor <- min(mat_cors, na.rm = TRUE)
max_cor <- max(mat_cors, na.rm = TRUE)
mid_cor <- mean(c(min_cor, max_cor))

# plot correlation matrix
ggplot(data = melted_cormat, aes(x = Var1, y = Var2, 
                                 label = format(round(cor, digits = 2), nsmall = 2), 
                                 fill = cor)) + 
  geom_tile(color = "white") + 
  scale_fill_gradient2(low = "gold", mid = "darkorange", high = "firebrick3", 
                       midpoint = mid_cor, limit = c(min_cor, max_cor), 
                       name = "Spearman\ncorrelation") + 
  geom_text(color = "black", size = 2.6) + 
  coord_fixed() + 
  ggtitle("Sample-to-sample rank correlations") + 
  theme_bw() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5), 
        axis.title.x = element_blank(), 
        axis.title.y = element_blank())

fn <- file.path(dir_plots, "ranks_heatmap_multipleSamples")
ggsave(paste0(fn, ".pdf"), width = 5, height = 4.25)
ggsave(paste0(fn, ".png"), width = 5, height = 4.25)

