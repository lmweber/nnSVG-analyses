####################################################################
# Simulations: plots - spatial coordinate masks and expression plots
# Lukas Weber, Mar 2023
####################################################################


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)


dir_sims <- here("outputs", "simulations")
dir_plots <- here("plots", "simulations")


# --------------------------------
# filenames for simulated datasets
# --------------------------------

sim_names_main <- c(
  "sim_largeBandwidth_fullExpr", 
  "sim_largeBandwidth_medExpr", 
  "sim_largeBandwidth_lowExpr", 
  "sim_medBandwidth_fullExpr", 
  "sim_medBandwidth_medExpr", 
  "sim_medBandwidth_lowExpr", 
  "sim_smallBandwidth_fullExpr", 
  "sim_smallBandwidth_medExpr", 
  "sim_smallBandwidth_lowExpr"
)

sim_names_shuffle <- c(
  "sim_shuffle00", 
  "sim_shuffle01", 
  "sim_shuffle02", 
  "sim_shuffle03", 
  "sim_shuffle04", 
  "sim_shuffle05", 
  "sim_shuffle06", 
  "sim_shuffle07", 
  "sim_shuffle08", 
  "sim_shuffle09", 
  "sim_shuffle10"
)


# -----------------------------------------------------------
# load spatial coordinates and masks; example gene expression
# -----------------------------------------------------------

# main simulations

res_list_main <- as.list(rep(NA, length(sim_names_main)))
names(res_list_main) <- sim_names_main

for (s in seq_along(sim_names_main)) {
  
  # load simulated dataset
  spe <- readRDS(here(dir_sims, "main", paste0("spe_", sim_names_main[s], ".rds")))
  
  # extract spatial coordinates and masks
  spatial_coords <- spatialCoords(spe)
  mask <- colData(spe)$mask
  stopifnot(nrow(spatial_coords) == nrow(mask))
  
  # expression of example expressed gene
  ix_example <- which(rowData(spe)$expressed)[1]
  expr_example <- logcounts(spe)[ix_example, ]
  
  res <- data.frame(
    spatial_coords, 
    mask, 
    sim_name = sim_names_main[s], 
    expression_strength = gsub("^.*Bandwidth_", "", gsub("Expr$", "", sim_names_main[s])), 
    expr_example = expr_example
  )
  
  res_list_main[[s]] <- res
}


# shuffled simulations

res_list_shuffle <- as.list(rep(NA, length(sim_names_shuffle)))
names(res_list_shuffle) <- sim_names_shuffle

for (s in seq_along(sim_names_shuffle)) {
  
  # load simulated dataset
  spe <- readRDS(here(dir_sims, "shuffle", paste0("spe_", sim_names_shuffle[s], ".rds")))
  
  # extract spatial coordinates and masks
  spatial_coords <- spatialCoords(spe)
  mask <- colData(spe)$mask
  stopifnot(nrow(spatial_coords) == nrow(mask))
  
  # expression of example expressed gene
  ix_example <- which(rowData(spe)$expressed)[1]
  expr_example <- logcounts(spe)[ix_example, ]
  
  res <- data.frame(
    spatial_coords, 
    mask, 
    sim_name = sim_names_shuffle[s], 
    expr_example = expr_example
  )
  
  res_list_shuffle[[s]] <- res
}


# -----------------------
# plots: main simulations
# -----------------------

df_plot_main <- do.call("rbind", res_list_main)
rownames(df_plot_main) <- NULL
df_plot_main$sim_name <- factor(df_plot_main$sim_name, 
                                levels = sim_names_main)
df_plot_main$expression_strength <- factor(df_plot_main$expression_strength, 
                                           levels = c("full", "med", "low"))
# arbitrary alpha values for best visualization in plots only
df_plot_main$alpha <- 1
df_plot_main$alpha[df_plot_main$mask & df_plot_main$expression_strength == "med"] <- 0.5
df_plot_main$alpha[df_plot_main$mask & df_plot_main$expression_strength == "low"] <- 0.2


# plot spatial coordinate masks

ggplot(df_plot_main, 
       aes(x = x, y = y, color = mask, alpha = alpha)) + 
  facet_wrap(~ sim_name, nrow = 3) + 
  geom_point(size = 0.2) + 
  coord_fixed() + 
  scale_color_manual(values = c("dodgerblue", "darkorange")) + 
  scale_alpha_continuous(range = range(df_plot_main$alpha), breaks = c(1, 0.5, 0.2), 
                         name = "expression\nstrength", labels = c("full", "medium", "low")) + 
  ggtitle("Simulated datasets: spatial coordinate masks") + 
  guides(color = guide_legend(override.aes = list(size = 2.5), order = 1), 
         alpha = guide_legend(override.aes = list(size = 2.5, color = "darkorange"))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

fn <- file.path(dir_plots, "simulations_masks_main")
ggsave(paste0(fn, ".pdf"), width = 7, height = 7)
ggsave(paste0(fn, ".png"), width = 7, height = 7)


# expression plots

ggplot(df_plot_main, 
       aes(x = x, y = y, color = expr_example)) + 
  facet_wrap(~ sim_name, nrow = 3) + 
  geom_point(size = 0.2) + 
  coord_fixed() + 
  scale_color_viridis_c(name = "logcounts") + 
  ggtitle("Simulated datasets: expression") + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

fn <- file.path(dir_plots, "simulations_expression_main")
ggsave(paste0(fn, ".pdf"), width = 7, height = 7)
ggsave(paste0(fn, ".png"), width = 7, height = 7)


# ---------------------------
# plots: shuffled simulations
# ---------------------------

df_plot_shuffle <- do.call("rbind", res_list_shuffle)
rownames(df_plot_shuffle) <- NULL

sim_labs <- paste0(as.numeric(
  gsub("sim_shuffle", "", sim_names_shuffle)) * 10, "% shuffled")

df_plot_shuffle$sim_name <- factor(df_plot_shuffle$sim_name, 
                                   levels = sim_names_shuffle, labels = sim_labs)


# plot spatial coordinate masks

ggplot(df_plot_shuffle, 
       aes(x = x, y = y, color = mask)) + 
  facet_wrap(~ sim_name, ncol = 4) + 
  geom_point(size = 0.01) + 
  coord_fixed() + 
  scale_color_manual(values = c("dodgerblue", "darkorange")) + 
  ggtitle("Simulated datasets: shuffled coordinates") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

fn <- file.path(dir_plots, "simulations_masks_shuffle")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5.6)
ggsave(paste0(fn, ".png"), width = 7, height = 5.6)


# expression plots

ggplot(df_plot_shuffle, 
       aes(x = x, y = y, color = expr_example)) + 
  facet_wrap(~ sim_name, ncol = 4) + 
  geom_point(size = 0.01) + 
  coord_fixed() + 
  scale_color_viridis_c(name = "logcounts") + 
  ggtitle("Simulated datasets: shuffled coordinates (expression)") + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme_bw() + 
  theme(panel.grid = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

fn <- file.path(dir_plots, "simulations_expression_shuffle")
ggsave(paste0(fn, ".pdf"), width = 7, height = 5.6)
ggsave(paste0(fn, ".png"), width = 7, height = 5.6)

