#####################################################
# Script to plot runtimes for scalability simulations
# Lukas Weber, Mar 2023
#####################################################

# comparisons between methods


library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------
# load results
# ------------

runtimes_nnSVG_DLPFC_twogenes <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_nnSVG_DLPFC_twogenes.rds"))
runtimes_SPARKX_DLPFC_twogenes <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_SPARKX_DLPFC_twogenes.rds"))
runtimes_SPARK_DLPFC_twogenes <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_SPARK_DLPFC_twogenes.rds"))
runtimes_SpatialDE_DLPFC_twogenes <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_SpatialDE_DLPFC_twogenes.rds"))

mat_runtimes_nnSVG <- do.call("rbind", runtimes_nnSVG_DLPFC_twogenes)
colnames(mat_runtimes_nnSVG) <- paste0("seed", 1:ncol(mat_runtimes_nnSVG))
df_runtimes_nnSVG <- cbind(
  n_spots = rownames(mat_runtimes_nnSVG), 
  method = "nnSVG", 
  as.data.frame(mat_runtimes_nnSVG))

mat_runtimes_SPARKX <- do.call("rbind", runtimes_SPARKX_DLPFC_twogenes)
colnames(mat_runtimes_SPARKX) <- paste0("seed", 1:ncol(mat_runtimes_SPARKX))
df_runtimes_SPARKX <- cbind(
  n_spots = rownames(mat_runtimes_SPARKX), 
  method = "SPARK-X", 
  as.data.frame(mat_runtimes_SPARKX))

mat_runtimes_SPARK <- do.call("rbind", runtimes_SPARK_DLPFC_twogenes)
colnames(mat_runtimes_SPARK) <- paste0("seed", 1:ncol(mat_runtimes_SPARK))
df_runtimes_SPARK <- cbind(
  n_spots = rownames(mat_runtimes_SPARK), 
  method = "SPARK", 
  as.data.frame(mat_runtimes_SPARK))

mat_runtimes_SpatialDE <- do.call("rbind", runtimes_SpatialDE_DLPFC_twogenes)
colnames(mat_runtimes_SpatialDE) <- paste0("seed", 1:ncol(mat_runtimes_SpatialDE))
df_runtimes_SpatialDE <- cbind(
  n_spots = rownames(mat_runtimes_SpatialDE), 
  method = "SpatialDE", 
  as.data.frame(mat_runtimes_SpatialDE))


# --------------
# generate plots
# --------------

df_nnSVG <- pivot_longer(df_runtimes_nnSVG, cols = starts_with("seed"), 
                         names_to = "seed", values_to = "runtime")
df_SPARKX <- pivot_longer(df_runtimes_SPARKX, cols = starts_with("seed"), 
                          names_to = "seed", values_to = "runtime")
df_SPARK <- pivot_longer(df_runtimes_SPARK, cols = starts_with("seed"), 
                         names_to = "seed", values_to = "runtime")
df_SpatialDE <- pivot_longer(df_runtimes_SpatialDE, cols = starts_with("seed"), 
                         names_to = "seed", values_to = "runtime")

df <- rbind(df_nnSVG, df_SPARKX, df_SPARK, df_SpatialDE)

df$n_spots <- as.numeric(df$n_spots)
df$method <- factor(df$method, levels = c("nnSVG", "SPARK-X", "SPARK", "SpatialDE"))

x_vals <- c(0, unique(df$n_spots))

pal <- c("blue3", "deepskyblue2", "sienna", "seagreen3")


ggplot(df, aes(x = n_spots, y = runtime, color = method, shape = method, 
               group = interaction(n_spots, method))) + 
  geom_point(stroke = 1.25, size = 2) + 
  scale_color_manual(values = pal) + 
  scale_shape_manual(values = c(4, 3, 1, 2)) + 
  scale_x_continuous(breaks = x_vals) + 
  scale_y_continuous(trans = scales::trans_new("cubrt", function(x) x^(1/3), function(x) x^3), 
                     limits = c(0, 500), 
                     breaks = c(1, 2, 5, c(1, 2, 5) * 10, c(1, 2, 5) * 100)) + 
  labs(x = "number of spots", 
       y = "runtime (sec) (cubic scale)") + 
  ggtitle("Scalability comparisons") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability_sims", "runtimes_nnSVG_DLPFC_comparisons.png"), 
       width = 5.25, height = 4)

