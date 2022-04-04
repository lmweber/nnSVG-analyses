#####################################################
# Script to plot runtimes for scalability simulations
# Lukas Weber, Feb 2022
#####################################################


library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------
# load results
# ------------

runtimes_nnSVG_DLPFC_singlegene <- readRDS(here(
  "outputs", "scalability", "runtimes_scalability_nnSVG_DLPFC_singlegene.rds"))
runtimes_nnSVG_mouseHPC_singlegene <- readRDS(here(
  "outputs", "scalability", "runtimes_scalability_nnSVG_mouseHPC_singlegene.rds"))

mat_runtimes_DLPFC <- do.call("rbind", runtimes_nnSVG_DLPFC_singlegene)
colnames(mat_runtimes_DLPFC) <- paste0("seed", 1:ncol(mat_runtimes_DLPFC))
df_runtimes_DLPFC <- cbind(
  n_spots = rownames(mat_runtimes_DLPFC), 
  dataset = "DLPFC", 
  as.data.frame(mat_runtimes_DLPFC))

mat_runtimes_mouseHPC <- do.call("rbind", runtimes_nnSVG_mouseHPC_singlegene)
colnames(mat_runtimes_mouseHPC) <- paste0("seed", 1:ncol(mat_runtimes_mouseHPC))
df_runtimes_mouseHPC <- cbind(
  n_spots = rownames(mat_runtimes_mouseHPC), 
  dataset = "mouseHPC", 
  as.data.frame(mat_runtimes_mouseHPC))


# --------------------------
# add linear and cubic lines
# --------------------------

# calculate linear slope
slope_DLPFC <- 
  (mean(mat_runtimes_DLPFC["1000", ]) - mean(mat_runtimes_DLPFC["500", ])) / (1000 - 500)
# calculate linear and cubic lines starting from first point
n_spots <- as.numeric(df_runtimes_DLPFC$n_spots)
lin <- mean(mat_runtimes_DLPFC[1, ]) + (slope_DLPFC * (n_spots - min(n_spots)))
cub <- lin + (lin - lin[1])^3
df_runtimes_DLPFC <- cbind(
  df_runtimes_DLPFC, 
  linear = lin, 
  cubic = cub
)

# calculate linear slope
slope_mouseHPC <- 
  (mean(mat_runtimes_mouseHPC["40000", ]) - mean(mat_runtimes_mouseHPC["20000", ])) / (40000 - 20000)
# calculate linear and cubic lines starting from first point
n_spots <- as.numeric(df_runtimes_mouseHPC$n_spots)
lin <- mean(mat_runtimes_mouseHPC[1, ]) + (slope_mouseHPC * (n_spots - min(n_spots)))
cub <- lin + (lin - lin[1])^3
df_runtimes_mouseHPC <- cbind(
  df_runtimes_mouseHPC, 
  linear = lin, 
  cubic = cub
)


# --------------
# generate plots
# --------------

# DLPFC dataset

df <- pivot_longer(df_runtimes_DLPFC, cols = starts_with("seed"), 
                   names_to = "seed", values_to = "runtime")
df$n_spots <- as.numeric(df$n_spots)
df$dataset <- as.factor(df$dataset)

x_vals <- c(0, unique(df$n_spots))

pal <- "purple3"

df_trends <- pivot_longer(df_runtimes_DLPFC[, c("n_spots", "linear", "cubic")], 
                          cols = c("linear", "cubic"), 
                          names_to = "trend", values_to = "trend_val")
df_trends$n_spots <- as.numeric(df_trends$n_spots)
df_trends$trend <- factor(df_trends$trend, levels = c("linear", "cubic"))

# seed for geom_jitter so no points missing
set.seed(6)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 100, size = 1.5, alpha = 0.75) + 
  geom_point(data = df_trends, aes(x = n_spots, y = trend_val, group = NULL), 
             color = "black", alpha = 0.5) + 
  geom_line(data = df_trends, aes(x = n_spots, y = trend_val, linetype = trend, group = NULL), 
            color = "black", alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  ylim(c(0, 8.6)) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_nnSVG_DLPFC_singlegene.png"), width = 6, height = 4.5)


# mouseHPC dataset

df <- pivot_longer(df_runtimes_mouseHPC, cols = starts_with("seed"), 
                   names_to = "seed", values_to = "runtime")
df$n_spots <- as.numeric(df$n_spots)
df$dataset <- as.factor(df$dataset)

x_vals <- c(0, unique(df$n_spots))

pal <- "red"

df_trends <- pivot_longer(df_runtimes_mouseHPC[, c("n_spots", "linear", "cubic")], 
                          cols = c("linear", "cubic"), 
                          names_to = "trend", values_to = "trend_val")
df_trends$n_spots <- as.numeric(df_trends$n_spots)
df_trends$trend <- factor(df_trends$trend, levels = c("linear", "cubic"))

# seed for geom_jitter so no points missing
set.seed(6)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(width = 2000, lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 1000, size = 1.5, alpha = 0.75) + 
  geom_point(data = df_trends, aes(x = n_spots, y = trend_val, group = NULL), 
             color = "black", alpha = 0.5) + 
  geom_line(data = df_trends, aes(x = n_spots, y = trend_val, linetype = trend, group = NULL), 
            color = "black", alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  ylim(c(0, 152)) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_nnSVG_mouseHPC_singlegene.png"), width = 6.25, height = 4.5)

