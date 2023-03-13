#####################################################
# Script to plot runtimes for scalability simulations
# Lukas Weber, Jan 2023
#####################################################

# with and without covariates


library(here)
library(dplyr)
library(tidyr)
library(ggplot2)


# ------------
# load results
# ------------

runtimes_nnSVG_mouseHPC_singlegene_withCovariates <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_nnSVG_mouseHPC_singlegene_withCovariates.rds"))
runtimes_nnSVG_mouseHPC_singlegene_noCovariates <- readRDS(here(
  "outputs", "scalability_sims", "runtimes_scalability_nnSVG_mouseHPC_singlegene_noCovariates.rds"))

mat_runtimes_mouseHPC_withCovariates <- do.call("rbind", runtimes_nnSVG_mouseHPC_singlegene_withCovariates)
colnames(mat_runtimes_mouseHPC_withCovariates) <- paste0("seed", 1:ncol(mat_runtimes_mouseHPC_withCovariates))
df_runtimes_mouseHPC_withCovariates <- cbind(
  n_spots = rownames(mat_runtimes_mouseHPC_withCovariates), 
  dataset = "mouseHPC_withCovariates", 
  as.data.frame(mat_runtimes_mouseHPC_withCovariates))

mat_runtimes_mouseHPC_noCovariates <- do.call("rbind", runtimes_nnSVG_mouseHPC_singlegene_noCovariates)
colnames(mat_runtimes_mouseHPC_noCovariates) <- paste0("seed", 1:ncol(mat_runtimes_mouseHPC_noCovariates))
df_runtimes_mouseHPC_noCovariates <- cbind(
  n_spots = rownames(mat_runtimes_mouseHPC_noCovariates), 
  dataset = "mouseHPC_noCovariates", 
  as.data.frame(mat_runtimes_mouseHPC_noCovariates))


# --------------------------------
# calculate linear and cubic lines
# --------------------------------

# calculate linear slope
slope_mouseHPC_withCovariates <- 
  (median(mat_runtimes_mouseHPC_withCovariates["15003", ]) - median(mat_runtimes_mouseHPC_withCovariates["500", ])) / (15003 - 500)
# calculate linear and cubic lines starting from first point
n_spots <- as.numeric(df_runtimes_mouseHPC_withCovariates$n_spots)
lin <- median(mat_runtimes_mouseHPC_withCovariates[1, ]) + (slope_mouseHPC_withCovariates * (n_spots - min(n_spots)))
cub <- lin + (lin - lin[1])^3
df_runtimes_mouseHPC_withCovariates <- cbind(
  df_runtimes_mouseHPC_withCovariates, 
  linear = lin, 
  cubic = cub
)

# calculate linear slope
slope_mouseHPC_noCovariates <- 
  (median(mat_runtimes_mouseHPC_noCovariates["15003", ]) - median(mat_runtimes_mouseHPC_noCovariates["500", ])) / (15003 - 500)
# calculate linear and cubic lines starting from first point
n_spots <- as.numeric(df_runtimes_mouseHPC_noCovariates$n_spots)
lin <- median(mat_runtimes_mouseHPC_noCovariates[1, ]) + (slope_mouseHPC_noCovariates * (n_spots - min(n_spots)))
cub <- lin + (lin - lin[1])^3
df_runtimes_mouseHPC_noCovariates <- cbind(
  df_runtimes_mouseHPC_noCovariates, 
  linear = lin, 
  cubic = cub
)


# --------------
# generate plots
# --------------

# mouseHPC dataset: with covariates

df <- pivot_longer(df_runtimes_mouseHPC_withCovariates, cols = starts_with("seed"), 
                   names_to = "seed", values_to = "runtime")
df$n_spots <- as.numeric(df$n_spots)
df$dataset <- as.factor(df$dataset)

x_vals <- c(0, unique(df$n_spots))

pal <- "firebrick3"

df_trends <- pivot_longer(df_runtimes_mouseHPC_withCovariates[, c("n_spots", "linear", "cubic")], 
                          cols = c("linear", "cubic"), 
                          names_to = "trend", values_to = "trend_val")
df_trends$n_spots <- as.numeric(df_trends$n_spots)
df_trends$trend <- factor(df_trends$trend, levels = c("linear", "cubic"))

# keep only linear trend
df_trends <- 
  filter(df_trends, trend == "linear") %>% 
  mutate(trend = droplevels(trend))


# seed for geom_jitter
set.seed(1)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(width = 2000, lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 1000, size = 1, alpha = 0.75) + 
  geom_point(data = df_trends, aes(x = n_spots, y = trend_val, group = NULL), 
             color = "black", alpha = 0.5) + 
  geom_line(data = df_trends, aes(x = n_spots, y = trend_val, linetype = trend, group = NULL), 
            color = "black", alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  scale_linetype_manual(values = "dashed") + 
  ylim(c(0, 100)) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG, single gene") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability_sims", "runtimes_nnSVG_mouseHPC_singlegene_withCovariates.png"), 
       width = 6, height = 4)


# mouseHPC dataset: no covariates

df <- pivot_longer(df_runtimes_mouseHPC_noCovariates, cols = starts_with("seed"), 
                   names_to = "seed", values_to = "runtime")
df$n_spots <- as.numeric(df$n_spots)
df$dataset <- as.factor(df$dataset)

x_vals <- c(0, unique(df$n_spots))

pal <- "firebrick3"
  
df_trends <- pivot_longer(df_runtimes_mouseHPC_noCovariates[, c("n_spots", "linear", "cubic")], 
                          cols = c("linear", "cubic"), 
                          names_to = "trend", values_to = "trend_val")
df_trends$n_spots <- as.numeric(df_trends$n_spots)
df_trends$trend <- factor(df_trends$trend, levels = c("linear", "cubic"))

# keep only linear trend
df_trends <- 
  filter(df_trends, trend == "linear") %>% 
  mutate(trend = droplevels(trend))


# seed for geom_jitter
set.seed(1)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(width = 2000, lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 1000, size = 1, alpha = 0.75) + 
  geom_point(data = df_trends, aes(x = n_spots, y = trend_val, group = NULL), 
             color = "black", alpha = 0.5) + 
  geom_line(data = df_trends, aes(x = n_spots, y = trend_val, linetype = trend, group = NULL), 
            color = "black", alpha = 0.5) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  scale_linetype_manual(values = "dashed") + 
  ylim(c(0, 100)) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG, single gene") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(size = 7, angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability_sims", "runtimes_nnSVG_mouseHPC_singlegene_noCovariates.png"), 
       width = 6, height = 4)

