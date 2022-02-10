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

# seed for geom_jitter so no points missing
set.seed(6)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 100, size = 1.5, alpha = 0.75) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  ylim(c(0, max(df$runtime))) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_nnSVG_DLPFC_singlegene.png"), width = 5.5, height = 4.5)


# mouseHPC dataset

df <- pivot_longer(df_runtimes_mouseHPC, cols = starts_with("seed"), 
                   names_to = "seed", values_to = "runtime")
df$n_spots <- as.numeric(df$n_spots)
df$dataset <- as.factor(df$dataset)

x_vals <- c(0, unique(df$n_spots))

pal <- "red"

# seed for geom_jitter so no points missing
set.seed(6)
ggplot(df, aes(x = n_spots, y = runtime, color = dataset, group = n_spots)) + 
  geom_boxplot(lwd = 0.75, outlier.shape = NA) + 
  geom_jitter(width = 100, size = 1.5, alpha = 0.75) + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  ylim(c(0, max(df$runtime))) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_nnSVG_mouseHPC_singlegene.png"), width = 5.5, height = 4.5)

