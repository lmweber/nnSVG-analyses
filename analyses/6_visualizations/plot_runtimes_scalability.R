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

runtimes_nnSVG_DLPFC_singlegene <- readRDS(here("outputs", "scalability", "runtimes_scalability_nnSVG_DLPFC_singlegene.rds"))
runtimes_nnSVG_mouseHPC_singlegene <- readRDS(here("outputs", "scalability", "runtimes_scalability_nnSVG_mouseHPC_singlegene.rds"))

x_vals <- sort(unique(as.numeric(c(0, 
  gsub("^n", "", names(runtimes_nnSVG_DLPFC_singlegene)), 
  gsub("^n", "", names(runtimes_nnSVG_mouseHPC_singlegene))))))
x_vals
x_nms <- paste0("n", x_vals)

df_runtimes <- data.frame(
  n_spots = x_vals, 
  DLPFC = NA, 
  mouseHPC = NA, 
  row.names = x_nms
)

df_runtimes[names(unlist(runtimes_nnSVG_DLPFC_singlegene)), "DLPFC"] <- unlist(runtimes_nnSVG_DLPFC_singlegene)
df_runtimes[names(unlist(runtimes_nnSVG_mouseHPC_singlegene)), "mouseHPC"] <- unlist(runtimes_nnSVG_mouseHPC_singlegene)

df_runtimes


# --------------
# generate plots
# --------------

# DLPFC dataset

df <- pivot_longer(na.omit(df_runtimes[, c("n_spots", "DLPFC")]), 
                   cols = "DLPFC", 
                   names_to = "dataset", values_to = "runtime")
df$dataset <- as.factor(df$dataset)

pal <- "purple3"

ggplot(df, aes(x = n_spots, y = runtime, color = dataset)) + 
  geom_line() + 
  geom_point() + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = c(0, df$n_spots)) + 
  ylim(c(0, max(df$runtime))) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG, single gene") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_DLPFC_singlegene.png"), width = 5.5, height = 4.5)


# mouseHPC dataset

df <- pivot_longer(na.omit(df_runtimes[, c("n_spots", "mouseHPC")]), 
                   cols = "mouseHPC", 
                   names_to = "dataset", values_to = "runtime")
df$dataset <- as.factor(df$dataset)

pal <- "red"

ggplot(df, aes(x = n_spots, y = runtime, color = dataset)) + 
  geom_line() + 
  geom_point() + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = c(0, df$n_spots)) + 
  ylim(c(0, max(df$runtime))) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: nnSVG, single gene") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_mouseHPC_singlegene.png"), width = 5.5, height = 4.5)

