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

df_runtimes <- data.frame(
  n_spots = as.numeric(gsub("^n", "", names(runtimes_nnSVG_DLPFC_singlegene))), 
  nnSVG = unname(unlist(runtimes_nnSVG_DLPFC_singlegene))
)


# --------------
# generate plots
# --------------

df <- pivot_longer(df_runtimes, cols = "nnSVG", 
                   names_to = "method", values_to = "runtime")

#pal <- unname(palette.colors(8, palette = "Okabe-Ito"))
pal <- "navy"

x_vals <- sort(unique(c(seq(0, 3000, by = 1000), df$n_spots)))

ggplot(df, aes(x = n_spots, y = runtime, color = method)) + 
  geom_line() + 
  geom_point() + 
  scale_color_manual(values = pal) + 
  scale_x_continuous(breaks = x_vals) + 
  ylim(c(0, max(df$runtime))) + 
  labs(x = "number of spots", 
       y = "runtime (sec)") + 
  ggtitle("Scalability: single gene") + 
  theme_bw() + 
  theme(panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

ggsave(here("plots", "scalability", "runtimes_singlegene.png"), width = 5.5, height = 4.5)

