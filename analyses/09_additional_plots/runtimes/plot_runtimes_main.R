#########################
# Script to plot runtimes
# Lukas Weber, May 2022
#########################

# runtimes for main results
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)


library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)


# directory to save plots
dir_plots <- here(file.path("plots", "runtimes"))


# ------------
# load results
# ------------

# runtimes
res_runtimes <- c(
  humanDLPFC_nnSVG = metadata(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds")))[["runtime"]], 
  mouseHPC_nnSVG = metadata(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG.rds")))[["runtime"]], 
  mouseOB_nnSVG = metadata(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG.rds")))[["runtime"]], 
  humanDLPFC_SPARKX = metadata(readRDS(here("outputs", "results", "spe_humanDLPFC_SPARKX.rds")))[["runtime"]], 
  mouseHPC_SPARKX = metadata(readRDS(here("outputs", "results", "spe_mouseHPC_SPARKX.rds")))[["runtime"]], 
  mouseOB_SPARKX = metadata(readRDS(here("outputs", "results", "spe_mouseOB_SPARKX.rds")))[["runtime"]]
)
res_runtimes

# numbers of spots
n_spots <- c(
  humanDLPFC = ncol(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  mouseHPC = ncol(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG.rds"))), 
  mouseOB = ncol(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG.rds")))
)
n_spots

# numbers of genes
n_genes <- c(
  humanDLPFC = nrow(readRDS(here("outputs", "results", "spe_humanDLPFC_nnSVG.rds"))), 
  mouseHPC = nrow(readRDS(here("outputs", "results", "spe_mouseHPC_nnSVG.rds"))), 
  mouseOB = nrow(readRDS(here("outputs", "results", "spe_mouseOB_nnSVG.rds")))
)
n_genes


# -------------
# plot runtimes
# -------------

df <- data.frame(
  dataset = gsub("_.*$", "", names(res_runtimes)), 
  method = gsub("^.*_", "", names(res_runtimes)), 
  runtime = unname(res_runtimes), 
  n_spots = rep(n_spots, times = 2)
  ) %>% 
  mutate(dataset = factor(dataset, levels = c("mouseOB", "humanDLPFC", "mouseHPC"))) %>% 
  mutate(method = gsub("SPARKX", "SPARK-X", method)) %>% 
  mutate(method = factor(method, levels = c("nnSVG", "SPARK-X")))

names(n_spots) <- paste0(names(n_spots), "\n(", n_spots, " spots)")


# seed for placement of text labels
set.seed(1)
ggplot(df, aes(x = n_spots, y = runtime, color = method, shape = dataset, 
               label = round(runtime))) + 
  geom_point(stroke = 1.25, size = 1.5) + 
  scale_shape_manual(values = c(1, 2, 3)) + 
  scale_color_manual(values = c("blue3", "deepskyblue2")) + 
  geom_text_repel(nudge_x = 100, size = 3.25, segment.color = NA, box.padding = 0.4, 
                  show.legend = FALSE) + 
  scale_x_continuous(breaks = n_spots) + 
  scale_y_log10() + 
  labs(x = "dataset", 
       y = "runtime (sec)") + 
  ggtitle("Runtimes") + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank())

fn <- file.path(dir_plots, "runtimes_main")
ggsave(paste0(fn, ".pdf"), width = 6, height = 4)
ggsave(paste0(fn, ".png"), width = 6, height = 4)

