#################################
# Script to calculate evaluations
# Lukas Weber, Feb 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(ggrepel)


# ------------
# load results
# ------------

# scalable methods: nnSVG, SPARK-X, HVGs

res_list <- list(
  mouseHPC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds"))), 
  mouseHPC_nnSVG_covariate = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_clusters.rds"))), 
  mouseHPC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX.rds"))), 
  mouseHPC_SPARKX_covariate = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_SlideSeqHippo_SPARKX_clusters.rds"))), 
  mouseHPC_HVGs = rowData(readRDS(here("outputs", "results", "HVGs", "spe_SlideSeqHippo_HVGs.rds")))
)

# add method names to all columns except gene IDs and gene names
colnames(res_list[["mouseHPC_nnSVG"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG"]]), "_nnSVG")[-1]
colnames(res_list[["mouseHPC_nnSVG_covariate"]])[-1] <- paste0(colnames(res_list[["mouseHPC_nnSVG_covariate"]]), "_nnSVG_covariate")[-1]
colnames(res_list[["mouseHPC_SPARKX"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX"]]), "_SPARKX")[-1]
colnames(res_list[["mouseHPC_SPARKX_covariate"]])[-1] <- paste0(colnames(res_list[["mouseHPC_SPARKX_covariate"]]), "_SPARKX_covariate")[-1]
colnames(res_list[["mouseHPC_HVGs"]])[-1] <- paste0(colnames(res_list[["mouseHPC_HVGs"]]), "_HVGs")[-1]

