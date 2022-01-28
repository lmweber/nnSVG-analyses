#################################
# Script to calculate evaluations
# Lukas Weber, Jan 2022
#################################

library(SpatialExperiment)
library(here)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(ggsci)


# ------------
# load results
# ------------

# DLPFC and mOB datasets

res_list <- list(
  DLPFC_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_DLPFC_nnSVG.rds"))), 
  DLPFC_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_DLPFC_SPARKX.rds"))), 
  mOB_nnSVG = rowData(readRDS(here("outputs", "results", "nnSVG", "spe_mOB_nnSVG.rds"))), 
  mOB_SPARKX = rowData(readRDS(here("outputs", "results", "SPARKX", "spe_mOB_SPARKX.rds")))
)

