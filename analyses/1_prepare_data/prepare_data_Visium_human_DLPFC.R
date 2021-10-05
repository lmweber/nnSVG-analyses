#################################
# Script to load and prepare data
# Lukas Weber, Oct 2021
#################################

# dataset: 10x Genomics Visium human dorsolateral prefrontal cortex (DLPFC)

# references:
# Maynard and Collado-Torres et al. (2021): https://www.nature.com/articles/s41593-020-00787-0
# spatialLIBD Bioconductor package: http://bioconductor.org/packages/spatialLIBD


# to run in interactive session on cluster:
# module load conda_R/4.1.x
# Rscript filename.R


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(dplyr)
library(here)


# ---------
# load data
# ---------

# load data object from STexampleData package

spe <- Visium_humanDLPFC()
dim(spe)
assayNames(spe)


# -------------------
# preprocessing steps
# -------------------

# using 'preprocessSVG()' function from nnSVG package

# preprocessing steps:
# - select spots over tissue
# - filter low-expressed genes
# - filter mitochondrial genes
# - calculate deviance residuals and/or logcounts

# set seed for reproducibility
set.seed(123)
spe <- preprocessSVG(spe, in_tissue = TRUE, 
                     filter_genes = 5, filter_mito = TRUE)

dim(spe)
assayNames(spe)
assays(spe)[["binomial_deviance_residuals"]][1:6, 1:6]
logcounts(spe)[1:6, 1:6]


# -----------
# save object
# -----------

file <- here("outputs", "SPE", "spe_DLPFC.rds")
saveRDS(spe, file = file)

