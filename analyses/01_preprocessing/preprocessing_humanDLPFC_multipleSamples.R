###################################
# Script to run preprocessing steps
# Lukas Weber, Feb 2023
###################################

# dataset: Visium human DLPFC (multiple samples)


library(SpatialExperiment)
library(spatialLIBD)
library(ExperimentHub)
library(nnSVG)
library(scater)
library(scran)
library(here)


# ---------
# load data
# ---------

# load dataset from spatialLIBD package as SpatialExperiment object

ehub <- ExperimentHub()

spe <- fetch_data(type = "spe", eh = ehub)

dim(spe)
table(colData(spe)$sample_id)


# -------------
# preprocessing
# -------------

# note: preprocessing including spot-level quality control (QC) and calculating
# logcounts has already been done in this dataset

# check object contains only spots over tissue
table(colData(spe)$in_tissue)

# check object contains logcounts
assayNames(spe)

# check if object contains mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_multipleSamples_preprocessed.rds")
saveRDS(spe, file = fn)

