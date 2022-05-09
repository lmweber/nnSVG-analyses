###################################
# Script to run preprocessing steps
# Lukas Weber, May 2022
###################################

# dataset: seqFISH mouse embryo


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scater)
library(scran)
library(scry)
library(here)


# ---------
# load data
# ---------

# load dataset as SpatialExperiment object from STexampleData package
spe <- seqFISH_mouseEmbryo()
dim(spe)


# -------------
# preprocessing
# -------------

# note: preprocessing has already been performed for this dataset
# gene filtering is not needed since this dataset contains a small set of targeted genes


# calculate deviance residuals using scry package
set.seed(123)
spe <- nullResiduals(
  spe, 
  assay = "counts", 
  fam = "binomial", 
  type = "deviance"
)

assayNames(spe)


# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_mouseEmbryo_preprocessed_noFilt.rds")
saveRDS(spe, file = fn)

