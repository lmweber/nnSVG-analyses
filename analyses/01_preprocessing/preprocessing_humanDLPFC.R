###################################
# Script to run preprocessing steps
# Lukas Weber, Mar 2022
###################################

# dataset: human DLPFC


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scry)
library(scran)
library(here)


# ---------
# load data
# ---------

# load dataset as SpatialExperiment object from STexampleData package
spe <- Visium_humanDLPFC()
dim(spe)


# -------------
# preprocessing
# -------------

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

dim(spe)

# filter low-expressed and mitochondrial genes
# using filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 3, 
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)

dim(spe)

# calculate log-transformed normalized counts using scran package
set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)

assayNames(spe)

# calculate deviance residuals using scry package
spe <- nullResiduals(
  spe, 
  assay = "counts", 
  fam = "binomial", 
  type = "deviance"
)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "SPE", "spe_humanDLPFC_preprocessed.rds")
saveRDS(spe, file = fn)

