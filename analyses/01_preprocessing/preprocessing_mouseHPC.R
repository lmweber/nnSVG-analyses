###################################
# Script to run preprocessing steps
# Lukas Weber, updated Mar 2022
###################################

# dataset: Slide-seqV2 mouse HPC


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
spe <- SlideSeqV2_mouseHPC()
dim(spe)


# -------------
# preprocessing
# -------------

# remove spots with NA cell type labels
spe <- spe[, !is.na(colData(spe)$celltype)]
dim(spe)

table(colData(spe)$celltype)


# filter low-expressed and mitochondrial genes
# using gene filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 1, 
  filter_genes_pcspots = 1, 
  filter_mito = TRUE
)

dim(spe)


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

fn <- here("outputs", "preprocessed", "spe_mouseHPC_preprocessed.rds")
saveRDS(spe, file = fn)

