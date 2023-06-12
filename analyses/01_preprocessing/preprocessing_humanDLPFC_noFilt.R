###################################
# Script to run preprocessing steps
# Lukas Weber, updated Jan 2023
###################################

# dataset: Visium human DLPFC, without filtering low-expressed genes


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
spe <- Visium_humanDLPFC()
dim(spe)


# -------------
# preprocessing
# -------------

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]

dim(spe)


# spot-level quality control (QC) using scater package

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 500
qc_detected <- colData(spe)$detected < 250
qc_mito <- colData(spe)$subsets_mito_percent > 30
qc_cell_count <- colData(spe)$cell_count > 12
# spots to discard
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
table(discard)
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]

dim(spe)


# filter mitochondrial genes
spe <- spe[!is_mito, ]

dim(spe)


# filter zero-expressed genes
is_zero <- rowSums(counts(spe)) == 0
table(is_zero)

spe <- spe[!is_zero, ]

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

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed_noFilt.rds")
saveRDS(spe, file = fn)

