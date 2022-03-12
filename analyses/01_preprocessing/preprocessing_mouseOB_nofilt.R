###################################
# Script to run preprocessing steps
# Lukas Weber, Mar 2022
###################################

# dataset: ST mouse OB, without filtering low-expressed genes


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scater)
library(scran)
library(here)


# ---------
# load data
# ---------

# load dataset as SpatialExperiment object from STexampleData package
spe <- ST_mouseOB()
dim(spe)


# -------------
# preprocessing
# -------------

# spot-level quality control (QC) using scater package

# identify mitochondrial genes
# note: mitochondrial genes have already been filtered out
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 500
qc_detected <- colData(spe)$detected < 250
# spots to discard
discard <- qc_lib_size | qc_detected
table(discard)
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]

dim(spe)


# filter zero-expressed genes
# note: no zero-expressed genes in this dataset
is_zero <- rowSums(counts(spe)) == 0
table(is_zero)


# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "SPE", "spe_mouseOB_preprocessed_nofilt.rds")
saveRDS(spe, file = fn)

