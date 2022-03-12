###################################
# Script to run preprocessing steps
# Lukas Weber, Mar 2022
###################################

# dataset: Slide-seqV2 mouse HPC, without filtering low-expressed genes


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
spe <- SlideSeqV2_mouseHPC()
dim(spe)


# -------------
# preprocessing
# -------------

# spot-level quality control (QC) using scater package

# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 20
qc_detected <- colData(spe)$detected < 15
qc_mito <- colData(spe)$subsets_mito_percent > 30
# spots to discard
discard <- qc_lib_size | qc_detected | qc_mito
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


# calculate log-transformed normalized counts using scran package
set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "SPE", "spe_mouseHPC_preprocessed_nofilt.rds")
saveRDS(spe, file = fn)

