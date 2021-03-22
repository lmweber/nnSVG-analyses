#########################
# SVGs spNNGP example
# Lukas Weber, March 2021
#########################

# ---------------------------------------------
# Load data, preprocessing, calculate logcounts
# ---------------------------------------------

# Code copied from current version of OSTA

# LOAD DATA

library(SpatialExperiment)
library(STexampleData)
spe <- load_data("Visium_humanDLPFC")

# QUALITY CONTROL (QC)

library(scater)
# subset to keep only spots over tissue
spe <- spe[, spatialData(spe)$in_tissue == 1]
# identify mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
# calculate per-spot QC metrics
spe <- addPerCellQC(spe, subsets = list(mito = is_mito))
# select QC thresholds
qc_lib_size <- colData(spe)$sum < 500
qc_detected <- colData(spe)$detected < 250
qc_mito <- colData(spe)$subsets_mito_percent > 30
qc_cell_count <- colData(spe)$cell_count > 12
# combined set of discarded spots
discard <- qc_lib_size | qc_detected | qc_mito | qc_cell_count
colData(spe)$discard <- discard
# filter low-quality spots
spe <- spe[, !colData(spe)$discard]

# NORMALIZATION

library(scran)
# quick clustering for pool-based size factors
set.seed(123)
qclus <- quickCluster(spe)
# calculate size factors
spe <- computeSumFactors(spe, cluster = qclus)
# calculate logcounts (log-transformed normalized counts)
spe <- logNormCounts(spe)

# FEATURE SELECTION

# remove mitochondrial genes
spe_nomito <- spe[!is_mito, ]
# fit mean-variance relationship
dec <- modelGeneVar(spe_nomito)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)


# -----------------
# Fit spNNGP models
# -----------------

# fitting a single model for one gene for now for testing purposes
# there are several tricks for speeding up a loop over all genes, including:
# - parallelization (one thread per gene)
# - re-use nearest neighbors info using argument 'neighbor.info'
# - filter out low-expressed genes and/or restrict to a set of top HVGs

# using logcounts for now
# but spNNGP can also fit discrete responses (counts)


# using implementation in spatzli package
# stores statistics, ranks, and runtimes in rowData of spe object

library(spatzli)

spe <- rankSVGsNNGP(spe, n_threads = 10)


