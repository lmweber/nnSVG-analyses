#########################
# SVGs BRISC example
# Lukas Weber, March 2021
#########################

# Similar to previous script using spNNGP but using BRISC instead


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
spe <- spe[!is_mito, ]
# fit mean-variance relationship
dec <- modelGeneVar(spe)
# select top HVGs
top_hvgs <- getTopHVGs(dec, prop = 0.1)


# ----------------
# Fit BRISC models
# ----------------

# fitting one model per gene
# BRISC requires using logcounts
# use default parameters (including x = NULL, which includes intercept by default)

library(BRISC)

y <- logcounts(spe)
dim(y)

coords <- spatialCoords(spe)
dim(coords)
head(coords)
# scale coordinates: required since our scale is too large
# to do: scale axes proportionally
coords <- apply(coords, 2, function(col) (col - min(col)) / (max(col) - min(col)))


# ------------------------------------------
# Example: fit BRISC model for a single gene
# ------------------------------------------

# extract responses for one gene
ix <- which(rowData(spe)$gene_name == "PCP4")
ix

# convert to dense vector
y_pcp4 <- y[ix, ]

head(y_pcp4)
stopifnot(length(y_pcp4) == nrow(coords))

# fit BRISC model for a single gene
# runtime: ~4 sec for one gene
runtime <- system.time({
  out_brisc <- BRISC_estimation(coords = coords, y = y_pcp4, x = NULL, n.neighbors = 15)
})

runtime


# some outputs: see ?BRISC_estimation for details

str(out_brisc, max.level = 1)

# estimated parameters
out_brisc$Theta

# tau squared parameter
out_brisc$Theta["tau.sq"]
# sigma squared parameter
out_brisc$Theta["sigma.sq"]

# fraction of spatial variance (FSV) using this parameterization (analogous to FSV in SpatialDE)
tau_sq <- out_brisc$Theta["tau.sq"]
sigma_sq <- out_brisc$Theta["sigma.sq"]
fsv <- sigma_sq / (sigma_sq + tau_sq)
fsv


# BRISC object
str(out_brisc$BRISC_Object)


# ------------------------------
# Bootstrap confidence intervals
# ------------------------------

# bootstrap with BRISC
runtime_bootstrap <- system.time({
  out_brisc_boot <- BRISC_bootstrap(out_brisc)
})

out_brisc_boot$confidence.interval

out_brisc_boot$boot.time

