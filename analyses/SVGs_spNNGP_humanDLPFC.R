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

library(spNNGP)

y <- logcounts(spe)
dim(y)

coords <- spatialCoords(spe)
dim(coords)
head(coords)

# parameters from example in spatial statistics course
n.samples <- 2000
starting <- list("phi" = 3/0.5, "sigma.sq" = 50, "tau.sq" = 1)
tuning <- list("phi" = 0.05, "sigma.sq" = 0.05, "tau.sq" = 0.05)
p <- 1
priors <- list("beta.Norm" = list(rep(0, p), diag(1000, p)), 
               "phi.Unif" = c(3/1, 3/0.1), "sigma.sq.IG" = c(2, 2), 
               "tau.sq.IG" = c(2, 0.1))


# -------------------------------------------
# Example: fit spNNGP model for a single gene
# -------------------------------------------

# extract responses for one gene
ix <- which(rowData(spe)$gene_name == "PCP4")
ix

# convert to dense vector
y_pcp4 <- y[ix, ]

stopifnot(length(y) == nrow(coords))

# fit spNNGP model for a single gene
# runtime: around 30 sec (for one gene - but scaling linearly in number of spots)
out_spnngp <- spNNGP(y_pcp4 ~ 1, coords = coords, starting = starting, method = "latent", n.neighbors = 5, 
                     tuning = tuning, priors = priors, cov.model = "exponential", 
                     n.samples = n.samples, return.neighbor.info = TRUE, n.omp.threads = 1)

# some outputs
# see documentation ?spNNGP for complete list

# matrix of posterior samples for spatial random effects (rows = spots, columns = posterior samples)
dim(out_spnngp$p.w.samples)
out_spnngp$p.w.samples[1:6, 1995:2000]

# nearest neighbor info that can be re-used in loop (may only be a small fraction of total runtime though)
str(out_spnngp$neighbor.info, max.level = 1)

# total runtime (for one gene)
out_spnngp$run.time


# outputs that can be used for ranking SVGs

# sum of absolute values of medians of posterior samples for spatial random effects
library(matrixStats)
med_spraneff <- rowMedians(out_spnngp$p.w.samples)
length(med_spraneff)
# sum of absolute values across spots
sum_abs_med_spraneff <- sum(abs(med_spraneff))
sum_abs_med_spraneff


