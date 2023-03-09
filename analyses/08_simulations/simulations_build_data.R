#######################################
# Simulations: build simulated datasets
# Lukas Weber, Mar 2023
#######################################

# this script builds simulated datasets based on empirical parameters from the
# Visium human DLPFC dataset


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scran)
library(scater)
library(here)


dir_sims <- here("outputs", "simulations")


# ---------
# load data
# ---------

spe <- Visium_humanDLPFC()

spe_full <- spe
dim(spe_full)

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)


# ---------------------------
# preprocessing and logcounts
# ---------------------------

# starting from spots over tissue

# filter low-expressed and mitochondrial genes
# using gene filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 3, 
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)
dim(spe)

# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
assayNames(spe)


# --------------------
# empirical parameters
# --------------------

# calculate empirical simulation parameters based on expression of gene MOBP in
# white matter in human DLPFC dataset

# identify gene and white matter
ix_gene <- which(rowData(spe)$gene_name == "MOBP")
is_gene <- rowData(spe)$gene_name == "MOBP"
is_wm <- colData(spe)$ground_truth == "WM"
is_wm[is.na(is_wm)] <- FALSE

# parameters: expressed vs. not expressed
par_meanLogcountsExpressed <- mean(logcounts(spe)[is_gene, is_wm])
par_varLogcountsExpressed <- var(logcounts(spe)[is_gene, is_wm])
par_meanLogcountsNotExpressed <- mean(logcounts(spe)[is_gene, !is_wm])
par_varLogcountsNotExpressed <- var(logcounts(spe)[is_gene, !is_wm])

# parameters: sparsity in expressed vs. not expressed regions
par_sparsityExpressedRegion <- mean(logcounts(spe)[is_gene, is_wm] == 0)
par_sparsityNotExpressedRegion <- mean(logcounts(spe)[is_gene, !is_wm] == 0)

# check
par_meanLogcountsExpressed  # 2.797321
par_varLogcountsExpressed  # 0.6211075
par_meanLogcountsNotExpressed  # 0.4786561
par_varLogcountsNotExpressed  # 0.5436924
par_sparsityExpressedRegion  # 0.01949318
par_sparsityNotExpressedRegion  # 0.6353167


# ---------------------
# additional parameters
# ---------------------

# additional simulation parameters

# radius of expressed regions
smallBandwidthRadius <- 0.05
mediumBandwidthRadius <- 0.125
largeBandwidthRadius <- 0.25

# number of non-expressed and expressed genes
n_genes <- 1000
n_genes_nonExpressed <- 900
n_genes_expressed <- 100

# strength of expression (relative to MOBP in WM region)
lowExpression <- 1/3
mediumExpression <- 2/3
fullExpression <- 1

# shuffling of coordinates for ablation simulations
propShufflePerIteration <- 0.1
numShuffleIterations <- 10


# -------------------
# spatial coordinates
# -------------------

# using spatial coordinates from Visium slide in human DLPFC dataset

coords <- spatialCoords(spe_full)
colnames(coords) <- c("x", "y")

# scale to range 0 to 1 in each dimension
coords[, 1] <- (coords[, 1] - min(coords[, 1])) / (max(coords[, 1]) - min(coords[, 1]))
coords[, 2] <- (coords[, 2] - min(coords[, 2])) / (max(coords[, 2]) - min(coords[, 2]))

rownames(coords) <- NULL

head(coords)
dim(coords)
summary(coords)

# number of spots
n_spots <- nrow(coords)


# ----------------------------------
# create masks for expressed regions
# ----------------------------------

# create masks to identify spatial coordinates within expressed regions

# large bandwidth

centersLargeBandwidth_x <- 0.5
centersLargeBandwidth_y <- 0.5
radiusLargeBandwidth <- largeBandwidthRadius

maskLargeBandwidth <- rep(FALSE, n_spots)

for (i in seq_len(n_spots)) {
  x <- coords[i, "x"]
  y <- coords[i, "y"]
  if (((x - centersLargeBandwidth_x)^2 + (y - centersLargeBandwidth_y)^2) < radiusLargeBandwidth^2) {
    maskLargeBandwidth[i] <- TRUE
  }
}


# medium bandwidth

centersMediumBandwidth_x <- c(0.275, 0.725, 0.275, 0.725)
centersMediumBandwidth_y <- c(0.275, 0.275, 0.725, 0.725)
radiusMediumBandwidth <- mediumBandwidthRadius

maskMediumBandwidth <- rep(FALSE, n_spots)

for (k in seq_along(centersMediumBandwidth_x)) {
  for (i in seq_len(n_spots)) {
    x <- coords[i, "x"]
    y <- coords[i, "y"]
    if (((x - centersMediumBandwidth_x[k])^2 + (y - centersMediumBandwidth_y[k])^2) < radiusMediumBandwidth^2) {
      maskMediumBandwidth[i] <- TRUE
    }
  }
}


# small bandwidth

centersSmallBandwidth_x <- c(0.175, 0.5, 0.825, 0.175, 0.5, 0.825, 0.175, 0.5, 0.825)
centersSmallBandwidth_y <- c(0.175, 0.175, 0.175, 0.5, 0.5, 0.5, 0.825, 0.825, 0.825)
radiusSmallBandwidth <- smallBandwidthRadius

maskSmallBandwidth <- rep(FALSE, n_spots)

for (k in seq_along(centersSmallBandwidth_x)) {
  for (i in seq_len(n_spots)) {
    x <- coords[i, "x"]
    y <- coords[i, "y"]
    if (((x - centersSmallBandwidth_x[k])^2 + (y - centersSmallBandwidth_y[k])^2) < radiusSmallBandwidth^2) {
      maskSmallBandwidth[i] <- TRUE
    }
  }
}


# -----------------------------------------------------------------------------
# functions to create log-expression values for noise genes and expressed genes
# -----------------------------------------------------------------------------

# note: using global arguments for parameters within functions
# note: set random seed before running functions

# noise genes

fn_createNoiseGene <- function() {
  # sample random values
  logexpr <- rnorm(n_spots, 
                   mean = par_meanLogcountsNotExpressed, 
                   sd = sqrt(par_varLogcountsNotExpressed))
  # truncate values to achieve sparsity
  q <- quantile(logexpr, par_sparsityNotExpressedRegion)
  logexpr[logexpr < q] <- 0
  # remove any negative values
  logexpr[logexpr < 0] <- 0
  # return values
  logexpr
}


# expressed genes

fn_createExpressedGene <- function() {
  # sample random values
  logexpr <- rnorm(n_spots, 
                   mean = par_meanLogcountsExpressed, 
                   sd = sqrt(par_varLogcountsExpressed))
  # truncate values to achieve sparsity
  q <- quantile(logexpr, par_sparsityExpressedRegion)
  logexpr[logexpr < q] <- 0
  # remove any negative values
  logexpr[logexpr < 0] <- 0
  # return values
  logexpr
}


# -------------------------
# create simulated datasets
# -------------------------

# gene names

gene_names <- sprintf("gene%04d", seq_len(n_genes))


# function to build simulated dataset

fn_buildSimulatedData <- function(mask, expressionStrength) {
  
  rowdata <- DataFrame(
    gene_name = gene_names, 
    expressed = c(rep(FALSE, n_genes_nonExpressed), rep(TRUE, n_genes_expressed)), 
    expression_strength = c(rep(0, n_genes_nonExpressed), rep(expressionStrength, n_genes_expressed))
  )
  coldata <- DataFrame(
    mask = mask
  )
  spatialcoords <- coords
  
  logcounts <- matrix(NA, nrow = n_genes, ncol = n_spots)
  set.seed(123)
  for (g in seq_len(n_genes)) {
    if (!rowdata$expressed[g]) {
      logcounts[g, ] <- fn_createNoiseGene()
    } else if (rowdata$expressed[g]) {
      logcounts[g, ] <- fn_createExpressedGene() * expressionStrength
    }
  }
  
  SpatialExperiment(
    rowData = rowdata, 
    colData = coldata, 
    spatialCoords = spatialcoords, 
    assays = list(logcounts = logcounts)
  )
}


# build datasets for combinations of parameters

set.seed(101)
spe_sim_largeBandwidth_fullExpr <- fn_buildSimulatedData(maskLargeBandwidth, fullExpression)
set.seed(102)
spe_sim_largeBandwidth_medExpr <- fn_buildSimulatedData(maskLargeBandwidth, mediumExpression)
set.seed(103)
spe_sim_largeBandwidth_lowExpr <- fn_buildSimulatedData(maskLargeBandwidth, lowExpression)

set.seed(104)
spe_sim_medBandwidth_fullExpr <- fn_buildSimulatedData(maskMediumBandwidth, fullExpression)
set.seed(105)
spe_sim_medBandwidth_medExpr <- fn_buildSimulatedData(maskMediumBandwidth, mediumExpression)
set.seed(106)
spe_sim_medBandwidth_lowExpr <- fn_buildSimulatedData(maskMediumBandwidth, lowExpression)

set.seed(107)
spe_sim_smallBandwidth_fullExpr <- fn_buildSimulatedData(maskSmallBandwidth, fullExpression)
set.seed(108)
spe_sim_smallBandwidth_medExpr <- fn_buildSimulatedData(maskSmallBandwidth, mediumExpression)
set.seed(109)
spe_sim_smallBandwidth_lowExpr <- fn_buildSimulatedData(maskSmallBandwidth, lowExpression)


# -----------------------
# save simulated datasets
# -----------------------

saveRDS(spe_sim_largeBandwidth_fullExpr, file.path(dir_sims, "spe_sim_largeBandwidth_fullExpr.rds"))
saveRDS(spe_sim_largeBandwidth_medExpr, file.path(dir_sims, "spe_sim_largeBandwidth_medExpr.rds"))
saveRDS(spe_sim_largeBandwidth_lowExpr, file.path(dir_sims, "spe_sim_largeBandwidth_lowExpr.rds"))

saveRDS(spe_sim_medBandwidth_fullExpr, file.path(dir_sims, "spe_sim_medBandwidth_fullExpr.rds"))
saveRDS(spe_sim_medBandwidth_medExpr, file.path(dir_sims, "spe_sim_medBandwidth_medExpr.rds"))
saveRDS(spe_sim_medBandwidth_lowExpr, file.path(dir_sims, "spe_sim_medBandwidth_lowExpr.rds"))

saveRDS(spe_sim_smallBandwidth_fullExpr, file.path(dir_sims, "spe_sim_smallBandwidth_fullExpr.rds"))
saveRDS(spe_sim_smallBandwidth_medExpr, file.path(dir_sims, "spe_sim_smallBandwidth_medExpr.rds"))
saveRDS(spe_sim_smallBandwidth_lowExpr, file.path(dir_sims, "spe_sim_smallBandwidth_lowExpr.rds"))

