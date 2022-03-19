##################################################
# Script for null simulations: preprocessing steps
# Lukas Weber, Mar 2022
##################################################

# dataset: Visium human DLPFC


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
set.seed(123)
qclus <- quickCluster(spe)
spe <- computeSumFactors(spe, cluster = qclus)
spe <- logNormCounts(spe)

assayNames(spe)


# ---------------------------
# permute spatial coordinates
# ---------------------------

# permute spatial coordinates to remove spatial correlation

set.seed(123)
ix_permute <- sample(seq_len(ncol(spe)))

stopifnot(length(ix_permute) == nrow(spatialCoords(spe)))

spatialCoords(spe) <- spatialCoords(spe)[ix_permute, ]
rownames(spatialCoords(spe)) <- rownames(colData(spe))


# -----------
# save object
# -----------

fn <- here("outputs", "null_sims", "spe_humanDLPFC_nullSim_preprocessed.rds")
saveRDS(spe, file = fn)

