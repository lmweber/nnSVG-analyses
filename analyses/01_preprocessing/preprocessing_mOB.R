###################################
# Script to run preprocessing steps
# Lukas Weber, Mar 2022
###################################

# dataset: ST mouse OB


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


# filter low-expressed and mitochondrial genes
# using gene filtering function from nnSVG package
# note: using higher filtering parameters for ST platform due to higher number
# of cells per spot (compared to Visium)
# note: mitochondrial genes have already been filtered out
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 5, 
  filter_genes_pcspots = 1, 
  filter_mito = FALSE
)

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

fn <- here("outputs", "SPE", "spe_mOB_preprocessed.rds")
saveRDS(spe, file = fn)

