###################################
# Script to run preprocessing steps
# Lukas Weber, Mar 2022
###################################

# dataset: Slide-seqV2 mouse HPC


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


# filter low-expressed and mitochondrial genes
# using gene filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 2, 
  filter_genes_pcspots = 0.3, 
  filter_mito = TRUE
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

fn <- here("outputs", "SPE", "spe_mouseHPC_preprocessed.rds")
saveRDS(spe, file = fn)

