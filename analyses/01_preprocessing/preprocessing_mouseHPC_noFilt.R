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

# remove spots with NA cell type labels
spe <- spe[, !is.na(colData(spe)$celltype)]
dim(spe)

table(colData(spe)$celltype)


# filter mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
rowData(spe)$gene_name[is_mito]

spe <- spe[!is_mito, ]

dim(spe)


# filter zero-expressed genes
is_zero <- rowSums(counts(spe)) == 0
table(is_zero)

spe <- spe[!is_zero, ]

dim(spe)


# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

assayNames(spe)


# -----------
# save object
# -----------

fn <- here("outputs", "preprocessed", "spe_mouseHPC_preprocessed_noFilt.rds")
saveRDS(spe, file = fn)

