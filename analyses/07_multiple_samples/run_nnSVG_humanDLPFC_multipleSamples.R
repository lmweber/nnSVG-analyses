#######################################
# Multiple samples: script to run nnSVG
# Lukas Weber, updated Mar 2023
#######################################

# dataset: Visium human DLPFC (multiple samples)
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)

# run nnSVG individually for each sample


# interactive session on compute cluster

# qrsh -pe local 10 -l mem_free=2G,h_vmem=3G,h_fsize=200G -now n
# cd /dcs04/hicks/data/lweber/nnSVG_analyses/nnSVG-analyses
# module load conda_R/4.2.x
# R


library(SpatialExperiment)
library(here)
library(nnSVG)
library(scran)
library(scater)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "multiple_samples", "spe_humanDLPFC_multipleSamples_preprocessed.rds")
spe <- readRDS(fn)

dim(spe)
table(colData(spe)$sample_id)

# 12 samples in this dataset
n_samples <- length(table(colData(spe)$sample_id))
n_samples


# --------------
# gene filtering
# --------------

# remove genes with zero expression
ix_zero_genes <- rowSums(counts(spe)) == 0
table(ix_zero_genes)

if (sum(ix_zero_genes) > 0) {
  spe <- spe[!ix_zero_genes, ]
}

dim(spe)

# remove spots with zero expression
ix_zero_spots <- colSums(counts(spe)) == 0
table(ix_zero_spots)

if (sum(ix_zero_spots) > 0) {
  spe <- spe[, !ix_zero_spots]
}

dim(spe)


# -------------------------
# run nnSVG for each sample
# -------------------------

# run nnSVG once per sample and store lists of top SVGs

# note: we perform additional gene filtering for nnSVG per sample, and 
# re-calculate logcounts after gene filtering for each sample

sample_ids <- unique(colData(spe)$sample_id)
sample_ids

res_list <- as.list(rep(NA, length(sample_ids)))
names(res_list) <- sample_ids

for (s in seq_along(sample_ids)) {
  
  # select sample
  ix <- colData(spe)$sample_id == sample_ids[s]
  spe_sub <- spe[, ix]
  
  dim(spe_sub)
  
  # run nnSVG filtering for mitochondrial genes and low-expressed genes
  spe_sub <- filter_genes(
    spe_sub, 
    filter_genes_ncounts = 3, 
    filter_genes_pcspots = 0.5, 
    filter_mito = TRUE
  )
  
  # remove any zeros introduced by filtering
  ix_zeros <- colSums(counts(spe_sub)) == 0
  if (sum(ix_zeros) > 0) {
    spe_sub <- spe_sub[, !ix_zeros]
  }
  
  dim(spe_sub)
  
  # re-calculate logcounts after filtering
  spe_sub <- computeLibraryFactors(spe_sub)
  spe_sub <- logNormCounts(spe_sub)
  
  # run nnSVG
  set.seed(123)
  spe_sub <- nnSVG(
    spe_sub, 
    n_threads = 10
  )
  
  # store results for this sample
  res_list[[s]] <- rowData(spe_sub)
}


# -----------
# save object
# -----------

file <- here("outputs", "multiple_samples", "res_humanDLPFC_nnSVG_multipleSamples.rds")
saveRDS(res_list, file = file)

