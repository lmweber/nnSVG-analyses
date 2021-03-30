################################################
# R script to calculate SVGs across whole sample
# using BRISC method
################################################


# ---------------------------------------------
# Load data, preprocessing, calculate logcounts
# ---------------------------------------------

# Code from current version of OSTA

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


# ----------------
# Fit BRISC models
# ----------------

# using implementation in spatzli package
# stores outputs and runtimes in rowData of spe object

# runtime: ~1.5 hours using 6 cores

library(spatzli)

# remove NA labeled spots
spe <- spe[, !is.na(colData(spe)$ground_truth)]
dim(spe)

# create model matrix from manually annotated layer labels
sum(is.na(colData(spe)$ground_truth))
x <- model.matrix(~ colData(spe)$ground_truth)
dim(x)

runtime <- system.time({
  spe <- rankSVGsBRISC(spe, x = x, n_threads = 6)
})

runtime


# -------------------------
# Store HVGs for comparison
# -------------------------

library(dplyr)

# store HVGs data frame and match to correct rows
colnames(dec) <- paste0("hvgs_", colnames(dec))
dec$gene_id <- rownames(dec)

all(rownames(dec) %in% rowData(spe)$gene_id)
rowdata <- left_join(as.data.frame(rowData(spe)), as.data.frame(dec), by = "gene_id")
rowdata <- DataFrame(rowdata)
rownames(rowdata) <- rownames(rowData(spe))

identical(rownames(spe), rownames(rowdata))
stopifnot(nrow(rowdata) == nrow(spe))

rowData(spe) <- rowdata

# calculate HVGs ranks
rank_hvgs <- seq_along(top_hvgs)
names(rank_hvgs) <- top_hvgs

rank_hvgs_all <- rep(NA, nrow(spe))
names(rank_hvgs_all) <- rownames(spe)
rank_hvgs_all[names(rank_hvgs)] <- rank_hvgs

stopifnot(length(rank_hvgs_all) == nrow(spe))

rowData(spe)$rank_hvgs <- unname(rank_hvgs_all)


# -----------
# Save object
# -----------

saveRDS(spe, file = "../outputs/spe_svgs_brisc_withinlayers.rds")


