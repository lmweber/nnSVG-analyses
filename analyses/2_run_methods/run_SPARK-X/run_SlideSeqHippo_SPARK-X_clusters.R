#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: SPARK-X (with covariates for clusters)
# dataset: Slide-seqV2 mouse hippocampus


library(SPARK)
library(SpatialExperiment)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

file <- here("outputs", "SPE", "spe_SlideSeqHippo.rds")
spe <- readRDS(file)

spe


# ----------
# run method
# ----------

# run method and save the following:
# - results
# - runtime
# - peak memory usage

# skip filtering since this was already done during preprocessing

# create model matrix of covariates for cell types
# remove NAs from cell type labels
spe <- spe[, !is.na(colData(spe)$celltype)]
X <- model.matrix(~ as.factor(colData(spe)$celltype))
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run SPARK-X with covariates
runtime <- system.time({
  sparkx_out <- sparkx(count_in = counts(spe), 
                       locus_in = spatialCoords(spe), 
                       X_in = NULL, 
                       numCores = 4, option = "mixture", verbose = TRUE)
})

# results for individual kernels
head(sparkx_out$stats)
head(sparkx_out$res_stest)
# results for combined kernels
head(sparkx_out$res_mtest)

# store results in SpatialExperiment object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_id))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "SPARK-X", "spe_SlideSeqHippo_SPARK-X_clusters.rds")
saveRDS(spe, file = file)

