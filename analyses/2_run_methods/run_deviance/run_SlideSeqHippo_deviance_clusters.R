#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: deviance (scry) (with covariates for clusters)
# dataset: Slide-seqV2 mouse hippocampus


library(scry)
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

# skip filtering since this was performed during preprocessing

# factor of cell type labels
# remove NAs from cell type labels
spe <- spe[, !is.na(colData(spe)$celltype)]
X <- as.factor(colData(spe)$celltype)
stopifnot(length(X) == ncol(spe))

# run deviance feature selection
runtime <- system.time({
  spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial", 
                                  batch = X)
})

# calculate ranks
rowData(spe)$rank <- rank(-1 * rowData(spe)$binomial_deviance, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "deviance", "spe_SlideSeqHippo_deviance_clusters.rds")
saveRDS(spe, file = file)

