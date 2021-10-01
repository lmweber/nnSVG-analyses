#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: nnSVG (with covariates for clusters)
# dataset: Spatial Transcriptomics (ST) mouse olfactory bulb (mOB)


library(nnSVG)
library(SpatialExperiment)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

file <- here("outputs", "SPE", "spe_mOB.rds")
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

# create model matrix of covariates for cell types using cluster labels
X <- model.matrix(~ colData(spe)$label)
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run nnSVG with covariates
runtime <- system.time({
  spe <- nnSVG(spe, x = X, 
               filter_genes = FALSE, filter_mito = FALSE, 
               n_threads = 4)
})

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mOB_nnSVG_clusters.rds")
saveRDS(spe, file = file)

