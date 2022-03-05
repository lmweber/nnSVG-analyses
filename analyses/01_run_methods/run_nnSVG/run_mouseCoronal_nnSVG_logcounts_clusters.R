#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: nnSVG (with covariates for clusters)
# dataset: 10x Genomics Visium mouse coronal brain section


library(nnSVG)
library(SpatialExperiment)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

file <- here("outputs", "SPE", "spe_mouseCoronal.rds")
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

# create model matrix of covariates for cell types using cluster labels
X <- model.matrix(~ colData(spe)$label)
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run nnSVG with covariates
runtime <- system.time({
  spe <- nnSVG(spe, x = X, 
               assay_name = "logcounts", 
               filter_genes = FALSE, filter_mito = FALSE, 
               n_threads = 10)
})

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "nnSVG", "spe_mouseCoronal_nnSVG_logcounts_clusters.rds")
saveRDS(spe, file = file)

