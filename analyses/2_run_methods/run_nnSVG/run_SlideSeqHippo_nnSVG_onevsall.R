#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: nnSVG (with covariates for clusters: one vs. all)
# dataset: Slide-seqV2 mouse hippocampus


library(nnSVG)
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

# CA3 vs. all other cell types
X <- model.matrix(~ as.factor(as.numeric(colData(spe)$celltype == "CA3")))
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run nnSVG with covariates
runtime <- system.time({
  spe <- nnSVG(spe, x = X, 
               assay_name = "binomial_deviance_residuals", 
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

file <- here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG_onevsall.rds")
saveRDS(spe, file = file)

