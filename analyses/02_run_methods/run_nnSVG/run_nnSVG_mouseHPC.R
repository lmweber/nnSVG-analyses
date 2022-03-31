#######################
# Script to run method
# Lukas Weber, Mar 2022
#######################

# method: nnSVG with covariates
# dataset: Slide-seqV2 mouse HPC


library(SpatialExperiment)
library(nnSVG)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "preprocessed", "spe_mouseHPC_preprocessed.rds")
spe <- readRDS(fn)

dim(spe)


# ----------
# run method
# ----------

# run method and save results, runtime, peak memory usage

# create model matrix for cell type labels
X <- model.matrix(~ colData(spe)$celltype)
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run nnSVG with covariates
set.seed(123)
runtime <- system.time({
  spe <- nnSVG(
    spe, 
    X = X, 
    assay_name = "logcounts", 
    n_neighbors = 10, 
    order = "AMMD", 
    n_threads = 10, 
    verbose = FALSE
  )
})

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mouseHPC_nnSVG.rds")
saveRDS(spe, file = file)

