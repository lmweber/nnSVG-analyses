###############################
# Script to run method
# Lukas Weber, updated Mar 2022
###############################

# method: nnSVG
# dataset: Visium human DLPFC, without filtering low-expressed genes


library(SpatialExperiment)
library(nnSVG)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "preprocessed", "spe_humanDLPFC_preprocessed_noFilt.rds")
spe <- readRDS(fn)

dim(spe)


# ----------
# run method
# ----------

# run method and save results, runtime, peak memory usage

# run nnSVG
set.seed(123)
runtime <- system.time({
  spe <- nnSVG(
    spe, 
    X = NULL, 
    assay_name = "logcounts", 
    n_neighbors = 15, 
    order = "Sum_coords", 
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

file <- here("outputs", "results", "spe_humanDLPFC_nnSVG_noFilt.rds")
saveRDS(spe, file = file)

