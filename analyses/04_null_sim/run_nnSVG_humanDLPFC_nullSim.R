#######################################
# Script to run method: null simulation
# Lukas Weber, Mar 2022
#######################################

# method: nnSVG
# dataset: Visium human DLPFC


library(SpatialExperiment)
library(nnSVG)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "null_sim", "spe_humanDLPFC_nullSim_preprocessed.rds")
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

file <- here("outputs", "null_sim", "spe_humanDLPFC_nnSVG_nullSim.rds")
saveRDS(spe, file = file)

