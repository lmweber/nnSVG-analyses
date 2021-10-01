#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: nnSVG
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

# run nnSVG
runtime <- system.time({
  spe <- nnSVG(spe, x = NULL, 
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

file <- here("outputs", "results", "spe_mOB_nnSVG.rds")
saveRDS(spe, file = file)

