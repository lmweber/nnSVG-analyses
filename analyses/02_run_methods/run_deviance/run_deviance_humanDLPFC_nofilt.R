#######################
# Script to run method
# Lukas Weber, Mar 2022
#######################

# method: deviance
# dataset: human DLPFC


library(SpatialExperiment)
library(scry)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "SPE", "spe_humanDLPFC_preprocessed_noFilt.rds")
spe <- readRDS(fn)

dim(spe)


# ----------
# run method
# ----------

# run method and save results, runtime, peak memory usage

# run deviance feature selection
set.seed(123)
runtime <- system.time({
  spe <- devianceFeatureSelection(spe, assay = "counts", fam = "binomial")
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

file <- here("outputs", "results", "deviance", "spe_humanDLPFC_deviance_noFilt.rds")
saveRDS(spe, file = file)

