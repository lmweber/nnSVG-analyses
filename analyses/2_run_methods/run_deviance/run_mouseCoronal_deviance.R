#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: deviance (scry)
# dataset: 10x Genomics Visium mouse coronal brain section


library(scry)
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

# skip filtering since this was performed during preprocessing

# run deviance feature selection
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

file <- here("outputs", "results", "deviance", "spe_mouseCoronal_deviance.rds")
saveRDS(spe, file = file)

