#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: HVGs
# dataset: Spatial Transcriptomics (ST) mouse olfactory bulb (mOB)


library(scran)
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

# run HVGs
runtime <- system.time({
  dec <- modelGeneVar(spe)
})

# store in object
stopifnot(all(rownames(dec) == rowData(spe)$gene_id))
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks according to 'bio' statistic
rowData(spe)$rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mOB_HVGs.rds")
saveRDS(spe, file = file)

