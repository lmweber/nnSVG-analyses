#######################
# Script to run method
# Lukas Weber, Oct 2021
#######################

# method: nnSVG
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

# run nnSVG
runtime <- system.time({
  spe <- nnSVG(spe, x = NULL, 
               assay_name = "binomial_deviance_residuals", 
               filter_genes = FALSE, filter_mito = FALSE, 
               n_threads = 10, on_disk = TRUE)
})

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "nnSVG", "spe_SlideSeqHippo_nnSVG.rds")
saveRDS(spe, file = file)

