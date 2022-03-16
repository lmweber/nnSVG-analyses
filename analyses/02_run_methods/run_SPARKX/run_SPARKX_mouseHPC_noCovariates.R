#######################
# Script to run method
# Lukas Weber, Mar 2022
#######################

# method: SPARK-X without covariates
# dataset: Slide-seqV2 mouse HPC


library(SpatialExperiment)
library(SPARK)
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

# run SPARK-X without covariates
set.seed(123)
runtime <- system.time({
  sparkx_out <- sparkx(
    count_in = counts(spe), 
    locus_in = spatialCoords(spe), 
    X_in = NULL, 
    numCores = 10, 
    option = "mixture", 
    verbose = FALSE
  )
})

# results for individual kernels
head(sparkx_out$stats)
head(sparkx_out$res_stest)
# results for combined kernels
head(sparkx_out$res_mtest)

# store results in SPE object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_name))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mouseHPC_SPARKX_noCovariates.rds")
saveRDS(spe, file = file)

