###############################
# Script to run method
# Lukas Weber, updated Mar 2022
###############################

# method: HVGs
# dataset: Slide-seqV2 mouse HPC, without filtering low-expressed genes


library(SpatialExperiment)
library(scran)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("outputs", "preprocessed", "spe_mouseHPC_preprocessed_noFilt.rds")
spe <- readRDS(fn)

dim(spe)


# ----------
# run method
# ----------

# run method and save results, runtime, peak memory usage

# run HVGs
set.seed(123)
runtime <- system.time({
  dec <- modelGeneVar(spe)
})

# store in object
stopifnot(all(rownames(dec) == rowData(spe)$gene_name))
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks
rowData(spe)$rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mouseHPC_HVGs_noFilt.rds")
saveRDS(spe, file = file)

