#######################
# Script to run method
# Lukas Weber, Mar 2022
#######################

# method: Moran's I
# dataset: Visium human DLPFC, without filtering low-expressed genes


library(SpatialExperiment)
library(Rfast2)
library(BiocParallel)
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

# calculate weights matrix
# note: not scaling weights matrix since this can cause numerical precision issues due to small values
d <- dist(spatialCoords(spe))
d_mat <- as.matrix(d)
w <- 1 / (d_mat ^ 2)
diag(w) <- 0

# run Moran's I
set.seed(123)
runtime <- system.time({
  res <- bplapply(seq_len(nrow(spe)), function(i) {
    y_i <- logcounts(spe)[i, ]
    out_i <- moranI(
      x = y_i, 
      w = w, 
      scaled = FALSE, 
      R = 0
    )
    unname(out_i[1])
  }, BPPARAM = MulticoreParam(workers = 10))
})

res <- unlist(res)

# store in object
stopifnot(length(res) == nrow(spe))
rowData(spe)$MoransI <- res

# calculate ranks
rowData(spe)$rank <- rank(-1 * rowData(spe)$MoransI, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_humanDLPFC_MoransI_noFilt.rds")
saveRDS(spe, file = file)

