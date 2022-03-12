#######################
# Script to run method
# Lukas Weber, Mar 2022
#######################

# method: SPARK-X with covariates
# dataset: Slide-seqV2 mouse HPC, without filtering low-expressed genes


library(SpatialExperiment)
library(SPARK)
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

# remove spots with NA cell type labels
spe <- spe[, !is.na(colData(spe)$celltype)]

# create model matrix for cell type CA3 vs. all other cell types
X <- model.matrix(~ as.factor(as.numeric(colData(spe)$celltype == "CA3")))
dim(X)
head(X)
stopifnot(nrow(X) == ncol(spe))

# run SPARK-X with covariates
set.seed(123)
runtime <- system.time({
  sparkx_out <- sparkx(
    count_in = counts(spe), 
    locus_in = spatialCoords(spe), 
    X_in = X, 
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
# note: some additional genes have been filtered out by SPARK-X
table(rowData(spe)$gene_name %in% rownames(sparkx_out$res_mtest))
# match rows and store NAs for missing genes
colnames(sparkx_out$res_mtest)
rowData(spe)$combinedPval <- NA
rowData(spe)$adjustedPval <- NA
rowData(spe)[rownames(sparkx_out$res_mtest), "combinedPval"] <- sparkx_out$res_mtest$combinedPval
rowData(spe)[rownames(sparkx_out$res_mtest), "adjustedPval"] <- sparkx_out$res_mtest$adjustedPval

# calculate ranks
rowData(spe)$rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("outputs", "results", "spe_mouseHPC_SPARKX_noFilt.rds")
saveRDS(spe, file = file)

