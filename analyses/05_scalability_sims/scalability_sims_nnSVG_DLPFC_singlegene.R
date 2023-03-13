####################################
# Script for scalability simulations
# Lukas Weber, Jan 2023
####################################

# method: nnSVG
# data set: human DLPFC
# filtering: with filtering of low-expressed genes (using nnSVG default filtering)

# single gene


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scran)
library(here)


# ---------
# load data
# ---------

spe <- Visium_humanDLPFC()

# keep only spots over tissue
spe <- spe[, colData(spe)$in_tissue == 1]
dim(spe)


# -----------
# subsampling
# -----------

n_all <- ncol(spe)
n_all

n <- c(200, 500, 1000, 2000, n_all)

ix <- as.list(rep(NA, length(n)))

for (i in seq_along(n)) {
  # seed for reproducibility
  set.seed(123)
  ix[[i]] <- sample(seq_len(n_all), n[i], replace = FALSE)
}


# -------------
# preprocessing
# -------------

# starting from spots over tissue

# filter low-expressed and mitochondrial genes
# using gene filtering function from nnSVG package
spe <- filter_genes(
  spe, 
  filter_genes_ncounts = 3, 
  filter_genes_pcspots = 0.5, 
  filter_mito = TRUE
)
dim(spe)

# calculate log-transformed normalized counts using scran package
# using library size normalization
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)
assayNames(spe)


# ------------------
# select single gene
# ------------------

ix_gene <- which(rowData(spe)$gene_name == "PCP4")
ix_gene

spe <- spe[ix_gene, ]
dim(spe)


# -----------------
# run nnSVG in loop
# -----------------

# run nnSVG in loop for each subsampled set of spots

runtimes <- as.list(rep(NA, length(n)))
names(runtimes) <- n

# repeat runs due to stochasticity
n_iters <- 10

for (i in seq_along(n)) {
  spe_sub <- spe[, ix[[i]]]
  runtimes_i <- rep(NA, n_iters)
  
  # seed for reproducibility
  # (note: set seed outside loop to allow variability across runs)
  set.seed(123)
  
  for (j in seq_len(n_iters)) {
    print(paste0("loop iteration i = ", i, ", n[i] = ", n[i], ", j = ", j))
    
    # run nnSVG (with one thread)
    runtime <- system.time({
      out <- nnSVG(
        spe_sub, 
        X = NULL, 
        assay_name = "logcounts", 
        n_neighbors = 10, 
        order = "AMMD", 
        n_threads = 1, 
        verbose = FALSE
      )
    })
    
    # 'elapsed' time is real human time
    runtimes_i[j] <- runtime[["elapsed"]]
  }
  
  runtimes[[i]] <- runtimes_i
}


# ------------
# save results
# ------------

file_runtimes <- here("outputs", "scalability_sims", "runtimes_scalability_nnSVG_DLPFC_singlegene.rds")
saveRDS(runtimes, file = file_runtimes)

