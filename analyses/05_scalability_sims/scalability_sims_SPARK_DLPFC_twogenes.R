####################################
# Script for scalability simulations
# Lukas Weber, Mar 2023
####################################

# method: SPARK
# data set: human DLPFC

# two genes (note SPARK / SPARK-X model fitting requires at least 2 genes)


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scater)
library(scran)
library(SPARK)
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


# ------------
# select genes
# ------------

ix_gene <- which(rowData(spe)$gene_name %in% c("MOBP", "PCP4"))
ix_gene

spe <- spe[ix_gene, ]
dim(spe)


# -----------------
# run SPARK in loop
# -----------------

# run SPARK in loop for each subsampled set of spots

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
    
    runtime <- system.time({
      
      # create SPARK object and run SPARK (with one thread)
      # using code from SPARK vignette
      
      # create SPARK object
      spark <- CreateSPARKObject(counts = as.matrix(counts(spe_sub)), 
                                 location = as.data.frame(spatialCoords(spe_sub)), 
                                 percentage = 0.1, 
                                 min_total_counts = 10)
      spark@lib_size <- apply(spark@counts, 2, sum)
      
      # estimate parameters under null
      spark <- spark.vc(spark, 
                        covariates = NULL, 
                        lib_size = spark@lib_size, 
                        num_core = 1, 
                        verbose = FALSE)
      
      # calculate p-values
      spark <- spark.test(spark, 
                          check_positive = TRUE, 
                          verbose = FALSE)
    })
    
    runtimes[[i]] <- runtimes_i
  }
}


# ------------
# save results
# ------------

file_runtimes <- here("outputs", "scalability_sims", "runtimes_scalability_SPARK_DLPFC_twogenes.rds")
saveRDS(runtimes, file = file_runtimes)

