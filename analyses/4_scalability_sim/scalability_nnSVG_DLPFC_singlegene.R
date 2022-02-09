####################################
# Script for scalability simulations
# Lukas Weber, Feb 2022
####################################

# nnSVG, human DLPFC dataset, single gene

# interactive cluster session
# qrsh -l mem_free=20G,h_vmem=22G,h_fsize=100G -now n
# module load conda_R/4.1.x


library(SpatialExperiment)
library(STexampleData)
library(nnSVG)
library(scry)
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

n <- c(100, 200, 500, 1000, 2000, n_all)

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

# filter low-expressed genes
filter_genes <- 1
n_spots <- ceiling(filter_genes / 100 * ncol(spe))
ix_remove <- rowSums(counts(spe) > 0) < n_spots
table(ix_remove)
message("removing ", sum(ix_remove), " out of ", nrow(spe), " genes due to low expression")
spe <- spe[!ix_remove, ]
dim(spe)

# filter mitochondrial genes
is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
table(is_mito)
message("removing ", sum(is_mito), " mitochondrial genes")
spe <- spe[!is_mito, ]
dim(spe)

# set seed for reproducibility
set.seed(123)

# calculate deviance residuals using scry package
spe <- nullResiduals(spe, assay = "counts", fam = "binomial", type = "deviance")
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
    
    # skip filtering since already performed above
    runtime <- system.time({
      out <- nnSVG(spe_sub, filter_genes = FALSE, filter_mito = FALSE, n_threads = 1)
    })
    
    # 'elapsed' time is real human time
    runtimes_i[j] <- runtime[["elapsed"]]
  }
  
  runtimes[[i]] <- runtimes_i
}


# ------------
# save results
# ------------

file_runtimes <- here("outputs", "scalability", "runtimes_scalability_nnSVG_DLPFC_singlegene.rds")
saveRDS(runtimes, file = file_runtimes)

