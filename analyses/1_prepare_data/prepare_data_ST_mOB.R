#################################
# Script to load and prepare data
# Lukas Weber, Oct 2021
#################################

# dataset: Spatial Transcriptomics (ST) mouse olfactory bulb (mOB)

# references:
# Stahl et al. (2016): https://www.science.org/doi/full/10.1126/science.aaf2403
# Miller et al. (2021): https://genome.cshlp.org/content/early/2021/05/25/gr.271288.120


library(nnSVG)
library(SpatialExperiment)
library(MERINGUE)
library(dplyr)
library(here)


# ---------
# load data
# ---------

# load data object from MERINGUE package

data(mOB)

names(mOB)
dim(mOB$counts)
mOB$counts[1000:1006, 1:6]
head(mOB$pos)
head(mOB$results)
str(mOB$annot)


# ------------------------
# create SpatialExperiment
# ------------------------

row_data <- data.frame(
  gene_name = rownames(mOB$counts)
)
res <- cbind(mOB$results, gene_name = rownames(mOB$results))
row_data <- DataFrame(left_join(row_data, res, by = "gene_name"))
head(row_data)

col_data <- data.frame(
  barcode_id = colnames(mOB$counts)
)
annot <- data.frame(
  barcode_id = names(mOB$annot), 
  layer_id = mOB$annot
)
col_data <- DataFrame(left_join(col_data, annot, by = "barcode_id"))
head(col_data)

spatial_coords <- as.matrix(mOB$pos[col_data$barcode_id, ])
head(spatial_coords)

stopifnot(nrow(spatial_coords) == nrow(col_data))
stopifnot(all(rownames(spatial_coords) == col_data$barcode_id))
stopifnot(all(colnames(mOB$counts) == col_data$barcode_id))

spe <- SpatialExperiment(
  assays = list(counts = mOB$counts), 
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords
)


# -------------------
# preprocessing steps
# -------------------

# using 'preprocessSVG()' function from nnSVG package

# preprocessing steps:
# - filter low-expressed genes
# - filter mitochondrial genes
# - normalization
# - calculate logcounts

# set seed for reproducibility
set.seed(123)
spe <- preprocessSVG(spe, in_tissue = FALSE, 
                     filter_genes = 20, filter_mito = TRUE)

dim(spe)
assayNames(spe)
logcounts(spe)[1:6, 1:6]


# ----------
# clustering
# ----------

# clustering to identify cell types

# using 'clusterSVG()' function from nnSVG package

# set seed for reproducibility
set.seed(123)
spe <- clusterSVG(spe, filter_mito = FALSE)

colData(spe)

# compare against ground truth labels provided with this dataset
table(truth = colData(spe)$layer_id, 
      clusters = colData(spe)$label)


# -----------
# save object
# -----------

file <- here("outputs", "SPE", "spe_mOB.rds")
saveRDS(spe, file = file)

