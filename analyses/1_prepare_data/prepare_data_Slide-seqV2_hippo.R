#################################
# Script to load and prepare data
# Lukas Weber, Oct 2021
#################################

# dataset: Slide-seqV2 mouse hippocampus

# references:
# Stickels et al. (2020): https://www.nature.com/articles/s41587-020-0739-1
# Cable et al. (2021): https://www.nature.com/articles/s41587-021-00830-w


# to run in interactive session on cluster:
# module load conda_R/4.1.x
# Rscript filename.R


library(SpatialExperiment)
library(nnSVG)
library(dplyr)
library(readr)
library(here)


# ---------
# load data
# ---------

# load data object from Stickels et al. (2020)

# load expression matrix
# runtime: several minutes
exprs <- read_table(here("..", "data", "Slide_seqV2_hippo", "raw_data", 
                         "Puck_200115_08.digital_expression.txt.gz"))
dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]

# get gene IDs from first column
gene_ids <- exprs$GENE
str(gene_ids)
length(gene_ids)

# get barcode IDs from column names (excluding gene IDs column)
bead_ids <- colnames(exprs)[-1]
str(bead_ids)
length(bead_ids)

# convert expression matrix to numeric matrix without gene IDs
exprs <- exprs[, -1]
exprs <- as.matrix(exprs)
stopifnot(nrow(exprs) == length(gene_ids))
rownames(exprs) <- gene_ids

dim(exprs)
format(object.size(exprs), units = "GB")
exprs[1:6, 1:6]


# load beads info
# note: contains mix of tab-delimited and comma-delimited
file <- file.path(here("..", "data", "Slide_seqV2_hippo", "raw_data"), 
                  "Puck_200115_08_bead_locations.csv")
bead_locations_colnames <- unlist(strsplit(readLines(file, n = 1), "\t"))
bead_locations <- read_csv(file, skip = 1, col_names = FALSE)
colnames(bead_locations) <- bead_locations_colnames
dim(bead_locations)
head(bead_locations)

stopifnot(nrow(bead_locations) == ncol(exprs))
stopifnot(all(bead_ids == bead_locations$barcodes))


# create SpatialExperiment
row_data <- DataFrame(
  gene_name = gene_ids
)

col_data <- DataFrame(
  barcode_id = bead_ids
)

spatial_data <- DataFrame(
  barcode_id = bead_ids
)

spatial_coords <- as.matrix(bead_locations[, c("xcoord", "ycoord")])
rownames(spatial_coords) <- bead_ids

exprs_sparse <- as(exprs, "dgCMatrix")

spe <- SpatialExperiment(
  rowData = row_data, 
  colData = col_data, 
  spatialCoords = spatial_coords, 
  spatialData = spatial_data, 
  assays = list(counts = exprs_sparse)
)

spe

format(object.size(spe), units = "MB")


# ---------------------
# load cell type labels
# ---------------------

# load predicted cell type labels from Cable et al. (2021), Figure 5A
# from data object provided by email by authors of Cable et al. (2021)

# load file containing RCTD outputs
out_rctd <- load(here("..", "data", "Slide_seqV2_hippo", "processed_data", 
                      "SlideseqHippo.RData"))

# extract top predicted cell type labels
labels_rctd <- results_df$first_type

# extract barcode IDs
barcodes_rctd <- puck@counts@Dimnames[[2]]

names(labels_rctd) <- barcodes_rctd

head(labels_rctd)
str(labels_rctd)

# store labels in SpatialExperiment object
colData(spe)$celltype <- NA
colData(spe)[names(labels_rctd), "celltype"] <- as.character(labels_rctd)


# -------------------
# preprocessing steps
# -------------------

# using 'preprocessSVG()' function from nnSVG package

# preprocessing steps:
# - filter low-expressed genes
# - filter mitochondrial genes
# - calculate deviance residuals and/or logcounts

spe <- preprocessSVG(spe, in_tissue = FALSE, 
                     filter_genes = 20, filter_mito = TRUE)

dim(spe)
assayNames(spe)
assays(spe)[["binomial_deviance_residuals"]][1:6, 1:6]
logcounts(spe)[1:6, 1:6]


# -----------
# save object
# -----------

file <- here("outputs", "SPE", "Slide_seqV2_hippo_SPE.rds")
saveRDS(spe, file = file)

