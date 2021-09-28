#' preprocessSVGs
#' 
#' Preprocessing steps to run SVGs methods
#' 
#' Convenience function to run several preprocessing steps to prepare data for
#' input to SVGs methods. Mainly intended for development work, but may also be
#' useful for others. The steps are intended for 10x Genomics Visium data. Code
#' is sourced from our online book \code{OSTA}, and makes use of the
#' \code{SpatialExperiment} structure.
#' 
#' 
#' @param spe \code{SpatialExperiment} Input data, assumed to be a
#'   \code{SpatialExperiment} object.
#' 
#' @param in_tissue \code{logical} Whether to keep only spots over tissue,
#'   identified by column \code{in_tissue} in \code{spatialData}. Default =
#'   TRUE.
#' 
#' @param filter_genes \code{integer} Filtering on UMI counts per gene. Genes
#'   with less than or equal to this number of total UMI counts across spots
#'   will be filtered out. Default = 10.
#' 
#' @param filter_mito \code{logical} Whether to filter out mitochondrial genes.
#'   Assumes \code{rowData} contains a column named \code{gene_name}. Default =
#'   TRUE.
#' 
#' @param seed \code{integer} Random seed for steps requiring random seed.
#'   Default = 123.
#' 
#' 
#' @return Returns a \code{SpatialExperiment} object that is ready to be
#'   provided as the input to SVGs methods.
#' 
#' 
#' @importFrom SpatialExperiment spatialData
#' @importFrom SingleCellExperiment counts
#' @importFrom SummarizedExperiment assayNames
#' @importFrom scran quickCluster computeSumFactors
#' @importFrom scuttle logNormCounts
#' @importFrom methods isClass
#' 
#' @export
#' 
#' @examples
#' paste0("to do")
#' 
preprocessSVGs <- function(spe, in_tissue = TRUE, 
                           filter_genes = 10, filter_mito = TRUE, 
                           seed = 123) {
  
  stopifnot(isClass(spe, "SpatialExperiment"))
  
  # spots over tissue
  
  if (in_tissue) {
    stopifnot("in_tissue" %in% colnames(spatialData(spe)))
    spe <- spe[, spatialData(spe)$in_tissue == 1]
  }
  
  # gene filtering
  
  if (filter_genes > 0 ) {
    stopifnot("counts" %in% assayNames(spe))
    sums <- rowSums(counts(spe))
    ix_remove <- sums <= filter_genes
    message("removing ", sum(ix_remove), " out of ", nrow(spe), " genes due to low counts")
    spe <- spe[!ix_remove, ]
  }
  
  if (filter_mito) {
    is_mito <- grepl("(^MT-)|(^mt-)", rowData(spe)$gene_name)
    message("removing ", sum(is_mito), " mitochondrial genes")
    spe <- spe[!is_mito, ]
  }
  
  # normalization and logcounts
  
  set.seed(seed)
  qclus <- quickCluster(spe)
  spe <- computeSumFactors(spe, cluster = qclus)
  spe <- logNormCounts(spe)
  
  # return object
  
  spe
}

