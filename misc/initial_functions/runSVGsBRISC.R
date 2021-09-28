#' runSVGsBRISC
#' 
#' Run method to identify spatially variable genes (SVGs) using BRISC.
#' 
#' Run method to identify spatially variable genes (SVGs) using BRISC
#' ("bootstrap for rapid inference on spatial covariances") methodology (Saha
#' and Datta 2018).
#' 
#' This function runs BRISC separately for each gene, using parallelization for
#' faster runtime using one core per BRISC run. The main outputs of interest are
#' the covariance parameter estimates stored in 'Theta' in the BRISC output
#' (sigma.sq, tau.sq, phi). We use these estimates to perform inference on the
#' 'sigma.sq' parameter, and to calculate an effect size estimate defined as the
#' proportion of spatial variance out of total variance, 'prop_sv' = 'sigma.sq /
#' (sigma.sq + tau.sq)'.
#' 
#' Significant SVGs can then be identified as those with a highly significant
#' p-value from 'sigma.sq' and large effect size 'prop_sv'.
#' 
#' Assumes the input object is a \code{SpatialExperiment} containing an assay
#' named \code{logcounts}, which has been filtered to exclude very low-expressed
#' genes, e.g. as prepared with \code{\link{preprocessSVGs}}.
#' 
#' 
#' @param spe \code{SpatialExperiment} Input object, assumed to be a
#'   \code{SpatialExperiment} containing an assay named \code{logcounts} and
#'   spatial coordinates accessible with \code{spatialCoords()}.
#' 
#' @param x \code{numeric matrix} Matrix of covariates, with number of rows
#'   (spots) matching the number of columns (spots) in \code{spe}. Default =
#'   NULL, which specifies an intercept-only model. See \code{BRISC}
#'   documentation for details.
#' 
#' @param lr_test \code{logical} Whether to return log likelihoods and calculate
#'   likelihood ratio tests compared to the null model without spatial terms. If
#'   TRUE, will calculate log likelihoods for the null models, calculate
#'   likelihood ratio tests using the asymptotic chi-squared distribution with 2
#'   degrees of freedom, and calculate adjusted p-values using the
#'   Benjamini-Hochberg method. Default = TRUE.
#' 
#' @param n_threads \code{integer} Number of threads for parallelization.
#'   Default = 1.
#' 
#' @param verbose \code{logical} Whether to display verbose output from
#'   \code{BRISC}. Default = FALSE.
#' 
#' 
#' @return Returns output values stored as new columns in \code{rowData} in the
#'   \code{spe} \code{SpatialExperiment} object.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames rowData 'rowData<-'
#' @importFrom BRISC BRISC_order BRISC_neighbor BRISC_estimation
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowMeans
#' 
#' @export
#' 
#' @examples
#' library(SpatialExperiment)
#' library(STexampleData)
#' library(spatzli)
#' 
#' spe <- Visium_humanDLPFC()
#' 
#' spe <- preprocessSVGs(spe)
#' 
#' # subset 1 gene
#' spe_1 <- spe[1, ]
#' system.time({
#'   spe_1 <- runSVGsBRISC(spe_1, verbose = TRUE)
#' })
#' 
#' # subset 100 genes and use parallelization
#' # spe_100 <- spe[1:100, ]
#' # spe_100 <- runSVGsBRISC(spe_100, n_threads = 4)
#' 
runSVGsBRISC <- function(spe, x = NULL, lr_test = TRUE, 
                         n_threads = 1, verbose = FALSE) {
  
  stopifnot("logcounts" %in% assayNames(spe))
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  y <- logcounts(spe)
  
  # ---------
  # run BRISC
  # ---------
  
  # scale coordinates proportionally
  coords <- spatialCoords(spe)
  range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # calculate ordering
  order_brisc <- BRISC_order(coords, order = "AMMD", verbose = verbose)
  
  # calculate nearest neighbors
  nn_brisc <- BRISC_neighbor(coords, n.neighbors = 15, n_omp = 1, 
                             search.type = "cb", ordering = order_brisc, 
                             verbose = verbose)
  
  # run BRISC using parallelization
  ix <- seq_len(nrow(y))
  out_brisc <- bplapply(ix, function(i) {
    # fit model (note: default if x is NULL is intercept-only model)
    y_i <- y[i, ]
    runtime <- system.time({
      out_i <- BRISC_estimation(coords = coords, y = y_i, x = x, 
                                cov.model = "exponential", 
                                ordering = order_brisc, neighbor = nn_brisc, 
                                verbose = verbose)
    })
    res_i <- c(
      out_i$Theta, 
      loglik = out_i$log_likelihood, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse output list into matrix
  mat_brisc <- do.call("rbind", out_brisc)
  
  # --------------------
  # calculate statistics
  # --------------------
  
  # mean logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    mean = rowMeans(y)
  )
  
  # variance of logcounts
  mat_brisc <- cbind(
    mat_brisc, 
    var = rowVars(as.matrix(y))
  )
  
  # spatial coefficient of variation
  mat_brisc <- cbind(
    mat_brisc, 
    spcov = sqrt(mat_brisc[, "sigma.sq"]) / mat_brisc[, "mean"]
  )
  
  # ratio of spatial to non-spatial variance
  mat_brisc <- cbind(
    mat_brisc, 
    ratio_sv = mat_brisc[, "sigma.sq"] / mat_brisc[, "tau.sq"]
  )
  
  # proportion of spatial variance (out of total variance)
  mat_brisc <- cbind(
    mat_brisc, 
    prop_sv = mat_brisc[, "sigma.sq"] / (mat_brisc[, "sigma.sq"] + mat_brisc[, "tau.sq"])
  )
  
  # ----------------------
  # likelihood ratio tests
  # ----------------------
  
  # calculate log likelihoods for non-spatial models
  if (lr_test) {
    loglik_lm <- sapply(seq_len(nrow(spe)), function(i) {
      y_i <- y[i, ]
      if (is.null(x)) {
        x <- rep(1, ncol(spe))
      }
      as.numeric(logLik(lm(y_i ~ x)))
    })
    
    mat_brisc <- cbind(
      mat_brisc, 
      loglik_lm = loglik_lm
    )
  }
  
  # calculate likelihood ratio test (Wilks' theorem, asymptotic chi-square with 
  # 2 degrees of freedom since 2 more parameters in full model)
  if (lr_test) {
    lr_stat = -2 * (mat_brisc[, "loglik_lm"] - mat_brisc[, "loglik"])
    pval <- 1 - pchisq(lr_stat, df = 2)
    padj <- p.adjust(pval, method = "BH")
    
    mat_brisc <- cbind(
      mat_brisc, 
      lr_stat = lr_stat, 
      pval = pval, 
      padj = padj
    )
  }
  
  # -------------------------------
  # return in rowData of spe object
  # -------------------------------
  
  stopifnot(nrow(spe) == nrow(mat_brisc))
  
  rowData(spe) <- cbind(rowData(spe), mat_brisc)
  
  spe
}

