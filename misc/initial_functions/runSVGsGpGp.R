#' runSVGsGpGp
#' 
#' Run method to identify spatially variable genes (SVGs) using GpGp.
#' 
#' Run method to identify spatially variable genes (SVGs) using GpGp (Guinness
#' 2018).
#' 
#' This function runs GpGp separately for each gene, using parallelization for
#' faster runtime using one core per GpGp run. The main outputs of interest are
#' the covariance parameter estimates stored in 'covparms' in the GpGp output
#' (for the 'exponential_isotropic' covariance function these are: variance,
#' range, nugget). Optionally, the 'range' parameter (parameterized as 'phi' in
#' other methods) can be fixed to a user-specified value, so that this parameter
#' is the same for all genes. We also use the log-likelihood ('loglik' in GpGp
#' output) to calculate likelihood ratio tests.
#' 
#' Note the parameterization used by GpGp: total variance = sigmasq + (tausq *
#' sigmasq); i.e. the nugget is 'tausq * sigmasq', not 'tausq' itself as in many
#' other methods.
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
#'   NULL, which specifies an intercept-only model. See \code{GpGp}
#'   documentation for details.
#' 
#' @param fix_param_range \code{numeric} Whether to use a fixed parameter value
#'   for 'range' covariance function parameter in 'exponential_isotropic'
#'   covariance function (corresponding to 'phi' in other parameterizations).
#'   For example, set to 0.5 to fix the parameter to this value. Default = NULL,
#'   which will estimate the parameter. See \code{GpGp} documentation for
#'   details.
#' 
#' @param n_neighbors \code{numeric} Number of nearest neighbors. See
#'   \code{GpGp} documentation for details.
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
#'   \code{GpGp}. Default = FALSE.
#' 
#' 
#' @return Returns output values stored as new columns in \code{rowData} in the
#'   \code{spe} \code{SpatialExperiment} object.
#' 
#' 
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SingleCellExperiment logcounts
#' @importFrom SummarizedExperiment assayNames rowData 'rowData<-'
#' @importFrom GpGp order_maxmin find_ordered_nn fit_model
#' @importFrom BiocParallel bplapply MulticoreParam
#' @importFrom Matrix rowMeans
#' @importFrom matrixStats rowVars
#' @importFrom stats lm logLik pchisq p.adjust
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
#'   spe_1 <- runSVGsGpGp(spe_1, verbose = TRUE)
#' })
#' 
#' # subset 100 genes and use parallelization
#' # spe_100 <- spe[1:100, ]
#' # spe_100 <- runSVGsGpGp(spe_100, n_threads = 4)
#' 
runSVGsGpGp <- function(spe, x = NULL, fix_param_range = NULL, n_neighbors = 15, 
                        lr_test = TRUE, n_threads = 1, verbose = FALSE) {
  
  stopifnot("logcounts" %in% assayNames(spe))
  
  if (!is.null(x)) stopifnot(nrow(x) == ncol(spe))
  
  y <- logcounts(spe)
  
  # --------
  # run GpGp
  # --------
  
  coords <- spatialCoords(spe)
  # scale coordinates proportionally
  # note: not needed when calculating ordering manually below
  #range_all <- max(apply(coords, 2, function(col) diff(range(col))))
  #coords <- apply(coords, 2, function(col) (col - min(col)) / range_all)
  
  # prepare for parallelized loop
  
  # calculate ordering only once
  ord <- order_maxmin(coords)
  # re-parameterize ordering as ranking (if needed)
  #ord_rank <- match(seq_len(nrow(coords)), ord)
  # calculate nearest neighbors only once
  nn <- find_ordered_nn(coords[ord, ], m = n_neighbors)
  
  # run GpGp using parallelization
  ix <- seq_len(nrow(y))
  out_gpgp <- bplapply(ix, function(i) {
    # fit model (note: default if x is NULL is intercept-only model)
    y_i <- y[i, ]
    runtime <- system.time({
      # note: using manual reordering, pre-calculated nearest neighbors, optionally fixed parameter
      if (is.null(fix_param_range)) {
        start_parms <- NULL
        fixed_parms <- NULL
      } else {
        start_parms <- c(0.1, fix_param_range, 0.1)
        fixed_parms = 2
      }
      out_i <- fit_model(y = y_i[ord], locs = coords[ord, ], X = x[ord, ], 
                         covfun_name = "exponential_isotropic", 
                         start_parms = start_parms, fixed_parms = fixed_parms, 
                         NNarray = nn, reorder = FALSE, m_seq = n_neighbors, 
                         silent = !verbose)
    })
    res_i <- c(
      sigmasq = out_i$covparms[1], 
      tausq_sigmasq = out_i$covparms[3], 
      loglik = out_i$loglik, 
      runtime = runtime[["elapsed"]]
    )
    res_i
  }, BPPARAM = MulticoreParam(workers = n_threads))
  
  # collapse output list into matrix
  mat_gpgp <- do.call("rbind", out_gpgp)
  
  # add column containing standard parameterization of nugget, i.e. tausq
  mat_gpgp <- cbind(
    mat_gpgp, 
    tausq = mat_gpgp[, "tausq_sigmasq"] * mat_gpgp[, "sigmasq"]
  )
  
  # --------------------
  # calculate statistics
  # --------------------
  
  # mean logcounts
  mat_gpgp <- cbind(
    mat_gpgp, 
    mean = rowMeans(y)
  )
  
  # variance of logcounts
  mat_gpgp <- cbind(
    mat_gpgp, 
    var = rowVars(as.matrix(y))
  )
  
  # spatial coefficient of variation
  mat_gpgp <- cbind(
    mat_gpgp, 
    spcov = sqrt(mat_gpgp[, "sigmasq"]) / mat_gpgp[, "mean"]
  )
  
  # ratio of spatial to non-spatial variance
  mat_gpgp <- cbind(
    mat_gpgp, 
    ratio_sv = mat_gpgp[, "sigmasq"] / mat_gpgp[, "tausq"]
  )
  
  # proportion of spatial variance (out of total variance)
  mat_gpgp <- cbind(
    mat_gpgp, 
    prop_sv = mat_gpgp[, "sigmasq"] / (mat_gpgp[, "sigmasq"] + mat_gpgp[, "tausq"])
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
    
    mat_gpgp <- cbind(
      mat_gpgp, 
      loglik_lm = loglik_lm
    )
  }
  
  # calculate likelihood ratio test (Wilks' theorem, asymptotic chi-square with 
  # 2 degrees of freedom since 2 more parameters in full model)
  if (lr_test) {
    lr_stat = -2 * (mat_gpgp[, "loglik_lm"] - mat_gpgp[, "loglik"])
    pval <- 1 - pchisq(lr_stat, df = 2)
    padj <- p.adjust(pval, method = "BH")
    
    mat_gpgp <- cbind(
      mat_gpgp, 
      lr_stat = lr_stat, 
      pval = pval, 
      padj = padj
    )
  }
  
  # -------------------------------
  # return in rowData of spe object
  # -------------------------------
  
  stopifnot(nrow(spe) == nrow(mat_gpgp))
  
  rowData(spe) <- cbind(rowData(spe), mat_gpgp)
  
  spe
}

