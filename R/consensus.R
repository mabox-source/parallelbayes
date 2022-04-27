#' The consensus Monte Carlo algorithm of Scott et al 2016
#'
#' Compute the weights for samples from the partial posterior distributions 
#' according to the algorithm of Scott et al 2016, resulting in a pooled, 
#' weighted sample that can be used to approximate expected values under the 
#' full data posterior distribution.
#'
#' When \code{theta} is a list, \code{consensus.weights} replicates on a single 
#' machine what would be performed on a cluster. Each element of \code{theta} 
#' would correspond to the samples from a single partial posterior distribution.
#'
#' This is implemented in a distributed manner using Spark when \code{theta} is 
#' a Spark table.
#'
#' Parameter samples \code{theta} should be sampled from the partial posterior 
#' distributions using "fractionated" priors as in Scott et al 2016. This is 
#' assumed to be the case.
#'
#' Output fields \code{cov_used.partial} and \code{cov_used.pooled} indicate 
#' whether we used the full covariance matrices of the partial or pooled 
#' samples, respectively, for weighting. In some cases the inverse covariance 
#' matrix is difficult to compute, in which case the covariances are ignored, 
#' as suggested by Scott et al 2016.
#'
#' @section References
#' \itemize{
#' \item{Scott, Steven L., Blocker, A.W., Bonassi, F.V., Chipman, H.A., George, E.I. and McCulloch, R.E., 2016. Bayes and big data: The consensus Monte Carlo algorithm. \emph{International Journal of Management Science and Engineering Management}, 11(2), pp.78-88.}
#' }
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param type an integer, either 1 or 2, specifying the weighting type to use. 
#' 1: use constant weighting. 2: weight samples using the sample covariance 
#' matrices of the partial posteriors.
#' @param return.pooled logical. If \code{TRUE}, pooled, weighted samples are 
#' returned. If Spark is used, this means returned to local memory (the calling 
#' environment) - in either case, weighted samples will be returned in a Spark 
#' table.
#' @param par.clust an optional cluster connection object from package 
#' \code{parallel}. Ignored if \code{theta} is a Spark table.
#' @param ncores an optional integer specifying the number of CPU cores to use 
#' (see \code{\link[parallel]{makeCluster}}). The default, 1, signifies that 
#' \code{parallel} will not be used.
#'
#' @return A list with elements: \code{w}, list of weighting matrices used; 
#' diagnostics \code{cov_used.partial} and \code{cov_used.pooled} (see 
#' details); a matrix \code{theta.pooled} of pooled, weighted samples if 
#' \code{return.pooled} is \code{TRUE}. If Spark is used, also field 
#' \code{theta.w}, a Spark table similar to \code{theta} of weighted samples.
#' 
#' @export
consensus.weights <- function(
  theta,
  type = 2,
  par.clust = NULL,
  ncores = 1
) {

  if (class(theta)[1] == "tbl_spark") {
    if (!require(sparklyr)) stop("sparklyr is required!")
    use_spark <- TRUE
  } else {
    if (!("list" %in% class(theta))) stop("theta must be a Spark table or a list!")
    use_spark <- FALSE
  }
  if (!use_spark && !is.null(par.clust) && class(par.clust)[1] == "SOCKcluster" && require(parallel)) {
    use_parallel <- TRUE
    new_cluster <- FALSE
  } else if (!use_spark && is.null(par.clust) && ncores > 1 && require(parallel)) {
    use_parallel <- TRUE
    new_cluster <- TRUE
    n_cores_available <- parallel::detectCores()
    if (ncores > n_cores_available) {
      ncores <- n_cores_available
      message(paste0("Using the maximum number of CPU cores available (", n_cores_available, ")"))
    }
    par.clust <- parallel::makeCluster(ncores)
  } else {
    if (!use_spark && ncores > 1) message("Package parallel not found, using ncores = 1.")
    ncores <- 1
    use_parallel <- FALSE
  }
}

