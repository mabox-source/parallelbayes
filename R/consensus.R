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
#' components of parameter vectors). They must also have the same number of 
#' rows (same number of samples from each partial posterior). Each list element 
#' corresponds to a single partial posterior.
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
#' @param cov.tol a number, the tolerance level for the reciprocal condition 
#' number of the covariance matrix. If less than this, the covariances are 
#' ignored. See details.
#'
#' @return A list with elements: \code{w}, list of weighting matrices used; 
#' diagnostics \code{cov_used.partial} and \code{cov_used.pooled} (see 
#' details); a matrix \code{theta.w.pooled} of pooled, weighted samples if 
#' \code{return.pooled} is \code{TRUE}. If Spark is used, also field 
#' \code{theta.w}, a Spark table similar to \code{theta} of weighted samples.
#' 
#' @export
consensus.weights <- function(
  theta,
  type = 2,
  return.pooled = FALSE,
  par.clust = NULL,
  ncores = 1,
  cov.tol = .Machine$double.eps
) {

  if (!(type %in% 1:2)) stop("type must be 1 or 2!")
  if (class(theta)[1] == "tbl_spark") {
    if (!require(sparklyr)) stop("sparklyr is required!")
    Hvec <- theta %>% select(1) %>% count()
    use_spark <- TRUE
  } else {
    if (!("list" %in% class(theta))) stop("theta must be a Spark table or a list!")
    Hvec <- sapply(theta, nrow)
    use_spark <- FALSE
  }
  if (any(Hvec != Hvec[1])) stop("Each set of samples must have the same number of realisations!")
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

  # First derive weights.

  ctx <- list(
    use_spark = use_spark,
    type = type,
    tol = cov.tol
  )
  if (use_spark) {
    spark.res <- spark_apply(theta, consensus.worker, context = ctx)
    local.res <- as.data.frame(spark.res)
    # Check for error flag.
    if (any(is.na(local.res[,2]))) stop("Insufficient useable samples!") 
    w <- split(local.res[,2], local.res[,1])
  } else if (use_parallel) {
    w <- parLapply(par.clust, theta, consensus.worker, context = ctx)
    if (any(sapply(w, function(wi) {is.na(wi[,1])}))) stop("Insufficient useable samples!")
  } else {
    w <- lapply(theta, consensus.worker, ctx)
    if (any(sapply(w, function(wi) {is.na(wi[,1])}))) stop("Insufficient useable samples!")
  }
  names(w) <- NULL
  cov_used.partial <- sapply(w, function(wi) {wi[1] > 0})
  w <- lapply(w, function(wi) {as.matrix(wi[,-1,drop = FALSE])})
  w <- lapply(w, function(wi) {dimnames(wi) <- NULL; wi})
  w.pooled <- Reduce("+", w)
  if (type == 2) {
    # Check reciprocal condition number. If it suggests we cannot invert w, use 
    # the variances only as suggested by Scott et al 2016.
    if (rcond(w.pooled) <= cov.tol) {
      w.pooled <- diag(diag(w.pooled))
      cov_used.pooled <- FALSE
    } else {cov_used.pooled <- TRUE}
  } else {cov_used.pooled <- FALSE}
  w.pooled <- solve(w.pooled)

  # Weight the samples.

  if (return.pooled) {
    if (use_spark) {
      theta.w.spark <- spark_apply(theta,
        function(df, context) {
          index <- df[1,1]
          n <- nrow(df) - 1
          theta.w <- matrix(NA, n, ncol(df))
          theta.w[,1] <- index
          theta.w[,2:ncol(theta.w)] <- as.matrix(df[,2:ncol(df),drop = FALSE]) %*% context$w[[index]]
          theta.w[,2:ncol(theta.w)] <- theta.w[,2:ncol(theta.w),drop = FALSE] %*% context$w.pooled
          theta.w
        },
        context = ctx
      )
      if (return.pooled) {
        theta.w.pooled <- as.data.frame(theta.w.spark)
        theta.w.pooled <- split(theta.w.pooled[,-1,drop = FALSE], theta.w.pooled[,1])
        theta.w.pooled <- lapply(theta.w.pooled, FUN = as.matrix)
        theta.w.pooled <- rowSums(abind::abind(theta.w.pooled, along = 3), dims = 2)
      }
    } else if (use_parallel) {
      theta.w.pooled <- rowSums(abind::abind(
        clusterMap(
          paar.clust,
          fun = function(th, wi) {th %*% wi},
          theta,
          w,
          SIMPLIFY = FALSE
        ),
        along = 3
      ), dims = 2) %*% w.pooled
    } else {
      theta.w.pooled <- rowSums(abind::abind(
        mapply(
          FUN = function(th, wi) {th %*% wi},
          theta,
          w,
          SIMPLIFY = FALSE
        ),
        along = 3
      ), dims = 2) %*% w.pooled
    }
    dimnames(theta.w.pooled) <- NULL
  }

  out <- list(
    w = w,
    cov_used.partial = cov_used.partial,
    cov_used.pooled = cov_used.pooled
  )
  if (use_spark) out$theta.w.spark <- theta.w.spark
  if (return.pooled) out$theta.w.pooled <- theta.w.pooled
  if (use_parallel && new_cluster) stopCluster(par.clust)

  out
}

#' Compute consensus algorithm weight matrices
#'
#' This is a worker function for computing weight matrices for samples from 
#' partial posteriors. If samples are held in local memory, should be applied 
#' iteratively using a loop, but can be applied in parallel using the 
#' \code{parallel} package or \code{sparklyr} if samples are distributed as a 
#' Spark table.
#'
#' \code{context} should be a list with the following fields:
#' \describe{
#' \item{\code{type}}{Integer, either 1 or 2, specifying the weighting type to 
#' use. 1: use constant weighting. 2: weight samples using the sample covariance 
#' matrices of the partial posteriors.}
#' \item{\code{use_spark}}{Logical. \code{TRUE} if a Spark cluster is available.}
#' \item{\code{tol}}{The tolerance level for the reciprocal condition number of 
#' the covariance matrix. If less than this, the covariances are ignored.}
#' }
#'
#' @section References
#' \itemize{
#' \item{Scott, Steven L., Blocker, A.W., Bonassi, F.V., Chipman, H.A., George, E.I. and McCulloch, R.E., 2016. Bayes and big data: The consensus Monte Carlo algorithm. \emph{International Journal of Management Science and Engineering Management}, 11(2), pp.78-88.}
#' }
#'
#' @param df a matrix or data.frame of samples. If using Spark the first column 
#' should contain integers identifying the partial posterior.
#' @param context a list of parameters common to all workers. See details.
#'
#' @return Consensus algorithm weighting matrix for the sample. If using Spark 
#' the first column will contain the index of the partial posterior from which 
#' the sample was drawn. The next column is a binary flag indicating whether 
#' (1) or not (0) the full covariance matrix was used. The number of remaining 
#' columns is equal to the number of rows. These columns constitute a square 
#' matrix: the weighting matrix for consensus algorithm \code{type}.
#' @export
consensus.worker <- function(df, context) {
  # On Spark the first column is the part label identifying the data shard.
  w.part <- matrix(NA, ncol(df), ncol(df) + 1)
  # Insufficient useable samples. Return NA matrix to signal error.
  if (nrow(df) < 2) return(w.part)
  # The shard label.
  if (context$use_spark) w.part[,1] <- df[1,1]
  # Weight by sample covariance matrix.
  if (context$type == 2) {
    if (context$use_spark) {
      v <- cov(df[,-1,drop = FALSE])
    } else {
      v <- cov(df)
    }
    # Check reciprocal condition number. If it suggests we cannot invert v, 
    # ignore the covariances as suggest by Scott et al 2016.
    is_well_conditioned <- rcond(v) > context$tol
    if (!is_well_conditioned) v <- diag(diag(v))
    v <- solve(v)

  # Equal weighting.
  } else {
    v <- diag(nrow(w.part))
    is_well_conditioned <- TRUE
  }
  w.part[,1 + context$use_spark * 1] <- is_well_conditioned * 1
  w.part[,(2 + context$use_spark * 1):ncol(w.part)] <- v
  w.part
}

#' Consensus algorithm jackknife bias correction
#'
#' Compute the jackknife estimate of the small sample bias that arises in the 
#' consensus Monte Carlo algorithm. The result can be subtracted from the 
#' consensus weighted sample to correct for this bias.
#'
#' @seealso \code{consensus.weights}
#'
#' @param theta a matrix containing the pooled samples from all partial 
#' posteriors, weighted using the consensus algorithm. See output field 
#' \code{theta.w.pooled} from \code{\link{consensus.weights}}.
#' @param subsamples a list of matrix, each consisting of (unweighted) samples 
#' of model parameters from a single partial posterior, chosen at random from 
#' all available samples from that partial posterior and representing a 
#' proportion \code{alpha} of all available samples. See argument \code{theta} 
#' to \code{\link{consensus.weights}}.
#' @param alpha a number in [0,1], the proportion of samples taken from all 
#' available from each partial posterior to form each element of 
#' \code{subsamples}.
#' @param type an integer, either 1 or 2, specifying the weighting type to use. 
#' 1: use constant weighting. 2: weight samples using the sample covariance 
#' matrices of the partial posteriors. Must be the same as used to compute the 
#' weighted samples in \code{theta}.
#' 
#' @return A list with two elements: matrix \code{theta} of bias-corrected 
#' pooled samples, similar to the argument of the same name, and \code{bias}, 
#' the bias correction vector.
#' @export
bias_correction <- function (
  theta,
  subsamples,
  alpha,
  type = 2
) {
  subsamples.w.pooled <- consensus.master(subsamples, type = type, return.pooled = TRUE)$theta.w.pooled

  bias <- (colMeans(subsamples.w.pooled) - colMeans(theta)) * alpha / (1 - alpha)

  theta <- sweep(theta, MARGIN = 2, STATS = bias, FUN = "-", check.margin = FALSE)

  return(list(theta = theta, bias = bias))
}

