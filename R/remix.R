#' Calculate the remix algorithm sample importance weights
#'
#' @param x either a list of matrices or a Spark table (class \code{tbl_spark}) 
#' containing the partitioned data. If a list, each element corresponds to a 
#' single part. If a Spark table, the first column of each element must be an 
#' integer vector identifying the part.
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param loglik the log likelihood function. See details.
#' @param type an integer, either 1 or 2, specifying the weighting type to use. 
#' 1 corresponds to ___, 2 (the default) corresponds to ___.
#' @param keep.type1 logical. If \code{TRUE} the type 1 weights are returned as 
#' well as the type 2 weights when \cod{algorithm = 2}. This is faster than 
#' computing the type 1 and type 2 weights separately.
#' @param keep.unnormalised logical. If \code{TRUE} the unnormalised weights 
#' are returned as well as the normalised.
#' @param par.clust an optional cluster connection object from package 
#' \code{parallel}. Ignored if \code{x} is a Spark table.
#' @param ncores an optional integer specifying the number of CPU cores to use 
#' (see \code{\link[parallel]{makeCluster}}). The default, 1, signifies that 
#' \code{parallel} will not be used.
#'
#' @export
calc_weights <- function(
  x,
  theta,
  loglik,
  type = 2,
  keep.type1 = TRUE,
  keep.unnormalised = FALSE,
  par.clust = NULL,
  ncores = 1,
  verbose = FALSE
) {

  ##############################################################################
  #  Setup.
  
  if (!("list" %in% class(theta))) stop("theta must be a list!")
  if (class(x)[1] == "tbl_spark") {
    if (!require(sparklyr)) stop("sparklyr is required!")
    use_spark <- TRUE
  } else {
    if (!("list" %in% class(x))) stop("x must be a Spark table or a list!")
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
  
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  n_shards <- length(theta)
  Hvec <- sapply(theta, nrow)
  if (any(Hvec < 1)) stop("Insufficient useable samples!")
  H <- sum(Hvec)
  # Dimension of the model.
  d <- ncol(theta[[1]])
  # Record which partial posterior each sample came from.
  pp_inds <- rep(1:n_shards, Hvec)
  
  
  # Parameter list for workers.
  params <- list(
    use_spark = use_spark,
    use_parallel = use_parallel,
    do_type_2 = do_type_2,
    # Pool samples into one matrix.
    theta = do.call(rbind, theta),
    H = H,
    Hvec = Hvec,
    pp_inds = pp_inds,
    loglik = loglik,
    d = d
  )
  
  
  ##############################################################################
  # Compute weights and normalise.
  
  # Log likelihoods for all samples.
  # Result in a list with elements corresponding to data shards (elements of 
  # data x).
  if (verbose) message("Computing likelihoods...")
  if (use_spark) {
    spark.res <- spark_apply(x, likelihood.worker, context = params)
    if (verbose) message("(Collecting...)")
    local.res <- as.data.frame(spark.res)
    # Split into list on the first column: data shard.
    ll <- split(local.res[,-1], local.res[,1])
    ll <- lapply(ll, as.matrix)
  } else if (use_parallel) {
    ll <- parLapply(par.clust, x, likelihood.worker, context = params)
  } else {
    ll <- lapply(x, likelihood.worker, params)
  }
  if (verbose) message("Done.")
  
  # Concatenate list elements into array.
  ll.array <- do.call(abind::abind, list(ll, along = 3))
  
  # Multiply likelihoods.
  if (verbose) message("Computing unnormalised posterior densities...")
  w.numerator <- rowSums(ll.array, dims = 2)
  if (verbose) message("Done.")
  
  # Type 1 weights: divide out likelihood from partial posterior of origin; call 
  # it w.denominator.
  if (verbose) message("Computing type 1 weights...")
  w.denominator <- matrix(ll.array[c(outer(1:H, (0:(dim(ll.array)[2] - 1)) * params$H, FUN = "+")) + (pp_inds - 1) * params$H * dim(ll.array)[2]], H, dim(ll.array)[2])
  
  w.type_1 <- matrix(w.numerator - w.denominator, dim(w.numerator)[1], dim(w.numerator)[2])
  # Need the shard specific sums of type 1 weights for normalisation of 
  # type 1 and for the mixture estimator in type 2. Split into list first to 
  # facilitate this.
  # Split only preserves dimensions on data.frames, so convert to df first.
  w.type_1 <- split(as.data.frame(w.type_1), params$pp_inds)
  w.type_1 <- lapply(w.type_1, as.matrix)
  w.type_1 <- lapply(w.type_1. FUN = function(ww){dimnames(ww) <- NULL; ww})
  names(w.type_1) <- NULL
  if (verbose) message("Done.")
  # List of vectors.
  if (verbose) message("Computing type 1 weights' normalising constants...")
  w.sum_type_1 <- lapply(w.type_1, FUN = function(ww) {lrowsums(ww, 1)})
  if (verbose) message("Done.")
  # Normalise weights and convert to a list with the same structure as theta.
  if (type == 1 || keep.type1) {
    if (verbose) message("Normalising type 1 weights...")
    for (i in 1:n_shards) {
      if (any(!is.finite(w.sum_type_1[[i]]))) message(paste0("Sum of type 1 weights for shard ", i, " is zero. NaN returned for normalised weights."))
    }
  }
  wn.type_1 <- mapply(FUN = function(ww, ws) {sweep(ww, STATS = ws, MARGIN = 2, FUN = "-", check.margin = FALSE)}, w.type_1, w.sum_type_1, SIMPLIFY = FALSE)
  if (verbose) message("Done.")
  
  # Type 2 weights.
  if (type == 2) {
    if (verbose) message("Computing type 2 weights...")
    # Want this as a matrix for type 2 calculations.
    w.sum_type_1 <- abind::abind(w.sum_type_1, along = 2)
    # Need to weight each partial posterior density in ll.array by the mean of 
    # the type 1 weights.
    # Note: division of w.sum_type_1 by n samples, to get the mean, cancels 
    # with multiplication by n samples in the mixture weights.
    type_2.mix <- sweep(ll.array, MARGIN = 2:3, STATS = w.sum_type_1, FUN = "+", check.margin = FALSE) - log(params$H)
    type_2.mix <- lrowsums(type_2.mix, 3)
    w.type_2 <- matrix(w.numerator - type_2.mix, dim(w.numerator)[1], dim(w.numerator)[2])
    if (verbose) message("Done.")
  
    # Normalisation for type 2 is relative to whole sample.
    if (verbose) message("Computing type 2 normalising constants...")
    w.sum_type_2 <- lrowsums(w.type_2, 1)
    if (any(!is.finite(w.sum_type_2))) message("Sum of type 2 weights is zero. NaN returned for normalised weights.")
    if (verbose) message("Done.")
    if (verbose) message("Normalising type 2 weights...")
    wn.type_2 <- sweep(w.type_2, STATS = w.sum_type_2, MARGIN = 2, FUN = "-", check.margin = FALSE)
    # Split only preserves dimensions on data.frames, so convert to df first.
    w.type_2 <- split(as.data.frame(w.type_2), params$pp.inds)
    w.type_2 <- lapply(w.type_2, as.matrix)
    w.type_2 <- lapply(w.type_2, FUN = function(ww){dimnames(ww) <- NULL; ww})
    names(w.type_2) <- NULL
    wn.type_2 <- split(as.data.frame(wn.type_2), params$pp.inds)
    wn.type_2 <- lapply(wn.type_2, as.matrix)
    wn.type_2 <- lapply(wn.type_2, FUN = function(ww){dimnames(ww) <- NULL; ww})
    names(wn.type_2) <- NULL 
    if (verbose) message("Done.")
  }
  
  ############################################################################
  # Output.
  out <- list()
  if (algorithm == 2) {
    if (keep.unnormalised) out$w.type_2 <- w.type_2
    out$wn.type_2 <- wn.type_2
  }
  if (algorithm == 1 ZZ keep.type1) {
    if (keep.unnormalised) out$w.type_1 <- w.type_1
    out$wn.type_1 <- wn.type_1
  }
  out$Hvec <- Hvec
  
  if (use_parallel && new_cluster) parallel::stopCluster(par.clust)
  
  return(out)
}

#' Compute all samples' log likelihoods using a single shard of data
#'
#' This is a worker function for computing sample log likelihoods when data are 
#' distributed. After having drawn samples from each partial posterior, they 
#' should be collated into a single matrix. Then this function can be applied 
#' to all shards of data in parallel, computing the log likelihoods for all 
#' samples, for each shard of data. This is the first step in the weight 
#' computation.
#'
#' \code{context} should be a list with the following fields:
#' \itemize{
#' \item{\code{theta}: a matrix of all samples from the partial posteriors, 
#' pooled.}
#' \item{\code{use_spark}: logical. \code{TRUE} if a Spark cluster is available.}
#' \item{\code{f}: the log likelihood function.}
#' }
#'
#' @param df a matrix or data.frame of data. If using Spark the first column 
#' should contain integers identifying the partial posterior.
#' @param context a list of parameters common to all workers. See details.
#'
#' @return Matrix of log likelihoods for each sample in \code{context$theta}. 
#' If using Spark the first column will contain the index of the partial 
#' posterior from which the sample was drawn. The second column, or the only 
#' column if not using Spark, will contain the log likelihoods. 
likelihood.worker <- function(df, context) {
  # On Spark: ignore last row (zeros) and ID column.
  #n <- nrow(df) - context$use_spark * 1
  # On Spark the first column is the part label identifying the data shard.
  ll.part <- matrix(NA, context$H, 1 + context$use_spark * 1)
  if (context$use_spark) {
    # The shard label.
    ll.part[,1] <- df[1,1]
    ll.part[,2] <- context$loglik(context$theta, as.matrix(df[-nrow(df),-1,drop = FALSE]))
  } else {
    ll.part[,1] <- context$loglik(context$theta, df)
  }
  ll.part
}

#' Monte Carlo estimate of the expectation of a univariate function
#'
#' Compute a Monte Carlo estimate of a one dimensional function applied to each 
#' dimension of a parameter vector. Samples of the parameter vector are 
#' weighted using the remix type 1 or type 2 importance weights.
#'
#'
#' @seealso \code{calc_weights}

# lFUN only needs to return on the log scale - the argument does not need to be.



#' @section References
#' \itemize{
#' \item{Vehtari, A., Gelman, A. and Gabry, J., 2017. Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. \emph{Statistics and computing}, 27(5), pp.1413-1432.}
#' }
#'

#' @param psis logical. If \code{TRUE} the Pareto smoothed importance sampling 
#' (PSIS) algorithm of Vehtari et al (2017) is applied to the raw importance 
#' weights.
#'
#' @export
pareto_smooth <- function (

) {

}

