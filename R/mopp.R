#' Calculate the MoPP algorithm sample importance weights
#'
#' Computes importance weights for the samples drawn from all partial posterior 
#' distributions to form a pooled, weighted sample that can be used to 
#' approximate expected values under the full data posterior distribution.
#'
#' Data \code{x} can be held in local memory or distributed on a cluster. 
#' Distributed computation is performed using Spark interfaced with the 
#' \code{sparklyr} package. In this case \code{x} is a Spark table, with each 
#' part containing a matrix of data (a data "shard"). If \code{x} are held in 
#' local memory it should be a list of matrices. In either case the matrices 
#' must have the same number of columns but can have different numbers of rows.
#'
#' Argument \code{loglik} is a list containing log likelihoods. It should have 
#' one element for each element of \code{x} (and \code{theta}), whose values 
#' are log likelihoods for all samples using just the data in the corresponding 
#' element of \code{x}. That is, each element of \code{loglik} contains log 
#' likelihoods for all samples in all elements of \code{theta}, so each element 
#' of \code{loglik} has the same number of rows, but computed using only a 
#' single data shard. This argument should be used when data shards are held by 
#' separate parties who are unable to communicate it to a single analyst, e.g. 
#' for privacy reasons. In this situation, all samples \code{theta} should be 
#' distributed amongst the data owners, who should compute the log likelihoods 
#' separately and send those to the analyst.
#'
#' Argument \code{loglik.fun} is a function that returns the log likelihood of 
#' parameter samples. It should have two arguments and the \code{...} argument. 
#' The first argument is a matrix of parameter samples in the same form as the 
#' elements of \code{theta}. The second argument is a matrix of data in the 
#' same form as the elements of \code{x}. The \code{...} should follow, and can 
#' be used in the definition of \code{loglik.fun} using the \code{list(...)} 
#' constuct. This allows additional named arguments to be supplied to the 
#' likelihood function (e.g. additional parameters). The returned value is a 
#' vector of log likelihoods with one value for each row of the first argument.
#'
#' Parameter samples \code{theta} should be sampled from the partial posterior 
#' distributions using proper priors rather than the "fractionated" priors of 
#' the consensus Monte Carlo algorithm of Scott et al 2016.
#'
#' @section Reusing output from \code{type == 1}:
#'
#' Some of the computational steps are common to the three versions of the 
#' algorithm. Therefore it can save computation to reuse the intermediate 
#' results when calling with \code{type == 2} or \code{type == 3} after an 
#' initial call with \code{type == 1}.
#'
#' If Laplace approximations were used in the initial run, the samples from 
#' those approximations should be supplied in elements of \code{theta}, 
#' appended after the samples from the partial posteriors. Arguments 
#' \code{laplace.type_1}, \code{laplace.type_2} and \code{laplace.type_3} 
#' should all be \code{FALSE}.
#'
#' @section References:
#' \itemize{
#' \item{Scott, Steven L., Blocker, A.W., Bonassi, F.V., Chipman, H.A., George, E.I. and McCulloch, R.E., 2016. Bayes and big data: The consensus Monte Carlo algorithm. \emph{International Journal of Management Science and Engineering Management}, 11(2), pp.78-88.}
#' }
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
#' @param loglik an optional list of numeric vectors, or single column 
#' matrices, whose values are log likelihoods. See details.
#' @param loglik.fun the log likelihood function. Optional if \code{loglik} is 
#' supplied. See details.
#' @param type an integer, either 1, 2 or 3, specifying the weighting type to 
#' use.
#' @param laplace.type_3.scale an optional matrix
#' @param laplace.type_3.dof an optional numeric. Must be greater than the dimension of the model plus 1.
#' @param w.type_1 an optional list of single column matrices containing the 
#' unnormalised MoPP type 1 weights, the same as the output field. Supply this 
#' to speed up computation of the type 2 weights if the type 1 weights have 
#' already been computated.
#' @param keep.type1 logical. If \code{TRUE} the type 1 weights are returned as 
#' well as the type 2/3 weights when \code{type = 2} or \code{type = 3}. This 
#' is faster than computing the different weights separately.
#' @param keep.unnormalised logical. If \code{TRUE} the unnormalised weights 
#' are returned as well as the normalised.
#' @param return.loglik logical. If \code{TRUE} the log likelihoods are 
#' returned as an element of the output list. The format is the same as the 
#' \code{loglik} argument, with elements as matrices.
#' @param par.clust an optional cluster connection object from package 
#' \code{parallel}. Ignored if \code{x} is a Spark table.
#' @param ncores an optional integer specifying the number of CPU cores to use 
#' (see \code{\link[parallel]{makeCluster}}). The default, 1, signifies that 
#' \code{parallel} will not be used.
#' @param ... additional parameters to be supplied to the log likelihood 
#' function. See details on \code{loglik.fun}.
#'
#' @return A list containing fields:
#' \item{Hvec}{A vector of the number of samples from each partial 
#' posterior.}
#' \item{wn.type_1, wn.type_2 or wn.type_3}{The normalised weighted of type 1, 
#' 2 or type 3, depending on the value of \code{type}.}
#' \item{wn.type_1}{Additionally returned if \code{keep.type1} is \code{TRUE}.}
#' \item{w.type_1}{The corresponding unnormalised weights, if 
#' \code{keep.unnormalised} is \code{TRUE} and \code{type = 1}}
#' \item{w.type_2}{The corresponding unnormalised weights, if 
#' \code{keep.unnormalised} is \code{TRUE} and \code{type = 2}}
#' \item{w.type_3}{The corresponding unnormalised weights, if 
#' \code{keep.unnormalised} is \code{TRUE} and \code{type = 3}}
#' \item{loglik}{log likelihoods, returned if \code{return.loglik} is 
#' \code{TRUE}}
#' \item{subsamples}{List of samples, like \code{theta}, after taking a 
#' subsample in the \code{type = 3} algorithm}
#' \item{subsample.inds}{List of sample indices indicating which were taken in 
#' the subsample in the \code{type = 3} algorithm}
#' \item{kl_hat}{Estimates of the KL divergence from the posterior distribution 
#' to each of the partial posteriors. Only returned if \code{type = 3}}
#' Weights and log likelihoods are returned as lists of matrices with 1 column 
#' and rows corresponding to sample. The list elements correspond to partial 
#' posterior (same as argument \code{theta}).
#'
#' @export
mopp.weights <- function(
  x,
  theta,
  loglik = NULL,
  loglik.fun = NULL,
  type = 1,
  subsample_size = NULL,
  laplace.type_1 = FALSE,
  laplace.type_2 = FALSE,
  laplace.type_3 = FALSE,
  laplace.type_2.sample_size = NULL,
  laplace.type_3.sample_size = NULL,
  laplace.type_3.scale = NULL,
  laplace.type_3.dof = NULL,
  params = NULL,
  w.type_1 = NULL,
  keep.type1 = TRUE,
  keep.unnormalised = FALSE,
  return.loglik = FALSE,
  par.clust = NULL,
  ncores = 1,
  verbose = FALSE,
  ...
) {

  ##############################################################################
  #  Setup.
  
  if (!("list" %in% class(theta))) stop("theta must be a list!")
  if (is.null(loglik) && is.null(loglik.fun)) stop("One of loglik or loglik.fun must be supplied!")
  if (!(type %in% 1:3)) stop("type must be 1, 2 or 3!")
  if (class(x)[1] == "tbl_spark") {
    if (!require(sparklyr)) stop("sparklyr is required!")
    use_spark <- TRUE
  } else {
    if (!("list" %in% class(x))) stop("x must be a Spark table or a list!")
    use_spark <- FALSE
  }
  if (!use_spark) par <- parallel.start(par.clust, ncores) 
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (is.null(params)) {
    params <- list(
      Hvec = sapply(theta, nrow),
      loglik.fun = loglik.fun,
      # Dimension of the model.
      d = ncol(theta[[1]]),
      n_shards = length(theta)
    )
    # Record which partial posterior each sample came from.
    params$pp.inds <- rep(1:length(params$Hvec), params$Hvec)
  }
  params$use_parallel <- par$valid
  # Parameter list for workers.
  ctx <- list(
    theta = do.call(rbind, theta),
    H = sum(params$Hvec),
    use_spark = use_spark,
    loglik.fun = loglik.fun,
    args = list(...)
  )
  if (any(params$Hvec < 1)) stop("Insufficient useable samples!")
  if (laplace.type_2 && is.null(laplace.type_2.sample_size)) laplace.type_2.sample_size <- max(params$Hvec)
  if (laplace.type_3 && is.null(laplace.type_3.sample_size)) laplace.type_3.sample_size <- max(params$Hvec)
  if (type == 3) {
    if (is.null(subsample_size)) {
      subsample_size <- min(params$Hvec)
      if (laplace.type_2 || laplace.type_3) subsample_size <- min(subsample_size, laplace.type_2.sample_size, laplace.type_3.sample_size)
    } else if (subsample_size > ctx$H) {
      stop("subsample_size is too large!")
    } else if (subsample_size > min(sapply(theta, nrow))) {
      warning("Using a subsample_size greater than the least number of samples from any partial posterior may result in a biased estimator.")
    }
  }
  

  # Sample from Laplace approximations.
# HVE NOT YET WORKED OUT HOW THIS WORKS WITH SPARK.
  if (laplace.type_1) {
    laplace.type_1.samples <- consensus.weights(
      theta,
      type = 2,
      return.pooled = TRUE,
      par.clust = par$par.clust,
    )$theta.w.pooled
    params$laplace.type_1.mean <- colMeans(laplace.type_1.samples)
    params$laplace.type_1.cov <- crossprod(sweep(laplace.type_1.samples, MARGIN = 2, STATS = params$laplace.type_1.mean, FUN = "-", check.margin = FALSE)) / (nrow(laplace.type_1.samples) - 1)
    ctx$theta <- rbind(ctx$theta, laplace.type_1.samples)
    ctx$H <- ctx$H + nrow(laplace.type_1.samples)
    params$pp.inds <- c(params$pp.inds, rep(max(params$pp.inds) + 1, nrow(laplace.type_1.samples)))
    params$Hvec <- c(params$Hvec, nrow(laplace.type_1.samples))
    # Records that Laplace type 1 was done.
    params$laplace.type_1 <- TRUE
  } else {
    laplace.type_1.samples <- matrix(NA, 0, params$d)
    if (is.null(params$laplace.type_1)) params$laplace.type_1 <- FALSE
  }
  if (laplace.type_2 || laplace.type_3) theta.pooled <- abind::abind(theta, along = 1)
  if (laplace.type_2) {
    params$laplace.type_2.mean <- colMeans(theta.pooled)
    # Sample covariance matrix.
    params$laplace.type_2.cov <- crossprod(sweep(theta.pooled, MARGIN = 2, STATS = params$laplace.type_2.mean, FUN = "-", check.margin = FALSE)) / (nrow(theta.pooled) - 1)
    # Sample from multivariate normal distribution.
    laplace.type_2.samples <- MASS::mvrnorm(laplace.type_2.sample_size, params$laplace.type_2.mean, params$laplace.type_2.cov)
    ctx$theta <- rbind(ctx$theta, laplace.type_2.samples)
    ctx$H <- ctx$H + laplace.type_2.sample_size
    params$pp.inds <- c(params$pp.inds, rep(max(params$pp.inds) + 1, laplace.type_2.sample_size))
    params$Hvec <- c(params$Hvec, laplace.type_2.sample_size)
    # Records that Laplace type 2 was done.
    params$laplace.type_2 <- TRUE
  } else {
    laplace.type_2.samples <- matrix(NA, 0, params$d)
    if (is.null(params$laplace.type_2)) params$laplace.type_2 <- FALSE
  }
  if (laplace.type_3) {
    if (is.null(laplace.type_3.scale)) laplace.type_3.scale <- matrix(2.5 ^ 2, params$d, params$d)
    if (is.null(laplace.type_3.dof)) laplace.type_3.dof <- params$d + 2
    params$laplace.type_3.mean <- colMeans(theta.pooled)
    S <- crossprod(abind::abind(lapply(theta, FUN = function(th) {sweep(th, MARGIN = 2, STATS = colMeans(th), FUN = "-", check.margin = FALSE)}), along = 1))
    # Sample from multivariate normal distribution.
    params$laplace.type_3.cov <- (S + laplace.type_3.scale) / (nrow(theta.pooled) + laplace.type_3.dof - params$d - 1)
    laplace.type_3.samples <- MASS::mvrnorm(laplace.type_3.sample_size, params$laplace.type_3.mean, params$laplace.type_3.cov)
    ctx$theta <- rbind(ctx$theta, laplace.type_3.samples)
    ctx$H <- ctx$H + laplace.type_3.sample_size
    params$pp.inds <- c(params$pp.inds, rep(max(params$pp.inds) + 1, laplace.type_3.sample_size))
    params$Hvec <- c(params$Hvec, laplace.type_3.sample_size)
    # Records that Laplace type 3 was done.
    params$laplace.type_3 <- TRUE
  } else {
    laplace.type_3.samples <- matrix(NA, 0, params$d)
    if (is.null(params$laplace.type_3)) params$laplace.type_3 <- FALSE
  }
  if (laplace.type_2 || laplace.type_3) rm(theta.pooled)
  
  
  ##############################################################################
  # Compute weights and normalise.
  
  # Log likelihoods for all samples.
  # Result in a list with elements corresponding to data shards (elements of 
  # data x).
  if (verbose) message("Computing likelihoods...")
  if (is.null(loglik)) {
    if (use_spark) {
      spark.res <- spark_apply(x, likelihood.worker, context = ctx)
      if (verbose) message("(Collecting...)")
      local.res <- as.data.frame(spark.res)
      # Split into list on the first column: data shard.
      loglik <- split(local.res[,-1], local.res[,1])
      loglik <- lapply(loglik, as.matrix)
    } else if (params$use_parallel) {
      loglik <- parallel::parLapply(par$par.clust, x, likelihood.worker, context = ctx)
    } else {
      loglik <- lapply(x, likelihood.worker, ctx)
    }
    if (verbose) message("Done.")
  } else {
    loglik <- lapply(loglik, FUN = as.matrix)
  }
  
  # Concatenate list elements into array.
  ll.array <- do.call(abind::abind, list(loglik, along = 3))
  
  # Multiply likelihoods.
  if (verbose) message("Computing unnormalised posterior densities...")
  w.numerator <- rowSums(ll.array, dims = 2)
  if (verbose) message("Done.")
  
  # Type 1 weights: divide out likelihood from partial posterior of origin; call 
  # it w.denominator.
  if (is.null(w.type_1)) {
    if (verbose) message("Computing type 1 weights...")
    # Weight denominator for non-Laplace samples.
    #w.denominator <- matrix(ll.array[c(outer(1:H, (0:(dim(ll.array)[2] - 1)) * H, FUN = "+")) + (params$pp.inds[1:H] - 1) * H * dim(ll.array)[2]], H, dim(ll.array)[2])
    H.nonlaplace <- sum(params$Hvec[1:params$n_shards])
    w.denominator <- matrix(ll.array[1:H.nonlaplace + (params$pp.inds[1:H.nonlaplace] - 1) * H.nonlaplace], H.nonlaplace, 1)
    # Weight denominator for Laplace samples.
    if (laplace.type_1) w.denominator <- c(w.denominator, mvtnorm::dmvnorm(laplace.type_1.samples, params$laplace.type_1.mean, params$laplace.type_1.cov, log = TRUE))
    if (laplace.type_2) w.denominator <- c(w.denominator, mvtnorm::dmvnorm(laplace.type_2.samples, params$laplace.type_2.mean, params$laplace.type_2.cov, log = TRUE))
    if (laplace.type_3) w.denominator <- c(w.denominator, mvtnorm::dmvnorm(laplace.type_3.samples, params$laplace.type_3.mean, params$laplace.type_3.cov, log = TRUE))
  
    w.type_1 <- matrix(w.numerator - w.denominator, dim(w.numerator)[1], dim(w.numerator)[2])
    # Need the shard specific sums of type 1 weights for normalisation of 
    # type 1 and for the mixture estimator in type 2. Split into list first to 
    # facilitate this.
    # Split only preserves dimensions on data.frames, so convert to df first.
    w.type_1 <- split(as.data.frame(w.type_1), params$pp.inds)
    w.type_1 <- lapply(w.type_1, as.matrix)
    w.type_1 <- lapply(w.type_1, FUN = function(ww){dimnames(ww) <- NULL; ww})
    names(w.type_1) <- NULL
    if (verbose) message("Done.")
  }
  norm <- mopp.normalise(
    w = w.type_1,
    type = 1,
    just_compute_constant = type %in% 2:3 && !keep.type1,
    par.clust = par$par.clust,
    verbose = verbose
  )
  w.sum_type_1 <- norm$w.sum
  wn.type_1 <- norm$wn
  
  # Type 2 weights.
  if (type == 2) {
    if (verbose) message("Computing type 2 weights...")
    # Want this as a matrix for type 2 calculations.
    w.sum_type_1 <- abind::abind(w.sum_type_1, along = 2)
    # Need to weight each partial posterior density in ll.array by the mean of 
    # the type 1 weights.
    # Note: division of w.sum_type_1 by n samples, to get the mean, cancels 
    # with multiplication by n samples in the mixture weights.
    type_2.mix <- sweep(ll.array, MARGIN = 2:3, STATS = w.sum_type_1, FUN = "+", check.margin = FALSE) - log(ctx$H)
    type_2.mix <- lrowsums(type_2.mix, 3, drop. = TRUE)
    w.type_2 <- matrix(w.numerator - type_2.mix, dim(w.numerator)[1], dim(w.numerator)[2])
    if (verbose) message("Done.")
  
    norm <- mopp.normalise(
      w = w.type_2,
      type = 2,
      pp.inds = params$pp.inds,
      just_compute_constant = FALSE,
      par.clust = par$par.clust,
      verbose = verbose
    )
    wn.type_2 <- norm$wn
    if (keep.unnormalised) {
      # Split only preserves dimensions on data.frames, so convert to df first.
      w.type_2 <- split(as.data.frame(w.type_2), params$pp.inds)
      w.type_2 <- lapply(w.type_2, as.matrix)
      w.type_2 <- lapply(w.type_2, FUN = function(ww){dimnames(ww) <- NULL; ww})
      names(w.type_2) <- NULL
    }

  } else if (type == 3) {
    if (verbose) message("Computing type 3 weights...")
    # Estimates of the KL divergences.
    # These are for the non-Laplace samples.
    kl_hat <- sapply(1:params$n_shards,
      FUN = function(j) {
        -sum(w.type_1[[j]]) / params$Hvec[j] +
        w.sum_type_1[[j]] - log(params$Hvec[j])
      }
    )
    # Add in the KL divergence estimates for the Laplace samples.
    if (params$laplace.type_1) {
      entropy <- 1 / 2 * determinant(2 * pi * exp(1) * params$laplace.type_1.cov, logarithm = TRUE)$modulus
      ind <- params$n_shards + 1
      kl_hat <- c(kl_hat, -sum(w.numerator[params$pp.inds == ind,,drop = FALSE]) / params$Hvec[ind] + w.sum_type_1[[ind]] - log(params$Hvec[ind]) - entropy)
    }
    if (params$laplace.type_2) {
      entropy <- 1 / 2 * determinant(2 * pi * exp(1) * params$laplace.type_2.cov, logarithm = TRUE)$modulus
      ind <- params$n_shards + params$laplace.type_1 + 1
      kl_hat <- c(kl_hat, -sum(w.numerator[params$pp.inds == ind,,drop = FALSE]) / params$Hvec[ind] + w.sum_type_1[[ind]] - log(params$Hvec[ind]) - entropy)
    }
    if (params$laplace.type_3) {
      entropy <- 1 / 2 * determinant(2 * pi * exp(1) * params$laplace.type_3.cov, logarithm = TRUE)$modulus
      ind <- params$n_shards + params$laplace.type_1 + params$laplace.type_2 + 1
      kl_hat <- c(kl_hat, -sum(w.numerator[params$pp.inds == ind,,drop = FALSE]) / params$Hvec[ind] + w.sum_type_1[[ind]] - log(params$Hvec[ind]) - entropy)
    }
    # Mixture component weights.
    q <- kl_hat
    # Non-positives could occur due to numerical error. Apply a correction.
    if (any(q <= 0)) {
      # Should only affect the Laplace approximations.
      if (any(which(q <= 0) <= params$n_shards)) stop("Non-positive KL divergence estimate detected for a partial posterior!")
      # Add 1 standard error of the minimum. Keep doing this until there are no 
      # non-positive.
      while (any(q <= 0)) {
        min.ind <- which.min(q)
        min.se <- sd(-w.numerator[params$pp.inds == ind,,drop = FALSE]) / sqrt(params$Hvec[ind])
        q <- q + min.se
      }
    }
    q <- 1 / q
    q <- q / sum(q)
    # Take subsample from theta.
    subsample.absolute_inds <- sample.int(ctx$H, size = subsample_size, prob = rep(q / params$Hvec, params$Hvec))
    split_inds <- findInterval(subsample.absolute_inds, cumsum(c(0, params$Hvec[1:(length(params$Hvec) - 1)])) + 1)
    subsample.inds <- split(subsample.absolute_inds, split_inds)
    if (laplace.type_1) theta <- c(theta, list(laplace.type_1.samples))
    if (laplace.type_2) theta <- c(theta, list(laplace.type_2.samples))
    if (laplace.type_3) theta <- c(theta, list(laplace.type_3.samples))
    if (!is.null(par.clust)) {
      subsample.inds <- parallel::clusterMap(
        par.clust,
        fun = function(inds, Hi) {inds - Hi},
        subsample.inds,
        cumsum(c(0, params$Hvec[1:(length(params$Hvec) - 1)])),
        SIMPLIFY = FALSE
      )
      subsamples <- parallel::clusterMap(
        par.clust,
        fun = function(th, inds) {th[inds,,drop = FALSE]},
        theta,
        subsample.inds,
        SIMPLIFY = FALSE
      )
    } else {
      subsample.inds <- mapply(
        FUN = function(inds, Hi) {inds - Hi},
        subsample.inds,
        cumsum(c(0, params$Hvec[1:(length(params$Hvec) - 1)])),
        SIMPLIFY = FALSE
      )
      subsamples <- mapply(
        FUN = function(th, inds) {th[inds,,drop = FALSE]},
        theta,
        subsample.inds,
        SIMPLIFY = FALSE
      )
    }
    # Want this as a matrix.
    w.sum_type_1 <- abind::abind(w.sum_type_1, along = 2)
    # Need to weight each partial posterior density in ll.array by the mean of 
    # the type 1 weights.
    # Note: division of w.sum_type_1 by n samples, to get the mean, cancels 
    # with multiplication by n samples in the mixture weights.
    type_3.mix <- sweep(sweep(ll.array[subsample.absolute_inds,,,drop = FALSE], MARGIN = 2:3, STATS = w.sum_type_1, FUN = "+", check.margin = FALSE), MARGIN = 3, STATS = log(q) - log(params$Hvec), FUN = "+", check.margin = FALSE)
    type_3.mix <- lrowsums(type_3.mix, 3, drop. = TRUE)
    w.type_3 <- matrix(w.numerator[subsample.absolute_inds,,drop = FALSE] - type_3.mix, subsample_size, dim(w.numerator)[2])
    if (verbose) message("Done.")
  
    norm <- mopp.normalise(
      w = w.type_3,
      type = 3,
      just_compute_constant = FALSE,
      par.clust = par$par.clust,
      verbose = verbose
    )
    wn.type_3 <- norm$wn
  }
  
  ############################################################################
  # Output.
  out <- list(params = params)
  if (type == 3) {
    out$subsamples <- subsamples
    out$subsample.inds <- subsample.inds
    if (keep.unnormalised) out$w.type_3 <- w.type_3
    out$wn.type_3 <- wn.type_3
    out$kl_hat <- kl_hat
  }
  if (type == 2) {
    if (keep.unnormalised) out$w.type_2 <- w.type_2
    out$wn.type_2 <- wn.type_2
  }
  if (type == 1 || keep.type1) {
    if (keep.unnormalised) out$w.type_1 <- w.type_1
    out$wn.type_1 <- wn.type_1
  }
  if (laplace.type_1) out$laplace.type_1.samples <- laplace.type_1.samples
  if (laplace.type_2) out$laplace.type_2.samples <- laplace.type_2.samples
  if (laplace.type_3) out$laplace.type_3.samples <- laplace.type_3.samples

  if (return.loglik) out$loglik <- loglik
  
  if (par$new) parallel::stopCluster(par$par.clust)
  
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
#' \describe{
#' \item{\code{theta}}{A matrix of all samples from the partial posteriors, 
#' pooled.}
#' \item{H}{Total number of samples (rows of \code{theta}).}
#' \item{\code{use_spark}}{Logical. \code{TRUE} if a Spark cluster is available.}
#' \item{\code{loglik.fun}}{The log likelihood function.}
#' \item{\code{args}}{A list of additional arguments to the log likelihood 
#' function. This list can be empty.}
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
    ll.part[,2] <- do.call(context$loglik.fun, c(list(context$theta, as.matrix(df[-nrow(df),-1,drop = FALSE])), context$args))
  } else {
    ll.part[,1] <- do.call(context$loglik.fun, c(list(context$theta, df), context$args))
  }
  ll.part
}

#' Monte Carlo estimate of the expectation of a univariate function
#'
#' Compute a Monte Carlo estimate of the expectation of a function of model 
#' parameters. Samples of the parameter vector are weighted using the MoPP type 
#' 1, 2 or 3 importance weights.
#'
#' \code{FUN} should take a matrix argument and return a matrix. The argument 
#' matrix should be a matrix of samples, just like the elements of 
#' \code{theta}. The returned matrix can have any number of columns but should 
#' have the same number of rows. The rows of both matrices correspond to 
#' samples.
#'
#' @seealso \code{mopp.quantile}, \code{mopp.weights}
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices (for \code{type = 1} or \code{type = 2}) or a 
#' matrix (\code{type = 3}), each containing normalised weights, one for each 
#' sample in \code{theta} and obtained using the \code{\link{mopp.weights}} 
#' function.
#' @param FUN a matrix-valued function that we wish to estimate the expected 
#' value of. See details.
#' @param type an integer, either 1, 2 or 3, specifying the weighting type used.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#'
#' @return A vector of weighted sample means of \code{FUN}, approximating the 
#' posterior expectation.
#' @export
mopp.mean <- function(
  theta,
  wn,
  FUN = identity,
  type = 2,
  Hvec = NULL
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (type != 3 && class(wn) != "list") stop("wn must be a list!")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (type != 3 && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (!(type %in% 1:3)) stop("type must be 1, 2 or 3!")
  if (any(sapply(theta, ncol) != ncol(theta[[1]]))) stop("All matrices in theta must have the same number of columns!")
  
  # In type 1 we need to add the mixture distribution weights (these are already 
  # implicit in the type 2 weight definition).
  if (type == 1) {
    if (is.null(Hvec)) {
      Hvec <- sapply(theta, nrow)
    } else if (length(Hvec) != length(theta)) {
      stop("There should be one element of Hvec for each element of theta!")
    }
    l_mixture_weight <- log(Hvec) - log(sum(Hvec))
    wn <- mapply(FUN = function(wi, m) {wi + m}, wn, l_mixture_weight, SIMPLIFY = FALSE)
  }
  
  # Collect samples and weights.
  theta <- abind::abind(theta, along = 1)
  if (type != 3) wn <- abind::abind(wn, along = 1)
  
  # Positivisation:
  # The expectation is estimated as a weighted sum. This can be computed 
  # without underflow if we take the log of FUN applied to theta. This requires 
  # us to temporarily multiply negative function values by -1.
  negatives <- FUN(theta) < 0
  
  # Apply FUN and take the weighted sum.
  d <- ncol(negatives)
  mu.hat <- rep(NA, d)
  for (j in 1:d) {
    if (any(negatives[,j])) {
      mu.hat[j] <- -exp(lrowsums(
          log(FUN(-theta[negatives[,j],,drop = FALSE])[,j]) + wn[negatives[,j]]
      ))
    } else {mu.hat[j] <- 0}
    if (any(!negatives[,j])) {
      mu.hat[j] <- mu.hat[j] + exp(lrowsums(
          log(FUN(theta[!negatives[,j],,drop = FALSE])[,j]) + wn[!negatives[,j]]
      ))
    }
  }
  
  return(mu.hat)
}

#' Monte Carlo estimate of posterior density
#'
#' Compute a KDE estimate of the posterior density evaluated at a number of 
#' values of a univariate random variable. The posterior is estimated with a 
#' Monte Carlo sample weighted by the MoPP type 1 or type 2 importance weights.
#'
#' The output of this function is the same as \code{mopp.mean} with a 
#' \code{FUN} argument being a Gaussian kernel function, evaluated over a range 
#' of values. This function is optimised to perform this without looping over 
#' target values.
#'
#' @seealso \code{mopp.mkde}, \code{mopp.mean}, \code{mopp.weights}
#'
#' @param x a vector of values for which the density estimate is to be computed.
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices (for \code{type = 1} or \code{type = 2}) or a 
#' matrix (\code{type = 3}), each containing normalised weights, one for each 
#' sample in \code{theta} and obtained using the \code{\link{mopp.weights}} 
#' function.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such 
#' that this is the standard deviation of the (Gaussian) smoothing kernel.
#' @param type an integer, either 0, 1, 2 or 3, specifying the weighting type 
#' used. Types 1, 2 and 3 refer to the MoPP weighting algorithms (see 
#' \code{mopp.weights}). Type 0 can be used to use uniform weights, which 
#' allows one to supply samples from another algorithm; in this case, \code{wn} 
#' is not required.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#' @param log. logical. If \code{TRUE}, density estimates are returned on the 
#' (natural) log scale.
#' @param mem.limit a positive value specifying the memory permitted for large 
#' arrays, in bytes. The function features a computational step that can 
#' consume a lot of memory. Increase this for a speed up or decrease it if 
#' memory is constrained.
#'
#' @return A vector the same length as \code{x} of density estimates.
#' @export
mopp.kde <- function(
  x,
  theta,
  wn,
  bw,
  type = 2,
  Hvec = NULL,
  log. = FALSE,
  mem.limit = 1024^3
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (!(type %in% 0:3)) stop("type must be 0, 1, 2 or 3!")
  if (!(type %in% c(0,3)) && class(wn) != "list") stop("wn must be a list!")
  if (!(type %in% c(0,3)) && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (type == 0) warning("Using uniform weights.")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (any(sapply(theta, ncol) != ncol(theta[[1]]))) stop("All matrices in theta must have the same number of columns!")
  if (bw <= 0) stop("bw must be positive!")

  if (is.null(Hvec)) {
    Hvec <- sapply(theta, nrow)
  } else if (length(Hvec) != length(theta)) {
    stop("There should be one element of Hvec for each element of theta!")
  }
  
  # In type 1 we need to add the mixture distribution weights (these are already 
  # implicit in the type 2 weight definition).
  if (type == 1) {
    l_mixture_weight <- log(Hvec) - log(sum(Hvec))
    wn <- mapply(FUN = function(wi, m) {wi + m}, wn, l_mixture_weight, SIMPLIFY = FALSE)
  }
  
  # Collect samples and weights.
  theta <- abind::abind(theta, along = 1)
  if (!(type %in% c(0, 3))) {
    wn <- unlist(wn)
  } else if (type == 3) {
    wn <- matrix(-log(H), H, 1)
  }
  
  n <- length(x)
  y <- rep(-Inf, n)
  buffer <- 1.5
  mem.req <- (n + sum(Hvec) * n) * 8 * buffer
  n_chunks <- min(ceiling(mem.req / mem.limit), n)
  nk <- ceiling(n / n_chunks)
  for (k in 1:n_chunks) {
    if ((k - 1) * nk + 1 > n) break
    chunk.inds <- ((k - 1) * nk + 1):min(k * nk, n)
    f <- outer(theta, x[chunk.inds], FUN = "-") / bw
    y[chunk.inds] <- lrowsums(-f ^ 2 / 2 + wn)
  }
  y <- y - (1 / 2) * log(2 * pi) - log(bw)

  if (!log.) y <- exp(y)

  return(y)
}

#' Monte Carlo estimate of multivariate posterior density
#'
#' Compute a KDE estimate of the posterior density evaluated at a number of 
#' values of a multivariate random variable. The posterior is estimated with a 
#' Monte Carlo sample weighted by the MoPP type 1 or type 2 importance weights.
#'
#' The output of this function is the same as \code{mopp.mean} with a 
#' \code{FUN} argument being a multivariate Gaussian kernel function, evaluated 
#' over a range of values. This function is optimised to perform this without 
#' looping over target values.
#'
#' @seealso \code{mopp.kde}, \code{mopp.mean}, \code{mopp.weights}
#'
#' @param x a vector of values for which the density estimate is to be computed.
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices (for \code{type = 1} or \code{type = 2}) or a 
#' matrix (\code{type = 3}), each containing normalised weights, one for each 
#' sample in \code{theta} and obtained using the \code{\link{mopp.weights}} 
#' function.
#' @param BW symmetric positive definite matrix, the smoothing bandwidth to be 
#' used, as a covariance matrix. The kernels are scaled such that this is the 
#' covariance matrix of the (multivariate Gaussian) smoothing kernel.
#' @param type an integer, either 0, 1, 2 or 3, specifying the weighting type 
#' used. Types 1, 2 and 3 refer to the MoPP weighting algorithms (see 
#' \code{mopp.weights}). Type 0 can be used to use uniform weights, which 
#' allows one to supply samples from another algorithm; in this case, \code{wn} 
#' is not required.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#' @param log. logical. If \code{TRUE}, density estimates are returned on the 
#' (natural) log scale.
#' @param mem.limit a positive value specifying the memory permitted for large 
#' arrays, in bytes. The function features a computational step that can 
#' consume a lot of memory. Increase this for a speed up or decrease it if 
#' memory is constrained.
#'
#' @return A vector the same length as \code{x} of density estimates.
#' @export
mopp.mkde <- function(
  x,
  theta,
  wn,
  BW,
  type = 2,
  Hvec = NULL,
  log. = FALSE,
  mem.limit = 1024^3
) {
  if (!(type %in% 0:3)) stop("type must be 0, 1, 2 or 3!")
  if (!(type %in% c(0,3)) && class(wn) != "list") stop("wn must be a list!")
  if (!(type %in% c(0,3)) && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (type == 0) warning("Using uniform weights.")
  if (class(theta) != "list") stop("theta must be a list!")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (any(sapply(theta, ncol) != ncol(theta[[1]]))) stop("All matrices in theta must have the same number of columns!")
  # Check BW is a symmetric positive definite matrix.
  if (class(BW) != "matrix") stop("BW must be a symmetric positive definite matrix!")
  if (!isSymmetric(BW) || any(eigen(BW, symmetric = TRUE)$values <= 0)) stop("BW must be a symmetric positive definite matrix!")

  if (is.null(Hvec)) {
    Hvec <- sapply(theta, nrow)
  } else if (length(Hvec) != length(theta)) {
    stop("There should be one element of Hvec for each element of theta!")
  }
  H <- as.numeric(sum(Hvec))
  d <- ncol(x)
  
  # In type 1 we need to add the mixture distribution weights (these are already 
  # implicit in the type 2 weight definition).
  if (type == 1) {
    l_mixture_weight <- log(Hvec) - log(H)
    wn <- mapply(FUN = function(wi, m) {wi + m}, wn, l_mixture_weight, SIMPLIFY = FALSE)
  }
  
  # Collect samples and weights.
  theta <- abind::abind(theta, along = 1)
  if (!(type %in% c(0, 3))) {
    wn <- unlist(wn)
  } else if (type == 3) {
    wn <- matrix(-log(H), H, 1)
  }

  #BW.eig <- eigen(BW)
  #BW.sqrt <- BW.eig$vectors %*% diag(sqrt(BW.eig$values)) %*% solve(BW.eig$vectors)
  #BW.prec <- solve(BW.sqrt)
  #BW.det <- determinant(BW)$modulus
  
  n <- nrow(x)
  y <- rep(-Inf, n)
  buffer <- 1.5
  mem.req <- (n + H * n + 2 * H * n * d) * 8 * buffer
  n_chunks <- min(ceiling(mem.req / mem.limit), n)
  nk <- ceiling(n / n_chunks)
  theta.rep <- array(rep(theta, each = nk), dim = c(nk, H, d))
  for (k in 1:n_chunks) {
    if ((k - 1) * nk + 1 > n) break
    chunk.inds <- ((k - 1) * nk + 1):min(k * nk, n)
    x.dist <- matrix(sweep(theta.rep[1:length(chunk.inds),,,drop = FALSE], MARGIN = c(1,3), STATS = x[chunk.inds,,drop = FALSE], FUN = "-", check.margin = FALSE), length(chunk.inds) * H, d)
    #x.dist <- t(BW.prec %*% t(matrix(sweep(theta.rep[1:length(chunk.inds),,,drop = FALSE], MARGIN = c(1,3), STATS = x[chunk.inds,,drop = FALSE], FUN = "-", check.margin = FALSE), length(chunk.inds) * H, d)))
    f <- sweep(
      matrix(mvtnorm::dmvnorm(x.dist, sigma = BW, log = TRUE, checkSymmetry = FALSE), length(chunk.inds), H),
      #matrix(mvtnorm::dmvnorm(x.dist, log = TRUE, checkSymmetry = FALSE), length(chunk.inds), H),
      MARGIN = 2,
      STATS = c(wn),
      FUN = "+",
      check.margin = FALSE
    )
    y[chunk.inds] <- lrowsums(f, 2)
  }
  #y <- y - 1 / 2 * BW.det

  if (!log.) y <- exp(y)

  return(y)
}

#' Monte Carlo estimate of a quantile of a univariate function
#'
#' Compute a Monte Carlo estimate of quantiles of a function of model 
#' parameters. Samples of the parameter vector are weighted using the MoPP type 
#' 1 or type 2 importance weights.
#'
#' Whilst \code{FUN} is matrix-valued, with possibly >1 columns, the estimated 
#' quantiles returned are quantiles of each dimension of the function 
#' separately, i.e. the marginals.
#'
#' \code{FUN} should take a matrix argument and return a matrix. The argument 
#' matrix should be a matrix of samples, just like the elements of 
#' \code{theta}. The returned matrix can have any number of columns but should 
#' have the same number of rows. The rows of both matrices correspond to 
#' samples.
#'
#' @seealso \code{mopp.mean}, \code{mopp.weights}
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param prob a probability, in [0, 1], corresponding to the quantile of 
#' \code{FUN} to be estimated.
#' @param wn a list of matrices (for \code{type = 1} or \code{type = 2}) or a 
#' matrix (\code{type = 3}), each containing normalised weights, one for each 
#' sample in \code{theta} and obtained using the \code{\link{mopp.weights}} 
#' function.
#' @param FUN a matrix-valued function that we wish to estimate the expected 
#' value of. See details.
#' @param type an integer, either 1, 2 or 3, specifying the weighting type used.
#' @param tol a number specifying the tolerance level for convergence of 
#' numerical root finding. See \code{\link{uniroot}} for details.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#'
#' @param A matrix with named rows matching the output of 
#' \code{\link{uniroot}}, the first row of which contains the weighted 
#' quantiles of \code{FUN}.
#' @export
mopp.quantile <- function(
  theta,
  prob,
  wn,
  FUN = identity,
  type = 2,
  tol = .Machine$double.eps ^ 0.5,
  Hvec = NULL
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (type != 3 && class(wn) != "list") stop("wn must be a list!")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (type != 3 && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (!(type %in% 1:3)) stop("type must be 1, 2 or 3!")
  if (any(sapply(theta, ncol) != ncol(theta[[1]]))) stop("All matrices in theta must have the same number of columns!")
  
  # In type 1 we need to add the mixture distribution weights (these are already 
  # implicit in the type 2 weight definition).
  if (type == 1) {
    if (is.null(Hvec)) {
      Hvec <- sapply(theta, nrow)
    } else if (length(Hvec) != length(theta)) {
      stop("There should be one element of Hvec for each element of theta!")
    }
    l_mixture_weight <- log(Hvec) - log(sum(Hvec))
    wn <- mapply(FUN = function(wi, m) {wi + m}, wn, l_mixture_weight, SIMPLIFY = FALSE)
  }
  
  # Collect samples and weights.
  theta <- abind::abind(theta, along = 1)
  if (type != 3) wn <- abind::abind(wn, along = 1)

  # Check dimension of function output.
  d <- ncol(FUN(theta[1,,drop = FALSE]))

  q.hat <- matrix(NA, 5, d, dimnames = list(c("root", "f.root", "iter", "init.it", "estim.prec")))
  
  mu.hat <- rep(NA, d)
  for (j in 1:d) {
    fs <- FUN(theta)[,j]
    a <- min(fs)
    fa <- lrowsums(wn[fs <= a]) - log(prob)
    if (fa > 0) {
      q.hat[,j] <- c(a, fa, 0, NA, 0)
      next
    }
    b <- max(fs)
    fb <- lrowsums(wn[fs <= b]) - log(prob)
    if (fb < 0) {
      q.hat[,j] <- c(b, fb, 0, NA, 0)
      next
    }
    q.hat[,j] <- unsplit(uniroot(
      function(z) {lrowsums(wn[fs <= z]) - log(prob)},
      lower = a,
      upper = b,
      tol = tol
    ), 1:5)
  }
  
  return(q.hat)
}

#' Normalise MoPP weights
#'
#' Performs the self-normalisation of the importance weights in MoPP, necessary 
#' for estimating posterior expectations with the weighted samples from partial 
#' posteriors.
#'
#' Weight normalisation is done differently for the two versions of MoPP, as 
#' specified by the \code{type} argument.
#'
#' The type 1 weights will be normalised relative to the partial posterior from 
#' which the corresponding samples were drawn. This means if there are $M$ 
#' partial posteriors, the sum of the normalised weights will be $M$.
#'
#' @param w either a single column matrix, if \code{type = 2} or 
#' \code{type = 3}, or a list of single column matrices, if \code{type = 1}, 
#' containing the unnormalised MoPP weights on the log scale.
#' @param type an integer, either 1, 2 or 3, specifying the weighting type used.
#' @param pp.inds integer vector, required for type 2 if \code{as.list} is 
#' \code{TRUE} (default), with one element for each row of \code{w}: the index 
#' of the partial posterior the corresponding sample is from.
#' @param just_compute_constant Logical. If \code{TRUE}, only the normalising 
#' constant(s) will be computed and returned - not the normalised weights.
#' @param as.list logical. If \code{TRUE} the normalised type 2 weights will be 
#' returned as a list, split according to \code{pp.inds}.
#' @param par.clust an optional cluster connection object from package 
#' \code{parallel}.
#' @return A list containing fields:
#' \item{wn}{A list of single column matrices, if \code{type = 1}, containing 
#' the normalised MoPP weights on the log scale. If \code{as.list} is 
#' \code{FALSE} and \code{type} is 2, these will be returned as a single 
#' matrix. If \code{type}, a single column matrix.}
#' \item{w.sum}{Either a numeric, if \code{type = 2} or \code{type = 3}, or a 
#' list of numerics, if \code{type = 1}, reporting the normalising constants 
#' used (log scale).}
#'
#' @export
mopp.normalise <- function(
  w,
  type,
  pp.inds = NULL,
  just_compute_constant = FALSE,
  as.list = TRUE,
  par.clust = NULL,
  verbose = FALSE
) {
  if (type == 1) {
    if (verbose) message("Computing type 1 weights' normalising constants...")
    if (!is.null(par.clust)) {
      w.sum <- parallel::parLapply(par.clust, w, fun = function(ww) {lrowsums(ww, 1)})
    } else {
      w.sum <- lapply(w, FUN = function(ww) {lrowsums(ww, 1)})
    }
    if (verbose) message("Done.")
    if (!just_compute_constant) {
      # Normalise weights and convert to a list with the same structure as theta.
      if (verbose) message("Normalising type 1 weights...")
      for (i in 1:length(w)) {
        if (any(!is.finite(w.sum[[i]]))) message(paste0("Sum of type 1 weights for shard ", i, " is zero. NaN returned for normalised weights."))
      }
      if (!is.null(par.clust)) {
        wn <- parallel::clusterMap(
            par.clust,
            fun = function(ww, ws) {sweep(ww, STATS = ws, MARGIN = 2, FUN = "-", check.margin = FALSE)},
            w,
            w.sum,
            SIMPLIFY = FALSE
          )
      } else {
        wn <- mapply(FUN = function(ww, ws) {sweep(ww, STATS = ws, MARGIN = 2, FUN = "-", check.margin = FALSE)}, w, w.sum, SIMPLIFY = FALSE)
      }
      if (verbose) message("Done.")
    } else {wn <- NULL}

  } else if (type %in% 2:3) {
    # Normalisation for type 2 is relative to whole sample.
    if (verbose) message(paste0("Computing type ", type, " normalising constants..."))
    w.sum <- lrowsums(w, 1)
    if (any(!is.finite(w.sum))) message(paste0("Sum of type ", type, " weights is zero. NaN returned for normalised weights."))
    if (verbose) message("Done.")
    if (!just_compute_constant) {
      if (verbose) message(paste("Normalising type ", type, " weights..."))
      wn <- sweep(w, STATS = w.sum, MARGIN = 2, FUN = "-", check.margin = FALSE)
      if (type == 2 && as.list) {
        # Split only preserves dimensions on data.frames, so convert to df first.
        wn <- split(as.data.frame(wn), pp.inds)
        wn <- lapply(wn, as.matrix)
        wn <- lapply(wn, FUN = function(ww){dimnames(ww) <- NULL; ww})
        names(wn) <- NULL
      }
      if (verbose) message("Done.")
    } else {wn <- NULL}

  } else {stop("type must be 1, 2 or 3!")}

  return(list(wn = wn, w.sum = w.sum))
}



#' Smooth the MoPP weights using the generalised Pareto distribution
#'
#' @section References:
#' \itemize{
#' \item{Vehtari, A., Simpson, D., Gelman, A. Yao, Y. and Gabry, J., 2015. Pareto smoothed importance sampling. \emph{arXiv preprint arXiv: 1507.02646.}}
#' \item{Vehtari, A., Gelman, A. and Gabry, J., 2017. Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. \emph{Statistics and computing}, 27(5), pp.1413-1432.}
#' }
#'
#' @param psis logical. If \code{TRUE} the Pareto smoothed importance sampling 
#' (PSIS) algorithm of Vehtari et al (2015) is applied to the raw importance 
#' weights.
#'
pareto_smooth <- function (

) {

}

