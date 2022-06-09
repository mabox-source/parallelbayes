#' Calculate the remix algorithm sample importance weights
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
#' @section References
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
#' @param type an integer, either 1 or 2, specifying the weighting type to use. 
#' @param w.type_1 an optional list of single column matrices containing the 
#' unnormalised remix type 1 weights, the same as the output field. Supply this 
#' to speed up computation of the type 2 weights if the type 1 weights have 
#' already been computated.
#' @param keep.type1 logical. If \code{TRUE} the type 1 weights are returned as 
#' well as the type 2 weights when \code{type = 2}. This is faster than 
#' computing the type 1 and type 2 weights separately.
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
#' \item{Hvec}{a vector of the number of samples from each partial 
#' posterior.}
#' \item{wn.type_1 or wn.type_2}{the normalised weighted of type 1 or type 2, 
#' depending on the value of \code{type}.}
#' \item{wn.type_1}{additionally returned if \code{keep.type1} is \code{TRUE}.}
#' \item{w.type_1 and/or w.type_2}{the corresponding unnormalised weights, if 
#' \code{keep.unnormalised} is \code{TRUE}
#' \item{loglik}{log likelihoods, returned if \code{return.loglik} is 
#' \code{TRUE}}
#' Weights and log likelihoods are returned as lists of matrices with 1 column 
#' and rows corresponding to sample. The list elements correspond to partial 
#' posterior (same as argument \code{theta}).
#'
#' @export
remix.weights <- function(
  x,
  theta,
  loglik = NULL,
  loglik.fun,
  type = 2,
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
  if (!(type %in% 1:2)) stop("type must be 1 or 2!")
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
  pp.inds <- rep(1:n_shards, Hvec)
  
  
  # Parameter list for workers.
  ctx <- list(
    theta = do.call(rbind, theta),
    H = H,
    use_spark = use_spark,
    loglik.fun = loglik.fun,
    args = list(...)
  )
  # Other parameters.
  params <- list(
    use_spark = use_spark,
    use_parallel = use_parallel,
    # Pool samples into one matrix.
    theta = do.call(rbind, theta),
    Hvec = Hvec,
    pp.inds = pp.inds,
    loglik.fun = loglik.fun,
    d = d
  )
  
  
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
    } else if (use_parallel) {
      loglik <- parLapply(par.clust, x, likelihood.worker, context = ctx)
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
    w.denominator <- matrix(ll.array[c(outer(1:H, (0:(dim(ll.array)[2] - 1)) * H, FUN = "+")) + (pp.inds - 1) * H * dim(ll.array)[2]], H, dim(ll.array)[2])
  
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
    wn.type_1 <- mapply(FUN = function(ww, ws) {sweep(ww, STATS = ws, MARGIN = 2, FUN = "-", check.margin = FALSE)}, w.type_1, w.sum_type_1, SIMPLIFY = FALSE)
    if (verbose) message("Done.")
  }
  
  # Type 2 weights.
  if (type == 2) {
    if (verbose) message("Computing type 2 weights...")
    # Want this as a matrix for type 2 calculations.
    w.sum_type_1 <- abind::abind(w.sum_type_1, along = 2)
    # Need to weight each partial posterior density in ll.array by the mean of 
    # the type 1 weights.
    # Note: division of w.sum_type_1 by n samples, to get the mean, cancels 
    # with multiplication by n samples in the mixture weights.
    type_2.mix <- sweep(ll.array, MARGIN = 2:3, STATS = w.sum_type_1, FUN = "+", check.margin = FALSE) - log(H)
    type_2.mix <- lrowsums(type_2.mix, 3, drop. = TRUE)
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
  if (type == 2) {
    if (keep.unnormalised) out$w.type_2 <- w.type_2
    out$wn.type_2 <- wn.type_2
  }
  if (type == 1 || keep.type1) {
    if (keep.unnormalised) out$w.type_1 <- w.type_1
    out$wn.type_1 <- wn.type_1
  }
  out$Hvec <- Hvec
  if (return.loglik) out$loglik <- loglik
  
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
#' parameters. Samples of the parameter vector are weighted using the remix 
#' type 1 or type 2 importance weights.
#'
#' \code{FUN} should take a matrix argument and return a matrix. The argument 
#' matrix should be a matrix of samples, just like the elements of 
#' \code{theta}. The returned matrix can have any number of columns but should 
#' have the same number of rows. The rows of both matrices correspond to 
#' samples.
#'
#' @seealso \code{remix.quantile}, \code{remix.weights}
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices, each containing normalised weights, one for 
#' each sample in \code{theta} and obtained using the 
#' \code{\link{remix.weights}} function.
#' @param FUN a matrix-valued function that we wish to estimate the expected 
#' value of. See details.
#' @param type an integer, either 1 or 2, specifying the weighting type used.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#'
#' @return A vector of weighted sample means of \code{FUN}, approximating the 
#' posterior expectation.
#' @export
remix.mean <- function(
  theta,
  wn,
  FUN = identity,
  type = 2,
  Hvec = NULL
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (class(wn) != "list") stop("wn must be a list!")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (!(type %in% 1:2)) stop("type must be 1 or 2!")
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
  wn <- abind::abind(wn, along = 1)
  
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
#' Monte Carlo sample weighted by the remix type 1 or type 2 importance 
#' weights.
#'
#' The output of this function is the same as \code{remix.mean} with a 
#' \code{FUN} argument being a Gaussian kernel function, evaluated over a range 
#' of values. This function is optimised to perform this without looping over 
#' target values.
#'
#' @seealso \code{remix.mkde}, \code{remix.mean}, \code{remix.weights}
#'
#' @param x a vector of values for which the density estimate is to be computed.
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices, each containing normalised weights, one for 
#' each sample in \code{theta} and obtained using the 
#' \code{\link{remix.weights}} function.
#' @param bw the smoothing bandwidth to be used. The kernels are scaled such 
#' that this is the standard deviation of the (Gaussian) smoothing kernel.
#' @param type an integer, either 0, 1 or 2, specifying the weighting type 
#' used. Types 1 and 2 refer to the remix weighting algorithms (see 
#' \code{remix.weights}). Type 0 can be used to use uniform weights, which 
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
remix.kde <- function(
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
  if (!(type %in% 0:2)) stop("type must be 0, 1 or 2!")
  if (type == 0) warning("Using uniform weights.")
  if (type != 0 && class(wn) != "list") stop("wn must be a list!")
  if (type != 0 && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
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
  if (type != 0) {
    wn <- unlist(wn)
  } else {
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
#' Monte Carlo sample weighted by the remix type 1 or type 2 importance 
#' weights.
#'
#' The output of this function is the same as \code{remix.mean} with a 
#' \code{FUN} argument being a multivariate Gaussian kernel function, evaluated 
#' over a range of values. This function is optimised to perform this without 
#' looping over target values.
#'
#' @seealso \code{remix.kde}, \code{remix.mean}, \code{remix.weights}
#'
#' @param x a vector of values for which the density estimate is to be computed.
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param wn a list of matrices, each containing normalised weights, one for 
#' each sample in \code{theta} and obtained using the 
#' \code{\link{remix.weights}} function.
#' @param BW symmetric positive definite matrix, the smoothing bandwidth to be 
#' used, as a covariance matrix. The kernels are scaled such that this is the 
#' covariance matrix of the (multivariate Gaussian) smoothing kernel.
#' @param type an integer, either 0, 1 or 2, specifying the weighting type 
#' used. Types 1 and 2 refer to the remix weighting algorithms (see 
#' \code{remix.weights}). Type 0 can be used to use uniform weights, which 
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
remix.mkde <- function(
  x,
  theta,
  wn,
  BW,
  type = 2,
  Hvec = NULL,
  log. = FALSE,
  mem.limit = 1024^3
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (!(type %in% 0:2)) stop("type must be 0, 1 or 2!")
  if (type == 0) warning("Using uniform weights.")
  if (type != 0 && class(wn) != "list") stop("wn must be a list!")
  if (type != 0 && any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
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
  if (type != 0) {
    wn <- unlist(wn)
  } else {
    wn <- matrix(-log(H), H, 1)
  }
  
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
    f <- sweep(
      matrix(mvtnorm::dmvnorm(x.dist, sigma = BW, log = TRUE, checkSymmetry = FALSE), length(chunk.inds), H),
      MARGIN = 2,
      STATS = c(wn),
      FUN = "+",
      check.margin = FALSE
    )
    y[chunk.inds] <- lrowsums(f, 2)
  }

  if (!log.) y <- exp(y)
  
  return(y)
}

#' Monte Carlo estimate of a quantile of a univariate function
#'
#' Compute a Monte Carlo estimate of quantiles of a function of model 
#' parameters. Samples of the parameter vector are weighted using the remix 
#' type 1 or type 2 importance weights.
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
#' @seealso \code{remix.mean}, \code{remix.weights}
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). Each list element corresponds to a single 
#' partial posterior.
#' @param prob a probability, in [0, 1], corresponding to the quantile of 
#' \code{FUN} to be estimated.
#' @param wn a list of matrices, each containing normalised weights, one for 
#' each sample in \code{theta} and obtained using the 
#' \code{\link{remix.weights}} function.
#' @param FUN a matrix-valued function that we wish to estimate the expected 
#' value of. See details.
#' @param type an integer, either 1 or 2, specifying the weighting type used.
#' @param tol a number specifying the tolerance level for convergence of 
#' numerical root finding. See \code{\link{uniroot}} for details.
#' @param Hvec an optional vector specifying the number of samples taken from 
#' each partial posterior.
#'
#' @param A matrix with named rows matching the output of 
#' \code{\link{uniroot}}, the first row of which contains the weighted 
#' quantiles of \code{FUN}.
#' @export
remix.quantile <- function(
  theta,
  prob,
  wn,
  FUN = identity,
  type = 2,
  tol = .Machine$double.eps ^ 0.5,
  Hvec = NULL
) {
  if (class(theta) != "list") stop("theta must be a list!")
  if (class(wn) != "list") stop("wn must be a list!")
  if (any(sapply(theta, class) != "matrix")) stop("theta must be a list of matrices!")
  if (any(sapply(wn, class) != "matrix")) stop("wn must be a list of matrices!")
  if (!(type %in% 1:2)) stop("type must be 1 or 2!")
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
  wn <- abind::abind(wn, along = 1)

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



#' Smooth the remix weights using the generalised Pareto distribution
#'
#' @section References
#' \itemize{
#' \item{Vehtari, A., Simpson, D., Gelman, A. Yao, Y. and Gabry, J., 2015. Pareto smoothed importance sampling. \emph{arXiv preprint arXiv: 1507.02646.}}
#' \item{Vehtari, A., Gelman, A. and Gabry, J., 2017. Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. \emph{Statistics and computing}, 27(5), pp.1413-1432.}
#' }
#'
#' @param psis logical. If \code{TRUE} the Pareto smoothed importance sampling 
#' (PSIS) algorithm of Vehtari et al (2015) is applied to the raw importance 
#' weights.
#'
#' @export
pareto_smooth <- function (

) {

}

