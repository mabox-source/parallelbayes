#' Single-dimension log-sums of an array-like object
#'
#' Given a vector, matrix or array of numeric values in the log domain, compute 
#' the sum of values on the linear domain across each row/column/higher 
#' dimension without risk of underflow.
#'
#' It is often useful in statistical computations to store probabilities on a 
#' logarithmic scale. However, when a sum of probabilities is required, 
#' performing the inverse log transformation before adding can result in 
#' numerical underflow, i.e. nonzero values are rounded to zero.
#'
#' This function computes these sums accurately by taking precaution against 
#' underflow, and does so in a vectorised manner for multidimensional objects.
#'
#' @param x vector, matrix of array of numbers to be summed.
#' @param dim integer specifying which dimension of \code{x} to take the sum 
#' over: \code{dim = 1} to sum values stored in columns, \code{dim = 2} to sum 
#' values stored in rows, etc. Ignored if \code{x} is a vector.
#' @param na.rm logical, irgnore NAs?
#' @param drop. logical, remove any dimensions of size 1 after summing?
#' @return An array-like object containing the log-sums of \code{x} taken over 
#' dimension \code{dim}. That is, if \code{y} is a vector stored across 
#' dimension \code{dim} of \code{x}, the value \code{log(sum(exp(y)))} is 
#' computed. See the examples for more explanation.
#' @seealso \code{rowsum}, \code{rowSums}
#' @examples
#' x <- log(array(1:24, dim = c(3,4,2)))
#' # This
#' lrowsums(x, 1)
#' # is equivalent to
#' log(apply(exp(x), 2:3, sum))
#' # And this
#' lrowsums(x, 3)
#' # is equivalent to
#' log(apply(exp(x), 1:2, sum))
#' @export
lrowsums <- function (x, dim = 1, na.rm = FALSE, drop. = FALSE) {

    if (!na.rm && anyNA(x)) stop("NA values detected in x!")
    if (any(is.infinite(x) & x > 0)) warning("Inf values detected in x!")
    if (na.rm) x[is.na(x)] <- -Inf

    if (class(x)[1] == "array") {
	a <- dim(x)
	if (dim > length(x)) stop("Invalid dim!")
	M <- length(a)
	D <- setdiff(1:M, dim)
	max_x <- apply(x, D, max)
	s <- log(colSums(aperm(exp(sweep(x, D, max_x, "-", check.margin = FALSE)), c(dim, D)), dims = 1)) + max_x
	s[max_x == -Inf] <- -Inf
	a[dim] <- 1
	dim(s) <- a
    } else if (class(x)[1] == "matrix") {
    # Sum down columns.
	if (dim == 1) {
	    max_x <- apply(x, 2, max)
	    s <- log(.rowSums(exp(t(x) - max_x), ncol(x), nrow(x))) + max_x
	} else if (dim == 2) {
	    max_x <- apply(x, 1, max)
	    s <- log(.rowSums(exp(x - max_x), nrow(x), ncol(x))) + max_x
	} else {
	    stop("Invalid dim!")
	}
	s[max_x == -Inf] <- -Inf
    } else {
	max_x <- max(x)
	if (max_x != -Inf) {
	    s <- log(sum(exp(x - max_x))) + max_x
	} else {
	    s <- -Inf
	}
    }
    if (drop.) s <- drop(s)
    return(s)
}

#' Partition matrix-like data into \code{r} parts
#'
#' Performs one of the following partitions: sequential (the first \code{n_i} 
#' elements go into part \code{i}, etc), unbalanced random (random allocation 
#' into parts which are allowed to have different size) or balanced random 
#' (random allocation into parts of the same size).
#'
#' If \code{x} is a matrix the partition is of sets of rows.
#'
#' Alternatively an integer scalar \code{x} can be supplied, in which case the 
#' partition is of the integers \code{1:x}.
#'
#' @param x a vector, matrix or data.frame of data to be partitioned.
#' @param part a vector the same length as \code{x} of part labels that can be 
#' used to specify how to partition \code{x}.
#' @param r an integer, the required number of parts.
#' @param part_weights a vector of length \code{r} specifying the relative 
#' sizes of the parts or probability of part membership in the unbalanced 
#' random allocation.
#' @param random logical. If \code{TRUE}, elements of \code{x} are assigned to 
#' parts at random.
#' @param balanced logical. If \code{TRUE}, and if \code{random} is 
#' \code{TRUE}, the parts will have the same number of elements (or as close as 
#' possible). \code{part_weights} will be ignored.
#' @return A list of nonintersecting subsets of \code{x}.
#' @export
partition <- function(x,
    r = 2,
    part = NULL,
    part_weights = rep(1, r),
    random = TRUE,
    balanced = FALSE
) {
    orig_class <- class(x)[1]
    orig_names <- unique(c(names(x), colnames(x)))
    x <- as.data.frame(x)
    n <- nrow(x)
    # If a scalar is supplied, partition the integers 1 to x.
    if (n == 1) {
	n <- unlist(x)
	x <- as.data.frame(1:n)
    }
    if (is.null(part)) {
	if (random && !balanced) {
	    part <- factor(sample(1:r, n, replace = TRUE, prob = part_weights), levels = 1:r)
	} else {
	    part_count <- floor(n * part_weights / sum(part_weights))
	    residue <- n - sum(part_count)
	    # Add residue to randomly chosen parts with probabilities given by the 
	    # part weights.
	    if (residue > 0) part_count <- part_count + c(rmultinom(1, residue, part_weights))
	    # Balanced random partition: just take random permutation of x first.
	    if (random) {
		x <- x[sample(1:n, n, replace = FALSE),,drop = FALSE]
	    }
	    part <- factor(rep(1:r, part_count), levels = 1:r)
	}
    }
    p <- split(x, part, drop = FALSE)

    if (orig_class == "matrix") {
	p <- lapply(p, FUN = as.matrix)
	p <- lapply(p, FUN = function(part) {colnames(part) = orig_names; part})
    } else if (orig_class != "data.frame") {
	p <- lapply(p, FUN = unlist)
	p <- lapply(p, FUN = function(part) {names(part) = NULL; part})
    }

    return(p)
}

#' Gibbs Sampling in a Logistic Regression Model
#'
#' Samples from the posterior distribution of model parameters (coefficients) 
#' in a logistic regression model with the logit link and using supplied data. 
#' This Gibbs sampling algorithm is due to Polson et al 2013.
#'
#' The model is specified via argument \code{params} (also returned as a field 
#' of the output list). \code{params} is a list with the following fields:
#'
#' \describe{
#' \item{\code{mu}}{Location parameter vector in multivariate normal prior for 
#' coefficients.}
#' \item{\code{Sigma}}{Covariance matrix in multivariate normal prior for 
#' coefficients.}
#' }
#'
#' @section References:
#' \itemize{
#' \item{Polson, N.G., Scott, J.G. and Windle, J., 2013. Bayesian inference for logistic model using P\out{&oacture;}lya-Gamma latent variables. \emph{Journal of the American statistical Association, 108} \bold{(504)}, pp.1339-1349.}
#' }
#'
#' @param x matrix of independent variables. Each column is interpreted as a 
#' predictor variable.
#' @param y integer or logical vector, the dependent variable.
#' @param weights an optional vector of sample weights to use in a binomial 
#' logistic regression where the dependent variable \code{y} is an integer 
#' vector. In this case \code{weights} are the corresponding case weights.
#' @param offset optional vector of offsets.
#' @param H an integer specifying the number of MCMC iterations and therefore 
#' of samples to draw.
#' @param params list of model parameters. See details.
#' @return A list containing matrix of samples, \code{samples}, and list 
#' \code{params}. 
#' @export
logistic.sampler <- function(
  y,
  x,
  weights = rep(1, nrow(x)),
  offset = NULL,
  H = 1000,
  params = list()
) {

  n <- nrow(x)
  m <- ncol(x)

  ##############################################################################
  # Setup model.
  if (is.null(params$mu)) params$mu <- rep(0, m)
  if (is.null(params$Sigma)) params$Sigma <- diag(m) * 2.5 ^ 2
  full_cond.mean_part <- crossprod(x, y - weights / 2)
  prior_prec <- solve(params$Sigma)
  samples <- list(
    theta = matrix(NA, H, m)
  )
  # Sample first value from prior distribution.
  samples$theta[1,] <- MASS::mvrnorm(1, params$mu, params$Sigma)

  ##############################################################################
  # Gibbs sampling.

  pb <- txtProgressBar(style = 3)
  for (h in 2:H) {

    ############################################################################
    # Sample Polya-gamma augmentation variables. See Polson et al 2013.
    if (is.null(offset)) {
      eta <- x %*% samples$theta[h - 1,]
    } else {
      eta <- x %*% samples$theta[h - 1,] + offset
    }
    aug <- BayesLogit::rpg(num = n, z = abs(eta), h = as.numeric(weights))

    ############################################################################
    # Sample coefficients from full conditional distribution.
    # Need covariance and mean.
    full_cond.cov <- crossprod(x, aug * x)
    full_cond.cov <- solve(full_cond.cov + prior_prec)
    if (is.null(offset)) {
      full_cond.mean <- full_cond.cov %*% (full_cond.mean_part + prior_prec %*% params$mu)
    } else {
      full_cond.mean <- full_cond.cov %*% (full_cond.mean_part - crossprod(x, aug * offset) + prior_prec %*% params$mu)
    }
    samples$theta[h,] <- MASS::mvrnorm(1, full_cond.mean, full_cond.cov)

    setTxtProgressBar(pb, (h - 1) / (H - 1))
  }
  close(pb)

  return(list(
    samples = samples,
    params = params
  ))
}
#' Verify an existing or initialise a new parallel cluster
#'
#' This is a wrapper of \code{parallel::detectCores} and 
#' \code{parallel::makeCluster}. 
#'
#' @param par.clust an optional cluster connection object from package 
#' \code{parallel}.
#' @param ncores an optional integer specifying the number of CPU cores to use 
#' (see \code{\link[parallel]{makeCluster}}).
#' @param forking logical. If \code{TRUE}, forking functions 
#' \code{\link[parallel]{mclapply}}, \code{\link[parallel]{mcmapply}} will be 
#' used and no connection object will be defined.
#'
#' @return A list containing fields:
#' \item{par.clust}{Cluster connection object from package \code{parallel}.}
#' \item{valid}{Logical. \code{TRUE} if \code{par.clust} is working.}
#' \item{new}{Logical. \code{TRUE} if \code{par.clust} was created by this 
#' function.}
#' \item{ncores}{Integer. The number of cores used by \code{par.clust}.}
#'
#' @export
parallel.start <- function(par.clust = NULL, ncores = 1, forking = FALSE) {
  if (!is.null(par.clust) && class(par.clust)[1] == "SOCKcluster" && require(parallel)) {
    valid <- TRUE
    new <- FALSE
  } else if (is.null(par.clust) && !forking && ncores > 1 && require(parallel)) {
    new <- TRUE
    n_cores_available <- parallel::detectCores()
    if (ncores > n_cores_available) {
      ncores <- n_cores_available
      message(paste0("Using the maximum number of CPU cores available (", n_cores_available, ")"))
    }
    par.clust <- parallel::makeCluster(ncores)
    valid <- TRUE
  } else {
    if (!forking && ncores > 1) message("Package parallel not found, using ncores = 1.")
    par.clust <- NULL
    valid <- FALSE
    new <- FALSE
  }

  return(list(par.clust = par.clust, valid = valid, new = new, ncores = ncores))
}
#' Multivariate normal density
#'
#' Probability density function for the multivariate normal distribution with 
#' mean \code{mean} and covariance matrix \code{sigma}.
#'
#' Although this functionality already exists in package mvtnorm's 
#' \code{\link[mvtnorm]{dmvnorm}}, that function was found to underflow in some 
#' parameterisations, even with \code{log = TRUE}. This does not underflow in 
#' those examples. This one might even be a bit faster.
#'
#' The precision matrix \code{prec} can be supplied as an argument as well as 
#' the covariance matrix to avoid unnecessary computation when this is already 
#' available. If not supplied it is computed from \code{sigma} using 
#' \code{\link{solve}}.
#'
#' @seealso \code{\link[mvtnorm]{dmvnorm}}
#'
#' @param x vector or matrix of quantiles. If \code{x} is a matrix, each row is 
#' taken to be a quantile.
#' @param mean mean vector, default is \code{rep(0, length = ncol(x))}.
#' @param sigma covariance matrix, default is ‘diag(ncol(x))’.
#' @param prec precision matrix, the matrix inverse of \code{sigma}. Optional. 
#' @param log logical; if \code{TRUE}, densities \code{d} are given as 
#' \code{log(d)}.
#'
#' @return A numeric vector containing the densities of \code{x} under the 
#' specified multivariate normal distribution.
#'
#' @examples
#' library(parallelbayes)
#' library(mvtnorm)
#' # These parameters cause underflow with dmvnorm:
#' m <- c(2.6, 6.3, -1.5, 1.6, 1, -0.8, 3.1, -1.6)
#' S <- matrix(c(20.1, 4.8, 5.6, 1.4, -1.5, 1.1, -1.1, 0.6, 4.8, 11.4, -1.6, 3.1, 0.8, -1.6, -2.9, 0.1, 5.6, -1.6, 3.5, 0.1, 0.1, -0.7, 1.1, -0.4, 1.4, 3.1, 0.1, 5.5, 2, 0.2, -0.8, 1.2, -1.5, 0.8, 0.1, 2, 10.6, 0.8, -3.1, -0.4, 1.1, -1.6, -0.7, 0.2, 0.8, 7.2, 1.2, 0.9, -1.1, -2.9, 1.1, -0.8, -3.1, 1.2, 3.5, 2.4, 0.6, 0.1, -0.4, 1.2, -0.4, 0.9, 2.4, 6.4), nrow = 8, ncol = 8)
#' dmvnorm(rep(0, 8), mean = m, sigma = S, log = TRUE)
#'
#' # The result should be:
#' dmnorm(rep(0, 8), mean = m, sigma = S, log = TRUE)
#'
#' @export
dmnorm <- function(
  x,
  mean = rep(0, ncol(x)),
  sigma = diag(ncol(x)),
  prec = NULL,
  log = FALSE
) {
  if (class(x)[1] != "matrix") x <- matrix(x, 1, length(x))
  if (is.null(prec)) prec <- solve(sigma)
  delta <- sweep(x, MARGIN = 2, STATS = mean, FUN = "-", check.margin = FALSE)
  dens <- -1 / 2 * determinant(2 * pi * sigma, logarithm = TRUE)$modulus - 1 / 2 * .rowSums(delta * t(tcrossprod(prec, delta)), nrow(x), ncol(x))
  if (!log) dens <- exp(dens)
  attributes(dens) <- NULL
  dens
}
