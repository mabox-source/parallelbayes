#' The density product estimator algorithm of Neiswanger et al 2013
#'
#' The NDPE and SPDE algorithms use kernel density estimation to smooth samples 
#' from partial posterior distributions and samples from a mixture distribution 
#' constructed using partial posterior samples which approximates the full data 
#' posterior distribution.
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
#' @section References:
#' \itemize{
#' \item{Neiswanger, W., Wang, C. and Xing, E., 2013. Asymptotically exact, embarrassingly parallel MCMC. \emph{arXiv preprint arXiv: 1311.4780.}}
#' \item{Scott, Steven L., Blocker, A.W., Bonassi, F.V., Chipman, H.A., George, E.I. and McCulloch, R.E., 2016. Bayes and big data: The consensus Monte Carlo algorithm. \emph{International Journal of Management Science and Engineering Management}, 11(2), pp.78-88.}
#' }
#'
#' @param theta a list of matrices, each containing samples of model parameters 
#' from partial posterior distributions. Each matrix should have the same 
#' number of columns, which correspond to model parameters (including 
#' components of parameter vectors). They must also have the same number of 
#' rows (same number of samples from each partial posterior). Each list element 
#' corresponds to a single partial posterior.
#' @param type a character vector specifying the algorithm to use: "ndpe" for 
#' the nonparametric density product estimator or "sdpe" for the semiparametric 
#' density product estimator.
#' @param recursive logical. If \code{TRUE}, the algorithm is iterated 
#' recursively on subsets of available partial posterior samples. This may 
#' improve the acceptance rate of the independent Metropolis within Gibbs 
#' sampler used by the algorithm, but will take longer.
#' @param subset.size an integer, the number of subsets to take at a time in 
#' the recursive subsetting version of the algorithm used when \code{recursive} 
#' is \code{TRUE}. Otherwise ignored.
#' @param cov.tol a number, the tolerance level for the reciprocal condition 
#' number of the covariance matrix. If less than this, the covariances are 
#' ignored. See details.
#'
#' @return  A list containing fields:
#' \item{theta.pooled}{A matrix of pooled samples approximating the posterior 
#' distribution.}
#' \item{acceptance_rates}{Either a list of matrices, if \code{recursive} is 
#' \code{TRUE}, or a single matrix of Metropolis acceptance rates. One row for 
#' each partition of the data parts into subsets and one column for each 
#' \code{subset.size}. If a list, one element for each recursive step.}
#' \item{cov_used.partial}{Diagnostics (see Details).}
#' \item{cov_used.pooled}{Diagnostics (see Details).}
#'
#' @export
dpe.master <- function(
  theta,
  type = "sdpe",
  recursive = FALSE,
  subset.size = 2,
  cov.tol = .Machine$double.eps
) {

  if (is.na(charmatch(tolower(type), c("ndpe", "sdpe")))) stop("type must match one of (ndpe, sdpe)!")
  use_sdpe <- charmatch(tolower(type), c("ndpe", "sdpe")) == 2
  M <- length(theta)
  H <- nrow(theta[[1]])
  d <- ncol(theta[[1]])
  if (any(sapply(theta, nrow) != H)) stop("Each set of samples (matrices in theta) must have the same number of realisations (rows)!")

  if (use_sdpe) {
    worker.mean <- lapply(theta, FUN = function(th) {.colMeans(th, H, d)})
    worker.cov <- lapply(theta, FUN = function(th) {cov(th)})
    # Check reciprocal condition number. If it suggests we cannot invert the 
    # covariance matrix, use the variances only as suggested by Scott et al 
    # 2016.
    cov_used.partial <- sapply(worker.cov, rcond) > cov.tol
    if (any(!cov_used.partial)) worker.cov[!cov_used.partial] <- lapply(worker.cov[!cov_used.partial], FUN = function(w) {diag(diag(w))})
    worker.prec <- lapply(worker.cov, solve)
    cov_used.pooled <- TRUE
  }

  # Get a list of all the indices of partial posteriors to combine in each 
  # iteration of the algorithm.
  if (recursive) {
    n_parts <- ceiling(M / subset.size)
    subsets <- partition(M, n_parts, balanced = TRUE)
  } else {
    n_parts <- 1
    subset.size <- M
    subsets <- list(1:M)
  }

  theta.pooled <- list()
  component.accepted <- array(FALSE, dim = c(H, subset.size, n_parts))

  # Iterate over subset combinations (parts of the partition subsets).
  for (k in 1:n_parts) {
    # Might be a singleton: nothing to combine with in this iteration.
    if (length(subsets[[k]]) == 1) {
      theta.pooled[[k]] <- theta[[subsets[[k]]]]
      next
    }

    theta.pooled[[k]] <- matrix(NA, nrow = H, ncol = d)

    # KDE bandwidths.
    bw.vec <- (1:H) ^ (-1 / (d + 4))

    if (use_sdpe) {
      consensus.prec <- rowSums(abind::abind(worker.prec[subsets[[k]]], along = 3), dims = 2)
      # Check reciprocal condition number. If it suggests we cannot invert the 
      # covariance matrix, use the variances only as suggested by Scott et al 
      # 2016.
      if (rcond(consensus.prec) <= cov.tol) {
        w.pooled <- diag(diag(consensus.prec))
        cov_used.pooled <- FALSE
      } else {cov_used.pooled <- TRUE}
      consensus.cov <- solve(consensus.prec)
      consensus.mean <- consensus.cov %*% rowSums(abind::abind(
        mapply(
          FUN = function(a, b) {b %*% a},
          worker.mean[subsets[[k]]],
          worker.prec[subsets[[k]]],
          SIMPLIFY = FALSE
        ),
        along = 3
      ), dims = 2)
    }

    # Sample initial component index vector.
    component.inds <- ceiling(runif(subset.size) * H)
    component <- abind::abind(mapply(function(th, ind) {th[ind,]}, theta[subsets[[k]]], component.inds, SIMPLIFY = FALSE), along = 2)
    component.cov <- calc_covariance(d, bw.vec[1], subset.size, consensus.prec, use_sdpe)
    component.mean <- calc_mean(component, d, subset.size, component.cov, bw.vec[1], consensus.prec, consensus.mean, use_sdpe)
    component.w <- calc_weight(component, component.mean, bw.vec[1], d, subset.size, worker.mean[subsets[[k]]], worker.cov[subsets[[k]]], consensus.prec, consensus.mean, use_sdpe)

    # MCMC iterations.
    for (h in 1:H) {
      # Same for all components.
      component.cov <- calc_covariance(d, bw.vec[h], subset.size, consensus.prec, use_sdpe)

      # Sample new component vector.
      for (i in 1:subset.size) {
        proposal.inds <- component.inds
        proposal.inds[i] <- ceiling(runif(1) * H)
        proposal <- abind::abind(mapply(function(th, ind) {th[ind,]}, theta[subsets[[k]]], proposal.inds, SIMPLIFY = FALSE), along = 2)
        proposal.mean <- calc_mean(proposal, d, subset.size, component.cov, bw.vec[h], consensus.prec, consensus.mean, use_sdpe)
        proposal.w <- calc_weight(proposal, proposal.mean, bw.vec[h], d, subset.size, worker.mean[subsets[[k]]], worker.cov[subsets[[k]]], consensus.prec, consensus.mean, use_sdpe)

        # Check if proposal accepted.
        u <- log(runif(1))
        if (u < proposal.w - component.w) {
          component.inds <- proposal.inds
          component.mean <- proposal.mean
          component.w <- proposal.w
          component.accepted[h,i,k] <- TRUE
        }
      }

      # Sample from current component.
      theta.pooled[[k]][h,] <- MASS::mvrnorm(1, component.mean, component.cov)
    }
  }

  acceptance_rates <- colMeans(component.accepted, dims = 1)
  # If recursive: make these lists, because the second dimension of the 
  # matrices can be of variable size.
  if (recursive) acceptance_rates <- list(acceptance_rates)
  if (recursive && use_sdpe) cov_used.partial <- list(cov_used.partial)

  # Recursive call.
  if (recursive && M > 1) {
    out <- dpe.master(
      theta.pooled,
      type = type,
      recursive = recursive,
      subset.size = subset.size,
      cov.tol = cov.tol
    )
    theta.pooled[[1]] <- out$theta.pooled
    acceptance_rates <- c(acceptance_rates, out$acceptance_rates)
    if (use_sdpe) {
      cov_used.partial <- c(cov_used.partial, out$cov_used.partial)
      cov_used.pooled <- c(cov_used.pooled, out$cov_used.pooled)
    }
  }

  out <- list(
    theta.pooled = theta.pooled[[1]],
    acceptance_rates = acceptance_rates
  )
  if (use_sdpe) {
    out$cov_used.partial = cov_used.partial
    out$cov_used.pooled = cov_used.pooled
  }

  out
}
calc_ndpe_weight <- function(theta, theta.mean, bw) {
  w.ndpe <- -1 / 2 / bw ^ 2 * sum(sweep(theta, MARGIN = 1, STATS = theta.mean, FUN = "-", check.margin = FALSE) ^ 2)
}
calc_covariance <- function(D, bw, K, consensus.prec, use_sdpe) {
  if (use_sdpe) {
    solve(K / bw ^ 2 * diag(D) + consensus.prec)
  } else {
    bw ^ 2 / K * diag(D)
  }
}
calc_mean <- function(theta, D, K, S, bw, consensus.prec, consensus.mean, use_sdpe) {
  # theta: D by K matrix of parameter samples (dimension D) from each 
  # partial posterior sampler (K).
  # S: the component covariance found using calc_covariance().
  theta.mean <- .rowMeans(theta, D, K)
  if (use_sdpe) {
    S %*% (K / bw ^ 2 * theta.mean + consensus.prec %*% consensus.mean)
  } else {
    theta.mean
  }
}
calc_weight <- function(theta, component.mean, bw, D, K, worker.mean, worker.cov, consensus.cov, consensus.mean, use_sdpe) {
  # theta: D by K matrix of parameter samples (dimension D) from each 
  # partial posterior sampler (K).
  if (use_sdpe) {
    theta.mean <- .rowMeans(theta, D, K)
    w.ndpe <- calc_ndpe_weight(theta, theta.mean, bw)
    dens.mean <- mvtnorm::dmvnorm(theta.mean, consensus.mean, consensus.cov + bw ^ 2 / K * diag(D), log = TRUE)
    denom <- 0
    for (k in 1:K) {
      denom <- denom + mvtnorm::dmvnorm(theta[,k], worker.mean[[k]], worker.cov[[k]], log = TRUE)
    }
    w.ndpe + dens.mean - denom
  } else {
    calc_ndpe_weight(theta, component.mean, bw)
  }
}

