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
#' @section References
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
#' @return
#' @export
dpe.master <- function(
  theta,
  type = "sdpe",
  recursive = FALSE,
  subset.size = 2,
  cov.tol = .Machine$double.eps
) {

if (is.na(charmatch(tolower(type), c("ndpe", "sdpe")))) stop("type must match one of (ndpe, sdpe)!")
}
