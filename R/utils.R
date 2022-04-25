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

    if (class(x) == "array") {
	a <- dim(x)
	if (dim > length(x)) stop("Invalid dim!")
	M <- length(a)
	D <- setdiff(1:M, dim)
	max_x <- apply(x, D, max)
	s <- log(colSums(aperm(exp(sweep(x, D, max_x, "-", check.margin = FALSE)), c(dim, D)), dims = 1)) + max_x
	s[max_x == -Inf] <- -Inf
	a[dim] <- 1
	dim(s) <- a
    } else if (class(x) == "matrix") {
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
