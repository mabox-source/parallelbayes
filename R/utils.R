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
    orig_class <- class(x)
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
