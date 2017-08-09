# Parallel version of quantile, mean and median

.dimensions <- function(x) {
   if (is.null(dim(x))) {
       return(length(x))
   }
   return(dim(x))
}

.shiftList <- function(c0, lst) {
  if (length(lst) == 0) return(list(c0))
  ll <- length(lst)
  out <- list()
  out[[1]] <- c0
  for (i in 1:ll) {
    out[[i+1]] <- lst[[i]]
  }
  out
}





#' Parallel quantile, median, mean
#' 
#' 
#' 
#' @aliases pquantile pquantile.fromList pmedian pmean pmedian pmean
#' @param prob A number or vector of probabilities at which to calculate the
#' quantile. See \code{\link{quantile}}.
#' @param \dots Numeric arguments. All arguments must have the same dimensions.
#' See details.
#' @return A vector or array containing the quantiles, medians, or means. The
#' dimensions are determined by the dimensions of the input arguments and
#' whether the \code{prob} input is scalar or a vector.  If any of the input
#' variables have \code{dimnames}, the first non-NULL dimnames are copied into
#' the output. ======= pmin pmean.fromList pminWhich.fromList Parallel
#' quantile, median, mean
#' 
#' Calculation of ``parallel'' quantiles, minima, maxima, medians, and means,
#' across given arguments or across lists
#' 
#' pquantile(prob, ...) pquantile.fromList(dataList, prob) pmedian(...)
#' pmean(..., weights = NULL) pmin(...) pmean.fromList(dataList, weights =
#' NULL) pminWhich.fromList(dataList)
#' 
#' \item{prob}{ A single probability at which to calculate the quantile. See
#' \code{\link{quantile}}.  } \item{dataList}{A list of numeric vectors or
#' arrays, all of the same length and dimensions, over which to calculate
#' ``parallel'' quantiles.} \item{weights}{Optional vector of the same length
#' as \code{dataList}, giving the weights to be used in the weighted mean. If
#' not given, unit weights will be used.} \item{list()}{ Numeric arguments. All
#' arguments must have the same dimensions. See details. }
#' 
#' Given numeric arguments, say x,y,z, of equal dimensions (and length), the
#' \code{pquantile} calculates and returns the quantile of the first components
#' of x,y,z, then the second components, etc. Similarly, \code{pmedian} and
#' \code{pmean} calculate the median and mean, respectively. The funtion
#' \code{pquantile.fromList} is identical to \code{pquantile} except that the
#' argument \code{dataList} replaces the ... in holding the numeric vectors
#' over which to calculate the quantiles.
#' 
#' \item{pquantile, pquantile.fromList}{A vector or array containing
#' quantiles.} \item{pmean, pmean.fromList}{A vector or array containing means.
#' } \item{pmin, pmedian}{A vector or array containing minima, medians, and
#' maxima, respectively.} \item{pminWhich.fromList}{A list with two components:
#' \code{min} gives the minima, \code{which} gives the indices of the elements
#' that are the minima.}
#' 
#' Dimensions are copied from dimensions of the input arguments. If any of the
#' input variables have \code{dimnames}, the first non-NULL dimnames are copied
#' into the output. >>>>>>> 80351c0... version 1.60
#' 
#' Calculation of ``parallel'' quantiles, medians, and means, across given
#' arguments.
#' 
#' <<<<<<< HEAD
#' 
#' Given the argumens, say x,y,z, of equal dimensions, the \code{pquantile}
#' calculates and returns the quantile of the first components of x,y,z, then
#' the second components, etc. Similarly, \code{pmedian} and \code{pmean}
#' calculate the median and mean, respectively. =======
#' 
#' \code{\link{quantile}}, \code{\link{median}}, \code{\link{mean}} for the
#' underlying statistics. >>>>>>> 80351c0... version 1.60
#' 
#' # Generate 2 simple matrices a <- matrix(c(1:12), 3, 4) b <- a + 1 c <- a +
#' 2
#' 
#' # Set the colnames on matrix a colnames(a) <- paste0("col_", c(1:4))
#' 
#' # Example use pquantile(prob = 0.5, a, b, c) <<<<<<< HEAD pquantile(prob =
#' c(0, 0.5, 1), a, b, c) ======= >>>>>>> 80351c0... version 1.60
#' 
#' pmean(a, b, c) pmedian(a, b, c)
#' 
#' \code{\link{pmin}} and \code{\link{pmax}} for analogous functions for
#' minimum and maximum,
#' 
#' \code{\link{quantile}}, \code{\link{median}}, \code{\link{mean}} for the
#' underlying statistics.
#' 
#' Peter Langfelder and Steve Horvath
#' 
#' misc
pquantile <- function(prob, ...) {
   pars <- list(...)
   nPars <- length(pars)
   dn <- NULL
   for (p in 1:nPars) {
       if (mode(pars[[p]]) != "numeric")
          stop("Argument number", p, " is not numeric.")
       if (p == 1) {
          dims <- .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dims))) {
             stop("Argument dimensions are not consistent.")
          }
       }
       if (prod(dims) == 0) {
           stop("Argument has zero dimension.")
       }
       if (is.null(dn)) {
           dn <- dimnames(pars[[p]])
       }
       pars[[p]] <- as.numeric(pars[[p]])
   }
   x <- as.matrix(as.data.frame(pars))
   if (any(is.na(x))) {
      warning("The input contains missing data that will be removed.")
   }
   q <- apply(x, 1, quantile, prob = prob, na.rm = TRUE)
   rnq <- rownames(q)
   if (length(dims) > 1) {
       if (length(prob) ==  1) {
           dim(q) <- dims
       } else {
           dim(q) <- c(length(prob), dims)
       }

   }
   if (!is.null(dn)) {
       if (length(prob) == 1) {
           dimnames(q) <- dn
       } else {
           dimnames(q) <- .shiftList(rnq, dn)
       }
   }
   q
}

#' @export
#' @rdname pquantile
pmedian <- function(...) {
    pquantile(0.5, ...)
}

#' @export
#' @rdname pquantile
pmean <- function(...) {
   pars <- list(...)
   nPars <- length(pars)
   dn <- NULL
   for (p in 1:nPars) {
       if (mode(pars[[p]]) != "numeric") {
          stop(paste("Argument number", p, " is not numeric."))
       }
       if (p == 1) {
          dim <- .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dim))) {
             stop("Argument dimensions are not consistent.")
          }
       }
       if (prod(dim) == 0) {
           stop(paste("Argument has zero dimension."))
       }
       if (is.null(dn)) {
           dn <- dimnames(pars[[p]])
       }
       pars[[p]] <- as.numeric(pars[[p]])
   }
   x <- as.matrix(as.data.frame(pars))
   if (any(is.na(x))) {
      warning("The input contains missing data that will be removed.")
   }
   q <- rowMeans(x, na.rm = TRUE)
   if (length(dim) > 1) {
       dim(q) <- dim
   }
   if (!is.null(dn)) {
       dimnames(q) <- dn
   }
   q
}

