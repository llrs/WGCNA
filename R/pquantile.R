# Parallel version of quantile, mean and median

.dimensions = function(x) {
   if (is.null(dim(x))) {
       return(length(x))
   }
   return(dim(x))
}

.shiftList = function(c0, lst) {
  if (length(lst) == 0) return(list(c0))
  ll = length(lst)
  out = list()
  out[[1]] = c0
  for (i in 1:ll)
    out[[i+1]] = lst[[i]]
  out
}



#' Parallel quantile, median, mean
#'
#' Calculation of ``parallel'' quantiles, medians, and means, across given
#' arguments.
#'
#' Given the argumens, say x,y,z, of equal dimensions, the \code{pquantile}
#' calculates and returns the quantile of the first components of x,y,z, then
#' the second components, etc. Similarly, \code{pmedian} and \code{pmean}
#' calculate the median and mean, respectively.
#'
#' @aliases pquantile pmedian pmean
#' @param prob A number or vector of probabilities at which to calculate the
#' quantile. See \code{\link{quantile}}.
#' @param \dots Numeric arguments. All arguments must have the same dimensions.
#' See details.
#' @return A vector or array containing the quantiles, medians, or means. The
#' dimensions are determined by the dimensions of the input arguments and
#' whether the \code{prob} input is scalar or a vector.  If any of the input
#' variables have \code{dimnames}, the first non-NULL dimnames are copied into
#' the output.
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{pmin}} and \code{\link{pmax}} for analogous functions
#' for minimum and maximum,
#'
#' \code{\link{quantile}}, \code{\link{median}}, \code{\link{mean}} for the
#' underlying statistics.
#' @keywords misc
#' @examples
#'
#'
#' # Generate 2 simple matrices
#' a = matrix(c(1:12), 3, 4)
#' b = a+ 1
#' c = a + 2
#'
#' # Set the colnames on matrix a
#'
#' colnames(a) = paste0("col_", c(1:4))
#'
#' # Example use
#'
#' pquantile(prob = 0.5, a, b, c)
#' pquantile(prob = c(0, 0.5, 1), a,b, c)
#'
#' pmean(a,b,c)
#' pmedian(a,b,c)
#'
#'
pquantile = function(prob, ...) {
   pars = list(...)
   nPars = length(pars)
   dn = NULL
   for (p in 1:nPars) {
       if (mode(pars[[p]]) != "numeric")
          stop("Argument number", p, " is not numeric.")
       if (p == 1) {
          dims = .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dims))) {
             stop("Argument dimensions are not consistent.")
          }
       }
       if (prod(dims) == 0) {
           stop("Argument has zero dimension.")
       }
       if (is.null(dn)) {
           dn = dimnames(pars[[p]])
       }
       pars[[p]] = as.numeric(pars[[p]])
   }
   x = as.matrix(as.data.frame(pars))
   if (any(is.na(x))) {
      warning("The input contains missing data that will be removed.")
   }
   q = apply(x, 1, quantile, prob = prob, na.rm = TRUE)
   rnq = rownames(q)
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

pmedian = function(...) { pquantile(0.5, ...)}

pmean = function(...) {
   pars = list(...)
   nPars = length(pars)
   dn = NULL
   for (p in 1:nPars) {
       if (mode(pars[[p]]) != "numeric") {
          stop(paste("Argument number", p, " is not numeric."))
       }
       if (p == 1) {
          dim = .dimensions(pars[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(pars[[p]]), dim))) {
             stop("Argument dimensions are not consistent.")
          }
       }
       if (prod(dim) == 0) {
           stop(paste("Argument has zero dimension."))
       }
       if (is.null(dn)) {
           dn = dimnames(pars[[p]])
       }
       pars[[p]] = as.numeric(pars[[p]])
   }
   x = as.matrix(as.data.frame(pars))
   if (any(is.na(x))) {
      warning("The input contains missing data that will be removed.")
   }
   q = rowMeans(x, na.rm = TRUE)
   if (length(dim) > 1) {
       dim(q) = dim
   }
   if (!is.null(dn)) {
       dimnames(q) = dn
   }
   q
}

