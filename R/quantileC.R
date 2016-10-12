# This function calls the C++ implementation of column quantile.



#' Fast colunm- and row-wise quantile of a matrix.
#' 
#' Fast calculation of column- and row-wise quantiles of a matrix at a single
#' probability. Implemented via compiled code, it is much faster than the
#' equivalent \code{apply(data, 2, quantile, prob = p)}.
#' 
#' At present, only one quantile type is implemented, namely the default type 7
#' used by R.
#' 
#' @aliases colQuantileC rowQuantileC
#' @param data a numerical matrix column-wise quantiles are desired. Missing
#' values are removed.
#' @param p a single probability at which the quantile is to be calculated.
#' @return A vector of length equal the number of columns (for
#' \code{colQuantileC}) or rows (for \code{rowQuantileC}) in \code{data}
#' containing the column- or row-wise quantiles.
#' @author Peter Langfelder
#' @seealso \code{\link[stats]{quantile}}
#' @keywords misc
colQuantileC = function(data, p)
{  
  data = as.matrix(data)

  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!")
  ncol = ncol(data)
  nrow = nrow(data)
  quantiles = rep(0, ncol)

  p = as.numeric(as.character(p))
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!")
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."))

  res = .C("quantileC", data = as.double(data), nrow = as.integer(nrow), ncol = as.integer(ncol), 
           p = as.double(p), quantiles = as.double(quantiles), NAOK = TRUE)

  res$quantiles
}

rowQuantileC = function(data, p)
{  
  data = as.matrix(data)

  #if (sum(is.na(data))>0) 
  #  stop("Missing values are not handled correctly yet. Sorry!")
  ncol = ncol(data)
  nrow = nrow(data)
  quantiles = rep(0, nrow)

  p = as.numeric(as.character(p))
  if (length(p) > 1)
    stop("This function only calculates one quantile at a time, for now. Sorry!")
  if ( (p<0) || (p>1) ) 
    stop(paste("Probability", p, "is out of the allowed range between 0 and 1."))

  res = .C("rowQuantileC", data = as.double(data), nrow = as.integer(nrow), ncol = as.integer(ncol), 
           p = as.double(p), quantiles = as.double(quantiles), NAOK = TRUE)

  res$quantiles
}
