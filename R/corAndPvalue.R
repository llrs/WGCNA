# Functions to calculate correlation and corresponding p-values. Geared towards cor p-values of large
# matrices where varying numbers of missing data make the number of observations vary for each pair of
# columns.


#' Calculation of correlations and associated p-values
#' 
#' A faster, one-step calculation of Student correlation p-values for multiple
#' correlations, properly taking into account the actual number of
#' observations.
#' 
#' The function calculates correlations of a matrix or of two matrices and the
#' corresponding Student p-values. The output is not as full-featured as
#' \code{\link{cor.test}}, but can work with matrices as input.
#' 
#' @param x a vector or a matrix
#' @param y a vector or a matrix. If \code{NULL}, the correlation of columns of
#' \code{x} will be calculated.
#' @param use determines handling of missing data. See \code{\link{cor}} for
#' details.
#' @param alternative specifies the alternative hypothesis and must be (a
#' unique abbreviation of) one of \code{"two.sided"}, \code{"greater"} or
#' \code{"less"}.  the initial letter.  \code{"greater"} corresponds to
#' positive association, \code{"less"} to negative association.
#' @param \dots other arguments to the function \code{\link{cor}}.
#' @return A list with the following components, each a matrix: \item{cor}{the
#' calculated correlations} \item{p}{the Student p-values corresponding to the
#' calculated correlations} \item{Z}{Fisher transforms of the calculated
#' correlations} \item{t}{Student t statistics of the calculated correlations}
#' \item{nObs}{Numbers of observations for the correlation, p-values etc.}
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{cor}} for calculation of correlations only;
#' 
#' \code{\link{cor.test}} for another function for significance test of
#' correlations
#' @references Peter Langfelder, Steve Horvath (2012) Fast R Functions for
#' Robust Correlations and Hierarchical Clustering.  Journal of Statistical
#' Software, 46(11), 1-17.  \url{http://www.jstatsoft.org/v46/i11/}
#' @keywords stats
#' @examples
#' 
#' # generate random data with non-zero correlation
#' set.seed(1);
#' a = rnorm(100);
#' b = rnorm(100) + a;
#' x = cbind(a, b);
#' # Call the function and display all results
#' corAndPvalue(x)
#' # Set some components to NA
#' x[c(1:4), 1] = NA
#' corAndPvalue(x)
#' # Note that changed number of observations.
#' 
corAndPvalue = function(x, y = NULL, 
                        use = "pairwise.complete.obs", 
                        alternative = c("two.sided", "less", "greater"), 
                        ...)
{ 
  ia = match.arg(alternative)


#' Fast calculations of Pearson correlation.
#' 
#' These functions implements a faster calculation of Pearson correlation.
#' 
#' The speedup against the R's standard \code{\link[stats]{cor}} function will
#' be substantial particularly if the input matrix only contains a small number
#' of missing data. If there are no missing data, or the missing data are
#' numerous, the speedup will be smaller but still present.
#' 
#' The fast calculations are currently implemented only for
#' \code{method="pearson"} and \code{use} either \code{"all.obs"} or
#' \code{"pairwise.complete.obs"}.  The \code{corFast} function is a wrapper
#' that calls the function \code{cor}. If the combination of \code{method} and
#' \code{use} is implemented by the fast calculations, the fast code is
#' executed; otherwise, R's own correlation \code{\link[stats]{cor}} is
#' executed.
#' 
#' The argument \code{quick} specifies the precision of handling of missing
#' data. Zero will cause all calculations to be executed precisely, which may
#' be significantly slower than calculations without missing data.
#' Progressively higher values will speed up the calculations but introduce
#' progressively larger errors. Without missing data, all column means and
#' variances can be pre-calculated before the covariances are calculated. When
#' missing data are present, exact calculations require the column means and
#' variances to be calculated for each covariance. The approximate calculation
#' uses the pre-calculated mean and variance and simply ignores missing data in
#' the covariance calculation. If the number of missing data is high, the
#' pre-calculated means and variances may be very different from the actual
#' ones, thus potentially introducing large errors.  The \code{quick} value
#' times the number of rows specifies the maximum difference in the number of
#' missing entries for mean and variance calculations on the one hand and
#' covariance on the other hand that will be tolerated before a recalculation
#' is triggered. The hope is that if only a few missing data are treated
#' approximately, the error introduced will be small but the potential speedup
#' can be significant.
#' 
#' @aliases cor1 corFast cor
#' @param x a numeric vector or a matrix. If \code{y} is null, \code{x} must be
#' a matrix.
#' @param y a numeric vector or a matrix. If not given, correlations of columns
#' of \code{x} will be calculated.
#' @param use a character string specifying the handling of missing data. The
#' fast calculations currently support \code{"all.obs"} and
#' \code{"pairwise.complete.obs"}; for other options, see R's standard
#' correlation function \code{\link[stats]{cor}}.  Abbreviations are allowed.
#' @param method a character string specifying the method to be used. Fast
#' calculations are currently available only for \code{"pearson"}.
#' @param quick real number between 0 and 1 that controls the precision of
#' handling of missing data in the calculation of correlations. See details.
#' @param cosine logical: calculate cosine correlation? Only valid for
#' \code{method="pearson"}. Cosine correlation is similar to Pearson
#' correlation but the mean subtraction is not performed. The result is the
#' cosine of the angle(s) between (the columns of) \code{x} and \code{y}.
#' @param cosineX logical: use the cosine calculation for \code{x}? This
#' setting does not affect \code{y} and can be used to give a hybrid
#' cosine-standard correlation.
#' @param cosineY logical: use the cosine calculation for \code{y}? This
#' setting does not affect \code{x} and can be used to give a hybrid
#' cosine-standard correlation.
#' @param drop logical: should the result be turned into a vector if it is
#' effectively one-dimensional?
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose Controls the level of verbosity. Values above zero will cause
#' a small amount of diagnostic messages to be printed.
#' @param indent Indentation of printed diagnostic messages. Each unit above
#' zero adds two spaces.
#' @return The matrix of the Pearson correlations of the columns of \code{x}
#' with columns of \code{y} if \code{y} is given, and the correlations of the
#' columns of \code{x} if \code{y} is not given.
#' @note The implementation uses the BLAS library matrix multiplication
#' function for the most expensive step of the calculation. Using a tuned,
#' architecture-specific BLAS may significantly improve the performance of this
#' function.
#' 
#' The values returned by the corFast function may differ from the values
#' returned by R's function \code{\link[stats]{cor}} by rounding errors on the
#' order of 1e-15.
#' @author Peter Langfelder
#' @seealso R's standard Pearson correlation function \code{\link{cor}}.
#' @references Peter Langfelder, Steve Horvath (2012) Fast R Functions for
#' Robust Correlations and Hierarchical Clustering.  Journal of Statistical
#' Software, 46(11), 1-17.  \url{http://www.jstatsoft.org/v46/i11/}
#' @keywords misc
#' @examples
#' 
#' 
#' ## Test the speedup compared to standard function cor
#' 
#' # Generate a random matrix with 200 rows and 1000 columns
#' 
#' set.seed(10)
#' nrow = 100;
#' ncol = 500;
#' data = matrix(rnorm(nrow*ncol), nrow, ncol);
#' 
#' ## First test: no missing data
#' 
#' system.time( {corStd = stats::cor(data)} );
#' 
#' system.time( {corFast = cor(data)} );
#' 
#' all.equal(corStd, corFast)
#' 
#' # Here R's standard correlation performs very well.
#' 
#' # We now add a few missing entries.
#' 
#' data[sample(nrow, 10), 1] = NA;
#' 
#' # And test the correlations again...
#' 
#' system.time( {corStd = stats::cor(data, use ='p')} );
#' 
#' system.time( {corFast = cor(data, use = 'p')} );
#' 
#' all.equal(corStd, corFast)
#' 
#' # Here the R's standard correlation slows down considerably
#' # while corFast still retains it speed. Choosing
#' # higher ncol above will make the difference more pronounced.
#' 
#' 
  cor = cor(x, y, use = use, ...)
  x = as.matrix(x)
  finMat = !is.na(x)
  if (is.null(y))
  {
    np = t(finMat) %*% finMat
  } else {
    y = as.matrix(y)
    np = t(finMat) %*% (!is.na(y))
  }
  Z = 0.5 * log( (1+cor)/(1-cor) ) * sqrt(np-2)
  if (ia=="two.sided")
  {
    T = sqrt(np - 2) * abs(cor)/sqrt(1 - cor^2)
    p = 2*pt(T, np - 2, lower.tail = FALSE)
  } else if (ia=="less")
  {
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = TRUE)
  } else if (ia=="greater")
  {
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = FALSE)
  }

  list(cor = cor, p = p, Z = Z, t = T, nObs = np)
}



#' Calculation of biweight midcorrelations and associated p-values
#' 
#' A faster, one-step calculation of Student correlation p-values for multiple
#' biweight midcorrelations, properly taking into account the actual number of
#' observations.
#' 
#' The function calculates the biweight midcorrelations of a matrix or of two
#' matrices and the corresponding Student p-values. The output is not as
#' full-featured as \code{\link{cor.test}}, but can work with matrices as
#' input.
#' 
#' @param x a vector or a matrix
#' @param y a vector or a matrix. If \code{NULL}, the correlation of columns of
#' \code{x} will be calculated.
#' @param use determines handling of missing data. See \code{\link{bicor}} for
#' details.
#' @param alternative specifies the alternative hypothesis and must be (a
#' unique abbreviation of) one of \code{"two.sided"}, \code{"greater"} or
#' \code{"less"}. the initial letter.  \code{"greater"} corresponds to positive
#' association, \code{"less"} to negative association.
#' @param \dots other arguments to the function \code{\link{bicor}}.
#' @return A list with the following components, each a marix: \item{bicor}{the
#' calculated correlations} \item{p}{the Student p-values corresponding to the
#' calculated correlations} \item{Z}{Fisher transform of the calculated
#' correlations} \item{t}{Student t statistics of the calculated correlations}
#' \item{nObs}{Numbers of observations for the correlation, p-values etc.}
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{bicor}} for calculation of correlations only;
#' 
#' \code{\link{cor.test}} for another function for significance test of
#' correlations
#' @references Peter Langfelder, Steve Horvath (2012) Fast R Functions for
#' Robust Correlations and Hierarchical Clustering.  Journal of Statistical
#' Software, 46(11), 1-17.  \url{http://www.jstatsoft.org/v46/i11/}
#' @keywords stats
#' @examples
#' 
#' # generate random data with non-zero correlation
#' set.seed(1);
#' a = rnorm(100);
#' b = rnorm(100) + a;
#' x = cbind(a, b);
#' # Call the function and display all results
#' bicorAndPvalue(x)
#' # Set some components to NA
#' x[c(1:4), 1] = NA
#' corAndPvalue(x)
#' # Note that changed number of observations.
#' 
bicorAndPvalue = function(x, y = NULL, use = "pairwise.complete.obs", 
                          alternative = c("two.sided", "less", "greater"), 
                          ...)
{
  ia = match.arg(alternative)
  cor = bicor(x, y, use = use, ...)
  x = as.matrix(x)
  finMat = !is.na(x)
  if (is.null(y))
  {
    np = t(finMat) %*% finMat
  } else {
    y = as.matrix(y)
    np = t(finMat) %*% (!is.na(y))
  }
  Z = 0.5 * log( (1+cor)/(1-cor) ) * sqrt(np-2)
  if (ia=="two.sided")
  {
    T = sqrt(np - 2) * abs(cor)/sqrt(1 - cor^2)
    p = 2*pt(T, np - 2, lower.tail = FALSE)
  } else if (ia=="less")
  {
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = TRUE)
  } else if (ia=="greater")
  {
    #Z = 0.5 * log( (1+cor)/(1-cor) ) * sqrt(np-3)
    T = sqrt(np - 2) * cor/sqrt(1 - cor^2)
    p = pt(T, np - 2, lower.tail = FALSE)
  }

  list(bicor = cor, p = p, Z = Z, t = T, nObs = np)
}


