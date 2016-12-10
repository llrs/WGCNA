# slight re-definition of the bicor function
# bicor ####
#' Biweight Midcorrelation
#'
#' Calculate biweight midcorrelation efficiently for matrices.
#'
#' This function implements biweight midcorrelation calculation (see
#' references). If \code{y} is not supplied, midcorrelation of columns of
#' \code{x} will be calculated; otherwise, the midcorrelation between columns
#' of \code{x} and \code{y} will be calculated. Thus, \code{bicor(x)} is
#' equivalent to \code{bicor(x,x)} but is more efficient.
#'
#' The options \code{robustX}, \code{robustY} allow the user to revert the
#' calculation to standard correlation calculation. This is important, for
#' example, if any of the variables is binary (or, more generally, discrete) as
#' in such cases the robust methods produce meaningless results.  If both
#' \code{robustX}, \code{robustY} are set to \code{FALSE}, the function
#' calculates the standard Pearson correlation (but is slower than the function
#' \code{\link{cor}}).
#'
#' The argument \code{quick} specifies the precision of handling of missing
#' data in the correlation calculations. Value \code{quick = 0} will cause all
#' calculations to be executed accurately, which may be significantly slower
#' than calculations without missing data. Progressively higher values will
#' speed up the calculations but introduce progressively larger errors. Without
#' missing data, all column meadians and median absolute deviations (MADs) can
#' be pre-calculated before the covariances are calculated. When missing data
#' are present, exact calculations require the column medians and MADs to be
#' calculated for each covariance. The approximate calculation uses the
#' pre-calculated median and MAD and simply ignores missing data in the
#' covariance calculation. If the number of missing data is high, the
#' pre-calculated medians and MADs may be very different from the actual ones,
#' thus potentially introducing large errors. The \code{quick} value times the
#' number of rows specifies the maximum difference in the number of missing
#' entries for median and MAD calculations on the one hand and covariance on
#' the other hand that will be tolerated before a recalculation is triggered.
#' The hope is that if only a few missing data are treated approximately, the
#' error introduced will be small but the potential speedup can be significant.
#'
#' The choice \code{"all"} for \code{pearsonFallback} is not fully implemented
#' in the sense that there are rare but possible cases in which the calculation
#' is equivalent to \code{"individual"}. This may happen if the \code{use}
#' option is set to \code{"pairwise.complete.obs"} and the missing data are
#' arranged such that each individual mad is non-zero, but when two columns are
#' analyzed together, the missing data from both columns may make a mad zero.
#' In such a case, the calculation is treated as Pearson, but other columns
#' will be treated as bicor.
#'
#' @param x a vector or matrix-like numeric object
#' @param y a vector or matrix-like numeric object
#' @param robustX use robust calculation for \code{x}?
#' @param robustY use robust calculation for \code{y}?
#' @param use specifies handling of \code{NA}s. One of (unique abbreviations
#' of) "all.obs", "pairwise.complete.obs".
#' @param maxPOutliers specifies the maximum percentile of data that can be
#' considered outliers on either side of the median separately. For each side
#' of the median, if higher percentile than \code{maxPOutliers} is considered
#' an outlier by the weight function based on \code{9*mad(x)}, the width of the
#' weight function is increased such that the percentile of outliers on that
#' side of the median equals \code{maxPOutliers}. Using \code{maxPOutliers=1}
#' will effectively disable all weight function broadening; using
#' \code{maxPOutliers=0} will give results that are quite similar (but not
#' equal to) Pearson correlation.
#' @param quick real number between 0 and 1 that controls the handling of
#' missing data in the calculation of correlations. See details.
#' @param pearsonFallback Specifies whether the bicor calculation should revert
#' to Pearson when median absolute deviation (mad) is zero. Recongnized values
#' are (abbreviations of) \code{"none", "individual", "all"}. If set to
#' \code{"none"}, zero mad will result in \code{NA} for the corresponding
#' correlation.  If set to \code{"individual"}, Pearson calculation will be
#' used only for columns that have zero mad. If set to \code{"all"}, the
#' presence of a single zero mad will cause the whole variable to be treated in
#' Pearson correlation manner (as if the corresponding \code{robust} option was
#' set to \code{FALSE}).
#' @param cosine logical: calculate cosine biweight midcorrelation?  Cosine
#' bicorrelation is similar to standard bicorrelation but the median
#' subtraction is not performed.
#' @param cosineX logical: use the cosine calculation for \code{x}? This
#' setting does not affect \code{y} and can be used to give a hybrid
#' cosine-standard bicorrelation.
#' @param cosineY logical: use the cosine calculation for \code{y}? This
#' setting does not affect \code{x} and can be used to give a hybrid
#' cosine-standard bicorrelation.
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose if non-zero, the underlying C function will print some
#' diagnostics.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of biweight midcorrelations. Dimnames on the result are set
#' appropriately.
#' @author Peter Langfelder
#' @references
#'
#' Peter Langfelder, Steve Horvath (2012) Fast R Functions for Robust
#' Correlations and Hierarchical Clustering. Journal of Statistical Software,
#' 46(11), 1-17. \url{http://www.jstatsoft.org/v46/i11/}
#'
#' "Dealing with Outliers in Bivariate Data: Robust Correlation", Rich
#' Herrington, http://www.unt.edu/benchmarks/archives/2001/december01/rss.htm
#'
#' "Introduction to Robust Estimation and Hypothesis Testing", Rand Wilcox,
#' Academic Press, 1997.
#'
#' "Data Analysis and Regression: A Second Course in Statistics", Mosteller and
#' Tukey, Addison-Wesley, 1977, pp. 203-209.
#' @keywords robust
bicor <- function(x, y = NULL, robustX = TRUE, robustY = TRUE, use = 'all.obs',
                  maxPOutliers = 1, quick = 0, pearsonFallback = "individual",
                  cosine = FALSE, cosineX = cosine, cosineY = cosine,
                  nThreads = 0, verbose = 0, indent = 0) {
    Cerrors <- c("Memory allocation error")
    nKnownErrors <- length(Cerrors)
    na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
    if (is.na(na.method)) {
        stop("Unrecognized parameter 'use'. Recognized values are \n",
             "'all.obs', 'pairwise.complete.obs'")
    }
    if (na.method == 1) {
        if (sum(is.na(x)) > 0) {
            stop("Missing values present in input variable 'x'. Consider ",
                 "using use = 'pairwise.complete.obs'.")
        }
        if (!is.null(y)) {
            if (sum(is.na(y)) > 0) {
                stop("Missing values present in input variable 'y'. Consider ",
                     "using use = 'pairwise.complete.obs'.")
            }
        }
    }

    fallback = pmatch(pearsonFallback, .pearsonFallbacks)
    if (is.na(na.method)) {
        .pearsonFallbacks <- c("none", "individual", "all")
        stop("Unrecognized 'pearsonFallback'. Recognized values are (unique ",
             "abbreviations of)\n",
             paste(.pearsonFallbacks, collapse = ", "))
    }
    if (quick < 0) {
        stop("quick must be non-negative.")
    }
    if (nThreads < 0) {
        stop("nThreads must be non-negative.")
    }
    if (is.null(nThreads) || (nThreads==0)) {
        nThreads <- .useNThreads()
    }
    x <- as.matrix(x)
    if (prod(dim(x)) == 0) {
        stop("'x' has a zero dimension.")
    }
    nNA <- 0
    err <- 0
    warnX <- 0
    warnY <- 0
    if (is.null(y)) {
        if (!robustX) {
            bi <- cor(x, use = use)
        } else {
            bi <- matrix(0, ncol(x), ncol(x))
            res <- .C("bicor1Fast", x = as.double(x), nrow = as.integer(nrow(x)), ncol = as.integer(ncol(x)),
                     maxPOutliers = as.double(maxPOutliers),
                     quick = as.double(quick),
                     fallback = as.integer(fallback),
                     cosine = as.integer(cosineX),
                     res = as.double(bi), nNA = as.integer(nNA),
                     err = as.integer(err),
                     warn = as.integer(warnX), nThreads = as.integer(nThreads),
                     verbose = as.integer(verbose), indent = as.integer(indent),
                     NAOK = TRUE, PACKAGE = "WGCNA")
        }
        dim(res$res) <- dim(bi)
        if (!is.null(dimnames(x)[[2]])) {
            dimnames(res$res) <- list(dimnames(x)[[2]],  dimnames(x)[[2]] )
        }
        .zeroMADWarnings <- c("Some results will be NA.",
                              paste("Pearson correlation was used for",
                                    "individual columns with zero (or missing)",
                                    "MAD."),
                              "Pearson correlation was used for entire variable.")
        if (res$warn > 0) {
            # For now have only one warning
            warning("bicor: zero MAD in variable 'x'. ",
                    .zeroMADWarnings[fallback])
        }
    } else {
        y <- as.matrix(y)
        if (prod(dim(y)) == 0) {
            stop("'y' has a zero dimension.")
        }
        if (nrow(x) != nrow(y)) {
            stop("'x' and 'y' have incompatible dimensions (unequal numbers",
                 " of rows).")
        }
        bi <- matrix(0, ncol(x), ncol(y))
        res <- .C("bicorFast", x = as.double(x), nrow = as.integer(nrow(x)),
                 ncolx = as.integer(ncol(x)), y = as.double(y),
                 ncoly = as.integer(ncol(y)),
                 robustX = as.integer(robustX), robustY = as.integer(robustY),
                 maxPOutliers = as.double(maxPOutliers),
                 quick = as.double(quick),
                 fallback = as.integer(fallback),
                 cosineX = as.integer(cosineX),
                 cosineY = as.integer(cosineY),
                 res = as.double(bi), nNA = as.integer(nNA),
                 err = as.integer(err),
                 warnX = as.integer(warnX),
                 warnY = as.integer(warnY),
                 nThreads = as.integer(nThreads),
                 verbose = as.integer(verbose), indent = as.integer(indent),
                 NAOK = TRUE,
                 PACKAGE = "WGCNA")
        dim(res$res) <- dim(bi)
        if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]])) {
            dimnames(res$res) = list(dimnames(x)[[2]], dimnames(y)[[2]])
        }
        if (res$warnX > 0) {
            warning("bicor: zero MAD in variable 'x'.",
                    .zeroMADWarnings[fallback])
        }
        if (res$warnY > 0){
            warning("bicor: zero MAD in variable 'y'.",
                    .zeroMADWarnings[fallback])
        }
    }
    if (res$err > 0) {
        stop("An error occurred in compiled code. Error code is ", err)
    }
    if (res$nNA > 0) {
        warning("Missing values generated in calculation of bicor.\nLikely ",
                "cause: too many missing entries, zero median absolute ",
                "deviation, or zero variance.")
    }
    res$res
}

# cor ####
# Code to call my implementation of correlation
# For less than 100 correlations, use stats::cor since that is usually faster, particularly when no missing
# data are present, likely due to the complicated threading I do in the WGCNA correlations.
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
#' ## Test the speedup compared to standard function cor
#' # Generate a random matrix with 200 rows and 1000 columns
#'
#' set.seed(10)
#' nrow = 100
#' ncol = 500
#' data = matrix(rnorm(nrow*ncol), nrow, ncol)
#'
#' ## First test: no missing data
#' system.time( {corStd = stats::cor(data)} )
#' system.time( {corFast = cor(data)} )
#'
#' all.equal(corStd, corFast)
#'
#' # Here R's standard correlation performs very well.
#' # We now add a few missing entries.
#'
#' data[sample(nrow, 10), 1] <- NA
#'
#' # And test the correlations again...
#'
#' system.time( {corStd = stats::cor(data, use ='p')} )
#' system.time( {corFast = cor(data, use = 'p')} )
#' all.equal(corStd, corFast)
#'
#' # Here the R's standard correlation slows down considerably
#' # while corFast still retains it speed. Choosing
#' # higher ncol above will make the difference more pronounced.
#' @export
cor <- function(x, y = NULL, use = "all.obs",
                method = c("pearson", "kendall", "spearman"),
                quick = 0, cosine = FALSE, cosineX = cosine, cosineY = cosine,
                drop = FALSE, nThreads = 0, verbose = 0, indent = 0) {
    na.method <- pmatch(use, c("all.obs", "complete.obs",
                               "pairwise.complete.obs", "everything",
                               "na.or.complete"), nomatch = 0)
    method <- match.arg(method)

    x <- as.matrix(x)
    nx <- ncol(x)
    if (!is.null(y)) {
        y <- as.matrix(y)
        ny <- ncol(y)
    } else {
        ny <- nx
    }

    if ((method == "pearson") && ( (na.method == 1) || (na.method == 3) )) {
        Cerrors <- c("Memory allocation error")
        nKnownErrors <- length(Cerrors)
        na.method <- pmatch(use, c("all.obs", "pairwise.complete.obs"))
        if (is.na(na.method))
            stop("Unrecognized parameter 'use'. Recognized values are \n",
                       "'all.obs', 'pairwise.complete.obs'")
        if (na.method == 1) {
            if (sum(is.na(x)) > 0) {
                stop("Missing values present in input variable 'x'. Consider ",
                     "using use = 'pairwise.complete.obs'.")
            }
            if (!is.null(y)) {
                if (sum(is.na(y)) > 0) {
                    stop("Missing values present in input variable 'y'. ",
                         "Consider using use = 'pairwise.complete.obs'.")}
            }
        }

        if (quick < 0) {
            stop("quick must be non-negative.")
        }
        if (nThreads < 0) {
            stop("nThreads must be non-negative.")
        }
        if (is.null(nThreads) || (nThreads == 0)) {
            nThreads = .useNThreads()
        }

        if (prod(dim(x)) == 0 ) {
            stop("'x' has a zero dimension.")
        }
        nNA <- as.integer(0)
        err <- as.integer(0)
        cosine <- as.integer(cosine)
        nThreads <- as.integer(nThreads)
        verbose <- as.integer(verbose)
        indent <- as.integer(indent)
        if (is.null(y)) {
            res <- .Call("cor1Fast_call", x,
                         quick, cosine,
                         nNA, err, nThreads,
                         verbose, indent, package = "WGCNA")
            if (!is.null(dimnames(x)[[2]])) {
                dimnames(res) <- list(dimnames(x)[[2]],  dimnames(x)[[2]] )
            }
        } else {
            y <- as.matrix(y)
            if (prod(dim(y)) == 0) {
                stop("'y' has a zero dimension.")
            }
            if (nrow(x) != nrow(y)) {
                stop("'x' and 'y' have incompatible dimensions (unequal ",
                     "numbers of rows).")
            }
            bi <- matrix(0, ncol(x), ncol(y))
            res <- .C("corFast", x = as.double(x), nrow = as.integer(nrow(x)),
                     ncolx = as.integer(ncol(x)),
                     y = as.double(y), ncoly = as.integer(ncol(y)),
                     quick = as.double(quick),
                     cosineX = as.integer(cosineX),
                     cosineY = as.integer(cosineY),
                     res = as.double(bi), nNA = as.integer(nNA),
                     err = as.integer(err),
                     nThreads = as.integer(nThreads),
                     verbose = as.integer(verbose), indent = as.integer(indent), NAOK = TRUE,
                     PACKAGE = "WGCNA")
            res <- res$res
            dim(res) <- dim(bi)
            if (!is.null(dimnames(x)[[2]]) || !is.null(dimnames(y)[[2]])) {
                dimnames(res) = list(dimnames(x)[[2]], dimnames(y)[[2]])
            }
        }
        if (err > 0) {
            if (err > nKnownErrors) {
                stop("An error occurred in compiled code. Error code is ", err)
            } else {
                stop(Cerrors[err], " occurred in compiled code. ")
            }
        }
        if (nNA > 0) {
            warning("Missing values generated in calculation of cor. ",
                    "Likely cause: too many missing entries or zero variance.")
        }
        if (drop) {
            res[, , drop = TRUE]
        } else {
            res
        }
    } else {
        stats::cor(x,y, use, method)
    }
}

# Wrappers ####
# Wrappers for compatibility with older scripts
#' @rdname cor
#' @export
cor1 <- function(x, use = "all.obs", verbose = 0, indent = 0) {
    cor(x, use = use, verbose = verbose, indent = indent)
}
#' @rdname cor
#' @export
corFast <- function(x, y = NULL, use = "all.obs",
                   quick = 0, nThreads = 0, verbose = 0, indent = 0) {
    cor(x,y, use, method = "pearson", quick, nThreads, verbose, indent)
}
