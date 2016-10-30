
# corPvalueFisher ####
#' Fisher's asymptotic p-value for correlation
#'
#' Calculates Fisher's asymptotic p-value for given correlations.
#'
#'
#' @param cor A vector of correlation values whose corresponding p-values are
#' to be calculated
#' @param nSamples Number of samples from which the correlations were
#' calculated
#' @param twoSided logical: should the calculated p-values be two sided?
#' @return A vector of p-values of the same length as the input correlations.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
corPvalueFisher <- function(cor, nSamples, twoSided = TRUE) {
    if (sum(abs(cor) > 1, na.rm = TRUE) > 0)
        stop("Some entries in 'cor' are out of normal range  - 1 to 1.")
    if (twoSided) {
        z = abs(0.5 * log((1 + cor) / (1 - cor)) * sqrt(nSamples - 3))
        2 * pnorm(-z)
    } else {
        # return a small p-value for positive correlations
        z = -0.5 * log((1 + cor) / (1 - cor)) * sqrt(nSamples - 3)
        pnorm(-z)
    }
}

# corPvalueStudent ####
#' Student asymptotic p-value for correlation
#'
#' Calculates Student asymptotic p-value for given correlations.
#'
#'
#' @param cor A vector of correlation values whose corresponding p-values are
#' to be calculated
#' @param nSamples Number of samples from which the correlations were
#' calculated
#' @return A vector of p-values of the same length as the input correlations.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
corPvalueStudent <- function(cor, nSamples) {
    T = sqrt(nSamples - 2) * cor / sqrt(1 - cor ^ 2)
    2 * pt(abs(T), nSamples - 2, lower.tail = FALSE)
}


# rankPvalue ####
#' Estimate the p-value for ranking consistently high (or low) on multiple
#' lists
#'
#' The function rankPvalue calculates the p-value for observing that an object
#' (corresponding to a row of the input data frame \code{datS}) has a
#' consistently high ranking (or low ranking) according to multiple ordinal
#' scores (corresponding to the columns of the input data frame \code{datS}).
#'
#' The function calculates asymptotic p-values (and optionally q-values) for
#' testing the null hypothesis that the values in the columns of datS are
#' independent. This allows us to find objects (rows) with consistently high
#' (or low) values across the columns.
#'
#' Example: Imagine you have 5 vectors of Z statistics corresponding to the
#' columns of datS. Further assume that a gene has ranks 1,1,1,1,20 in the 5
#' lists. It seems very significant that the gene ranks number 1 in 4 out of
#' the 5 lists. The function rankPvalue can be used to calculate a p-value for
#' this occurrence.
#'
#' The function uses the central limit theorem to calculate asymptotic p-values
#' for two types of test statistics that measure consistently high or low
#' ordinal values. The first method (referred to as percentile rank method)
#' leads to accurate estimates of p-values if datS has at least 4 columns but
#' it can be overly conservative.  The percentile rank method replaces each
#' column datS by the ranked version rank(datS[,i]) (referred to ask low
#' ranking) and by rank(-datS[,i]) (referred to as high ranking). Low ranking
#' and high ranking allow one to find consistently small values or consistently
#' large values of datS, respectively.  All ranks are divided by the maximum
#' rank so that the result lies in the unit interval [0,1]. In the following,
#' we refer to rank/max(rank) as percentile rank. For a given object
#' (corresponding to a row of datS) the observed percentile rank follows
#' approximately a uniform distribution under the null hypothesis. The test
#' statistic is defined as the sum of the percentile ranks (across the columns
#' of datS). Under the null hypothesis that there is no relationship between
#' the rankings of the columns of datS, this (row sum) test statistic follows a
#' distribution that is given by the convolution of random uniform
#' distributions. Under the null hypothesis, the individual percentile ranks
#' are independent and one can invoke the central limit theorem to argue that
#' the row sum test statistic follows asymptotically a normal distribution.  It
#' is well-known that the speed of convergence to the normal distribution is
#' extremely fast in case of identically distributed uniform distributions.
#' Even when datS has only 4 columns, the difference between the normal
#' approximation and the exact distribution is negligible in practice (Killmann
#' et al 2001). In summary, we use the central limit theorem to argue that the
#' sum of the percentile ranks follows a normal distribution whose mean and
#' variance can be calculated using the fact that the mean value of a uniform
#' random variable (on the unit interval) equals 0.5 and its variance equals
#' 1/12.
#'
#' The second method for calculating p-values is referred to as scale method.
#' It is often more powerful but its asymptotic p-value can only be trusted if
#' either datS has a lot of columns or if the ordinal scores (columns of datS)
#' follow an approximate normal distribution.  The scale method scales (or
#' standardizes) each ordinal variable (column of datS) so that it has mean 0
#' and variance 1. Under the null hypothesis of independence, the row sum
#' follows approximately a normal distribution if the assumptions of the
#' central limit theorem are met. In practice, we find that the second approach
#' is often more powerful but it makes more distributional assumptions (if datS
#' has few columns).
#'
#' @param datS a data frame whose rows represent objects that will be ranked.
#' Each column of \code{datS} represents an ordinal variable (which can take on
#' negative values). The columns correspond to (possibly signed) object
#' significance measures, e.g., statistics (such as Z statistics), ranks, or
#' correlations.
#' @param columnweights allows the user to input a vector of non-negative
#' numbers reflecting weights for the different columns of \code{datZ}. If it
#' is set to \code{NULL} then all weights are equal.
#' @param na.last controls the treatment of missing values (NAs) in the rank
#' function. If \code{TRUE}, missing values in the data are put last (i.e. they
#' get the highest rank values). If \code{FALSE}, they are put first; if
#' \code{NA}, they are removed; if \code{"keep"} they are kept with rank NA.
#' See \code{\link{rank}} for more details.
#' @param ties.method represents the ties method used in the rank function for
#' the percentile rank method. See \code{\link{rank}} for more details.
#' @param calculateQvalue logical: should q-values be calculated? If set to
#' TRUE then the function calculates corresponding q-values (local false
#' discovery rates) using the qvalue package, see Storey JD and Tibshirani R.
#' (2003). This option assumes that qvalue package has been installed.
#' @param pValueMethod determines which method is used for calculating
#' p-values. By default it is set to "all", i.e. both methods are used. If it
#' is set to "rank" then only the percentile rank method is used. If it set to
#' "scale" then only the scale method will be used.
#' @return
#'
#' A list whose actual content depends on which p-value methods is selected,
#' and whether q0values are calculated. The following inner components are
#' calculated, organized in outer components \code{datoutrank} and
#' \code{datoutscale},:
#'
#' \item{pValueExtremeRank}{This is the minimum between pValueLowRank and
#' pValueHighRank, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowRank}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the rank method.}
#'
#' \item{pValueHighRank}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the rank method.}
#'
#' \item{pValueExtremeScale}{This is the minimum between pValueLowScale and
#' pValueHighScale, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowScale}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the Scale method.}
#'
#' \item{pValueHighScale}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the Scale method.}
#'
#' \item{qValueExtremeRank}{local false discovery rate (q-value) corresponding
#' to the p-value pValueExtremeRank}
#'
#' \item{qValueLowRank}{local false discovery rate (q-value) corresponding to
#' the p-value pValueLowRank}
#'
#' \item{qValueHighRank}{local false discovery rate (q-value) corresponding to
#' the p-value pValueHighRank}
#'
#' \item{qValueExtremeScale}{local false discovery rate (q-value) corresponding
#' to the p-value pValueExtremeScale}
#'
#' \item{qValueLowScale}{local false discovery rate (q-value) corresponding to
#' the p-value pValueLowScale}
#'
#' \item{qValueHighScale}{local false discovery rate (q-value) corresponding to
#' the p-value pValueHighScale}
#' @author Steve Horvath
#' @seealso \code{\link{rank}}, \code{\link{qvalue}}
#' @references
#' Killmann F, VonCollani E (2001) A Note on the Convolution of the
#' Uniform and Related Distributions and Their Use in Quality Control. Economic
#' Quality Control Vol 16 (2001), No. 1, 17-41.ISSN 0940-5151
#'
#' Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide
#' experiments. Proceedings of the National Academy of Sciences, 100:
#' 9440-9445.
#' @keywords misc
rankPvalue <- function(datS,
                       columnweights = NULL,
                       na.last = "keep",
                       ties.method = "average",
                       calculateQvalue = TRUE,
                       pValueMethod = "all") {
    no.rows = dim(datS)[[1]]
    no.cols = dim(datS)[[2]]
    if (!is.null(columnweights) & no.cols  != length(columnweights))
        stop(
            "The number of components of the vector columnweights is unequal",
            "to the number of columns of datS. Hint: consider transposing ",
            "datS."
        )

    if (!is.null(columnweights)) {
        if (min(columnweights, na.rm = TRUE) < 0)
            stop(
                "At least one component of columnweights is negative, which ",
                "makes no sense. The entries should be positive numbers"
            )
        if (sum(is.na(columnweights)) > 0)
            stop(
                "At least one component of columnweights is missing, which ",
                "makes no sense. The entries should be positive numbers"
            )
        if (sum(columnweights) != 1) {
            # warning("The entries of columnweights do not sum to 1.
            # Therefore, they will divided by the sum. Then the resulting
            # weights sum to 1.")
            columnweights = columnweights / sum(columnweights)
        }
    }

    if (pValueMethod  != "scale") {
        percentilerank1 <- function(x) {
            R1 = rank(x, ties.method = ties.method, na.last = na.last)
            (R1 - .5) / max(R1, na.rm = TRUE)
        }

        datrankslow = apply(datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankslow = t(t(datrankslow) * columnweights)
        }
        datSpresent = !is.na(datS) + 0
        if (!is.null(columnweights)) {
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expectedsum = rowSums(datSpresent, na.rm = TRUE)  *
            0.5
        varsum = rowSums(datSpresent ^ 2, na.rm = TRUE) * 1 / 12
        observed.sumPercentileslow = as.numeric(rowSums(datrankslow,
                                                        na.rm = TRUE))
        Zstatisticlow = (observed.sumPercentileslow - expectedsum) / sqrt(varsum)
        datrankshigh = apply(-datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankshigh = t(t(datrankshigh) * columnweights)
        }
        observed.sumPercentileshigh = as.numeric(rowSums(datrankshigh,
                                                         na.rm = TRUE))
        Zstatistichigh = (observed.sumPercentileshigh - expectedsum) / sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = pmin(pValueLow, pValueHigh)
        datoutrank = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme = pmin(qValueLow, qValueHigh)
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutrank = data.frame(datoutrank, datq)
            names(datoutrank) = paste0(names(datoutrank), "Rank")
        }
    }
    if (pValueMethod  != "rank") {
        datSpresent = !is.na(datS) + 0
        scaled.datS = scale(datS)
        if (!is.null(columnweights)) {
            scaled.datS = t(t(scaled.datS) * columnweights)
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expected.value = rep(0, no.rows)
        varsum = rowSums(datSpresent ^ 2) * 1
        observed.sumScaleddatS = as.numeric(rowSums(scaled.datS, na.rm = TRUE))
        Zstatisticlow = (observed.sumScaleddatS - expected.value) / sqrt(varsum)
        scaled.minusdatS = scale(-datS)
        if (!is.null(columnweights)) {
            scaled.minusdatS = t(t(scaled.minusdatS) * columnweights)
        }
        observed.sumScaledminusdatS = as.numeric(rowSums(scaled.minusdatS,
                                                         na.rm = TRUE))
        Zstatistichigh = (observed.sumScaledminusdatS - expected.value) /
            sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = 2 * pnorm(-abs(Zstatisticlow))
        datoutscale = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme[rest1] = qvalue(pValueExtreme[rest1])$qvalues
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutscale = data.frame(datoutscale, datq)
        }
        names(datoutscale) = paste0(names(datoutscale), "Scale")
    }
    if (pValueMethod == "rank") {
        datout = datoutrank
    }
    if (pValueMethod == "scale") {
        datout = datoutscale
    }
    if (pValueMethod  != "rank" & pValueMethod  != "scale")
        datout = data.frame(datoutrank, datoutscale)
    datout
} # End of function

