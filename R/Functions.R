# A collection of auxiliar functions

# collectGarbage ####
#' Iterative garbage collection.
#'
#' Performs garbage collection until free memory idicators show no change.
#'
#'
#' @return None.
#' @author Steve Horvath
#' @keywords utilities
collectGarbage <- function() {
    while (gc()[2, 4]  != gc()[2, 4] | gc()[1, 4]  != gc()[1, 4]) {}
}

# alignExpr ####
#' Align expression data with given vector
#'
#' Multiplies genes (columns) in given expression data such that their
#' correlation with given reference vector is non-negative.
#'
#' The function basically multiplies each column in \code{datExpr} by the sign
#' of its correlation with \code{y}. If \code{y} is not given, the first column
#' in \code{datExpr} will be used as the reference vector.
#'
#' @param datExpr expression data to be aligned. A data frame with columns
#' corresponding to genes and rows to samples.
#' @param y reference vector of length equal the number of samples (rows) in
#' \code{datExpr}
#' @return A data frame containing the aligned expression data, of the same
#' dimensions as the input data frame.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
alignExpr <- function(datExpr, y = NULL) {
    if (!is.null(y) & nrow(as.matrix(datExpr)) != length(y)) {
        stop("Incompatible number of samples in 'datExpr' and 'y'.")
    }
    if (is.null(y)){
        y <- as.numeric(datExpr[, 1])
    }
    sign1 <- sign(as.numeric(cor(y, datExpr, use = "p")))
    as.data.frame(scale(t(t(datExpr) * sign1)))
} # end of function alignExpr

# vectorizeMatrix ####
#' Turn a matrix into a vector of non-redundant components
#'
#' A convenient function to turn a matrix into a vector of non-redundant
#' components. If the matrix is non-symmetric, returns a vector containing all
#' entries of the matrix. If the matrix is symmetric, only returns the upper
#' triangle and optionally the diagonal.
#'
#'
#' @param M the matrix or data frame to be vectorized.
#' @param diag logical: should the diagonal be included in the output?
#' @return A vector containing the non-redundant entries of the input matrix.
#' @author Steve Horvath
#' @keywords misc
vectorizeMatrix <- function(M, diag = FALSE) {
    if (is.null(dim(M)) | length(dim(M)) != 2) {
        stop("The input of the vectorize function is not a matrix or data ",
             "frame.")
    }

    # now we check whether the matrix is symmetrical
    if (ncol(M) == nrow(M)) {
        M <- as.matrix(M)
        Mtranspose <- t(M)
        abs.difference <- max(abs(M - Mtranspose), na.rm = TRUE)
        if (abs.difference < 10 ^ (-14)) {
            out <- M[upper.tri(M, diag)]
        }
        else {
            out <- as.vector(M)
        }
    } else {
        out <- as.vector(M)
    }
    out
} # end

# metaZfunction ####
#' Meta-analysis Z statistic
#'
#' The function calculates a meta analysis Z statistic based on an input data
#' frame of Z statistics.
#'
#' For example, if datZ has 3 columns whose columns are labelled Z1,Z2,Z3 then
#' ZMeta= (Z1+Z2+Z3)/sqrt(3). Under the null hypothesis (where all Z statistics
#' follow a standard normal distribution and the Z statistics are independent),
#' ZMeta also follows a standard normal distribution.  To calculate a 2 sided
#' p-value, one an use the following code pvalue=2*pnorm(-abs(ZMeta) )
#'
#' @param datZ Matrix or data frame of Z statistics (assuming standard normal
#' distribution under the null hypothesis). Rows correspond to genes, columns
#' to independent data sets.
#' @param columnweights optional vector of non-negative numbers for weighing
#' the columns of datZ.
#' @return Vector of meta analysis Z statistic. Under the null hypothesis this
#' should follow a standard normal distribution.
#' @author Steve Horvath
#' @keywords misc
metaZfunction <- function(datZ, columnweights = NULL) {
    if (!is.null(columnweights))  {
        datZ <- t(t(datZ) *  columnweights)
    }
    datZpresent <- !is.na(datZ) + 0.0
    if (!is.null(columnweights))  {
        datZpresent <- t(t(datZpresent) *  columnweights)
    }
    sumZ <- as.numeric(rowSums(datZ, na.rm = TRUE))
    variance <- as.numeric(rowSums(datZpresent ^ 2))
    sumZ / sqrt(variance)
}

# prepComma ####
#' Prepend a comma to a non-empty string
#'
#' Utility function that prepends a comma before the input string if the string
#' is non-empty.
#'
#'
#' @param s Character string.
#' @return If \code{s} is non-empty, returns \code{paste(",", s)}, otherwise
#' returns s.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' prepComma("abc")
#' prepComma("")
#'
prepComma <- function(s) {
    if (s == "") {
        return(s)
    }
    paste(", ", s)
}

# prependZeros ####
#' Pad numbers with leading zeros to specified total width
#'
#' This function pads the specified numbers with zeros to a specified total
#' width.
#'
#'
#' @param x Vector of numbers to be padded.
#' @param width Width to pad the numbers to.
#' @return Character vector with the 0-padded numbers.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' prependZeros(1:10)
#' prependZeros(1:10, 4)
#'
prependZeros <- function(x, width = max(nchar(x))) {
    lengths <- nchar(x)
    if (width < max(lengths))
        stop("Some entries of 'x' are too long.")
    out <- as.character(x)
    n <- length(x)
    for (i in 1:n) {
        if (lengths[i] < width) {
            out[i] <- paste0(paste(rep("0", width - lengths[i]), collapse = ""),
                            x[i])
        }
    }
    out
}
