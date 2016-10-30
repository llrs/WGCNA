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
    while (gc()[2, 4]  != gc()[2, 4] | gc()[1, 4]  != gc()[1, 4]) {
    }
}

# clusterCoef ####
#' Clustering coefficient calculation
#'
#' This function calculates the clustering coefficients for all nodes in the
#' network given by the input adjacency matrix.
#'
#'
#' @param adjMat adjacency matrix
#' @return A vector of clustering coefficients for each node.
#' @author Steve Horvath
#' @keywords misc
clusterCoef <- function(adjMat) {
    checkAdjMat(adjMat)
    diag(adjMat) = 0
    nNodes = dim(adjMat)[[1]]
    computeLinksInNeighbors <- function(x, imatrix) {
        x %*% imatrix %*% x
    }
    nolinksNeighbors <- c(rep(-666, nNodes))
    total.edge <- c(rep(-666, nNodes))
    maxh1 = max(as.dist(adjMat))
    minh1 = min(as.dist(adjMat))
    if (maxh1 > 1 | minh1 < 0) {
        stop("The adjacency matrix contains entries that are larger than 1 or
                smaller than 0: max  = ", maxh1, ", min  = ", minh1)
    }
    nolinksNeighbors <- apply(adjMat, 1, computeLinksInNeighbors,
                              imatrix = adjMat)
    plainsum  <- apply(adjMat, 1, sum)
    squaresum <- apply(adjMat ^ 2, 1, sum)
    total.edge = plainsum ^ 2 - squaresum
    CChelp = rep(-666, nNodes)
    CChelp = ifelse(total.edge == 0, 0, nolinksNeighbors / total.edge)
    CChelp
} # end of function

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
    if (!is.null(y) & dim(as.matrix(datExpr))[[1]]  != length(y))
        stop("Incompatible number of samples in 'datExpr' and 'y'.")
    if (is.null(y))
        y <- as.numeric(datExpr[, 1])
    sign1 <- sign(as.numeric(cor(y, datExpr, use = "p")))
    as.data.frame(scale(t(t(datExpr) * sign1)))
} # end of function alignExpr


# this function can be used to rank the values in x. Ties are broken by the
# method first.
# This function does not appear to be used anywhere in these functions.
# rank1 <- function(x){
#    rank(x, ties.method = "first")
#}

# randIndex ####
#  Rand index calculation
# this function is used for computing the Rand index below...
#
.choosenew <- function(n, k) {
    n <- c(n)
    out1 <- rep(0, length(n))
    for (i in c(1:length(n))) {
        if (n[i] < k) {
            out1[i] <- 0
        }
        else {
            out1[i] <- choose(n[i], k)
        }
    }
    out1
}
#' Rand index of two partitions
#'
#' Computes the Rand index, a measure of the similarity between two
#' clusterings.
#'
#'
#' @param tab a matrix giving the cross-tabulation table of two clusterings.
#' @param adjust logical: should the "adjusted" version be computed?
#' @return the Rand index of the input table.
#' @author Steve Horvath
#' @references W. M. Rand (1971). "Objective criteria for the evaluation of
#' clustering methods". Journal of the American Statistical Association 66:
#' 846-850
#' @keywords misc
randIndex <- function(tab, adjust = TRUE) {
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    nn <- 0
    m <- nrow(tab)
    n <- ncol(tab)
    for (i in 1:m) {
        c <- 0
        for (j in 1:n) {
            a <- a + .choosenew(tab[i, j], 2)
            nj <- sum(tab[, j])
            c <- c + .choosenew(nj, 2)
        }
        ni <- sum(tab[i,])
        b <- b + .choosenew(ni, 2)
        nn <- nn + ni
    }
    if (adjust) {
        d <- .choosenew(nn, 2)
        adrand <- (a - (b * c) / d) / (0.5 * (b + c) - (b * c) / d)
        adrand
    } else {
        b <- b - a
        c <- c - a
        d <- .choosenew(nn, 2) - a - b - c
        rand <- (a + d) / (a + b + c + d)
        rand
    }
}

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
    if (is.null(dim(M)))
        stop("The input of the vectorize function is not a matrix or data frame.")
    if (length(dim(M)) != 2)
        stop("The input of the vectorize function is not a matrix or data frame.")
    # now we check whether the matrix is symmetrical
    if (dim(M)[[1]] == dim(M)[[2]]) {
        M = as.matrix(M)
        Mtranspose = t(M)
        abs.difference = max(abs(M - Mtranspose), na.rm = TRUE)
        if (abs.difference < 10 ^ (-14)) {
            out = M[upper.tri(M, diag)]
        }
        else
            out = as.vector(M)
    } else
        out = as.vector(M)
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
        datZ = t(t(datZ) *  columnweights)
    }
    datZpresent = !is.na(datZ) + 0.0
    if (!is.null(columnweights))  {
        datZpresent = t(t(datZpresent) *  columnweights)
    }
    sumZ = as.numeric(rowSums(datZ, na.rm = TRUE))
    variance = as.numeric(rowSums(datZpresent ^ 2))
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
    if (s == "")
        return (s)
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
    lengths = nchar(x)
    if (width < max(lengths))
        stop("Some entries of 'x' are too long.")
    out = as.character(x)
    n = length(x)
    for (i in 1:n) {
        if (lengths[i] < width) {
            out[i] <- paste0(paste(rep("0", width - lengths[i]), collapse = ""),
                            x[i])
        }
    }
    out
}
