
# normalizeLabels ####
#' Transform numerical labels into normal order.
#'
#' Transforms numerical labels into normal order, that is the largest group
#' will be labeled 1, next largest 2 etc. Label 0 is optionally preserved.
#'
#'
#' @param labels Numerical labels.
#' @param keepZero If \code{TRUE} (the default), labels 0 are preserved.
#' @return A vector of the same length as input, containing the normalized
#' labels.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
normalizeLabels <- function(labels, keepZero = TRUE) {
    if (keepZero) {
        NonZero = (labels != 0)
    } else {
        NonZero = rep(TRUE, length(labels))
    }
    f = as.numeric(factor(labels[NonZero]))
    t = table(labels[NonZero])
    # print(t)
    r = rank(-as.vector(t), ties.method = "first")
    norm_labs = rep(0, times = length(labels))
    norm_labs[NonZero] = r[f]
    norm_labs
}

# moduleNumber ####
#' Fixed-height cut of a dendrogram.
#'
#' Detects branches of on the input dendrogram by performing a fixed-height
#' cut.
#'
#' All contiguous branches below the height \code{cutHeight} that contain at
#' least \code{minSize} objects are assigned unique positive numerical labels
#' all unassigned objects are assigned label 0.
#'
#' @param dendro a hierarchical clustering dendorgram such as one returned by
#' \code{hclust}.
#' @param cutHeight Maximum joining heights that will be considered.
#' @param minSize Minimum cluster size.
#' @return A vector of numerical labels giving the assigment of each object.
#' @note The numerical labels may not be sequential. See
#' \code{\link{normalizeLabels}} for a way to put the labels into a standard
#' order.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso
#' \code{\link{hclust}}, \code{\link{cutree}},
#' \code{\link{normalizeLabels}}
#' @keywords cluster
moduleNumber <- function(dendro,
                         cutHeight = 0.9,
                         minSize = 50) {
    Branches = cutree(dendro, h = cutHeight)
    NOnBranches = table(Branches)
    TrueBranch = NOnBranches >= minSize
    Branches[!TrueBranch[Branches]] = 0

    Branches
}

# numbers2colors ####
#' Color representation for a numeric variable
#'
#' The function creates a color represenation for the given numeric input.
#'
#' Each column of \code{x} is processed individually, meaning that the color
#' palette is adjusted individually for each column of \code{x}.
#'
#' @param x a vector or matrix of numbers. Missing values are allowed and will
#' be assigned the color given in \code{naColor}. If a matrix, each column of
#' the matrix is processed separately and the return value will be a matrix of
#' colors.
#' @param signed logical: should \code{x} be considered signed? If \code{TRUE},
#' the default setting is to use to use a palette that starts with green for
#' the most negative values, continues with white for values around zero and
#' turns red for positive values. If \code{FALSE}, the default palette ranges
#' from white for minimum values to red for maximum values. If not given, the
#' behaviour is controlled by values in \code{x}: if there are both positive
#' and negative values, \code{signed} will be considered \code{TRUE}, otherwise
#' \code{FALSE}.
#' @param centered logical. If \code{TRUE} and \code{signed==TRUE}, numeric
#' value zero will correspond to the middle of the color palette. If
#' \code{FALSE} or \code{signed==FALSE}, the middle of the color palette will
#' correspond to the average of the minimum and maximum value. If neither
#' \code{signed} nor \code{centered} are given, \code{centered} will follow
#' \code{signed} (see above).
#' @param lim optional specification of limits, that is numeric values that
#' should correspond to the first and last entry of \code{colors}.
#' @param commonLim logical: should limits be calculated separately for each
#' column of x, or should the limits be the same for all columns? Only applies
#' if \code{lim} is \code{NULL}.
#' @param colors color palette to represent the given numbers.
#' @param naColor color to represent missing values in \code{x}.
#' @return A vector or matrix (of the same dimensions as \code{x}) of colors.
#' @author Peter Langfelder
#' @seealso \code{\link{labels2colors}} for color coding of ordinal labels.
#' @keywords misc
numbers2colors <- function(x,
                           signed = NULL,
                           centered = signed,
                           lim = NULL,
                           commonLim = FALSE,
                           colors = if (signed) {
                               blueWhiteRed(100)
                           } else {
                               blueWhiteRed(100)[51:100]
                           },
                           naColor = "grey") {
    x = as.matrix(x)
    if (!is.numeric(x))
        stop("'x' must be numeric. For a factor, please use as.numeric(x) in
             the call.")
    if (is.null(signed)) {
        if (any(x < 0, na.rm = TRUE) & any(x > 0, na.rm = TRUE))
        {
            signed = TRUE
        } else
            signed = FALSE
    }
    if (is.null(centered))
        centered = signed

    if (is.null(lim)) {
        if (signed & centered) {
            max = apply(abs(x), 2, max, na.rm = TRUE)
            lim = as.matrix(cbind(-max, max))
        } else {
            lim = as.matrix(cbind(
                apply(x, 2, min, na.rm = TRUE),
                apply(x, 2, max, na.rm = TRUE)
            ))
        }
        if (commonLim)
            lim = c(min(lim[, 1], na.rm = TRUE), max(lim[, 2], na.rm = TRUE))
    }
    if (is.null(dim(lim))) {
        if (length(lim) != 2)
            stop("'lim' must be a vector of length 2 or a matrix with 2 columns.")
        if (!is.numeric(lim))
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != 2)
            stop("'lim' must be finite.")
        lim = t(as.matrix(lim))
    } else {
        if (ncol(x) != nrow(lim))
            stop("Incompatible numbers of columns in 'x' and rows in 'lim'.")
        if (!is.numeric(lim))
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != length(lim))
            stop("'lim' must be finite.")
    }

    xMin = matrix(lim[, 1],
                  nrow = nrow(x),
                  ncol = ncol(x),
                  byrow = TRUE)
    xMax = matrix(lim[, 2],
                  nrow = nrow(x),
                  ncol = ncol(x),
                  byrow = TRUE)

    if (sum(xMin == xMax) > 0)
        warning("(some columns in) 'x' are constant. Their color will be the
                color of NA.")

    xx = x
    xx[is.na(xx)] = ((xMin + xMax)[is.na(xx)]) / 2
    if (sum(x < xMin, na.rm = TRUE) > 0)
    {
        warning("Some values of 'x' are below given minimum and will be
                truncated to the minimum.")
        x[xx < xMin] = xMin[xx < xMin]
    }

    if (sum(x > xMax, na.rm = TRUE) > 0)
    {
        warning("Some values of 'x' are above given maximum and will be
                truncated to the maximum.")
        x[xx > xMax] = xMax[xx > xMax]
    }

    mmEq = xMin == xMax

    nColors = length(colors)

    xCol = array(naColor, dim = dim(x))

    xInd = (x - xMin) / (xMax - xMin)
    xInd[xInd == 1] = 1 - 0.5 / nColors
    xCol[!mmEq] = colors[as.integer(xInd[!mmEq] * nColors) + 1]
    xCol[is.na(xCol)] = naColor

    xCol
}
