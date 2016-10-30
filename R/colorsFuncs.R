# Convert modules numbers to colors

#' Convert numerical labels to colors.
#'
#' Converts a vector or array of numerical labels into a corresponding vector
#' or array of colors corresponding to the labels.
#'
#' If \code{labels} is numeric, it is used directly as index to the standard
#' color sequence. If 0 is present among the labels and \code{zeroIsGrey=TRUE},
#' labels 0 are given grey color.
#'
#' If \code{labels} is not numeric, its columns are turned into factors and the
#' numeric representation of each factor is used to assign the corresponding
#' colors. In this case \code{commonColorCode} governs whether each column gets
#' its own color code, or whether the color code will be universal.
#'
#' The standard sequence start with well-distinguishable colors, and after
#' about 40 turns into a quasi-random sampling of all colors available in R
#' with the exception of all shades of grey (and gray).
#'
#' If the input \code{labels} have a dimension attribute, it is copied into the
#' output, meaning the dimensions of the returned value are the same as those
#' of the input \code{labels}.
#'
#' @param labels Vector or matrix of non-negative integer or other (such as
#' character) labels. See details.
#' @param zeroIsGrey If TRUE, labels 0 will be assigned color grey. Otherwise,
#' labels below 1 will trigger an error.
#' @param colorSeq Color sequence corresponding to labels. If not given, a
#' standard sequence will be used.
#' @param naColor Color that will encode missing values.
#' @param commonColorCode logical: if \code{labels} is a matrix, should each
#' column have its own colors?
#' @return A vector or array of character strings of the same length or
#' dimensions as \code{labels}.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords color
#' @examples
#'
#' labels = c(0:20);
#' labels2colors(labels);
#' labels = matrix(letters[1:9], 3,3);
#' labels2colors(labels)
#' # Note the difference when commonColorCode = FALSE
#' labels2colors(labels, commonColorCode = FALSE)
#'
labels2colors <- function(labels, zeroIsGrey = TRUE, colorSeq = NULL,
                          naColor = "grey",
                          commonColorCode = TRUE) {

    if (is.null(colorSeq)) colorSeq = standardColors()

    if (is.numeric(labels)) {
        if (zeroIsGrey) minLabel = 0 else minLabel = 1
        if (any(labels<0, na.rm = TRUE)) minLabel = min(c(labels), na.rm = TRUE)
        nLabels = labels
    } else {
        if (commonColorCode) {
            factors = factor(c(as.matrix(as.data.frame(labels))))
            nLabels = as.numeric(factors)
            dim(nLabels) = dim(labels)
        } else {
            labels = as.matrix(as.data.frame(labels))
            factors = list()
            for (c in 1:ncol(labels))
                factors[[c]] = factor(labels[, c])
            nLabels = sapply(factors, as.numeric)
        }
    }

    if (max(nLabels, na.rm = TRUE) > length(colorSeq)) {
        nRepeats = as.integer((max(labels) - 1)/length(colorSeq)) + 1
        warning(paste(
            "labels2colors: Number of labels exceeds number of avilable colors.",
            "Some colors will be repeated", nRepeats, "times."))
        extColorSeq = colorSeq
        for (rep in 1:nRepeats)
            extColorSeq = c(extColorSeq, paste0(colorSeq, ".", rep))
    } else {
        nRepeats = 1
        extColorSeq = colorSeq
    }
    colors = rep("grey", length(nLabels))
    fin = !is.na(nLabels)
    colors[!fin] = naColor
    finLabels = nLabels[fin]
    colors[fin][finLabels != 0] = extColorSeq[finLabels[finLabels != 0]]
    if (!is.null(dim(labels)))
        dim(colors) = dim(labels)
    colors
}

#' Colors this library uses for labeling modules.
#'
#' Returns the vector of color names in the order they are assigned by other
#' functions in this library.
#'
#'
#' @param n Number of colors requested. If \code{NULL}, all (approx. 450)
#' colors will be returned. Any other invalid argument such as less than one or
#' more than maximum (\code{length(standardColors())}) will trigger an error.
#' @return A vector of character color names of the requested length.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords color misc
#' @examples
#'
#' standardColors(10);
#'
standardColors <- function(n = NULL){
    if (is.null(n)) {
        return(.GlobalStandardColors)
    }
    if ((n > 0) && (n <= length(.GlobalStandardColors))) {
        return(.GlobalStandardColors[c(1:n)])
    } else {
        stop("Invalid number of standard colors requested.")
    }
}

#' Show colors used to label modules
#'
#' The function plots a barplot using colors that label modules.
#'
#' To see the first \code{n} colors, use argument \code{colors =
#' standardColors(n)}.
#'
#' @param colors colors to be displayed. Defaults to all colors available for
#' module labeling.
#' @return None.
#' @author Peter Langfelder
#' @seealso \code{\link{standardColors}}
#' @keywords misc
#' @examples
#'
#' displayColors(standardColors(10))
#'
displayColors <- function(colors = NULL) {
    if (is.null(colors))
        colors = standardColors()
    barplot(rep(1, length(colors)), col = colors, border = colors)
}
