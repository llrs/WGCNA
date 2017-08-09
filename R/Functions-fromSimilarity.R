# Functions to perform WGCNA from similarity input.



#' Construct a network from a matrix
#'
#' Constructs a network
#'
#' If \code{signed} is \code{FALSE}, the matrix \code{mat} is first converted
#' to its absolute value.
#'
#' This function then symmetrizes the matrix using the \code{symmetrizeMethod}
#' component-wise on \code{mat} and \code{t(mat)} (i.e., the transpose of
#' \code{mat}).
#'
#' In the next step, the symmetrized matrix is linearly scaled to the interval
#' [0,1] using either \code{min} and \code{max} (each either supplied or
#' determined from the matrix). Values outside of the [min, max] range are
#' truncated to \code{min} or \code{max}.
#'
#' Lastly, the adjacency is calculated by rasing the matrix to \code{power}.
#' The diagonal of the result is set to \code{diagEntry}. Note that most WGCNA
#' functions expect the diagonal of an adjacency matrix to be 1.
#'
#' @param mat matrix to be turned into a network. Must be square.
#' @param symmetrizeMethod method for symmetrizing the matrix. The method will
#' be applied to each component of mat and its transpose.
#' @param signed logical: should the resulting network be signed? Unsigned
#' networks are constructed from \code{abs(mat)}.
#' @param min minimum allowed value for \code{mat}. If \code{NULL}, the actual
#' attained minimum of \code{mat} will be used. Missing data are ignored.
#' Values below \code{min} are truncated to \code{min}.
#' @param max maximum allowed value for \code{mat}. If \code{NULL}, the actual
#' attained maximum of \code{mat} will be used. Missing data are ignored.
#' Values below \code{max} are truncated to \code{max}.
#' @param power the soft-thresholding power.
#' @param diagEntry the value of the entries on the diagonal in the result.
#' This is usally 1 but some applications may require a zero (or even NA)
#' diagonal.
#' @return The adjacency matrix that encodes the network.
#' @author Peter Langfelder
#' @seealso \code{adjacency} for calculation of a correlation network
#' (adjacency) from a numeric matrix such as expression data
#'
#' \code{adjacency.fromSimilarity} for simpler calculation of a network from a
#' symmetric similarity matrix.
#' @keywords misc
#' @examples
#'
#' mat <- matrix(rnorm(25), ncol = 5)
#' matrixToNetwork(mat, "max")
#'
#' @export matrixToNetwork
matrixToNetwork <- function(mat, symmetrizeMethod = c("average", "min", "max"),
                            signed = TRUE, min = NULL, max = NULL, power = 12,
                            diagEntry = 1) {
    sm <- match.arg(symmetrizeMethod, c("average", "min", "max"))
    if (is.na(sm)){
        stop("Unrecognized or non - unique 'symmetrizeMethod'.")
    }
    mat <- as.matrix(mat)

    nd <- length(dim(mat))
    if (nd != 2) {
        stop("'mat' appears to have incorrect type; must be a 2 - dimensional",
             "square matrix.")
    }
    if (ncol(mat) != nrow(mat)) {
        stop("'mat' must be a square matrix.")
    }
    if (!signed) {
        mat <- abs(mat)
    }
    if (sm == 1) {
        mat <- (mat + t(mat))/2
    } else if (sm == 2) {
        mat <- pmin(mat, t(mat), na.rm = TRUE)
    } else {
        mat <- pmax(mat, t(mat), na.rm = TRUE)
    }
    if (is.null(min)) {
        min <- min(mat, na.rm = TRUE)
    } else {
        mat[mat < min] <- min
    }

    if (is.null(max)) {
        max <- max(mat, na.rm = TRUE)
    } else {
        mat[mat > max] <- max
    }

    adj <- ((mat - min)/(max - min)) ^ power
    diag(adj) = diagEntry
    checkAdjMat(adj, max = max(diagEntry, adj, 1, na.rm = TRUE))
    adj
}

#' @name checkAdjMat
#' @rdname checkAdjMat
#' @export
checkSimilarity <- function(similarity, min =  - 1, max = 1) {
    checkAdjMat(similarity, min, max)
}
# Kept for backward compatibility
#' @name adjacency
#' @rdname adjacency
#' @examples
#' datExpr <- matrix(rnorm(100), 10, 10)
#' similarity <- similarity(datExpr)
#' adj <-  adjacency.fromSimilarity(similarity)
#' @export
adjacency.fromSimilarity <- function(similarity, ...) {
    checkSimilarity(similarity)
    adjacency(similarity, ...)
}

# Kept for backward compatibility


#' Analysis of scale free topology for hard-thresholding.
#'
#' Analysis of scale free topology for multiple hard thresholds. The aim is to
#' help the user pick an appropriate threshold for network construction.
#'
#' The function calculates unsigned networks by thresholding the correlation
#' matrix using thresholds given in \code{cutVector}. For each power the scale
#' free topology fit index is calculated and returned along with other
#' information on connectivity.
#'
#' @param similarity Matrix whose values are between -1 and 1
#' @param ... Other arguments from pickHardThreshold
#' @param data expression data in a matrix or data frame. Rows correspond to
#' samples and columns to genes.
#' @param dataIsExpr logical: should the data be interpreted as expression (or
#' other numeric) data, or as a similarity matrix of network nodes?
#' @param RsquaredCut desired minimum scale free topology fitting index
#' \eqn{R^2}.
#' @param cutVector a vector of hard threshold cuts for which the scale free
#' topology fit indices are to be calculated.
#' @param moreNetworkConcepts logical: should additional network concepts be
#' calculated? If \code{TRUE}, the function will calculate how the network
#' density, the network heterogeneity, and the network centralization depend on
#' the power. For the definition of these additional network concepts, see
#' Horvath and Dong (2008).  PloS Comp Biol.
#' @param removeFirst should the first bin be removed from the connectivity
#' histogram?
#' @param nBreaks number of bins in connectivity histograms
#' @param corFnc a character string giving the correlation function to be used
#' in adjacency calculation.
#' @param corOptions further options to the correlation function specified in
#' \code{corFnc}.
#' @return A list with the following components:
#'
#' \item{cutEstimate}{ estimate of an appropriate hard-thresholding cut: the
#' lowest cut for which the scale free topology fit \eqn{R^2} exceeds
#' \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} for all cuts,
#' \code{NA} is returned. }
#'
#' \item{fitIndices}{ a data frame containing the fit indices for scale free
#' topology. The columns contain the hard threshold, Student p-value for the
#' correlation threshold, adjusted \eqn{R^2} for the linear fit, the linear
#' coefficient, adjusted \eqn{R^2} for a more complicated fit models, mean
#' connectivity, median connectivity and maximum connectivity.  If input
#' \code{moreNetworkConcepts} is \code{TRUE}, 3 additional columns containing
#' network density, centralization, and heterogeneity.}
#' @author Steve Horvath
#' @seealso \code{\link{signumAdjacency}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @keywords misc
#' @export pickHardThreshold.fromSimilarity
#' @rdname pickHardThreshold
pickHardThreshold.fromSimilarity <- function(similarity, ...) {
    checkSimilarity(similarity)
    pickHardThreshold(similarity, dataIsExpr = FALSE, ...)
}



#' Analysis of scale free topology for soft-thresholding
#'
#' Analysis of scale free topology for multiple soft thresholding powers. The
#' aim is to help the user pick an appropriate soft-thresholding power for
#' network construction.
#'
#' The function calculates weighted networks either by interpreting \code{data}
#' directly as similarity, or first transforming it to similarity of the type
#' specified by \code{networkType}. The weighted networks are obtained by
#' raising the similarity to the powers given in \code{powerVector}.  For each
#' power the scale free topology fit index is calculated and returned along
#' with other information on connectivity.
#'
#' On systems with multiple cores or processors, the function pickSoftThreshold
#' takes advantage of parallel processing if the function
#' \code{\link{enableWGCNAThreads}} has been called to allow parallel
#' processing and set up the parallel calculation back-end.
#'
#' @param similarity Matrix whose values are between -1 and 1
#' @param ... Other arguments from pickSoftThreshold
#' @param data expression data in a matrix or data frame. Rows correspond to
#' samples and columns to genes.
#' @param dataIsExpr logical: should the data be interpreted as expression (or
#' other numeric) data, or as a similarity matrix of network nodes?
#' @param RsquaredCut desired minimum scale free topology fitting index
#' \eqn{R^2}.
#' @param powerVector a vector of soft thresholding powers for which the scale
#' free topology fit indices are to be calculated.
#' @param removeFirst should the first bin be removed from the connectivity
#' histogram?
#' @param nBreaks number of bins in connectivity histograms
#' @param blockSize block size into which the calculation of connectivity
#' should be broken up. If not given, a suitable value will be calculated using
#' function \code{blockSize} and printed if \code{verbose>0}. If R runs into
#' memory problems, decrease this value.
#' @param corFnc the correlation function to be used in adjacency calculation.
#' @param corOptions a list giving further options to the correlation function
#' specified in \code{corFnc}.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param moreNetworkConcepts logical: should additional network concepts be
#' calculated? If \code{TRUE}, the function will calculate how the network
#' density, the network heterogeneity, and the network centralization depend on
#' the power. For the definition of these additional network concepts, see
#' Horvath and Dong (2008). PloS Comp Biol.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{powerEstimate}{ estimate of an appropriate soft-thresholding power:
#' the lowest power for which the scale free topology fit \eqn{R^2} exceeds
#' \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} for all powers,
#' \code{NA} is returned. }
#'
#' \item{fitIndices}{ a data frame containing the fit indices for scale free
#' topology. The columns contain the soft-thresholding power, adjusted
#' \eqn{R^2} for the linear fit, the linear coefficient, adjusted \eqn{R^2} for
#' a more complicated fit models, mean connectivity, median connectivity and
#' maximum connectivity. If input \code{moreNetworkConcepts} is \code{TRUE}, 3
#' additional columns containing network density, centralization, and
#' heterogeneity.}
#' @author Steve Horvath and Peter Langfelder with improvements from Alexey
#' Segushichev
#' @seealso \code{\link{adjacency}}, \code{\link{softConnectivity}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @keywords misc misc
#' @export pickSoftThreshold.fromSimilarity
#' @rdname pickSoftThreshold
pickSoftThreshold.fromSimilarity <- function(similarity, ...) {
    checkSimilarity(similarity)
    pickSoftThreshold(similarity, dataIsExpr = FALSE, ...)

}
