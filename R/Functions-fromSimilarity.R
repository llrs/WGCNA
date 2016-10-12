# Functions to perform WGCNA from similarity input.
#' @rdname matrixToNetwork
#' @name matrixToNetwork
#' @title Construct a network from a matrix
#' @description
#' Constructs a network
#' @param mat matrix to be turned into a network. Must be square.
#' @param symmetrizeMethod method for symmetrizing the matrix. The method will
#' be applied to each component of mat and its transpose.
#' @param signed logical: should the resulting network be signed? Unsigned
#' networks are constructed from abs(mat).
#' @param min minimum allowed value for mat. If NULL, the actual attained
#' minimum of mat will be used. Missing data are ignored. Values below min are
#' truncated to min.
#' @param max maximum allowed value for mat. If NULL, the actual attained
#' maximum of mat will be used. Missing data are ignored. Values below max are
#' truncated to max.
#' @param power the soft-thresholding power.
#' @param the value of the entries on the diagonal in the result. This is usally
#' 1 but some applications may require a zero (or even NA) diagonal.
#' @details
#' If signed is FALSE, the matrix mat is first converted to its absolute value.
#'
#' This function then symmetrizes the matrix using the symmetrizeMethod
#' component-wise on mat and t(mat) (i.e., the transpose of mat).
#'
#' In the next step, the symmetrized matrix is linearly scaled to the interval
#' [0,1] using either min and max (each either supplied or determined from the
#' matrix). Values outside of the [min, max] range are truncated to min or max.
#'
#' Lastly, the adjacency is calculated by rasing the matrix to power. The
#' diagonal of the result is set to diagEntry. Note that most WGCNA functions
#' expect the diagonal of an adjacency matrix to be 1.
#' @return
#' The adjacency matrix that encodes the network
#' @author
#' Peter Langfelder
#' @seealso
#' \code{\link{adjacency}} for calculation of a correlation network (adjacency)
#' from a numeric matrix such as expression data, and adjacency.fromSimilarity
#' for simpler calculation of a network from a symmetric similarity matrix.
#' @export
matrixToNetwork <- function(mat,
                            symmetrizeMethod = c("average", "min", "max"),
                            signed = TRUE,
                            min = NULL,
                            max = NULL,
                            power = 12,
                            diagEntry = 1) {
    sm = match.arg(symmetrizeMethod)
    if (is.na(sm)){
        stop("Unrecognized or non - unique 'symmetrizeMethod'.")
    }
    mat = as.matrix(mat)

    nd = 0
    x = try({nd = dim(mat)})
    if ((class(x) == 'try - error') | (nd != 2))
        stop("'mat' appears to have incorrect type; must be a 2 - dimensional
             square matrix.")

    if (ncol(mat) != nrow(mat))
        stop("'mat' must be a square matrix.")

    if (!signed) mat = abs(mat)

    if (sm == 1) {
        mat = (mat + t(mat))/2
    } else if (sm == 2) {
        mat = pmin(mat, t(mat), na.rm = TRUE)
    } else
        mat = pmax(mat, t(mat), na.rm = TRUE)

  if (is.null(min)) {
    min = min(mat, na.rm = TRUE)
  } else
    mat[mat < min] = min

  if (is.null(max)) {
    max = max(mat, na.rm = TRUE)
  } else
    mat[mat > max] = max

adj = ((mat - min)/(max - min))^power

  diag(adj) = diagEntry

  adj
}
#' @name checkAdjMat
#' @rdname checkAdjMat
#' @export
checkSimilarity <- function(similarity, min =  - 1, max = 1)
{
  checkAdjMat(similarity, min, max)
}

#' @name adjacency
#' @rdname adjacency
#' @examples
#' adj <-  adjacency.fromSimilarity(similarity)
#' @export
adjacency.fromSimilarity <- function(similarity, type = "unsigned",
                                     power = if (type == "distance") 1 else 6) {
  checkSimilarity(similarity)
  adjacency(similarity, type = type, power = power, corFnc = "I",
            corOptions = "", distFnc = "I",
            distOptions = "")
}

pickHardThreshold.fromSimilarity <- function(similarity,
                                             RsquaredCut = 0.85,
                                             cutVector = seq(0.1, 0.9,
                                                             by = 0.05),
                                             moreNetworkConcepts = FALSE ,
                                             removeFirst = FALSE,
                                             nBreaks = 10) {
    checkSimilarity(similarity)
    pickHardThreshold(similarity, dataIsExpr = FALSE,
                      RsquaredCut =  RsquaredCut, cutVector = cutVector,
                      moreNetworkConcepts = moreNetworkConcepts,
                      removeFirst = removeFirst,
                      nBreaks = nBreaks, corFnc = "I", corOptions = "")
}

pickSoftThreshold.fromSimilarity <- function(similarity,
                                     RsquaredCut = 0.85,
                                     powerVector = c(seq(1, 10, by = 1),
                                                     seq(12, 20, by = 2)),
                                     removeFirst = FALSE, nBreaks = 10,
                                     blockSize = 1000,
                                     networkType = "unsigned",
                                     moreNetworkConcepts = FALSE,
                                     verbose = 0, indent = 0) {
  checkSimilarity(similarity)
  pickSoftThreshold(similarity, dataIsExpr = FALSE,
                    RsquaredCut =  RsquaredCut, powerVector = powerVector,
                    removeFirst = removeFirst, nBreaks = nBreaks,
                    blockSize = blockSize, networkType = networkType,
                    moreNetworkConcepts = moreNetworkConcepts,
                    verbose = verbose, indent = indent)

}
