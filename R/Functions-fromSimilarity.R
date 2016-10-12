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
#' @export matrixToNetwork
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
