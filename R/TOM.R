# GTOMdist ####
#' Generalized Topological Overlap Measure
#'
#' Generalized Topological Overlap Measure, taking into account interactions of
#' higher degree. This function computes a TOM dissimilarity
#'
#' WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains
#' binary entries.
#'
#' @param adjMat adjacency matrix with entries between [0, 1]
#' @param degree integer specifying the maximum degree to be calculated.
#' @return Matrix of the same dimension as the input \code{adjMat}.
#' @author Steve Horvath and Andy Yip
#' @references Yip A, Horvath S (2007) Gene network interconnectedness and the
#' generalized topological overlap measure. BMC Bioinformatics 2007, 8:22
#' @keywords misc
GTOMdist <- function(adjMat, degree = 1) {
    maxh1 = max(as.dist(adjMat))
    minh1 = min(as.dist(adjMat))
    if (degree != round(abs(degree)))
        stop("'degree' must be a positive integer.")
    if (maxh1 > 1 | minh1 < 0)
        stop(
            paste(
                "Entries of the adjacency matrix are not between 0 and 1: max  = ",
                maxh1,
                ", min  = ",
                minh1
            )
        )

    if (max(c(as.dist(abs(
        adjMat - t(adjMat)
    )))) > 0)
        stop("Given adjacency matrix is not symmetric.")

    B <- adjMat
    if (degree >= 2)
        for (i in 2:degree) {
            diag(B) <- diag(B) + 1
            # Calculates the number of paths with length at most degree connecting
            #  a pair
            B = B %*% adjMat
        }
    # this gives the degree - step reachability from a node to another
    B <- (B > 0)
    diag(B) <- 0 # exclude each node being its own neighbor
    # this gives the number of common degree-step-neighbor that a pair of nodes
    # share
    B <- B %*% B

    Nk <- diag(B)
    B <- B + adjMat # numerator
    diag(B) <- 1
    denomTOM = outer(Nk, Nk, FUN = "pmin") + 1 - adjMat
    diag(denomTOM) <- 1
    1 - B / denomTOM   # this turns the TOM matrix into a dissimilarity
}

# vectorTOM ####
# vectorTOM: calculate TOM of a vector (or a 'small' matrix) with expression
# data. If the number of columns in vect is small (or 1), number of columns in
# datExpr can be large.
#' Topological overlap for a subset of the whole set of genes
#'
#' This function calculates topological overlap of a small set of vectors with
#' respect to a whole data set.
#'
#' Topological overlap can be viewed as the normalized count of shared
#' neighbors encoded in an adjacency matrix. In this case, the adjacency matrix
#' is calculated between the columns of \code{vect} and \code{datExpr} and the
#' topological overlap of vectors in \code{vect} measures the number of shared
#' neighbors in \code{datExpr} that vectors of \code{vect} share.
#'
#' @param datExpr a data frame containing the expression data of the whole set,
#' with rows corresponding to samples and columns to genes.
#' @param vect a single vector or a matrix-like object containing vectors whose
#' topological overlap is to be calculated.
#' @param subtract1 logical: should calculation be corrected for
#' self-correlation? Set this to \code{TRUE} if \code{vect} contains a subset
#' of \code{datExpr}.
#' @param blockSize maximum block size for correlation calculations. Only
#' important if \code{vect} contains a large number of columns.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function.
#' @param networkType character string giving network type. Allowed values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, \code{"signed
#' hybrid"}. See \code{\link{adjacency}}.
#' @param power soft-thresholding power for network construction.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of dimensions \code{n*n}, where \code{n} is the number of
#' columns in \code{vect}.
#' @author Peter Langfelder
#' @seealso
#' \code{\link{TOMsimilarity}} for standard calculation of topological
#' overlap.
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
vectorTOM <-function(datExpr, vect, subtract1 = FALSE, blockSize = 2000,
                     corFnc = "cor", corOptions = "use = 'p'",
                     networkType = "unsigned", power = 6, verbose = 1,
                     indent = 0) {
    spaces = indentSpaces(indent)

    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType))
        stop(paste(
            "Unrecognized 'networkType'. Recognized values are",
            paste(.networkTypes, collapse = ", ")
        ))

    if (is.null(dim(vect)))
    {
        vect = as.matrix(vect)
        vectIsVector = TRUE
    } else
        vectIsVector = FALSE

    if (nrow(vect) != nrow(datExpr))
        stop("Input error: numbers of samples in 'vect' and
             'datExpr' must be the same.")

    if (ncol(vect) > blockSize)
      stop("Input error: number of columns in 'vect' is too large.",
           "If you are certain you want to try anyway, increase 'blockSize'",
           "to at least the number of columns in 'vect'.")

    corEval = parse(text = paste(corFnc,
                                 "(datExpr, vect ", prepComma(corOptions), ")"))
    corVE = eval(corEval)
    if (intType == 1) {
        corVE = abs(corVE)
    } else if (intType == 2) {
        corVE = (1 + corVE) / 2
    } else if (intType == 3) {
        corVE[corVE < 0] = 0
    } else
        stop("Unrecognized networkType argument. Recognized values are ",
             "'unsigned', 'signed', and 'signed hybrid'.")

    corVE = corVE ^ power

    subtract1 = as.numeric(subtract1)

    nVect = ncol(vect)
    nGenes = ncol(datExpr)
    TOM = matrix(nrow = nGenes, ncol = nVect)

    if (verbose > 0) {
        if (verbose > 1) {
            cat(paste(
                spaces,
                "Calculating TOM of a set of vectors with genes"
            ))
        }
        pind = initProgInd()
    }
    start = 1
    denomArr = array(0, dim = c(2, blockSize, nVect))
    while (start <= nGenes) {
        end = min(start + blockSize - 1, nGenes)
        blockInd = c(start:end)
        corEval = parse(text = paste(
            corFnc,
            "(datExpr[, blockInd], datExpr ",
            prepComma(corOptions),
            ")"
        ))
        corEE = eval(corEval)
        if (intType == 1) {
            corEE = abs(corEE)
        } else if (intType == 2) {
            corEE = (1 + corEE) / 2
        } else if (intType == 3) {
            corEE[corEE < 0] = 0
        }
        corEE = corEE ^ power
        num = corEE %*% corVE  - subtract1 * corVE[blockInd,]
        kV = apply(corVE, 2, sum, na.rm = TRUE) - subtract1
        kE = apply(corEE, 1, sum, na.rm = TRUE) - 1
        denomArr[1, 1:(end - start + 1),] = matrix(kV,
                                                   nrow = end - start + 1,
                                                   ncol = nVect,
                                                   byrow = TRUE)
        denomArr[2, 1:(end - start + 1),] = matrix(kE, nrow = end - start + 1,
                                                   ncol = nVect)
        denom = apply(denomArr[, 1:(end - start + 1),], c(2, 3), min) +
            1 - corVE[blockInd,]
        TOM[blockInd,] = num / denom
        if (verbose > 0)
            pind = updateProgInd(end / nGenes, pind)
        start = end + 1
        collectGarbage()
    }
    if (verbose > 0)
        printFlush(" ")

    TOM
}

# subsetTOM  ####
# subsetTOM: calculate TOM of a subset of vectors with respect to a full set of
# vectors.
#' Topological overlap for a subset of a whole set of genes
#'
#' This function calculates topological overlap of a subset of vectors with
#' respect to a whole data set.
#'
#' This function is designed to calculated topological overlaps of small
#' subsets of large expression data sets, for example in individual modules.
#'
#' @param datExpr a data frame containing the expression data of the whole set,
#' with rows corresponding to samples and columns to genes.
#' @param subset a single logical or numeric vector giving the indices of the
#' nodes for which the TOM is to be calculated.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function.
#' @param networkType character string giving network type. Allowed values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, \code{"signed
#' hybrid"}. See \code{\link{adjacency}}.
#' @param power soft-thresholding power for network construction.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of dimensions \code{n*n}, where \code{n} is the number of
#' entries selected by \code{block}.
#' @author Peter Langfelder
#' @seealso
#' \code{\link{TOMsimilarity}} for standard calculation of topological
#' overlap.
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
subsetTOM <- function(datExpr,
                      subset,
                      corFnc = "cor",
                      corOptions = "use = 'p'",
                      networkType = "unsigned",
                      power = 6,
                      verbose = 1,
                      indent = 0) {
    spaces = indentSpaces(indent)

    if (!is.null(dim(subset)))
        stop("'subset' must be a dimensionless vector.")

    if (is.null(dim(datExpr)))
        stop("'datExpr' must be a matrix or data frame.")
    if (length(dim(datExpr)) != 2)
        stop("'datExpr' must be two - dimensional.")

    nGenes = ncol(datExpr)

    if (is.logical(subset))
        subset = c(1:nGenes)[subset]

    nBlock = length(subset)

    if (any(!is.finite(subset)))
        stop("Entries of 'subset' must all be finite.")

    if (min(subset) < 1 | max(subset) > nGenes)
      stop("Some entries of 'subset' are out of range.",
           "\nNote: 'subset' must contain indices of the subset",
           "for which the TOM is calculated.")

    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType))
        stop("Unrecognized 'networkType'. Recognized values are",
            paste(.networkTypes, collapse = ", "))

    adj = adjacency(
        datExpr,
        subset,
        power = power,
        type = networkType,
        corFnc = corFnc,
        corOptions = corOptions
    )

    adj[is.na(adj)] = 0
    num = t(adj) %*% adj - adj[subset,]

    k = apply(adj, 2, sum)

    kMat = matrix(k, nBlock, nBlock)

    denom = pmin(kMat, t(kMat)) - adj[subset,]

    TOM = num / denom
    diag(TOM) = 1

    TOM
}
