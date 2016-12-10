# coClustering ####

#' Co-clustering measure of cluster preservation between two clusterings
#'
#' The function calculates the co-clustering statistics for each module in the
#' reference clustering.
#'
#' Co-clustering of cluster q in the reference clustering and cluster q' in the
#' test clustering measures the overlap of clusters q and q' by the number of
#' tuplets that can be chosen from the overlap of clusters q and q' relative to
#' the number of tuplets in cluster q. To arrive at a co-clustering measure for
#' cluster q, we sum the co-clustering of q and q' over all clusters q' in the
#' test clustering. A value close to 1 indicates high preservation of the
#' reference cluster in the test clustering, while a value close to zero
#' indicates a low preservation.
#'
#' @param clusters.ref Reference input clustering. A vector in which each
#' element gives the cluster label of an object.
#' @param clusters.test Test input clustering. Must be a vector of the same
#' size as \code{cluster.ref}.
#' @param tupletSize Co-clutering tuplet size.
#' @param unassignedLabel Optional specification of a clustering label that
#' denotes unassigned objects. Objects with this label are excluded from the
#' calculation.
#' @return A vector in which each component corresponds to a cluster in the
#' reference clustering. Entries give the co-clustering measure of cluster
#' preservation.
#' @author Peter Langfelder
#' @seealso \code{\link{modulePreservation}} for a large suite of module
#' preservation statistics \code{\link{coClustering.permutationTest}} for a
#' permutation test for co-clustering significance
#' @references For example, see Langfelder P, Luo R, Oldham MC, Horvath S
#' (2011) Is My Network Module Preserved and Reproducible? PLoS Comput Biol
#' 7(1): e1001057. Co-clustering is discussed in the Methods Supplement
#' (Supplementary text 1) of that article.
#' @keywords misc
#' @examples
#' # An example with random (unrelated) clusters:
#' set.seed(1)
#' nModules <- 10
#' nGenes <- 1000
#' cl1 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' cl2 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' coClustering(cl1, cl2)
#' coClustering(cl1, cl2, tupletSize = 3)
#' # For the same reference and test clustering:
#' coClustering(cl1, cl1)
coClustering <- function(clusters.ref, clusters.test, tupletSize = 2,
                        unassignedLabel = 0) {
  overlap <- table(clusters.test, clusters.ref)
  greyRow <- rownames(overlap) == as.character(unassignedLabel)
  greyCol <- colnames(overlap) == as.character(unassignedLabel)

  refModSizes <- table(clusters.ref)

  ccNumer <- apply(overlap[!greyRow, !greyCol, drop = FALSE], 2,
                   choose, tupletSize)
  ccDenom <- choose(refModSizes[!greyCol], tupletSize)

  apply(ccNumer, 2, sum)/ccDenom
}


#' Permutation test for co-clustering
#'
#' This function calculates permutation Z statistics that measure how different
#' the co-clustering of modules in a reference and test clusterings is from
#' random.
#'
#' This function performs a permutation test to determine whether observed
#' co-clustering statistics are significantly different from those expected by
#' chance. It returns the observed co-clustering as well as the permutation Z
#' statistic, calculated as \code{(observed - mean)/sd}, where \code{mean} and
#' \code{sd} are the mean and standard deviation of the co-clustering when the
#' test clustering is repeatedly randomly permuted.
#'
#' @param clusters.ref Reference input clustering. A vector in which each
#' element gives the cluster label of an object.
#' @param clusters.test Test input clustering. Must be a vector of the same
#' size as \code{cluster.ref}.
#' @param tupletSize Co-clutering tuplet size.
#' @param nPermutations Number of permutations to execute. Since the function
#' calculates parametric p-values, a relatively small number of permutations
#' (at least 50) should be sufficient.
#' @param unassignedLabel Optional specification of a clustering label that
#' denotes unassigned objects. Objects with this label are excluded from the
#' calculation.
#' @param randomSeed Random seed for initializing the random number generator.
#' If \code{NULL}, the generator is not initialized (useful for calling the
#' function sequentially). The default assures reproducibility.
#' @param verbose If non-zero, function will print out progress messages.
#' @param indent Indentation for progress messages. Each unit adds two spaces.
#' @return \item{observed }{the observed co-clustering measures for clusters in
#' \code{clusters.ref} } \item{Z}{permutation Z statics}
#' \item{permuted.mean}{means of the co-clustering measures when the test
#' clustering is permuted} \item{permuted.sd}{standard deviations of the
#' co-clustering measures when the test clustering is permuted}
#' \item{permuted.cc}{values of the co-clustering measure for each permutation
#' of the test clustering. A matrix of dimensions (number of
#' permutations)x(number of clusters in reference clustering). }
#' @author Peter Langfelder
#' @seealso \code{\link{coClustering}} for calculation of the "observed"
#' co-clustering measure \code{\link{modulePreservation}} for a large suite of
#' module preservation statistics
#' @references For example, see Langfelder P, Luo R, Oldham MC, Horvath S
#' (2011) Is My Network Module Preserved and Reproducible? PLoS Comput Biol
#' 7(1): e1001057. Co-clustering is discussed in the Methods Supplement
#' (Supplementary text 1) of that article.
#' @keywords misc
#' @examples
#'
#'   set.seed(1)
#'   nModules = 5
#'   nGenes = 100
#'   cl1 = sample(c(1:nModules), nGenes, replace = TRUE)
#'   cl2 = sample(c(1:nModules), nGenes, replace = TRUE)
#'
#'   cc = coClustering(cl1, cl2)
#'
#'   # Choose a low number of permutations to make the example fast
#'   ccPerm = coClustering.permutationTest(cl1, cl2, nPermutations = 20, verbose = 1)
#'
#'   ccPerm$observed
#'   ccPerm$Z
#'
#'   # Combine cl1 and cl2 to obtain clustering that is somewhat similar to cl1:
#'
#'   cl3 = cl2
#'   from1 = sample(c(TRUE, FALSE), nGenes, replace = TRUE)
#'   cl3[from1] = cl1[from1]
#'
#'   ccPerm = coClustering.permutationTest(cl1, cl3, nPermutations = 20, verbose = 1)
#'
#'   # observed co-clustering is higher than before:
#'   ccPerm$observed
#'
#'   # Note the high preservation Z statistics:
#'   ccPerm$Z
#'
coClustering.permutationTest = function(clusters.ref, clusters.test,
                                        tupletSize=2, nPermutations = 100, unassignedLabel=0,
                                        randomSeed = 12345, verbose = 0, indent = 0)
{
  spaces = indentSpaces(indent)
  if (!is.null(randomSeed)) set.seed(randomSeed)
  observed = coClustering(clusters.ref, clusters.test, tupletSize, unassignedLabel)

  nModules = length(observed)

  permValues = matrix(NA, nPermutations, nModules)

  if (verbose > 0)
    pind = initProgInd(paste0(spaces, "Running permutations: "), " done")
  for (p in 1:nPermutations)
  {
    ctPerm = sample(clusters.test)
    permValues[p, ] = as.numeric(coClustering(clusters.ref, ctPerm,
                                        tupletSize, unassignedLabel))
    if (verbose > 0) pind = updateProgInd(p/nPermutations, pind)
  }
  if (verbose > 0) printFlush("")
  means = colMeans(permValues)
  sds = apply(permValues, 2, sd, na.rm = TRUE)
  list(observed = observed, Z = (observed-means)/sds, permuted.mean = means, permuted.sd = sds,
       permuted.cc = permValues)
}

# randIndex ####
#  Rand index calculation
# this function is used for computing the Rand index below...
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
    diag(adjMat) <- 0
    nNodes <- ncol(adjMat)
    computeLinksInNeighbors <- function(x, imatrix) {
        x %*% imatrix %*% x
    }
    nolinksNeighbors <- c(rep(-666, nNodes))
    total.edge <- c(rep(-666, nNodes))
    maxh1 <- max(as.dist(adjMat))
    minh1 <- min(as.dist(adjMat))
    if (maxh1 > 1 | minh1 < 0) {
        stop("The adjacency matrix contains entries that are larger than 1 or
                smaller than 0: max  = ", maxh1, ", min  = ", minh1)
    }
    nolinksNeighbors <- apply(adjMat, 1, computeLinksInNeighbors,
                              imatrix = adjMat)
    plainsum  <- apply(adjMat, 1, sum)
    squaresum <- apply(adjMat ^ 2, 1, sum)
    total.edge <- plainsum ^ 2 - squaresum
    CChelp <- rep(-666, nNodes)
    CChelp <- ifelse(total.edge == 0, 0, nolinksNeighbors / total.edge)
    CChelp
} # end of function
