#' @rdname softConnectivity
#' @rdname softConnectivity
#' @export
softConnectivity.fromSimilarity <- function(similarity, type = "unsigned",
                                            power = if (type ==  "signed") 15 else 6,
                                            blockSize = 1500, verbose = 2, indent = 0) {
    checkSimilarity(similarity)
    softConnectivity(similarity, corFnc = "I", corOptions = "",
                     type = type, power = power,
                     blockSize = blockSize, verbose = verbose, indent = indent)
}

# This function is useful for speeding up the connectivity calculation.
# The idea is to partition the adjacency matrix into consecutive baches of a
# given size. In principle, the larger the block size the faster is the
# calculation. But smaller blockSizes require less memory...
# Input: gene expression data set where  * rows *  correspond to microarray samples
# and columns correspond to genes. If fewer than minNSamples contain gene
# expression information for a given gene, then its connectivity is set to
# missing.


#' Calculates connectivity of a weighted network.
#' 
#' Given expression data or a similarity, the function constructs the adjacency
#' matrix and for each node calculates its connectivity, that is the sum of the
#' adjacency to the other nodes.
#' 
#' 
#' @aliases softConnectivity softConnectivity.fromSimilarity
#' @param datExpr a data frame containing the expression data, with rows
#' corresponding to samples and columns to genes.
#' @param similarity a similarity matrix: a square symmetric matrix with
#' entries between -1 and 1.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function.
#' @param type network type. Allowed values are (unique abbreviations of)
#' \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}.
#' @param power soft thresholding power.
#' @param blockSize block size in which adjacency is to be calculated. Too low
#' (say below 100) may make the calculation inefficient, while too high may
#' cause R to run out of physical memory and slow down the computer. Should be
#' chosen such that an array of doubles of size (number of genes) * (block
#' size) fits into available physical memory.
#' @param minNSamples minimum number of samples available for the calculation
#' of adjacency for the adjacency to be considered valid.  If not given,
#' defaults to the greater of \code{..minNSamples} (currently 4) and number of
#' samples divided by 3.  If the number of samples falls below this threshold,
#' the connectivity of the corresponding gene will be returned as \code{NA}.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A vector with one entry per gene giving the connectivity of each
#' gene in the weighted network.
#' @author Steve Horvath
#' @seealso \code{\link{adjacency}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
#' @export softConnectivity
softConnectivity <- function(datExpr,
                             corFnc = "cor", corOptions = "use = 'p'",
                             type = "unsigned",
                             power = if (type == "signed") 15 else 6,
                             blockSize = 1500, minNSamples = NULL,
                             verbose = 2, indent = 0)
{
    spaces = indentSpaces(indent)
    nGenes = dim(datExpr)[[2]]

    if (blockSize * nGenes > .largestBlockSize) {
        blockSize = as.integer(.largestBlockSize/nGenes)}
    nSamples = dim(datExpr)[[1]]
    if (is.null(minNSamples))
    {
        minNSamples = max(..minNSamples, nSamples/3)
    }

    if (nGenes<..minNGenes | nSamples<minNSamples)
        stop(paste("Error: Something seems to be wrong. \n",
                   "   Make sure that the input data frame has genes as rows and array
                   samples as columns.\n",
                   "   Alternatively, there seem to be fewer than", ..minNGenes,
                   "genes or fewer than",
                   minNSamples, "samples."))
    if (nGenes<nSamples)
        printFlush("Warning: There are fewer genes than samples in the
                   function softConnectivity. Maybe you should transpose the data?")


    k = rep(NA, nGenes)
    start = 1
    if (verbose > 0) {
        printFlush(paste(spaces,
                         "softConnectivity: FYI: connecitivty of genes with less than",
                         ceiling(minNSamples),
                         "valid samples will be returned as NA."))
        cat(paste(spaces, "..calculating connectivities.."))
        pind = initProgInd()
    }
    while (start < nGenes)
    {
        end = min(start + blockSize - 1, nGenes)
        index1 = start:end
        ad1 = adjacency(datExpr, index1, power = power, type = type,
                        corFnc = corFnc, corOptions = corOptions)
        k[index1] = colSums(ad1, na.rm = TRUE) - 1
        # If fewer than minNSamples contain gene expression information for a given
        # gene, then we set its connectivity to 0.
        NoSamplesAvailable = colSums(!is.na(datExpr[, index1]))
        k[index1][NoSamplesAvailable< minNSamples] = NA
        if (verbose > 0) pind = updateProgInd(end/nGenes, pind)
        start = end + 1
    }
    if (verbose > 0) printFlush("")
    k
} # end of function
