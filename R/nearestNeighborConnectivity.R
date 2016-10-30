
# nearestNeighborConnectivity ####
#' Connectivity to a constant number of nearest neighbors
#'
#' Given expression data and basic network parameters, the function calculates
#' connectivity of each gene to a given number of nearest neighbors.
#'
#' Connectivity of gene \code{i} is the sum of adjacency strengths between gene
#' \code{i} and other genes; in this case we take the \code{nNeighbors} nodes
#' with the highest connection strength to gene \code{i}. The adjacency
#' strengths are calculated by correlating the given expression data using the
#' function supplied in \code{corFNC} and transforming them into adjacency
#' according to the given network \code{type} and \code{power}.
#'
#' @param datExpr a data frame containing expression data, with rows
#' corresponding to samples and columns to genes. Missing values are allowed
#' and will be ignored.
#' @param nNeighbors number of nearest neighbors to use.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param type a character string encoding network type. Recognized values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored upon
#' exit.
#' @param verbose integer controlling the level of verbosity. 0 means silent.
#' @param indent integer controlling indentation of output. Each unit above 0
#' adds two spaces.
#' @return A vector with one component for each gene containing the nearest
#' neighbor connectivity.
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}}, \code{\link{softConnectivity}}
#' @keywords misc
nearestNeighborConnectivity <- function(datExpr,
                                        nNeighbors = 50,
                                        power = 6,
                                        type = "unsigned",
                                        corFnc = "cor",
                                        corOptions = "use = 'p'",
                                        blockSize = 1000,
                                        sampleLinks = NULL,
                                        nLinks = 5000,
                                        setSeed = 38457,
                                        verbose = 1,
                                        indent = 0) {
    spaces = indentSpaces(indent)
    nGenes = dim(datExpr)[2]
    nSamples = dim(datExpr)[1]

    if (is.null(sampleLinks)) {
        sampleLinks = (nGenes > nLinks)
    }

    if (sampleLinks) {
        nLinks = min(nLinks, nGenes)
    } else {
        nLinks = nGenes
    }
    #printFlush(paste("blockSize  = ", blockSize))
    #printFlush(paste("nGenes  = ", nGenes))
    #printFlush(paste(".largestBlockSize  = ", .largestBlockSize))

    if (blockSize * nLinks > .largestBlockSize) {
        blockSize = as.integer(.largestBlockSize / nLinks)
    }
    intNetworkType = charmatch(type, .networkTypes)
    if (is.na(intNetworkType)) {
        stop("Unrecognized networkType argument. Recognized values are ",
             "(unique abbreviations of)",
             paste(.networkTypes, collapse = ", ")
        )
    }
    subtract = rep(1, nGenes)
    if (sampleLinks) {
        if (verbose > 0) {
            printFlush(
                paste(
                    spaces,
                    "nearestNeighborConnectivity: selecting sample pool
                    of size",
                    nLinks,
                    ".."
                )
            )
        }
        sd = apply(datExpr, 2, sd, na.rm = TRUE)
        order = order(-sd)
        saved = FALSE
        if (exists(".Random.seed")) {
            saved = TRUE
            savedSeed = .Random.seed
            if (is.numeric(setSeed))
                set.seed(setSeed)
        }
        samplePool = order[sample(x = nGenes, size = nLinks)]
        if (saved) {
            .Random.seed <<- savedSeed
        }
        poolExpr = datExpr[, samplePool]
        subtract[-samplePool] = 0
    }

    if (verbose > 0) {
        printFlush(
            paste(
                spaces,
                "nearestNeighborConnectivity: received",
                "dataset with nGenes  = ",
                as.character(nGenes)
            )
        )
        cat(
            paste(
                spaces,
                "..using nNeighbors  = ",
                nNeighbors,
                "and blockSize  = ",
                blockSize,
                "  "
            )
        )
        pind = initProgInd(trailStr = " done")
    }

    nearestNeighborConn = rep(0, nGenes)

    nBlocks = as.integer((nGenes - 1) / blockSize)
    SetRestrConn = NULL
    start = 1
    if (sampleLinks) {
        corEval = parse(text = paste(
            corFnc,
            "(poolExpr, datExpr[, blockIndex] ",
            prepComma(corOptions),
            ")"
        ))
    } else {
        corEval = parse(text = paste(
            corFnc,
            "(datExpr, datExpr[, blockIndex] ",
            prepComma(corOptions),
            ")"
        ))
    }

    while (start <= nGenes) {
        end = start + blockSize - 1
        if (end > nGenes) {
            end = nGenes
        }
        blockIndex = c(start:end)
        #if (verbose > 1) printFlush(paste(spaces, "..working on genes", start,
        #"through", end, "of", nGenes))
        c = eval(corEval)
        if (intNetworkType == 1) {
            c = abs(c)
        } else if (intNetworkType == 2) {
            c = (1 + c) / 2
        } else if (intNetworkType == 3) {
            c[c < 0] = 0
        } else {
            stop("Internal error: NetworkType has wrong value:",
                 intNetworkType, ". Sorry!"
            )
        }
        adj_mat = as.matrix(c ^ power)
        adj_mat[is.na(adj_mat)] = 0
        sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors + 1),])
        nearestNeighborConn[blockIndex] = apply(sortedAdj, 2, sum) - subtract[blockIndex]
        start = end + 1
        if (verbose > 0) {
            pind = updateProgInd(end / nGenes, pind)
        }
        collectGarbage()
    }
    if (verbose > 0) {
        printFlush(" ")
    }
    nearestNeighborConn
}

# nearestNeighbourConnectivityMS ####
#' Connectivity to a constant number of nearest neighbors across multiple data
#' sets
#'
#' Given expression data from several sets and basic network parameters, the
#' function calculates connectivity of each gene to a given number of nearest
#' neighbors in each set.
#'
#' Connectivity of gene \code{i} is the sum of adjacency strengths between gene
#' \code{i} and other genes; in this case we take the \code{nNeighbors} nodes
#' with the highest connection strength to gene \code{i}. The adjacency
#' strengths are calculated by correlating the given expression data using the
#' function supplied in \code{corFNC} and transforming them into adjacency
#' according to the given network \code{type} and \code{power}.
#'
#' @param multiExpr expression data in multi-set format. A vector of lists, one
#' list per set. In each list there must be a component named \code{data} whose
#' content is a matrix or dataframe or array of dimension 2 containing the
#' expression data. Rows correspond to samples and columns to genes (probes).
#' @param nNeighbors number of nearest neighbors to use.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param type a character string encoding network type. Recognized values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored after.
#' @param verbose integer controlling the level of verbosity. 0 means silent.
#' @param indent integer controlling indentation of output. Each unit above 0
#' adds two spaces.
#' @return A matrix in which columns correspond to sets and rows to genes; each
#' entry contains the nearest neighbor connectivity of the corresponding gene.
#' @author Peter Langfelder
#' @seealso
#' \code{\link{adjacency}}, \code{\link{softConnectivity}},
#' \code{\link{nearestNeighborConnectivity}}
#' @keywords misc
nearestNeighborConnectivityMS <- function(multiExpr,
                                          nNeighbors = 50,
                                          power = 6,
                                          type = "unsigned",
                                          corFnc = "cor",
                                          corOptions = "use = 'p'",
                                          blockSize = 1000,
                                          sampleLinks = NULL,
                                          nLinks = 5000,
                                          setSeed = 36492,
                                          verbose = 1,
                                          indent = 0) {
    spaces = indentSpaces(indent)
    setsize = checkSets(multiExpr)
    nGenes = setsize$nGenes
    nSamples = setsize$nSamples
    nSets = setsize$nSets

    if (is.null(sampleLinks)) {
        sampleLinks = (nGenes > nLinks)
    }

    if (sampleLinks)
        nLinks = min(nLinks, nGenes)
    else
        nLinks = nGenes

    #printFlush(paste("blockSize  = ", blockSize))
    #printFlush(paste("nGenes  = ", nGenes))
    #printFlush(paste(".largestBlockSize  = ", .largestBlockSize))

    if (blockSize * nLinks > .largestBlockSize)
        blockSize = as.integer(.largestBlockSize / nLinks)

    if (length(power) == 1)
    {
        power = rep(power, nSets)
    } else if (length(power) != nSets)
        stop("Invalid arguments: length of 'power' must equal number sets in
             'multiExpr'")

    intNetworkType = charmatch(type, .networkTypes)
    if (is.na(intNetworkType))
        stop(
            paste(
                "Unrecognized networkType argument. Recognized values are
                (unique abbreviations of)",
                paste(.networkTypes, collapse = ", ")
            )
        )

    subtract = rep(1, nGenes)
    if (sampleLinks) {
        if (verbose > 0)
            printFlush(
                paste(
                    spaces,
                    "nearestNeighborConnectivityMS: selecting
                    sample pool of size",
                    nLinks,
                    ".."
                )
            )
        sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE)
        order = order(-sd)
        saved = FALSE
        if (exists(".Random.seed"))
        {
            saved = TRUE
            savedSeed = .Random.seed
            if (is.numeric(setSeed))
                set.seed(setSeed)
        }
        samplePool = order[sample(x = nGenes, size = nLinks)]
        if (saved) {
            .Random.seed <<- savedSeed
        }
        subtract[-samplePool] = 0
    }

    if (verbose > 0)
        printFlush(
            paste(
                spaces,
                "nearestNeighborConnectivityMS: received",
                nSets,
                "datasets with nGenes  = ",
                as.character(nGenes)
            )
        )
    if (verbose > 0)
        printFlush(paste(spaces, "  Using nNeighbors  = ",
                         nNeighbors))

    nearestNeighborConn = matrix(nrow = nGenes, ncol = nSets)

    if (sampleLinks) {
        corEval = parse(
            text = paste(
                corFnc,
                "(multiExpr[[set]]$data[, samplePool],
                multiExpr[[set]]$data[, blockIndex] ",
                prepComma(corOptions),
                ")"
            )
        )
    } else {
        corEval = parse(
            text = paste(
                corFnc,
                "(multiExpr[[set]]$data,
                multiExpr[[set]]$data[, blockIndex] ",
                prepComma(corOptions),
                ")"
            )
        )
    }


    for (set in 1:nSets) {
        if (verbose > 0) {
            cat(paste(spaces, "  Working on set", set))
            pind = initProgInd(trailStr = " done")
        }
        nBlocks = as.integer((nGenes - 1) / blockSize)
        SetRestrConn = NULL
        start = 1
        while (start <= nGenes) {
            end = start + blockSize - 1
            if (end > nGenes)
                end = nGenes
            blockIndex = c(start:end)
            #if (verbose > 1) printFlush(paste(spaces, " .. working on genes",
            #start, "through", end, "of", nGenes))
            c = eval(corEval)
            if (intNetworkType == 1) {
                c = abs(c)
            } else if (intNetworkType == 2) {
                c = (1 + c) / 2
            } else if (intNetworkType == 3) {
                c[c < 0] = 0
            } else {
                stop(
                    "Internal error: intNetworkType has wrong value:",
                    intNetworkType,
                    ". Sorry!"
                )
            }
            adj_mat = as.matrix(c ^ power[set])
            adj_mat[is.na(adj_mat)] = 0
            sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors + 1),])
            nearestNeighborConn[blockIndex, set] = apply(sortedAdj, 2,
                                                         sum) -
                subtract[blockIndex]
            collectGarbage()
            start = end + 1
            if (verbose > 0)
                pind = updateProgInd(end / nGenes, pind)
            collectGarbage()
        }
        if (verbose > 0)
            printFlush(" ")
    }
    nearestNeighborConn
}
