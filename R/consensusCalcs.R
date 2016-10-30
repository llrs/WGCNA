
# consensusOrderMEs ####
.clustOrder <- function(distM,
                        greyLast = TRUE,
                        greyName = paste0(moduleColor.getMEprefix(), "grey")) {
    distM = as.matrix(distM)
    distNames = dimnames(distM)[[1]]
    greyInd = match(greyName, distNames)
    if (greyLast && !is.na(greyInd)) {
        clusterMEs = (greyName != distNames)
        if (sum(clusterMEs) > 1) {
            h = fastcluster::hclust(as.dist(distM[clusterMEs, clusterMEs]),
                                    method = "average")
            order = h$order
            if (sum(order >= greyInd) > 0) {
                order[order >= greyInd] = order[order >= greyInd] + 1
            }
            order = c(order, greyInd)
        } else if (ncol(distM) > 1) {
            if (greyInd == 1) {
                order = c(2, 1)
            } else {
                order = c(1, 2)
            }
        } else {
            order = 1
        }
    } else {
        if (length(distM) > 1) {
            h = fastcluster::hclust(as.dist(distM), method = "average")
            order = h$order
        } else {
            order = 1
        }
    }
    order
}

#' Put close eigenvectors next to each other in several sets.
#'
#' Reorder given (eigen-)vectors such that similar ones (as measured by
#' correlation) are next to each other. This is a multi-set version of
#' \code{\link{orderMEs}}; the dissimilarity used can be of consensus type (for
#' each pair of eigenvectors the consensus dissimilarity is the maximum of
#' individual set dissimilarities over all sets) or of majority type (for each
#' pair of eigenvectors the consensus dissimilarity is the average of
#' individual set dissimilarities over all sets).
#'
#' Ordering module eigengenes is useful for plotting purposes. This function
#' calculates the consensus or majority dissimilarity of given eigengenes over
#' the sets specified by \code{useSets} (defaults to all sets). A hierarchical
#' dendrogram is calculated using the dissimilarity and the order given by the
#' dendrogram is used for the eigengenes in all other sets.
#'
#' @param MEs Module eigengenes of several sets in a multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, with each list corresponding to
#' one dataset and the module eigengenes in the component \code{data}, that is
#' \code{MEs[[set]]$data[sample, module]} is the expression of the eigengene of
#' module \code{module} in sample \code{sample} in dataset \code{set}. The
#' number of samples can be different between the sets, but the modules must be
#' the same.
#' @param useAbs Controls whether vector similarity should be given by absolute
#' value of correlation or plain correlation.
#' @param useSets Allows the user to specify for which sets the eigengene
#' ordering is to be performed.
#' @param greyLast Normally the color grey is reserved for unassigned genes;
#' hence the grey module is not a proper module and it is conventional to put
#' it last. If this is not desired, set the parameter to \code{FALSE}.
#' @param greyName Name of the grey module eigengene.
#' @param method A character string giving the method to be used calculating
#' the consensus dissimilarity. Allowed values are (abbreviations of)
#' \code{"consensus"} and \code{"majority"}. The consensus dissimilarity is
#' calculated as the maximum of given set dissimilarities for
#' \code{"consensus"} and as the average for \code{"majority"}.
#' @return A vector of lists of the same type as \code{MEs} containing the
#' re-ordered eigengenes.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso
#' \code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}},
#' \code{\link{orderMEs}}
#' @keywords misc
consensusOrderMEs <- function(MEs,
                              useAbs = FALSE,
                              useSets = NULL,
                              greyLast = TRUE,
                              greyName = paste0(moduleColor.getMEprefix(),
                                                "grey"),
                              method = "consensus") {
    # Debugging code:
    #printFlush("consensusOrderMEs:")
    #size = checkSets(MEs)
    #print(size)
    # end debuging code
    Diss = consensusMEDissimilarity(MEs,
                                    useAbs = useAbs,
                                    useSets = useSets,
                                    method = method)
    order = .clustOrder(Diss, greyLast, greyName)
    #print(order)
    orderMEs(
        MEs,
        greyLast = greyLast,
        greyName = greyName,
        order = order,
        useSets = useSets
    )
}

# consensusMEDissimilarity ####
#' Consensus dissimilarity of module eigengenes.
#'
#' Calculates consensus dissimilarity \code{(1-cor)} of given module eigengenes
#' relaized in several sets.
#'
#' This function calculates the individual set dissimilarities of the given
#' eigengenes in each set, then takes the (parallel) maximum or average over
#' all sets. For details on the structure of imput data, see
#' \code{\link{checkSets}}.
#'
#' @param MEs Module eigengenes of the same modules in several sets.
#' @param useAbs Controls whether absolute value of correlation should be used
#' instead of correlation in the calculation of dissimilarity.
#' @param useSets If the consensus is to include only a selection of the given
#' sets, this vector (or scalar in the case of a single set) can be used to
#' specify the selection. If \code{NULL}, all sets will be used.
#' @param method A character string giving the method to use. Allowed values
#' are (abbreviations of) \code{"consensus"} and \code{"majority"}. The
#' consensus dissimilarity is calculated as the minimum of given set
#' dissimilarities for \code{"consensus"} and as the average for
#' \code{"majority"}.
#' @return A dataframe containing the matrix of dissimilarities, with
#' \code{names} and \code{rownames} set appropriately.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{checkSets}}
#' @keywords misc
consensusMEDissimilarity <- function(MEs, useAbs = FALSE, useSets = NULL,
                                     method = "consensus") {
    methods = c("consensus", "majority")
    m = charmatch(method, methods)
    if (is.na(m))
        stop("Unrecognized method given. Recognized values are",
             paste(methods, collapse  = ", "))

    nSets = length(MEs)
    MEDiss = vector(mode = "list", length = nSets)
    if (is.null(useSets))
        useSets = c(1:nSets)
    for (set in useSets)
    {
        if (useAbs)
        {
            diss = 1 - abs(cor(MEs[[set]]$data, use = "p"))
        } else
        {
            diss = 1 - cor(MEs[[set]]$data, use = "p")
        }
        MEDiss[[set]] = list(Diss = diss)
    }

    for (set in useSets)
        if (set == useSets[1])
        {
            ConsDiss = MEDiss[[set]]$Diss
        } else {
            if (m == 1) {
                ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss)
            } else {
                ConsDiss = ConsDiss + MEDiss[[set]]$Diss
            }
        }

    if (m == 2)
        ConsDiss = ConsDiss / nSets

    ConsDiss = as.data.frame(ConsDiss)
    names(ConsDiss) = names(MEs[[useSets[1]]]$data)
    rownames(ConsDiss) = names(MEs[[useSets[1]]]$data)

    ConsDiss
}

# Quantile normalization
# normalize each column such that (column) quantiles are the same
# The final value for each quantile is the 'summaryType' of the corresponding
# quantiles across the columns

.equalizeQuantiles <- function(data, summaryType = c("median", "mean")) {
    summaryType = match.arg(summaryType)
    data.sorted = apply(data, 2, sort)

    if (summaryType == "median") {
        refSample = rowMedians(data.sorted, na.rm = TRUE)
    } else if (summaryType == "mean") {
        refSample = rowMeans(data.sorted, na.rm = TRUE)
    }
    ranks = round(colRanks(data, ties.method = "average", preserveShape = TRUE))
    out = refSample [ranks]
    dim(out) = dim(data)
    dimnames(out) = dimnames(data)

    out
}

.turnVectorIntoDist <- function(x, size, Diag, Upper) {
    attr(x, "Size") = size
    attr(x, "Diag") = FALSE
    attr(x, "Upper") = FALSE
    class(x) = c("dist", class(x))
    x
}

.turnDistVectorIntoMatrix <- function(x, size, Diag, Upper, diagValue) {
    mat = as.matrix(.turnVectorIntoDist(x, size, Diag, Upper))
    if (!Diag) {
        diag(mat) = diagValue
    }
    mat
}

# This function calculates consensus dissimilarity of module eigengenes

.consensusMEDissimilarity <- function(multiMEs,
                                      useSets = NULL,
                                      corFnc = cor,
                                      corOptions = list(use = 'p'),
                                      equalizeQuantiles = FALSE,
                                      quantileSummary = "mean",
                                      consensusQuantile = 0,
                                      useAbs = FALSE,
                                      greyMEname = "ME0") {
    nSets <- checkSets(multiMEs)$nSets
    init <- multiMEs[[1]]$data
    useMEs <- c(1:ncol(init))[names(init) != greyMEname]
    useNames <- names(init)[useMEs]
    nUseMEs <- length(useMEs)
    #  if (nUseMEs<2)
    #    stop("Something is wrong: there are two or more proper modules,
    #    but less than two proper",
    #         "eigengenes. Please check that the grey color label and module
    #         eigengene label",
    #         "are correct.")

    if (is.null(useSets)) {
        useSets <- c(1:nSets)
    }
    nUseSets <- length(useSets)
    MEDiss <- array(NA, dim = c(nUseMEs, nUseMEs, nUseSets))
    for (set in useSets) {
        corOptions$x <- multiMEs[[set]]$data[, useMEs]
        if (useAbs) {
            diss <- 1 - abs(do.call(corFnc, corOptions))
        } else {
            diss <- 1 - do.call(corFnc, corOptions)
        }
        MEDiss[, , set] <- diss
    }

    if (equalizeQuantiles) {
        distMat <- apply(MEDiss, 3, function(x) {
            as.numeric(as.dist(x))
        })
        dim(distMat) <- c(nUseMEs * (nUseMEs - 1) / 2, nUseSets)
        normalized <- .equalizeQuantiles(distMat, summaryType = quantileSummary)
        MEDiss <- apply(normalized, 2, .turnDistVectorIntoMatrix,
                        size = nUseMEs, Diag = FALSE, Upper = FALSE,
                        diagValue = 0)
    }

    ConsDiss <- apply(MEDiss, c(1:2), quantile, probs = 1 - consensusQuantile,
                      names = FALSE, na.rm = TRUE)
    colnames(ConsDiss) <- useNames
    rownames(ConsDiss) <- useNames
    ConsDiss
}
