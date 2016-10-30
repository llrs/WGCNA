
# CorrelationPreservation ####
#' Preservation of eigengene correlations
#'
#' Calculates a summary measure of preservation of eigengene correlations
#' across data sets
#'
#' The function calculates the preservation of correlation of each eigengene
#' with all other eigengenes (optionally except the 'grey' eigengene) in all
#' pairs of sets.
#'
#' @param multiME consensus module eigengenes in a multi-set format. A vector
#' of lists with one list corresponding to each set. Each list must contain a
#' component \code{data} that is a data frame whose columns are consensus
#' module eigengenes.
#' @param setLabels names to be used for the sets represented in
#' \code{multiME}.
#' @param excludeGrey logical: exclude the 'grey' eigengene from preservation
#' measure?
#' @param greyLabel module label corresponding to the 'grey' module. Usually
#' this will be the character string \code{"grey"} if the labels are colors,
#' and the number 0 if the labels are numeric.
#' @return A data frame whose rows correspond to consensus module eigengenes
#' given in the input \code{multiME}, and columns correspond to all possible
#' set comparisons. The two sets compared in each column are indicated in the
#' column name.
#' @author Peter Langfelder
#' @seealso
#' \code{\link{multiSetMEs}} and module\code{\link{checkSets}} in
#' package moduleColor for more on eigengenes and the multi-set format
#' @references
#' Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
correlationPreservation <- function(multiME, setLabels, excludeGrey = TRUE,
                                    greyLabel = "grey") {
    nSets = length(multiME)
    if (nSets != length(setLabels))
        stop("The lengths of multiME and setLabels must equal.")
    if (nSets <= 1)
        stop("Something is wrong with argument multiME: its length is 0 or 1")
    Names = names(multiME[[1]]$data)
    if (excludeGrey) {
        Use = substring(Names, 3) != greyLabel
    } else {
        Use = rep(TRUE, times = length(Names))
    }
    No.Mods = ncol(multiME[[1]]$data[, Use])
    CP = matrix(0, nrow = No.Mods, ncol = nSets * (nSets - 1) / 2)
    diag(CP) = 1
    CPInd = 1
    CPNames = NULL
    for (i in 1:(nSets - 1))
        for (j in (i + 1):nSets) {
            corME1 = cor(multiME[[i]]$data[, Use], use = "p")
            corME2 = cor(multiME[[j]]$data[, Use], use = "p")
            d = 1 - abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2)) ^
                                 2))
            CP[, CPInd] = apply(d, 1, sum) - 1
            CPNames = c(CPNames,
                        paste(setLabels[i], "::", setLabels[j], collapse = ""))
            CPInd = CPInd + 1
        }
    CPx = as.data.frame(CP)
    names(CPx) = CPNames
    rownames(CPx) = Names[Use]
    CPx
}

# setCorrelationPreservation ####
#' Summary correlation preservation measure
#'
#' Given consensus eigengenes, the function calculates the average correlation
#' preservation pair-wise for all pairs of sets.
#'
#' For each pair of sets, the function calculates the average preservation of
#' correlation among the eigengenes. Two preservation measures are available,
#' the abosolute preservation (high if the two correlations are similar and low
#' if they are different), and the hyperbolically scaled preservation, which
#' de-emphasizes preservation of low correlation values.
#'
#' @param multiME consensus module eigengenes in a multi-set format. A vector
#' of lists with one list corresponding to each set. Each list must contain a
#' component \code{data} that is a data frame whose columns are consensus
#' module eigengenes.
#' @param setLabels names to be used for the sets represented in
#' \code{multiME}.
#' @param excludeGrey logical: exclude the 'grey' eigengene from preservation
#' measure?
#' @param greyLabel module label corresponding to the 'grey' module. Usually
#' this will be the character string \code{"grey"} if the labels are colors,
#' and the number 0 if the labels are numeric.
#' @param method character string giving the correlation preservation measure
#' to use. Recognized values are (unique abbreviations of) \code{"absolute"},
#' \code{"hyperbolic"}.
#' @return A data frame with each row and column corresponding to a set given
#' in \code{multiME}, containing the pairwise average correlation preservation
#' values. Names and rownames are set to entries of \code{setLabels}.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{multiSetMEs}} for module eigengene calculation;
#'
#' \code{\link{plotEigengeneNetworks}} for eigengene network visualization.
#' @references
#' Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
setCorrelationPreservation <- function(multiME, setLabels, excludeGrey = TRUE,
                                       greyLabel = "grey",
                                       method = "absolute") {
    m = charmatch(method, c("absolute", "hyperbolic"))
    if (is.na(m)) {
        stop("Unrecognized method given. Recognized methods are absolute,
             hyperbolic.")
    }
    nSets = length(multiME)
    if (nSets != length(setLabels))
        stop("The lengths of multiME and setLabels must equal.")
    if (nSets <= 1)
        stop("Something is wrong with argument multiME: its length is 0 or 1")
    Names = names(multiME[[1]]$data)
    if (excludeGrey) {
        Use = substring(Names, 3) != greyLabel
    } else {
        Use = rep(TRUE, times = length(Names))
    }
    No.Mods = ncol(multiME[[1]]$data[, Use])
    SCP = matrix(0, nrow = nSets, ncol = nSets)
    diag(SCP) = 0
    for (i in 1:(nSets - 1))
        for (j in (i + 1):nSets) {
            corME1 = cor(multiME[[i]]$data[, Use], use = "p")
            corME2 = cor(multiME[[j]]$data[, Use], use = "p")
            if (m == 1) {
                d = 1 - abs(corME1 - corME2) / 2
            } else {
                d = 1 - abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2)) ^ 2))
            }
            SCP[i, j] = sum(d[upper.tri(d)]) / sum(upper.tri(d))
            SCP[j, i] = SCP[i, j]
        }
    SCPx = as.data.frame(SCP)
    names(SCPx) = setLabels
    rownames(SCPx) = setLabels
    SCPx
    }

# preservationNetworkConnectivity ####
#' Network preservation calculations
#'
#' This function calculates several measures of gene network preservation.
#' Given gene expression data in several individual data sets, it calculates
#' the individual adjacency matrices, forms the preservation network and
#' finally forms several summary measures of adjacency preservation for each
#' node (gene) in the network.
#'
#' The preservation network is formed from adjacencies of compared sets. For
#' 'complete' preservations, all given sets are compared at once; for
#' 'pairwise' preservations, the sets are compared in pairs. Unweighted
#' preservations are simple mean preservations for each node; their weighted
#' counterparts are weighted averages in which a preservation of adjacencies
#' \eqn{A^{(1)}_{ij}}{A[i,j; 1]} and \eqn{A^{(2)}_{ij}}{A[i,j; 2]} of nodes
#' \eqn{i,j} between sets 1 and 2 is weighted by \eqn{[ (A^{(1)}_{ij} +
#' A^{(2)}_{ij} )/2]^weightPower}{ ( (A[i,j; 1]+A[i,j; 2])/2)^weightPower}. The
#' hyperbolic preservation is based on \eqn{tanh[( max -
#' min)/(max+min)^2]}{tanh[( max - min)/(max+min)^2]}, where \eqn{max}{max} and
#' \eqn{min}{min} are the componentwise maximum and minimum of the compared
#' adjacencies, respectively.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param useSets optional specification of sets to be used for the
#' preservation calculation. Defaults to using all sets.
#' @param useGenes optional specification of genes to be used for the
#' preservation calculation. Defaults to all genes.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param networkType a character string encoding network type. Recognized
#' values are (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored upon
#' exit.
#' @param weightPower power with which higher adjacencies will be weighted in
#' weighted means
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{pairwise}{ a matrix with rows corresponding to genes and columns to
#' unique pairs of given sets, giving the pairwise preservation of the
#' adjacencies connecting the gene to all other genes.}
#'
#' \item{complete}{ a vector with one entry for each input gene containing the
#' complete mean preservation of the adjacencies connecting the gene to all
#' other genes.}
#'
#' \item{pairwiseWeighted}{ a matrix with rows corresponding to genes and
#' columns to unique pairs of given sets, giving the pairwise weighted
#' preservation of the adjacencies connecting the gene to all other genes.}
#'
#' \item{completeWeighted}{ a vector with one entry for each input gene
#' containing the complete weighted mean preservation of the adjacencies
#' connecting the gene to all other genes.}
#'
#' \item{pairwiseHyperbolic}{ a matrix with rows corresponding to genes and
#' columns to unique pairs of given sets, giving the pairwise hyperbolic
#' preservation of the adjacencies connecting the gene to all other genes.}
#'
#' \item{completeHyperbolic}{ a vector with one entry for each input gene
#' containing the complete mean hyperbolic preservation of the adjacencies
#' connecting the gene to all other genes.}
#'
#' \item{pairwiseWeightedHyperbolic}{ a matrix with rows corresponding to genes
#' and columns to unique pairs of given sets, giving the pairwise weighted
#' hyperbolic preservation of the adjacencies connecting the gene to all other
#' genes.}
#'
#' \item{completeWeightedHyperbolic}{ a vector with one entry for each input
#' gene containing the complete weighted hyperbolic mean preservation of the
#' adjacencies connecting the gene to all other genes.}
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}} for calculation of adjacency;
#' @references
#' Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
preservationNetworkConnectivity <- function(multiExpr,
                                            useSets = NULL,
                                            useGenes = NULL,
                                            corFnc = "cor",
                                            corOptions = "use = 'p'",
                                            networkType = "unsigned",
                                            power = 6,
                                            sampleLinks = NULL,
                                            nLinks = 5000,
                                            blockSize = 1000,
                                            setSeed = 12345,
                                            weightPower = 2,
                                            verbose = 2,
                                            indent = 0) {
    spaces = indentSpaces(indent)

    size = checkSets(multiExpr)
    nGenes = size$nGenes
    nSets = size$nSets
    if (!is.null(useSets) || !is.null(useGenes)) {
        if (is.null(useSets))
            useSets = c(1:nSets)
        if (is.null(useGenes))
            useGenes = c(1:nGenes)
        useExpr = vector(mode = "list", length = length(useSets))
        for (set in 1:length(useSets))
            useExpr[[set]] = list(data = multiExpr[[useSets[set]]]$data[, useGenes])
        multiExpr = useExpr
        rm(useExpr)
        collectGarbage()
    }
    size = checkSets(multiExpr)
    nGenes = size$nGenes
    nSets = size$nSets

    if (is.null(sampleLinks)) {
        sampleLinks = (nGenes > nLinks)
    }

    if (sampleLinks)
        nLinks = min(nLinks, nGenes)
    else
        nLinks = nGenes

    if (blockSize * nLinks > .largestBlockSize)
        blockSize = as.integer(.largestBlockSize / nLinks)

    intNetworkType = charmatch(networkType, .networkTypes)
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
                    "preservationNetworkConnectivity:
                    selecting sample pool of size",
                    nLinks,
                    ".."
                )
            )
        sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE)
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
        subtract[-samplePool] = 0
    }

    nPairComps = nSets * (nSets  - 1) / 2

    allPres = rep(NA, nGenes)
    allPresW = rep(NA, nGenes)
    allPresH = rep(NA, nGenes)
    allPresWH = rep(NA, nGenes)

    pairPres = matrix(nGenes, nPairComps)
    pairPresW = matrix(nGenes, nPairComps)
    pairPresH = matrix(nGenes, nPairComps)
    pairPresWH = matrix(nGenes, nPairComps)

    compNames = NULL
    for (set1 in 1:(nSets - 1))
        for (set2 in (set1 + 1):nSets)
            compNames = c(compNames, paste(set1, "vs", set2))

    dimnames(pairPres) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresW) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresH) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresWH) = list(names(multiExpr[[1]]$data), compNames)

    if (verbose > 0)
    {
        pind = initProgInd(trailStr = " done")
    }

    nBlocks = as.integer((nGenes - 1) / blockSize)
    SetRestrConn = NULL
    start = 1
    if (sampleLinks)
    {
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

    while (start <= nGenes) {
        end = start + blockSize - 1
        if (end > nGenes)
            end = nGenes
        blockIndex = c(start:end)
        nBlockGenes = end - start + 1
        blockAdj = array(0, dim = c(nSets, nLinks, nBlockGenes))
        #if (verbose > 1) printFlush(paste(spaces, "..working on genes", start,
        #"through", end, "of", nGenes))
        for (set in 1:nSets) {
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
            adj_mat = as.matrix(c ^ power)
            if (sum(is.na(adj_mat)) > 0)
                stop(
                    "NA values present in adjacency - this function cannot
                    handle them yet. Sorry!"
                )
            adj_mat[is.na(adj_mat)] = 0
            blockAdj[set, ,] = adj_mat
        }
        min = matrix(0, nLinks, nBlockGenes)
        which = matrix(0, nLinks, nBlockGenes)
        res = .C(
            "minWhichMin",
            as.double(blockAdj),
            as.integer(nSets),
            as.integer(nLinks * nBlockGenes),
            min = as.double(min),
            as.double(which)
        )
        min[,] = res$min
        max = matrix(0, nLinks, nBlockGenes)
        res = .C(
            "minWhichMin",
            as.double(-blockAdj),
            as.integer(nSets),
            as.integer(nLinks * nBlockGenes),
            min = as.double(min),
            as.double(which)
        )
        max[,] = -res$min
        rm(res)
        diff = max - min
        allPres[blockIndex] = (apply(1 - diff, 2, sum) - subtract[blockIndex]) /
            (nLinks - subtract[blockIndex])
        weight = ((max + min) / 2) ^ weightPower
        allPresW[blockIndex] = (apply((1 - diff) * weight, 2, sum) -
                                    subtract[blockIndex]) /
            (apply(weight, 2, sum) - subtract[blockIndex])
        hyp = 1 - tanh(diff / (max + min) ^ 2)
        allPresH[blockIndex] = (apply(hyp, 2, sum) - subtract[blockIndex]) /
            (nLinks - subtract[blockIndex])
        allPresWH[blockIndex] = (apply(hyp * weight, 2, sum) -
                                     subtract[blockIndex]) /
            (apply(weight, 2, sum) - subtract[blockIndex])

        compNames = NULL
        compInd = 1
        for (set1 in 1:(nSets - 1))
            for (set2 in (set1 + 1):nSets) {
                diff = abs(blockAdj[set1, ,] - blockAdj[set2, ,])
                compNames = c(compNames, paste(set1, "vs", set2))
                pairPres[blockIndex, compInd] = (apply(1 - diff, 2, sum) - subtract[blockIndex]) /
                    (nLinks - subtract[blockIndex])
                weight = ((blockAdj[set1, ,] + blockAdj[set2, ,]) / 2) ^
                    weightPower
                pairPresW[blockIndex, compInd] = (apply((1 - diff) * weight, 2, sum) - subtract[blockIndex]) /
                    (apply(weight, 2, sum) - subtract[blockIndex])
                hyp = 1 - tanh(diff / (blockAdj[set1, ,] + blockAdj[set2, ,]) ^
                                   2)
                pairPresH[blockIndex, compInd] = (apply(hyp, 2, sum) - subtract[blockIndex]) /
                    (nLinks - subtract[blockIndex])
                pairPresWH[blockIndex, compInd] = (apply(hyp * weight, 2, sum) - subtract[blockIndex]) /
                    (apply(weight, 2, sum) - subtract[blockIndex])
                compInd = compInd + 1
            }

        start = end + 1
        if (verbose > 0)
            pind = updateProgInd(end / nGenes, pind)
        collectGarbage()
    }
    if (verbose > 0)
        printFlush(" ")
    list(
        pairwise = pairPres,
        complete = allPres,
        pairwiseWeighted = pairPresW,
        completeWeighted = allPresW,
        pairwiseHyperbolic = pairPresH,
        completeHyperbolic = allPresH,
        pairwiseWeightedHyperbolic = pairPresWH,
        completeWeightedHyperbolic = allPresWH
    )
}
