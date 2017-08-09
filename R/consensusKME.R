# consensusKME ####
.interleave <- function(matrices, nameBase = names(matrices), sep = ".",
                        baseFirst = TRUE) {
    # Drop null entries in the list
    keep = sapply(matrices, function(x)
        ! is.null(x))
    nameBase = nameBase[keep]
    matrices = matrices[keep]

    nMats = length(matrices)
    nCols = ncol(matrices[[1]])

    dims = lapply(matrices, dim)

    if (baseFirst){
        for (m in 1:nMats) {
            colnames(matrices[[m]]) <- paste0(nameBase[m], sep,
                                              colnames(matrices[[m]]))
        }
    } else {
        for (m in 1:nMats)
            colnames(matrices[[m]]) <- paste0(colnames(matrices[[m]]),
                                              sep, nameBase[m])
    }

    out = as.data.frame(lapply(1:nCols,
                               function(index, matrices)
                                   as.data.frame(lapply(matrices,
                                                        function(x, i)
                                                            x[, i, drop = FALSE], index)),
                               matrices))

    rownames(out) = rownames(matrices[[1]])
    out
}



#' Calculate consensus kME (eigengene-based connectivities) across multiple
#' data sets.
#' 
#' Calculate consensus kME (eigengene-based connectivities) across multiple
#' data sets, typically following a consensus module analysis.
#' 
#' The function \code{corAndPvalueFnc} is currently is expected to accept
#' arguments \code{x} (gene expression profiles), \code{y} (eigengene
#' expression profiles), and \code{alternative} with possibilities at least
#' \code{"greater", "two.sided"}.  Any additional arguments can be passed via
#' \code{corOptions}.
#' 
#' The function \code{corAndPvalueFnc} should return a list which at the least
#' contains (1) a matrix of associations of genes and eigengenes (this
#' component should have the name given by \code{corComponent}), and (2) a
#' matrix of the corresponding p-values, named "p" or "p.value". Other
#' components are optional but for full functionality should include (3)
#' \code{nObs} giving the number of observations for each association (which is
#' the number of samples less number of missing data - this can in principle
#' vary from association to association), and (4) \code{Z} giving a Z static
#' for each observation. If these are missing, \code{nObs} is calculated in the
#' main function, and calculations using the Z statistic are skipped.
#' 
#' @param multiExpr Expression (or other numeric) data in a multi-set format. A
#' vector of lists; in each list there must be a component named `data' whose
#' content is a matrix or dataframe or array of dimension 2.
#' @param moduleLabels Module labels: one label for each gene in
#' \code{multiExpr}.
#' @param multiEigengenes Optional eigengenes of modules specified in
#' \code{moduleLabels}. If not given, will be calculated from \code{multiExpr}.
#' @param consensusQuantile Quantile for the consensus calculation. Should be a
#' number between 0 (minimum) and 1.
#' @param signed logical: should the network be considered signed? In signed
#' networks (\code{TRUE}), negative kME values are not considered significant
#' and the corresponding p-values will be one-sided. In unsigned networks
#' (\code{FALSE}), negative kME values are considered significant and the
#' corresponding p-values will be two-sided.
#' @param useModules Optional specification of module labels to which the
#' analysis should be restricted. This could be useful if there are many
#' modules, most of which are not interesting. Note that the "grey" module
#' cannot be used with \code{useModules}.
#' @param metaAnalysisWeights Optional specification of meta-analysis weights
#' for each input set. If given, must be a numeric vector of length equal the
#' number of input data sets (i.e., \code{length(multiExpr)}). These weights
#' will be used in addition to constant weights and weights proportional to
#' number of samples (observations) in each set.
#' @param corAndPvalueFnc Function that calculates associations between
#' expression profiles and eigengenes. See details.
#' @param corOptions List giving additional arguments to function
#' \code{corAndPvalueFnc}. See details.
#' @param corComponent Name of the component of output of
#' \code{corAndPvalueFnc} that contains the actual correlation.
#' @param getQvalues logical: should q-values (estimates of FDR) be calculated?
#' @param useRankPvalue Logical: should the \code{\link{rankPvalue}} function
#' be used to obtain alternative meta-analysis statistics?
#' @param rankPvalueOptions Additional options for function
#' \code{\link{rankPvalue}}. These include \code{na.last} (default
#' \code{"keep"}), \code{ties.method} (default \code{"average"}),
#' \code{calculateQvalue} (default copied from input \code{getQvalues}), and
#' \code{pValueMethod} (default \code{"scale"}). See the help file for
#' \code{\link{rankPvalue}} for full details.
#' @param setNames names for the input sets. If not given, will be taken from
#' \code{names(multiExpr)}. If those are \code{NULL} as well, the names will be
#' \code{"Set_1", "Set_2", ...}.
#' @param excludeGrey logical: should the grey module be excluded from the kME
#' tables? Since the grey module is typically not a real module, it makes
#' little sense to report kME values for it.
#' @param greyLabel label that labels the grey module.
#' @return Data frame with the following components (for easier readability the
#' order here is not the same as in the actual output): \item{ID}{Gene ID,
#' taken from the column names of the first input data set}
#' 
#' \item{consensus.kME.1, consensus.kME.2, ...}{Consensus kME (that is, the
#' requested quantile of the kMEs in the individual data sets)in each module
#' for each gene across the input data sets. The module labels (here 1, 2,
#' etc.) correspond to those in \code{moduleLabels}.}
#' 
#' \item{weightedAverage.equalWeights.kME1, weightedAverage.equalWeights.kME2,
#' }{ Average kME in each module for each gene across the input data sets.
#' }\item{...}{ Average kME in each module for each gene across the input data
#' sets. }
#' 
#' \item{weightedAverage.RootDoFWeights.kME1, }{ Weighted average kME in each
#' module for each gene across the input data sets. The weight of each data set
#' is proportional to the square root of the number of samples in the set.
#' }\item{weightedAverage.RootDoFWeights.kME2, ...}{ Weighted average kME in
#' each module for each gene across the input data sets. The weight of each
#' data set is proportional to the square root of the number of samples in the
#' set. }
#' 
#' \item{weightedAverage.DoFWeights.kME1, weightedAverage.DoFWeights.kME2, }{
#' Weighted average kME in each module for each gene across the input data
#' sets. The weight of each data set is proportional to number of samples in
#' the set. }\item{...}{ Weighted average kME in each module for each gene
#' across the input data sets. The weight of each data set is proportional to
#' number of samples in the set. }
#' 
#' \item{weightedAverage.userWeights.kME1, weightedAverage.userWeights.kME2, }{
#' (Only present if input \code{metaAnalysisWeights} is non-NULL.) Weighted
#' average kME in each module for each gene across the input data sets. The
#' weight of each data set is given in \code{metaAnalysisWeights}.}\item{...}{
#' (Only present if input \code{metaAnalysisWeights} is non-NULL.) Weighted
#' average kME in each module for each gene across the input data sets. The
#' weight of each data set is given in \code{metaAnalysisWeights}.}
#' 
#' \item{meta.Z.equalWeights.kME1, meta.Z.equalWeights.kME2, ...}{Meta-analysis
#' Z statistic for kME in each module, obtained by weighing the Z scores in
#' each set equally. Only returned if the function \code{corAndPvalueFnc}
#' returns the Z statistics corresponding to the correlations.}
#' 
#' \item{meta.Z.RootDoFWeights.kME1, meta.Z.RootDoFWeights.kME2, ...}{
#' Meta-analysis Z statistic for kME in each module, obtained by weighing the Z
#' scores in each set by the square root of the number of samples. Only
#' returned if the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the correlations.}
#' 
#' \item{meta.Z.DoFWeights.kME1, meta.Z.DoFWeights.kME2, ...}{Meta-analysis Z
#' statistic for kME in each module, obtained by weighing the Z scores in each
#' set by the number of samples. Only returned if the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations.}
#' 
#' \item{meta.Z.userWeights.kME1, meta.Z.userWeights.kME2, ...}{Meta-analysis Z
#' statistic for kME in each module, obtained by weighing the Z scores in each
#' set by \code{metaAnalysisWeights}.  Only returned if
#' \code{metaAnalysisWeights} is non-NULL and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations.}
#' 
#' \item{meta.p.equalWeights.kME1, meta.p.equalWeights.kME2, ...}{ p-values
#' obtained from the equal-weight meta-analysis Z statistics. Only returned if
#' the function \code{corAndPvalueFnc} returns the Z statistics corresponding
#' to the correlations. }
#' 
#' \item{meta.p.RootDoFWeights.kME1, meta.p.RootDoFWeights.kME2, ...}{ p-values
#' obtained from the meta-analysis Z statistics with weights proportional to
#' the square root of the number of samples. Only returned if the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations. }
#' 
#' \item{meta.p.DoFWeights.kME1, meta.p.DoFWeights.kME2, ...}{ p-values
#' obtained from the degree-of-freedom weight meta-analysis Z statistics. Only
#' returned if the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the correlations. }
#' 
#' \item{meta.p.userWeights.kME1, meta.p.userWeights.kME2, ...}{ p-values
#' obtained from the user-supplied weight meta-analysis Z statistics. Only
#' returned if \code{metaAnalysisWeights} is non-NULL and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations. }
#' 
#' \item{meta.q.equalWeights.kME1, meta.q.equalWeights.kME2, ...}{ q-values
#' obtained from the equal-weight meta-analysis p-values. Only present if
#' \code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc}
#' returns the Z statistics corresponding to the kME values.}
#' 
#' \item{meta.q.RootDoFWeights.kME1, meta.q.RootDoFWeights.kME2, ...}{ q-values
#' obtained from the meta-analysis p-values with weights proportional to the
#' square root of the number of samples. Only present if \code{getQvalues} is
#' \code{TRUE} and the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the kME values.}
#' 
#' \item{meta.q.DoFWeights.kME1, meta.q.DoFWeights.kME2, ...}{ q-values
#' obtained from the degree-of-freedom weight meta-analysis p-values. Only
#' present if \code{getQvalues} is \code{TRUE} and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the kME
#' values.}
#' 
#' \item{meta.q.userWeights.kME1, meta.q.userWeights.kME2, ...}{ q-values
#' obtained from the user-specified weight meta-analysis p-values. Only present
#' if \code{metaAnalysisWeights} is non-NULL, \code{getQvalues} is \code{TRUE}
#' and the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the kME values.}
#' 
#' The next set of columns contain the results of function
#' \code{\link{rankPvalue}} and are only present if input \code{useRankPvalue}
#' is \code{TRUE}. Some columns may be missing depending on the options
#' specified in \code{rankPvalueOptions}. We explicitly list columns that are
#' based on weighing each set equally; names of these columns carry the suffix
#' \code{.equalWeights}
#' 
#' \item{pValueExtremeRank.ME1.equalWeights, }{ This is the minimum between
#' pValueLowRank and pValueHighRank, i.e. min(pValueLow,
#' pValueHigh)}\item{pValueExtremeRank.ME2.equalWeights, ...}{ This is the
#' minimum between pValueLowRank and pValueHighRank, i.e. min(pValueLow,
#' pValueHigh)}
#' 
#' \item{pValueLowRank.ME1.equalWeights, pValueLowRank.ME2.equalWeights, ...}{
#' Asymptotic p-value for observing a consistently low value across the columns
#' of datS based on the rank method.}
#' 
#' \item{pValueHighRank.ME1.equalWeights, pValueHighRank.ME2.equalWeights, }{
#' Asymptotic p-value for observing a consistently low value across the columns
#' of datS based on the rank method.}\item{...}{ Asymptotic p-value for
#' observing a consistently low value across the columns of datS based on the
#' rank method.}
#' 
#' \item{pValueExtremeScale.ME1.equalWeights, }{ This is the minimum between
#' pValueLowScale and pValueHighScale, i.e. min(pValueLow,
#' pValueHigh)}\item{pValueExtremeScale.ME2.equalWeights, ...}{ This is the
#' minimum between pValueLowScale and pValueHighScale, i.e. min(pValueLow,
#' pValueHigh)}
#' 
#' \item{pValueLowScale.ME1.equalWeights, pValueLowScale.ME2.equalWeights, }{
#' Asymptotic p-value for observing a consistently low value across the columns
#' of datS based on the Scale method.}\item{...}{ Asymptotic p-value for
#' observing a consistently low value across the columns of datS based on the
#' Scale method.}
#' 
#' \item{pValueHighScale.ME1.equalWeights, pValueHighScale.ME2.equalWeights, }{
#' Asymptotic p-value for observing a consistently low value across the columns
#' of datS based on the Scale method.}\item{...}{ Asymptotic p-value for
#' observing a consistently low value across the columns of datS based on the
#' Scale method.}
#' 
#' \item{qValueExtremeRank.ME1.equalWeights, }{ local false discovery rate
#' (q-value) corresponding to the p-value
#' pValueExtremeRank}\item{qValueExtremeRank.ME2.equalWeights, ...}{ local
#' false discovery rate (q-value) corresponding to the p-value
#' pValueExtremeRank}
#' 
#' \item{qValueLowRank.ME1.equalWeights, qValueLowRank.ME2.equalWeights, ...}{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueLowRank}
#' 
#' \item{qValueHighRank.ME1.equalWeights, lueHighRank.ME2.equalWeights, ...}{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueHighRank}
#' 
#' \item{qValueExtremeScale.ME1.equalWeights, }{ local false discovery rate
#' (q-value) corresponding to the p-value
#' pValueExtremeScale}\item{qValueExtremeScale.ME2.equalWeights, ...}{ local
#' false discovery rate (q-value) corresponding to the p-value
#' pValueExtremeScale}
#' 
#' \item{qValueLowScale.ME1.equalWeights, qValueLowScale.ME2.equalWeights, }{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueLowScale}\item{...}{ local false discovery rate (q-value)
#' corresponding to the p-value pValueLowScale}
#' 
#' \item{qValueHighScale.ME1.equalWeights,qValueHighScale.ME2.equalWeights, }{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueHighScale}\item{...}{ local false discovery rate (q-value)
#' corresponding to the p-value pValueHighScale}
#' 
#' \item{...}{Analogous columns corresponding to weighing individual sets by
#' the square root of the number of samples, by number of samples, and by user
#' weights (if given). The corresponding column name suffixes are
#' \code{.RootDoFWeights}, \code{.DoFWeights}, and \code{.userWeights}.}
#' 
#' The following set of columns summarize kME in individual input data sets.
#' 
#' \item{kME1.Set_1, kME1.Set_2, ..., kME2.Set_1, kME2.Set_2, ...}{ kME values
#' for each gene in each module in each given data set. }
#' 
#' \item{p.kME1.Set_1, p.kME1.Set_2, ..., p.kME2.Set_1, p.kME2.Set_2, ...}{
#' p-values corresponding to kME values for each gene in each module in each
#' given data set. }
#' 
#' \item{q.kME1.Set_1, q.kME1.Set_2, ..., q.kME2.Set_1, q.kME2.Set_2, ...}{
#' q-values corresponding to kME values for each gene in each module in each
#' given data set. Only returned if \code{getQvalues} is \code{TRUE}. }
#' 
#' \item{Z.kME1.Set_1, Z.kME1.Set_2, ..., Z.kME2.Set_1, Z.kME2.Set_2, ...}{ Z
#' statistics corresponding to kME values for each gene in each module in each
#' given data set. Only present if the function \code{corAndPvalueFnc} returns
#' the Z statistics corresponding to the kME values. }
#' @author Peter Langfelder
#' @seealso \link{signedKME} for eigengene based connectivity in a single data
#' set. \link{corAndPvalue}, \link{bicorAndPvalue} for two alternatives for
#' calculating correlations and the corresponding p-values and Z scores. Both
#' can be used with this function.
#' @references Langfelder P, Horvath S., WGCNA: an R package for weighted
#' correlation network analysis. BMC Bioinformatics. 2008 Dec 29; 9:559.
#' @keywords misc
consensusKME <-function(multiExpr,
                        moduleLabels,
                        multiEigengenes = NULL,
                        consensusQuantile = 0,
                        signed = TRUE,
                        useModules = NULL,
                        metaAnalysisWeights = NULL,
                        corAndPvalueFnc = corAndPvalue,
                        corOptions = list(),
                        corComponent = "cor",
                        getQvalues = FALSE,
                        useRankPvalue = TRUE,
                        rankPvalueOptions = list(calculateQvalue = getQvalues,
                                                 pValueMethod = "scale"),
                        setNames = NULL,
                        excludeGrey = TRUE,
                        greyLabel = if (is.numeric(moduleLabels))
                            0
                        else
                            "grey") {
    corAndPvalueFnc = match.fun(corAndPvalueFnc)

    size = checkSets(multiExpr)
    nSets = size$nSets
    nGenes = size$nGenes
    nSamples = size$nSamples

    if (!is.null(metaAnalysisWeights))
        if (length(metaAnalysisWeights) != nSets)
            stop("Length of 'metaAnalysisWeights' must equal number of",
                 " input sets.")

    if (!is.null(useModules))
    {
        if (greyLabel %in% useModules)
            stop(
                paste(
                    "Grey module (or module 0) cannot be used with ",
                    "'useModules'.\n Use 'excludeGrey = FALSE' to obtain ",
                    "results for the grey module as well."
                )
            )
        keep = moduleLabels %in% useModules
        if (sum(keep) == 0)
            stop("Incorrectly specified 'useModules': no such module(s).")
        moduleLabels [!keep] = greyLabel
    }

    if (is.null(multiEigengenes))
        multiEigengenes = multiSetMEs(
            multiExpr,
            universalColors = moduleLabels,
            verbose = 0,
            excludeGrey = excludeGrey,
            grey = greyLabel
        )

    modLevels = substring(colnames(multiEigengenes[[1]]$data), 3)
    nModules = length(modLevels)

    kME = p = Z = nObs = array(NA, dim = c(nGenes, nModules, nSets))

    corOptions$alternative = c("two.sided", "greater")[signed + 1]

    haveZs = FALSE
    for (set in 1:nSets) {
        corOptions$x = multiExpr[[set]]$data
        corOptions$y = multiEigengenes[[set]]$data
        cp = do.call(corAndPvalueFnc, args = corOptions)
        corComp = grep(corComponent, names(cp))
        pComp = match("p", names(cp))
        if (is.na(pComp))
            pComp = match("p.value", names(cp))
        if (is.na(pComp))
            stop("Function `corAndPvalueFnc' did not return a p-value.")
        kME[, , set] = cp[[corComp]]
        p[, , set] = cp[[pComp]]
        if (!is.null(cp$Z)) {
            Z[, , set] = cp$Z
            haveZs = TRUE
        }
        if (!is.null(cp$nObs)) {
            nObs[, , set] = cp$nObs
        } else
            nObs[, , set] = t(is.na(multiExpr[[set]]$data)) %*% (!is.na(multiEigengenes[[set]]$data))
    }

    if (getQvalues)
    {
        q = apply(p, c(2:3), qvalue.restricted)
    } else
        q = NULL
    # not neccessary since weighted average also contains it
    # kME.average = rowMeans(kME, dims = 2)

    powers = c(0, 0.5, 1)
    nPowers = length(powers)
    nWeights = nPowers+!is.null(metaAnalysisWeights)
    weightNames = c("equalWeights",
                    "RootDoFWeights",
                    "DoFWeights",
                    "userWeights")[1:nWeights]
    kME.weightedAverage = array(NA, dim = c(nGenes, nWeights, nModules))
    for (m in 1:nWeights) {
        if (m <= nPowers) {
            weights = nObs ^ powers[m]
        } else
            weights = array(rep(metaAnalysisWeights, rep(nGenes * nModules,
                                                         nSets)),
                            dim = c(nGenes, nModules, nSets))
        kME.weightedAverage[, m,] = rowSums(kME * weights, na.rm = TRUE,
                                            dims = 2) /
            rowSums(weights, dims = 2, na.rm = TRUE)
    }

    dim(kME.weightedAverage) = c(nGenes * nWeights, nModules)

    if (any(is.na(kME))) {
        kME.consensus.1 = apply(kME,
                                c(1, 2),
                                quantile,
                                prob = consensusQuantile,
                                na.rm = TRUE)
        kME.consensus.2 = apply(kME,
                                c(1, 2),
                                quantile,
                                prob = 1 - consensusQuantile,
                                na.rm = TRUE)
        kME.median = apply(kME, c(1, 2), median, na.rm = TRUE)
    } else {
        kME.consensus.1 = matrix(colQuantileC(t(
            matrix(kME, nGenes * nModules, nSets)
        ),
        p = consensusQuantile),
        nGenes,
        nModules)
        kME.consensus.2 = matrix(colQuantileC(t(
            matrix(kME, nGenes * nModules, nSets)
        ),
        p = 1 - consensusQuantile),
        nGenes,
        nModules)
        kME.median = matrix(colQuantileC(t(
            matrix(kME, nGenes * nModules, nSets)
        ), p = 0.5),
        nGenes, nModules)
    }
    kME.consensus = ifelse(kME.median > 0, kME.consensus.1, kME.consensus.2)

    kME.consensus[kME.consensus * kME.median < 0] = 0

    # Prepare identifiers for the variables (genes)
    if (is.null(colnames(multiExpr[[1]]$data)))
    {
        ID = paste0("Variable.", 1:nGenes)
    } else
        ID = colnames(multiExpr[[1]]$data)

    # Get meta - Z, - p, - q values
    if (haveZs)
    {
        Z.kME.meta = p.kME.meta = array(0, dim = c(nGenes, nWeights, nModules))
        if (getQvalues)
            q.kME.meta = array(0, dim = c(nGenes, nWeights, nModules))
        for (m in 1:nWeights) {
            if (m <= nPowers) {
                weights = nObs ^ powers[m]
            } else
                weights = array(rep(metaAnalysisWeights,
                                    rep(nGenes * nModules, nSets)),
                                dim = c(nGenes, nModules, nSets))

            Z1 = rowSums(Z * weights, na.rm = TRUE, dims = 2) / sqrt(rowSums(weights ^
                                                                                 2, na.rm = TRUE, dims = 2))
            if (signed) {
                p1 = pnorm(Z1, lower.tail = FALSE)
            } else
                p1 = 2 * pnorm(abs(Z1), lower.tail = FALSE)
            Z.kME.meta[, m,] = Z1
            p.kME.meta[, m,] = p1
            if (getQvalues) {
                q1 = apply(p1, 2, qvalue.restricted)
                q.kME.meta[, m,] = q1
            }
        }
        dim(Z.kME.meta) = dim(p.kME.meta) = c(nGenes *  nWeights, nModules)
        if (getQvalues) {
            dim(q.kME.meta) = c(nGenes * nWeights, nModules)
        } else
            q.kME.meta = NULL
    } else {
        Z.kME.meta = p.kME.meta = q.kME.meta = NULL
    }

    # Call rankPvalue

    if (useRankPvalue) {
        for (mod in 1:nModules)
            for (m in 1:nWeights)
            {
                if (m <= nPowers) {
                    weights = nObs[, mod,] ^ powers[m]
                } else
                    weights = matrix(metaAnalysisWeights, nGenes, nSets,
                                     byrow = TRUE)
                # rankPvalue requires a vector of weights... so compress the weights to a vector.
                # Output a warning if the compression loses information.
                nDifferent = apply(weights, 2, function(x) {
                    length(unique(x))
                })
                if (any(nDifferent) > 1)
                    printFlush(
                        paste(
                            "Warning in consensusKME: rankPvalue requires compressed weights.\n",
                            "Some weights may not be entirely accurate."
                        )
                    )
                cw = colMeans(weights, na.rm = TRUE)
                rankPvalueOptions$columnweights = cw / sum(cw)

                rankPvalueOptions$datS = kME[, mod,]
                rp1 = do.call(rankPvalue, rankPvalueOptions)
                colnames(rp1) = paste0(colnames(rp1),
                                       ".ME",
                                       modLevels[mod],
                                       ".",
                                       weightNames[m])
                if (mod == 1 && m == 1) {
                    rp = rp1
                } else
                    rp = cbind(rp, rp1)
            }
    }

    # Format the output... this will entail some rearranging of the individual set results.
    if (is.null(setNames))
        setNames = names(multiExpr)

    if (is.null(setNames))
        setNames = paste0("Set_", c(1:nSets))

    if (!haveZs)
        Z = NULL

    keep = c(TRUE, TRUE, getQvalues, haveZs)
    varNames = c("kME", "p.kME", "q.kME", "Z.kME")[keep]
    nVars = sum(keep)

    dimnames(kME) = list(colnames(multiExpr),
                         paste0("k", colnames(multiEigengenes)),
                         setNames)

    dimnames(p) = list(colnames(multiExpr),
                       paste0("p.k", colnames(multiEigengenes)),
                       setNames)

    if (getQvalues)
        dimnames(q) = list(colnames(multiExpr),
                           paste0("q.k", colnames(multiEigengenes)),
                           setNames)

    if (haveZs)
        dimnames(Z) = list(colnames(multiExpr),
                           paste0("Z.k", colnames(multiEigengenes)),
                           setNames)


    varList = list(
        kME = kME,
        p = p,
        q = if (getQvalues)
            q
        else
            NULL,
        Z = if (haveZs)
            Z
        else
            NULL
    )
    varList.interleaved = lapply(varList,
                                 function(arr) {
                                     if (!is.null(dim(arr))) {
                                         split = lapply(1:dim(arr)[3], function(i)
                                             arr[, , i])
                                         .interleave(split, nameBase = setNames,
                                                     baseFirst = FALSE)
                                     } else
                                         NULL
                                 })

    # the following seems to choke on larger data sets, at least in R 3.2.1
    # combined = array(c (kME, p, q, Z), dim = c(nGenes, nModules, nSets, nVars))
    # recast = matrix(c(cast(melt(combined), X1~X4~X3~X2)), nGenes,
    # nSets * nModules * nVars)

    # ... so I will replace it with more cumbersome but hopefully workable code.

    recast = .interleave(varList.interleaved,
                         nameBase = rep("", 4),
                         sep = "")

    combinedMeta.0 = rbind(kME.consensus,
                           kME.weightedAverage,
                           Z.kME.meta,
                           p.kME.meta,
                           q.kME.meta)

    combinedMeta = matrix(combinedMeta.0, nGenes,
                          (1 + nWeights + (2 * haveZs + haveZs * getQvalues) *
                               nWeights) * nModules)
    metaNames = c(
        "consensus.kME",
        paste0("weightedAverage.", weightNames, ".kME"),
        paste0("meta.Z.", weightNames, ".kME"),
        paste0("meta.p.", weightNames, ".kME"),
        paste0("meta.q.", weightNames, ".kME")
    )[c(
        rep(TRUE, nWeights + 1),
        rep(haveZs, nWeights),
        rep(haveZs, nWeights),
        rep(haveZs && getQvalues, nWeights)
    )]
    nMetaVars = length(metaNames)
    colnames(combinedMeta) = paste0 (rep(metaNames, nModules),
                                     rep(modLevels, rep(nMetaVars, nModules)))

    if (useRankPvalue) {
        out = data.frame(ID = ID, combinedMeta, rp, recast)
    } else
        out = data.frame(ID = ID, combinedMeta, recast)

    out
}

# Meta - analysis
.isBinary <- function(multiTrait) {
    bin = TRUE
    for (set in 1:length(multiTrait))
        if (length(sort(unique(multiTrait[[set]]$data))) > 2)
            bin = FALSE
        bin
}
