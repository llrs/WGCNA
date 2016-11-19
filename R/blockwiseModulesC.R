#  TOM similarity via a call to a compiled code.


#' Topological overlap matrix
#'
#' Calculation of the topological overlap matrix from given expression data.
#'
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. NAs are allowed, but not too many.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pairwise.complete.obs} option.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param power soft-thresholding power for netwoek construction.
#' @param TOMType one of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If
#' \code{"none"}, adjacency will be used for clustering. If \code{"unsigned"},
#' the standard TOM will be used (more generally, TOM function will receive the
#' adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
#' correlations between neighbors.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are \code{"min"} giving the standard TOM described in
#' Zhang and Horvath (2005), and \code{"mean"} in which the \code{min} function
#' in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
#' better results but at this time should be considered experimental.
#' @param maxPOutliers only used for \code{corType=="bicor"}. Specifies the
#' maximum percentile of data that can be considered outliers on either side of
#' the median separately. For each side of the median, if higher percentile
#' than \code{maxPOutliers} is considered an outlier by the weight function
#' based on \code{9*mad(x)}, the width of the weight function is increased such
#' that the percentile of outliers on that side of the median equals
#' \code{maxPOutliers}. Using \code{maxPOutliers=1} will effectively disable
#' all weight function broadening; using \code{maxPOutliers=0} will give
#' results that are quite similar (but not equal to) Pearson correlation.
#' @param quickCor real number between 0 and 1 that controls the handling of
#' missing data in the calculation of correlations. See details.
#' @param pearsonFallback Specifies whether the bicor calculation, if used,
#' should revert to Pearson when median absolute deviation (mad) is zero.
#' Recongnized values are (abbreviations of) \code{"none", "individual",
#' "all"}. If set to \code{"none"}, zero mad will result in \code{NA} for the
#' corresponding correlation.  If set to \code{"individual"}, Pearson
#' calculation will be used only for columns that have zero mad.  If set to
#' \code{"all"}, the presence of a single zero mad will cause the whole
#' variable to be treated in Pearson correlation manner (as if the
#' corresponding \code{robust} option was set to \code{FALSE}). Has no effect
#' for Pearson correlation. See \code{\link{bicor}}.
#' @param cosineCorrelation logical: should the cosine version of the
#' correlation calculation be used? The cosine calculation differs from the
#' standard one in that it does not subtract the mean.
#' @param replaceMissingAdjacencies logical: should missing values in the
#' calculation of adjacency be replaced by 0?
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix holding the topological overlap.
#' @author Peter Langfelder
#' @seealso \code{\link{TOMsimilarity}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
#' @export
TOMsimilarityFromExpr = function(datExpr, corType = "pearson",
                                 networkType = "unsigned",
                                 power = 6, TOMType = "signed", TOMDenom = "min",
                                 maxPOutliers = 1,
                                 quickCor = 0,
                                 pearsonFallback = "individual",
                                 cosineCorrelation = FALSE,
                                 replaceMissingAdjacencies = FALSE,
                                 nThreads = 0,
                                 verbose = 1, indent = 0)
{
  corTypeC = as.integer(pmatch(corType, .corTypes)-1)
  if (is.na(corTypeC))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  TOMTypeC = as.integer(pmatch(TOMType, .TOMTypes)-1)
  if (is.na(TOMTypeC)) {
    stop(paste("Invalid 'TOMType'. Recognized values are",
               paste(.TOMTypes, collapse = ", ")))
  }

  TOMDenomC = as.integer(pmatch(TOMDenom, .TOMDenoms)-1)
  if (is.na(TOMDenomC)) {
    stop(paste("Invalid 'TOMDenom'. Recognized values are",
               paste(.TOMDenoms, collapse = ", ")))
  }
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) {
      stop("maxPOutliers must be between 0 and 1.")
  }
  if (quickCor < 0) {
      stop("quickCor must be positive.")
  }
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) {
      stop("maxPOutliers must be between 0 and 1.")
  }

  fallback = as.integer(pmatch(pearsonFallback, .pearsonFallbacks))
  if (is.na(fallback))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))

  if (nThreads < 0) stop("nThreads must be positive.")
  if (is.null(nThreads) || (nThreads==0)) nThreads = .useNThreads()

  if ( (power<1) | (power>30) ) stop("power must be between 1 and 30.")

  networkTypeC = as.integer(charmatch(networkType, .networkTypes)-1)
  if (is.na(networkTypeC))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))
  dimEx = dim(datExpr)
  if (length(dimEx)!=2) stop("datExpr has incorrect dimensions.")
  nGenes = dimEx[2]
  nSamples = dimEx[1]
  warn = as.integer(0)

  datExpr = as.matrix(datExpr)

  tom = .Call("tomSimilarity_call", datExpr,
        as.integer(corTypeC), as.integer(networkTypeC), as.double(power),
        as.integer(TOMTypeC), as.integer(TOMDenomC),
        as.double(maxPOutliers), as.double(quickCor),
        as.integer(fallback), as.integer(cosineCorrelation),
        as.integer(replaceMissingAdjacencies),
        warn,
        as.integer(nThreads), as.integer(verbose), as.integer(indent), PACKAGE = "WGCNA")

  diag(tom) = 1
  return (tom)
}



#
# TOMsimilarity (from adjacency)


#' Topological overlap matrix similarity and dissimilarity
#'
#' Calculation of the topological overlap matrix from a given adjacency matrix.
#'
#'
#' The functions perform basically the same calculations of topological
#' overlap. \code{TOMdist} turns the overlap (which is a measure of similarity)
#' into a measure of dissimilarity by subtracting it from 1.
#'
#' Basic checks on the adjacency matrix are performed and missing entries are
#' replaced by zeros.  If \code{TOMType = "unsigned"}, entries of the adjacency
#' matrix are required to lie between 0 and 1; for \code{TOMType = "signed"}
#' they can be between -1 and 1. In both cases the resulting TOM entries, as
#' well as the corresponding dissimilarities, lie between 0 and 1.
#'
#' The underlying C code assumes that the diagonal of the adjacency matrix
#' equals 1. If this is not the case, the diagonal of the input is set to 1
#' before the calculation begins.
#'
#' @aliases TOMsimilarity TOMdist
#' @param adjMat adjacency matrix, that is a square, symmetric matrix with
#' entries between 0 and 1 (negative values are allowed if
#' \code{TOMType=="signed"}).
#' @param TOMType a character string specifying TOM type to be calculated. One
#' of \code{"unsigned"}, \code{"signed"}. If \code{"unsigned"}, the standard
#' TOM will be used (more generally, TOM function will receive the adjacency as
#' input). If \code{"signed"}, TOM will keep track of the sign of the adjacency
#' between neighbors.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are \code{"min"} giving the standard TOM described in
#' Zhang and Horvath (2005), and \code{"mean"} in which the \code{min} function
#' in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
#' better results but at this time should be considered experimental.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix holding the topological overlap.
#' @author Peter Langfelder
#' @seealso \code{\link{TOMsimilarityFromExpr}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
TOMsimilarity = function(adjMat, TOMType = "unsigned", TOMDenom = "min",
                         verbose = 1, indent = 0) {
  TOMTypeC = pmatch(TOMType, .TOMTypes)-1
  if (is.na(TOMTypeC)) {
    stop(paste("Invalid 'TOMType'. Recognized values are",
               paste(.TOMTypes, collapse = ", ")))
  }
  if (TOMTypeC == 0) {
    stop("'TOMType' cannot be 'none' for this function.")
  }
  TOMDenomC = pmatch(TOMDenom, .TOMDenoms)-1
  if (is.na(TOMDenomC)) {
    stop(paste("Invalid 'TOMDenom'. Recognized values are",
               paste(.TOMDenoms, collapse = ", ")))
  }
  checkAdjMat(adjMat, min = if (TOMTypeC==2) -1 else 0, max = 1)

  if (sum(is.na(adjMat))>0) {
      adjMat[is.na(adjMat)] = 0
  }
  if (any(diag(adjMat)!=1)) {
      diag(adjMat) = 1
  }

  nGenes = dim(adjMat)[1]

  tom = matrix(0, nGenes, nGenes)

  tomResult = .C("tomSimilarityFromAdj", as.double(as.matrix(adjMat)),
                 as.integer(nGenes), as.integer(TOMTypeC),
                 as.integer(TOMDenomC), tom = as.double(tom),
                 as.integer(verbose), as.integer(indent), PACKAGE = "WGCNA")

  tom[,] = tomResult$tom
  diag(tom) = 1
  rm(tomResult)
  collectGarbage()
  tom
}

#
# TOMdist (from adjacency)
#' @rdname TOMsimilarity
#' @name TOMsimilarity
#' @export
TOMdist = function(adjMat, TOMType = "unsigned", TOMDenom = "min", verbose = 1,
                   indent = 0) {
  1-TOMsimilarity(adjMat, TOMType, TOMDenom, verbose, indent)
}


#==========================================================================================================
#
# blockwiseModules
#
#==========================================================================================================
# Function to calculate modules and eigengenes from all genes.



#' Automatic network construction and module detection
#'
#' This function performs automatic network construction and module detection
#' on large expression datasets in a block-wise manner.
#'
#' Before module detection starts, genes and samples are optionally checked for
#' the presence of \code{NA}s. Genes and/or samples that have too many
#' \code{NA}s are flagged as bad and removed from the analysis; bad genes will
#' be automatically labeled as unassigned, while the returned eigengenes will
#' have \code{NA} entries for all bad samples.
#'
#' If \code{blocks} is not given and the number of genes exceeds
#' \code{maxBlockSize}, genes are pre-clustered into blocks using the function
#' \code{\link{projectiveKMeans}}; otherwise all genes are treated in a single
#' block.
#'
#' For each block of genes, the network is constructed and (if requested)
#' topological overlap is calculated. If requested, the topological overlaps
#' are returned as part of the return value list. Genes are then clustered
#' using average linkage hierarchical clustering and modules are identified in
#' the resulting dendrogram by the Dynamic Hybrid tree cut. Found modules are
#' trimmed of genes whose correlation with module eigengene (KME) is less than
#' \code{minKMEtoStay}. Modules in which fewer than \code{minCoreKMESize} genes
#' have KME higher than \code{minCoreKME} are disbanded, i.e., their
#' constituent genes are pronounced unassigned.
#'
#' After all blocks have been processed, the function checks whether there are
#' genes whose KME in the module they assigned is lower than KME to another
#' module. If p-values of the higher correlations are smaller than those of the
#' native module by the factor \code{reassignThresholdPS}, the gene is
#' re-assigned to the closer module.
#'
#' In the last step, modules whose eigengenes are highly correlated are merged.
#' This is achieved by clustering module eigengenes using the dissimilarity
#' given by one minus their correlation, cutting the dendrogram at the height
#' \code{mergeCutHeight} and merging all modules on each branch. The process is
#' iterated until no modules are merged. See \code{\link{mergeCloseModules}}
#' for more details on module merging.
#'
#' The argument \code{quick} specifies the precision of handling of missing
#' data in the correlation calculations. Zero will cause all calculations to be
#' executed precisely, which may be significantly slower than calculations
#' without missing data. Progressively higher values will speed up the
#' calculations but introduce progressively larger errors. Without missing
#' data, all column means and variances can be pre-calculated before the
#' covariances are calculated. When missing data are present, exact
#' calculations require the column means and variances to be calculated for
#' each covariance. The approximate calculation uses the pre-calculated mean
#' and variance and simply ignores missing data in the covariance calculation.
#' If the number of missing data is high, the pre-calculated means and
#' variances may be very different from the actual ones, thus potentially
#' introducing large errors.  The \code{quick} value times the number of rows
#' specifies the maximum difference in the number of missing entries for mean
#' and variance calculations on the one hand and covariance on the other hand
#' that will be tolerated before a recalculation is triggered. The hope is that
#' if only a few missing data are treated approximately, the error introduced
#' will be small but the potential speedup can be significant.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. NAs are allowed, but not too many.
#' @param checkMissingData logical: should data be checked for excessive
#' numbers of missing entries in genes and samples, and for genes with zero
#' variance? See details.
#' @param blocks optional specification of blocks in which hierarchical
#' clustering and module detection should be performed. If given, must be a
#' numeric vector with one entry per column (gene) of \code{exprData} giving
#' the number of the block to which the corresponding gene belongs.
#' @param maxBlockSize integer giving maximum block size for module detection.
#' Ignored if \code{blocks} above is non-NULL. Otherwise, if the number of
#' genes in \code{datExpr} exceeds \code{maxBlockSize}, genes will be
#' pre-clustered into blocks whose size should not exceed \code{maxBlockSize}.
#' @param blockSizePenaltyPower number specifying how strongly blocks should be
#' penalized for exceeding the maximum size. Set to a lrge number or \code{Inf}
#' if not exceeding maximum block size is very important.
#' @param nPreclusteringCenters number of centers for pre-clustering. Larger
#' numbers typically results in better but slower pre-clustering.
#' @param randomSeed integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit. If \code{NULL} is given, the function will not save and
#' restore the seed.
#' @param loadTOM logical: should Topological Overlap Matrices be loaded from
#' previously saved files (\code{TRUE}) or calculated (\code{FALSE})? It may be
#' useful to load previously saved TOM matrices if these have been calculated
#' previously, since TOM calculation is often the most computationally
#' expensive part of network construction and module identification. See
#' \code{saveTOMs} and \code{saveTOMFileBase} below for when and how TOM files
#' are saved, and what the file names are. If \code{loadTOM} is \code{TRUE} but
#' the files cannot be found, or do not contain the correct TOM data, TOM will
#' be recalculated.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pairwise.complete.obs} option.
#' @param maxPOutliers only used for \code{corType=="bicor"}. Specifies the
#' maximum percentile of data that can be considered outliers on either side of
#' the median separately. For each side of the median, if higher percentile
#' than \code{maxPOutliers} is considered an outlier by the weight function
#' based on \code{9*mad(x)}, the width of the weight function is increased such
#' that the percentile of outliers on that side of the median equals
#' \code{maxPOutliers}. Using \code{maxPOutliers=1} will effectively disable
#' all weight function broadening; using \code{maxPOutliers=0} will give
#' results that are quite similar (but not equal to) Pearson correlation.
#' @param quickCor real number between 0 and 1 that controls the handling of
#' missing data in the calculation of correlations. See details.
#' @param pearsonFallback Specifies whether the bicor calculation, if used,
#' should revert to Pearson when median absolute deviation (mad) is zero.
#' Recongnized values are (abbreviations of) \code{"none", "individual",
#' "all"}. If set to \code{"none"}, zero mad will result in \code{NA} for the
#' corresponding correlation.  If set to \code{"individual"}, Pearson
#' calculation will be used only for columns that have zero mad.  If set to
#' \code{"all"}, the presence of a single zero mad will cause the whole
#' variable to be treated in Pearson correlation manner (as if the
#' corresponding \code{robust} option was set to \code{FALSE}). Has no effect
#' for Pearson correlation.  See \code{\link{bicor}}.
#' @param cosineCorrelation logical: should the cosine version of the
#' correlation calculation be used? The cosine calculation differs from the
#' standard one in that it does not subtract the mean.
#' @param power soft-thresholding power for network construction.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param replaceMissingAdjacencies logical: should missing values in the
#' calculation of adjacency be replaced by 0?
#' @param TOMType one of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If
#' \code{"none"}, adjacency will be used for clustering. If \code{"unsigned"},
#' the standard TOM will be used (more generally, TOM function will receive the
#' adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
#' correlations between neighbors.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are \code{"min"} giving the standard TOM described in
#' Zhang and Horvath (2005), and \code{"mean"} in which the \code{min} function
#' in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
#' better results but at this time should be considered experimental.
#' @param getTOMs deprecated, please use saveTOMs below.
#' @param saveTOMs logical: should the consensus topological overlap matrices
#' for each block be saved and returned?
#' @param saveTOMFileBase character string containing the file name base for
#' files containing the consensus topological overlaps. The full file names
#' have \code{"block.1.RData"}, \code{"block.2.RData"} etc. appended. These
#' files are standard R data files and can be loaded using the
#' \code{\link{load}} function.
#' @param deepSplit integer value between 0 and 4. Provides a simplified
#' control over how sensitive module detection should be to module splitting,
#' with 0 least and 4 most sensitive. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param detectCutHeight dendrogram cut height for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minModuleSize minimum module size for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxCoreScatter maximum scatter of the core for a branch to be a
#' cluster, given as the fraction of \code{cutHeight} relative to the 5th
#' percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minGap minimum cluster gap given as the fraction of the difference
#' between \code{cutHeight} and the 5th percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxAbsCoreScatter maximum scatter of the core for a branch to be a
#' cluster given as absolute heights. If given, overrides
#' \code{maxCoreScatter}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for
#' more details.
#' @param minAbsGap minimum cluster gap given as absolute height difference. If
#' given, overrides \code{minGap}. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minSplitHeight Minimum split height given as the fraction of the
#' difference between \code{cutHeight} and the 5th percentile of joining
#' heights. Branches merging below this height will automatically be merged.
#' Defaults to zero but is used only if \code{minAbsSplitHeight} below is
#' \code{NULL}.
#' @param minAbsSplitHeight Minimum split height given as an absolute height.
#' Branches merging below this height will automatically be merged. If not
#' given (default), will be determined from \code{minSplitHeight} above.
#' @param useBranchEigennodeDissim Logical: should branch eigennode (eigengene)
#' dissimilarity be considered when merging branches in Dynamic Tree Cut?
#' @param minBranchEigennodeDissim Minimum consensus branch eigennode
#' (eigengene) dissimilarity for branches to be considerd separate. The branch
#' eigennode dissimilarity in individual sets is simly 1-correlation of the
#' eigennodes; the consensus is defined as quantile with probability
#' \code{consensusQuantile}.
#' @param stabilityLabels Optional matrix of cluster labels that are to be used
#' for calculating branch dissimilarity based on split stability. The number of
#' rows must equal the number of genes in \code{multiExpr}; the number of
#' columns (clusterings) is arbitrary. See
#' \code{\link{branchSplitFromStabilityLabels}} for details.
#' @param minStabilityDissim Minimum stability dissimilarity criterion for two
#' branches to be considered separate. Should be a number between 0
#' (essentially no dissimilarity required) and 1 (perfect dissimilarity or
#' distinguishability based on \code{stabilityLabels}). See
#' \code{\link{branchSplitFromStabilityLabels}} for details.
#' @param pamStage logical.  If TRUE, the second (PAM-like) stage of module
#' detection will be performed.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param pamRespectsDendro Logical, only used when \code{pamStage} is
#' \code{TRUE}.  If \code{TRUE}, the PAM stage will respect the dendrogram in
#' the sense an object can be PAM-assigned only to clusters that lie below it
#' on the branch that the object is merged into.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minCoreKME a number between 0 and 1. If a detected module does not
#' have at least \code{minModuleKMESize} genes with eigengene connectivity at
#' least \code{minCoreKME}, the module is disbanded (its genes are unlabeled
#' and returned to the pool of genes waiting for mofule detection).
#' @param minCoreKMESize see \code{minCoreKME} above.
#' @param minKMEtoStay genes whose eigengene connectivity to their module
#' eigengene is lower than \code{minKMEtoStay} are removed from the module.
#' @param reassignThreshold p-value ratio threshold for reassigning genes
#' between modules. See Details.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param impute logical: should imputation be used for module eigengene
#' calculation? See \code{\link{moduleEigengenes}} for more details.
#' @param trapErrors logical: should errors in calculations be trapped?
#' @param numericLabels logical: should the returned modules be labeled by
#' colors (\code{FALSE}), or by numbers (\code{TRUE})?
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows).  If
#' zero, the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @param \dots Other arguments.
#' @return
#' A list with the following components:
#'
#' \item{colors }{ a vector of color or numeric module labels for all genes.}
#'
#' \item{unmergedColors }{ a vector of color or numeric module labels for all
#' genes before module merging.}
#'
#' \item{MEs }{ a data frame containing module eigengenes of the found modules
#' (given by \code{colors}).}
#'
#' \item{goodSamples}{numeric vector giving indices of good samples, that is
#' samples that do not have too many missing entries. }
#'
#' \item{goodGenes}{ numeric vector giving indices of good genes, that is genes
#' that do not have too many missing entries.}
#'
#' \item{dendrograms}{ a list whose components conatain hierarchical clustering
#' dendrograms of genes in each block. }
#'
#' \item{TOMFiles}{ if \code{saveTOMs==TRUE}, a vector of character strings,
#' one string per block, giving the file names of files (relative to current
#' directory) in which blockwise topological overlaps were saved. }
#'
#' \item{blockGenes}{ a list whose components give the indices of genes in each
#' block. }
#'
#' \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector
#' of length equal number of genes giving the block label for each gene. Note
#' that block labels are not necessarilly sorted in the order in which the
#' blocks were processed (since we do not require this for the input
#' \code{blocks}). See \code{blockOrder} below. }
#'
#' \item{blockOrder}{ a vector giving the order in which blocks were processed
#' and in which \code{blockGenes} above is returned. For example,
#' \code{blockOrder[1]} contains the label of the first-processed block. }
#'
#' \item{MEsOK}{logical indicating whether the module eigengenes were
#' calculated without errors. }
#' @note If the input dataset has a large number of genes, consider carefully
#' the \code{maxBlockSize} as it significantly affects the memory footprint
#' (and whether the function will fail with a memory allocation error). From a
#' theoretical point of view it is advantageous to use blocks as large as
#' possible; on the other hand, using smaller blocks is substantially faster
#' and often the only way to work with large numbers of genes. As a rough
#' guide, it is unlikely a standard desktop computer with 4GB memory or less
#' will be able to work with blocks larger than 8000 genes.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{goodSamplesGenes}} for basic quality control and filtering
#'
#' \code{\link{adjacency}}, \code{\link{TOMsimilarity}} for network
#' construction
#'
#' \code{\link[stats]{hclust}} for hierarchical clustering
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in
#' hierarchical clustering dendrograms
#'
#' \code{\link{mergeCloseModules}} for merging of close modules.
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
blockwiseModules <- function(
  # Input data
  datExpr,
  # Data checking options
  checkMissingData = TRUE,
  # Options for splitting data into blocks
  blocks = NULL,
  maxBlockSize = 5000,
  blockSizePenaltyPower = 5,
  nPreclusteringCenters = as.integer(min(ncol(datExpr)/20, 100*ncol(datExpr)/maxBlockSize)),
  randomSeed = 12345,
  # load TOM from previously saved file?
  loadTOM = FALSE,
  # Network construction arguments: correlation options
  corType = "pearson",
  maxPOutliers = 1,
  quickCor = 0,
  pearsonFallback = "individual",
  cosineCorrelation = FALSE,
  # Adjacency function options
  power = 6,
  networkType = "unsigned",
  replaceMissingAdjacencies = FALSE,
  # Topological overlap options
  TOMType = "signed",
  TOMDenom = "min",
  # Saving or returning TOM
  getTOMs = NULL,
  saveTOMs = FALSE,
  saveTOMFileBase = "blockwiseTOM",
  # Basic tree cut options
  deepSplit = 2,
  detectCutHeight = 0.995,
  minModuleSize = min(20, ncol(datExpr)/2 ),
  # Advanced tree cut options
  maxCoreScatter = NULL, minGap = NULL,
  maxAbsCoreScatter = NULL, minAbsGap = NULL,
  minSplitHeight = NULL, minAbsSplitHeight = NULL,
  useBranchEigennodeDissim = FALSE,
  minBranchEigennodeDissim = mergeCutHeight,
  stabilityLabels = NULL,
  minStabilityDissim = NULL,
  pamStage = TRUE, pamRespectsDendro = TRUE,
  # Gene reassignment, module trimming, and module "significance" criteria
  reassignThreshold = 1e-6,
  minCoreKME = 0.5,
  minCoreKMESize = minModuleSize/3,
  minKMEtoStay = 0.3,
  # Module merging options
  mergeCutHeight = 0.15,
  impute = TRUE,
  trapErrors = FALSE,
  # Output options
  numericLabels = FALSE,
  # Options controlling behaviour
  nThreads = 0,
  verbose = 0, indent = 0,
  ...) {
  spaces = indentSpaces(indent)

  if (verbose>0)
     printFlush(paste(spaces, "Calculating module eigengenes block-wise from all genes"))

  seedSaved = FALSE
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       seedSaved = TRUE
       savedSeed = .Random.seed
    }
    set.seed(randomSeed)
  }

  intCorType = pmatch(corType, .corTypes)
  if (is.na(intCorType))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  intTOMType = pmatch(TOMType, .TOMTypes)
  if (is.na(intTOMType))
    stop(paste("Invalid 'TOMType'. Recognized values are", paste(.TOMTypes, collapse = ", ")))

  TOMDenomC = pmatch(TOMDenom, .TOMDenoms)-1
  if (is.na(TOMDenomC))
    stop(paste("Invalid 'TOMDenom'. Recognized values are", paste(.TOMDenoms, collapse = ", ")))

  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.")
  if (quickCor < 0) stop("quickCor must be positive.")
  if (nThreads < 0) stop("nThreads must be positive.")
  if (is.null(nThreads) || (nThreads==0)) nThreads = .useNThreads()

  if ( (power<1) | (power>30) ) stop("power must be between 1 and 30.")
  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(fallback))
      stop(paste0("Unrecognized value '", pearsonFallback, "' of argument 'pearsonFallback'.",
           "Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))


  dimEx = dim(datExpr)
  if (length(dimEx)!=2) stop("datExpr has incorrect dimensions.")
  nGenes = dimEx[2]
  nSamples = dimEx[1]
  allLabels = rep(0, nGenes)
  AllMEs = NULL
  allLabelIndex = NULL

  if (maxBlockSize >= floor(sqrt(2^31)) )
    stop("'maxBlockSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (!is.null(blocks) && (length(blocks)!=nGenes))
    stop("Input error: the length of 'geneRank' does not equal the number of genes in given 'datExpr'.")

  if (!is.null(getTOMs))
    warning("getTOMs is deprecated, please use saveTOMs instead.")

  # Check data for genes and samples that have too many missing values

  if (checkMissingData)
  {
    gsg = goodSamplesGenes(datExpr, verbose = verbose - 1, indent = indent + 1)
    if (!gsg$allOK) datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
    nGGenes = sum(gsg$goodGenes)
    nGSamples = sum(gsg$goodSamples)
  } else {
    nGGenes = nGenes
    nGSamples = nSamples
    gsg = list(goodSamples = rep(TRUE, nSamples), goodGenes = rep(TRUE, nGenes), allOK = TRUE)
  }

  if (any(is.na(datExpr)))
  {
     datExpr.scaled.imputed = t(impute.knn(t(scale(datExpr)))$data)
  } else
     datExpr.scaled.imputed = scale(datExpr)

  corFnc = .corFnc[intCorType]
  corOptions = list(use = 'p')
  signed = networkType %in% c("signed", "signed hybrid")

  # Set up advanced tree cut methods

  otherArgs = list(...)

  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = list("branchEigengeneDissim")
    externalSplitOptions = list(list( corFnc = corFnc, corOptions = corOptions,
                                      signed = signed))
    nExternalBranchSplitFnc = 1
    externalSplitFncNeedsDistance = FALSE
    minExternalSplit = minBranchEigennodeDissim
  } else {
    branchSplitFnc = list()
    externalSplitOptions = list()
    externalSplitFncNeedsDistance = logical(0)
    nExternalBranchSplitFnc = 0
    minExternalSplit = numeric(0)
  }

  if (!is.null(stabilityLabels))
  {
    branchSplitFnc = c(branchSplitFnc, "branchSplitFromStabilityLabels")
    minExternalSplit = c(minExternalSplit, minStabilityDissim)
    externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, FALSE)
    print(dim(stabilityLabels))
    externalSplitOptions = c(externalSplitOptions, list(list(stabilityLabels = stabilityLabels)))
  }

  if ("useBranchSplit" %in% names(otherArgs))
  {
    if (otherArgs$useBranchSplit)
    {
      nExternalBranchSplitFnc = nExternalBranchSplitFnc + 1
      branchSplitFnc[[nExternalBranchSplitFnc]] = "branchSplit"
      externalSplitOptions[[nExternalBranchSplitFnc]] = list(discardProp = 0.08, minCentralProp = 0.75,
                       nConsideredPCs = 3, signed = signed, getDetails = FALSE)
      externalSplitFncNeedsDistance[nExternalBranchSplitFnc] = FALSE
      minExternalSplit[ nExternalBranchSplitFnc] = otherArgs$minBranchSplit
    }
  }

  # Split data into blocks if needed

  if (is.null(blocks))
  {
    if (nGGenes > maxBlockSize)
    {
      if (verbose>1) printFlush(paste(spaces, "....pre-clustering genes to determine blocks.."))
      clustering = projectiveKMeans(datExpr, preferredSize = maxBlockSize,
                                    checkData = FALSE,
                                    sizePenaltyPower = blockSizePenaltyPower,
                                    nCenters = nPreclusteringCenters,
                                    verbose = verbose-2, indent = indent + 1)
      gBlocks = .orderLabelsBySize(clustering$clusters)
      if (verbose > 2) { printFlush("Block sizes:"); print(table(gBlocks)); }
    } else
      gBlocks = rep(1, nGGenes)
    blocks = rep(NA, nGenes)
    blocks[gsg$goodGenes] = gBlocks
  } else {
    gBlocks = blocks[gsg$goodGenes]
  }

  blockLevels = as.numeric(levels(factor(gBlocks)))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  # Initialize various variables

  dendros = list()

  TOMFiles = rep("", nBlocks)
  blockGenes = list()

  maxUsedLabel = 0
  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))

    blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]]
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    selExpr = as.matrix(datExpr[, block])
    nBlockGenes = length(block)
    TOMFiles[blockNo] = paste0(saveTOMFileBase, "-block.", blockNo, ".RData")
    if (loadTOM)
    {
      if (verbose > 2)
        printFlush(paste(spaces, "  ..loading TOM for block", blockNo, "from file", TOMFiles[blockNo]))
      x = try(load(file =TOMFiles[blockNo]), silent = TRUE)
      if (x!="TOM")
      {
        loadTOM = FALSE
        printFlush(paste0("Loading of TOM in block ", blockNo, " failed:\n file ",
                          TOMFiles[blockNo],
                          "\n  either does not exist or does not contain the object 'TOM'.\n",
                          "  Will recalculate TOM."))
      } else if (!inherits(TOM, "dist")) {
         printFlush(paste0("TOM file ", TOMFiles[blockNo],
                           " does not contain object of the right type or size.\n",
                           " Will recalculate TOM."))
      } else {
        size.1 = attr(TOM, "Size")
        if (length(size.1)!=1 || size.1!=nBlockGenes)
        {
           printFlush(paste0("TOM file ", TOMFiles[blockNo],
                             " does not contain object of the right type or size.\n",
                            " Will recalculate TOM."))
           loadTOM = FALSE
        } else {
           tom = as.matrix(TOM)
           rm(TOM)
           collectGarbage()
        }
      }
    }

    if (!loadTOM)
    {
      # Calculate TOM by calling a custom C function:
      callVerb = max(0, verbose - 1); callInd = indent + 2
      CcorType = intCorType - 1
      CnetworkType = intNetworkType - 1
      CTOMType = intTOMType -1

      warn = as.integer(0)
      tom = .Call("tomSimilarity_call", selExpr,
          as.integer(CcorType), as.integer(CnetworkType), as.double(power), as.integer(CTOMType),
          as.integer(TOMDenomC),
          as.double(maxPOutliers),
          as.double(quickCor),
          as.integer(fallback),
          as.integer(cosineCorrelation),
          as.integer(replaceMissingAdjacencies),
          warn, as.integer(nThreads),
          as.integer(callVerb), as.integer(callInd), PACKAGE = "WGCNA")

      # FIXME: warn if necessary

      if (saveTOMs)
      {
        TOM = as.dist(tom)
        TOMFiles[blockNo] = paste0(saveTOMFileBase, "-block.", blockNo, ".RData")
        if (verbose > 2)
          printFlush(paste(spaces, "  ..saving TOM for block", blockNo, "into file", TOMFiles[blockNo]))
        save(TOM, file =TOMFiles[blockNo])
        rm (TOM)
        collectGarbage()
      }
    }
    dissTom = 1-tom
    rm(tom)
    collectGarbage()
    if (verbose>2) printFlush(paste(spaces, "....clustering.."))

    dendros[[blockNo]] = fastcluster::hclust(as.dist(dissTom), method = "average")

    if (verbose>2) printFlush(paste(spaces, "....detecting modules.."))
    datExpr.scaled.imputed.block = datExpr.scaled.imputed[, block]
    if (nExternalBranchSplitFnc > 0) for (extBSFnc in 1:nExternalBranchSplitFnc)
      externalSplitOptions[[extBSFnc]]$expr = datExpr.scaled.imputed.block

    collectGarbage()

    blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]],
                           deepSplit = deepSplit,
                           cutHeight = detectCutHeight, minClusterSize = minModuleSize,
                           method ="hybrid", distM = dissTom,
                           maxCoreScatter = maxCoreScatter, minGap = minGap,
                           maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                           minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,
                           externalBranchSplitFnc = branchSplitFnc,
                           minExternalSplit = minExternalSplit,
                           externalSplitOptions = externalSplitOptions,
                           externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                           assumeSimpleExternalSpecification = FALSE,

                           pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                           verbose = verbose-3, indent = indent + 2), silent = FALSE)
    collectGarbage()
    if (verbose > 8)
    {
      labels0 = blockLabels
      if (interactive())
        plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), dendroLabels = FALSE,
           main = paste("Block", blockNo),
           rowText = blockLabels, textPositions = 1, rowTextAlignment = "center")
      if (FALSE)
        plotDendroAndColors(dendros[[blockNo]], labels2colors(allLabels), dendroLabels = FALSE,
           main = paste("Block", blockNo))
    }
    if (class(blockLabels)=='try-error')
    {
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** cutreeDynamic returned the following error:\n",
                         spaces, blockLabels, spaces,
                         "Stopping the module detection here."))
      } else
        warning(paste("blockwiseModules: cutreeDynamic returned the following error:\n",
                      "      ", blockLabels, "---> Continuing with next block. "))
      next
    }
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
          printFlush(paste(spaces, "No modules detected in block", blockNo))
      }
      blockNo = blockNo + 1
      next
    }

    blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel
    maxUsedLabel = max(blockLabels)

    if (verbose>2) printFlush(paste(spaces, "....calculating module eigengenes.."))
    MEs = try(moduleEigengenes(selExpr[, blockLabels!=0], blockLabels[blockLabels!=0], impute = impute,
                           # subHubs = TRUE, trapErrors = FALSE,
                           verbose = verbose - 3, indent = indent + 2), silent = TRUE)
    if (class(MEs)=='try-error')
    {
      if (trapErrors)
      {
        if (verbose>0) {
          printFlush(paste(spaces, "*** moduleEigengenes failed with the following message:"))
          printFlush(paste(spaces, "       ", MEs))
          printFlush(paste(spaces, "    ---> Stopping module detection here."))
        } else
          warning(paste("blockwiseModules: moduleEigengenes failed with the following message:",
                        "\n     ", MEs, "---> Continuing with next block. "))
        next
      } else stop(MEs)
    }

    #propMEs = as.data.frame(MEs$eigengenes[, names(MEs$eigengenes)!="ME0"])

    propMEs = MEs$eigengenes
    blockLabelIndex = as.numeric(substring(names(propMEs), 3))

    deleteModules = NULL
    changedModules = NULL

    # find genes whose closest module eigengene has cor higher than minKMEtoJoin , record blockLabels and
    # remove them from the pool
    # This block has been removed for now because it conflicts with the assumption that the blocks cannot
    # be changed until they have been clustered.
    #unassGenes = c(c(1:nGGenes)[-block][allLabels[-block]==0], block[blockLabels==0])
    #if (length(unassGenes) > 0)
    #{
    #  corEval = parse(text = paste(.corFnc[intCorType], "(datExpr[, unassGenes], propMEs,",
    #                               .corOptions[intCorType], ")"))
    #  KME = eval(corEval)
    #  if (intNetworkType==1) KME = abs(KME)
    #  KMEmax = apply(KME, 1, max)
    #  ClosestModule = blockLabelIndex[apply(KME, 1, which.max)]
    #  assign = (KMEmax >= minKMEtoJoin )
    #  if (sum(assign>0))
    #  {
    #    allLabels[unassGenes[assign]] = ClosestModule[assign]
    #    changedModules = union(changedModules, ClosestModule[assign])
    #  }
    #}

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff.

    if (verbose>2) printFlush(paste(spaces, "....checking modules for statistical meaningfulness.."))
    for (mod in 1:ncol(propMEs))
    {
      modGenes = (blockLabels==blockLabelIndex[mod])
      corEval = parse(text = paste(corFnc, "(selExpr[, modGenes], propMEs[, mod]",
                                   prepComma(.corOptions[intCorType]), ")"))
      KME = as.vector(eval(corEval))
      if (intNetworkType==1) KME = abs(KME)
      if (sum(KME>minCoreKME) < minCoreKMESize)
      {
        blockLabels[modGenes] = 0
        deleteModules = c(deleteModules, mod)
        if (verbose>3)
          printFlush(paste0(spaces, "    ..deleting module ", mod, ": of ",
                            sum(modGenes),
                     " total genes in the module\n       only ",  sum(KME>minCoreKME),
                     " have the requisite high correlation with the eigengene."))
      } else if (sum(KME<minKMEtoStay)>0)
      {
        # Remove genes whose KME is too low:
        if (verbose > 2)
          printFlush(paste(spaces, "    ..removing", sum(KME<minKMEtoStay),
                           "genes from module", mod, "because their KME is too low."))
        blockLabels[modGenes][KME < minKMEtoStay] = 0
        if (sum(blockLabels[modGenes]>0) < minModuleSize)
        {
          deleteModules = c(deleteModules, mod)
          blockLabels[modGenes] = 0
          if (verbose>3)
            printFlush(paste0(spaces, "    ..deleting module ",blockLabelIndex[mod],
                     ": not enough genes in the module after removal of low KME genes."))
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod])
        }
      }
    }

    # Remove modules that are to be removed

    if (!is.null(deleteModules))
    {
       propMEs = propMEs[, -deleteModules, drop = FALSE]
       modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
       blockLabels[modGenes] = 0
       modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
       allLabels[modAllGenes] = 0
       blockLabelIndex = blockLabelIndex[-deleteModules]
    }

    # Check if any modules are left

    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
        printFlush(paste(spaces, "No significant modules detected in block", blockNo))
      }
      blockNo = blockNo + 1
      next
    }

    # Update allMEs and allLabels

    if (is.null(AllMEs))
    {
      AllMEs = propMEs
    } else
      AllMEs = cbind(AllMEs, propMEs)

    allLabelIndex = c(allLabelIndex, blockLabelIndex)

    assigned = block[blockLabels!=0]
    allLabels[gsg$goodGenes][assigned] = blockLabels[blockLabels!=0]

    rm(dissTom)
    collectGarbage()

    #if (blockNo < nBlocks)
    #{
    #  leftoverBlockGenes = block[allLabels[block]==0]
    #  nLeftBG = length(leftoverBlockGenes)
    #  blocksLeft = c((blockNo+1):nBlocks)
    #  blockSizes = as.vector(table(blocks))[blocksLeft]
    #  blocksOpen = blocksLeft[blockSizes < maxBlockSize]
    #  nBlocksOpen = length(blocksOpen)
    #  if ((nLeftBG>0) && (nBlocksOpen>0))
    #  {
    #    openSizes = blockSizes[blocksOpen]
    #    centers = matrix(0, nGSamples, nBlocksOpen)
    #    for (cen in 1:nBlocksOpen)
    #      centers[, cen] = svd(datExpr[, blocks==blocksOpen[cen]], nu=1, nv=0)$u[,1]
    #    dst = geneCenterDist(datExpr[, block], centers)
    #    rowMat = matrix(c(1:nrow(dst)), nrow(dst), ncol(dst))
    #    colMat = matrix(c(1:ncol(dst)), nrow(dst), ncol(dst), byrow = TRUE)
    #    dstOrder = order(as.vector(dst))
    #    while ((nLeftBG>0) && (nBlocksOpen>0))
    #    {
    #      gene = colMat[dstOrder[1]]
    #      rowMatV = as.vector(rowMat)
    #      colMatV = as.vector(colMat)


    blockNo = blockNo + 1
  }

  # Check whether any of the already assigned genes should be re-assigned

  deleteModules = NULL
  goodLabels = allLabels[gsg$goodGenes]
  reassignIndex = rep(FALSE, length(goodLabels))
  if (sum(goodLabels!=0) > 0)
  {
     propLabels = goodLabels[goodLabels!=0]
     assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels!=0]
     corEval = parse(text = paste(corFnc, "(datExpr[, goodLabels!=0], AllMEs",
                                  prepComma(.corOptions[intCorType]), ")"))
     KME = eval(corEval)
     if (intNetworkType == 1) KME = abs(KME)
     nMods = ncol(AllMEs)
     for (mod in 1:nMods)
     {
       modGenes = c(1:length(propLabels))[propLabels==allLabelIndex[mod]]
       KMEmodule = KME[modGenes, mod]
       KMEbest = apply(KME[modGenes, , drop = FALSE], 1, max)
       candidates = (KMEmodule < KMEbest)
       candidates[!is.finite(candidates)] = FALSE

       if (FALSE)
       {
         modDiss = dissTom[goodLabels==allLabelIndex[mod], goodLabels==allLabelIndex[mod]]
         mod.k = colSums(modDiss)
         boxplot(mod.k~candidates)
       }

       if (sum(candidates) > 0)
       {
         pModule = corPvalueFisher(KMEmodule[candidates], nSamples)
         whichBest = apply(KME[modGenes[candidates], , drop = FALSE], 1, which.max)
         pBest = corPvalueFisher(KMEbest[candidates], nSamples)
         reassign = ifelse(is.finite(pBest/pModule), (pBest/pModule < reassignThreshold), FALSE)
         if (sum(reassign)>0)
         {
           if (verbose > 2)
             printFlush(paste(spaces, " ..reassigning", sum(reassign),
                                      "genes from module", mod, "to modules with higher KME."))
           allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
           changedModules = union(changedModules, whichBest[reassign])
           if (sum(modGenes)-sum(reassign) < minModuleSize)
           {
             deleteModules = c(deleteModules, mod)
           } else
             changedModules = union(changedModules, mod)
         }
       }
     }
  }

  # Remove modules that are to be removed

  if (!is.null(deleteModules))
  {
     AllMEs = AllMEs[, -deleteModules, drop = FALSE]
     genes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
     allLabels[genes] = 0
     allLabelIndex = allLabelIndex[-deleteModules]
     goodLabels = allLabels[gsg$goodGenes]
  }

  if (verbose>1) printFlush(paste(spaces, "..merging modules that are too close.."))
  if (numericLabels) {
    colors = allLabels
  } else {
    colors = labels2colors(allLabels)
  }
  mergedAllColors = colors
  MEsOK = TRUE
  mergedMods = try(mergeCloseModules(datExpr, colors[gsg$goodGenes], cutHeight = mergeCutHeight,
                                 relabel = TRUE, # trapErrors = FALSE,
                                 impute = impute,
                                 verbose = verbose-2, indent = indent + 2), silent = TRUE)
  if (class(mergedMods)=='try-error')
  {
    warning(paste("blockwiseModules: mergeCloseModules failed with the following error message:\n    ",
                  mergedMods, "\n--> returning unmerged colors.\n"))
    MEs = try(moduleEigengenes(datExpr, colors[gsg$goodGenes], # subHubs = TRUE, trapErrors = FALSE,
                               impute = impute,
                               verbose = verbose-3, indent = indent+3), silent = TRUE)
    if (class(MEs) == 'try-error')
    {
      if (!trapErrors) stop(MEs)
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** moduleEigengenes failed with the following error message:"))
        printFlush(paste(spaces, "     ", MEs))
        printFlush(paste(spaces, "*** returning no module eigengenes.\n"))
      } else
        warning(paste("blockwiseModules: moduleEigengenes failed with the following error message:\n    ",
                      MEs, "\n--> returning no module eigengenes.\n"))
      allSampleMEs = NULL
      MEsOK = FALSE
    } else {
      if (sum(!MEs$validMEs)>0) {
        colors[gsg$goodGenes] = MEs$validColors
        MEs = MEs$eigengenes[, MEs$validMEs]
      } else MEs = MEs$eigengenes
      allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, ncol = ncol(MEs)))
      allSampleMEs[gsg$goodSamples, ] = MEs[,]
      names(allSampleMEs) = names(MEs)
    }
  } else {
    mergedAllColors[gsg$goodGenes] = mergedMods$colors
    allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples, ncol = ncol(mergedMods$newMEs)))
    allSampleMEs[gsg$goodSamples, ] = mergedMods$newMEs[,]
    names(allSampleMEs) = names(mergedMods$newMEs)
  }

  if (seedSaved) {
      .Random.seed <<- savedSeed
  }

  if (!saveTOMs) {
      TOMFiles = NULL
  }

  list(colors = mergedAllColors,
       unmergedColors = colors,
       MEs = allSampleMEs,
       goodSamples = gsg$goodSamples,
       goodGenes = gsg$goodGenes,
       dendrograms = dendros,
       TOMFiles = TOMFiles,
       blockGenes = blockGenes,
       blocks = blocks,
       MEsOK = MEsOK)
}

#==================================================================================
#
# Helper functions
#
#==================================================================================

# order labels by size

.orderLabelsBySize = function(labels, exclude = NULL)
{
  levels.0 = sort(unique(labels))
  levels = levels.0[ !levels.0 %in% exclude]
  levels.excl = levels.0 [levels.0 %in% exclude]
  rearrange = labels %in% levels
  tab = table(labels [ rearrange ])
  rank = rank(-tab, ties.method = "first")

  oldOrder = c(levels.excl, names(tab))
  newOrder = c(levels.excl, names(tab)[rank])
  if (is.numeric(labels)) newOrder = as.numeric(newOrder)

  newOrder[ match(labels, oldOrder) ]
}

#======================================================================================================
#
# Re-cut trees for blockwiseModules
#
#======================================================================================================



#' Repeat blockwise module detection from pre-calculated data
#'
#' Given consensus networks constructed for example using
#' \code{\link{blockwiseModules}}, this function (re-)detects modules in them
#' by branch cutting of the corresponding dendrograms. If repeated branch cuts
#' of the same gene network dendrograms are desired, this function can save
#' substantial time by re-using already calculated networks and dendrograms.
#'
#'
#' For details on blockwise module detection, see
#' \code{\link{blockwiseModules}}. This function implements the module
#' detection subset of the functionality of \code{\link{blockwiseModules}}
#' network construction and clustering must be performed in advance. The
#' primary use of this function is to experiment with module detection settings
#' without having to re-execute long network and clustering calculations whose
#' results are not affected by the cutting parameters.
#'
#' This function takes as input the networks and dendrograms that are produced
#' by \code{\link{blockwiseModules}}.  Working block by block, modules are
#' identified in the dendrogram by the Dynamic Hybrid Tree Cut algorithm. Found
#' modules are trimmed of genes whose correlation with module eigengene (KME)
#' is less than \code{minKMEtoStay}. Modules in which fewer than
#' \code{minCoreKMESize} genes have KME higher than \code{minCoreKME} are
#' disbanded, i.e., their constituent genes are pronounced unassigned.
#'
#' After all blocks have been processed, the function checks whether there are
#' genes whose KME in the module they assigned is lower than KME to another
#' module. If p-values of the higher correlations are smaller than those of the
#' native module by the factor \code{reassignThresholdPS}, the gene is
#' re-assigned to the closer module.
#'
#' In the last step, modules whose eigengenes are highly correlated are merged.
#' This is achieved by clustering module eigengenes using the dissimilarity
#' given by one minus their correlation, cutting the dendrogram at the height
#' \code{mergeCutHeight} and merging all modules on each branch. The process is
#' iterated until no modules are merged. See \code{\link{mergeCloseModules}}
#' for more details on module merging.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. NAs are allowed, but not too many.
#' @param goodSamples a logical vector specifying which samples are considered
#' "good" for the analysis. See \code{\link{goodSamplesGenes}}.
#' @param goodGenes a logical vector with length equal number of genes in
#' \code{multiExpr} that specifies which genes are considered "good" for the
#' analysis. See \code{\link{goodSamplesGenes}}.
#' @param blocks specification of blocks in which hierarchical clustering and
#' module detection should be performed. A numeric vector with one entry per
#' gene of \code{multiExpr} giving the number of the block to which the
#' corresponding gene belongs.
#' @param TOMFiles a vector of character strings specifying file names in which
#' the block-wise topological overlaps are saved.
#' @param dendrograms a list of length equal the number of blocks, in which
#' each component is a hierarchical clustering dendrograms of the genes that
#' belong to the block.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pariwise.complete.obs} option.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param deepSplit integer value between 0 and 4. Provides a simplified
#' control over how sensitive module detection should be to module splitting,
#' with 0 least and 4 most sensitive. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param detectCutHeight dendrogram cut height for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minModuleSize minimum module size for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxCoreScatter maximum scatter of the core for a branch to be a
#' cluster, given as the fraction of \code{cutHeight} relative to the 5th
#' percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minGap minimum cluster gap given as the fraction of the difference
#' between \code{cutHeight} and the 5th percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxAbsCoreScatter maximum scatter of the core for a branch to be a
#' cluster given as absolute heights. If given, overrides
#' \code{maxCoreScatter}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for
#' more details.
#' @param minAbsGap minimum cluster gap given as absolute height difference. If
#' given, overrides \code{minGap}. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minSplitHeight Minimum split height given as the fraction of the
#' difference between \code{cutHeight} and the 5th percentile of joining
#' heights. Branches merging below this height will automatically be merged.
#' Defaults to zero but is used only if \code{minAbsSplitHeight} below is
#' \code{NULL}.
#' @param minAbsSplitHeight Minimum split height given as an absolute height.
#' Branches merging below this height will automatically be merged. If not
#' given (default), will be determined from \code{minSplitHeight} above.
#' @param useBranchEigennodeDissim Logical: should branch eigennode (eigengene)
#' dissimilarity be considered when merging branches in Dynamic Tree Cut?
#' @param minBranchEigennodeDissim Minimum consensus branch eigennode
#' (eigengene) dissimilarity for branches to be considerd separate. The branch
#' eigennode dissimilarity in individual sets is simly 1-correlation of the
#' eigennodes; the consensus is defined as quantile with probability
#' \code{consensusQuantile}.
#' @param pamStage logical.  If TRUE, the second (PAM-like) stage of module
#' detection will be performed.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param pamRespectsDendro Logical, only used when \code{pamStage} is
#' \code{TRUE}. If \code{TRUE}, the PAM stage will respect the dendrogram in
#' the sense an object can be PAM-assigned only to clusters that lie below it
#' on the branch that the object is merged into.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minCoreKME a number between 0 and 1. If a detected module does not
#' have at least \code{minModuleKMESize} genes with eigengene connectivity at
#' least \code{minCoreKME}, the module is disbanded (its genes are unlabeled
#' and returned to the pool of genes waiting for mofule detection).
#' @param minCoreKMESize see \code{minCoreKME} above.
#' @param minKMEtoStay genes whose eigengene connectivity to their module
#' eigengene is lower than \code{minKMEtoStay} are removed from the module.
#' @param reassignThreshold p-value ratio threshold for reassigning genes
#' between modules. See Details.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param impute logical: should imputation be used for module eigengene
#' calculation? See \code{\link{moduleEigengenes}} for more details.
#' @param trapErrors logical: should errors in calculations be trapped?
#' @param numericLabels logical: should the returned modules be labeled by
#' colors (\code{FALSE}), or by numbers (\code{TRUE})?
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @param ... Other arguments.
#' @return A list with the following components:
#'
#' \item{colors }{ a vector of color or numeric module labels for all genes.}
#'
#' \item{unmergedColors }{ a vector of color or numeric module labels for all
#' genes before module merging.}
#'
#' \item{MEs }{ a data frame containing module eigengenes of the found modules
#' (given by \code{colors}).}
#'
#' \item{MEsOK}{logical indicating whether the module eigengenes were
#' calculated without errors. }
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{blockwiseModules}} for full module calculation
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in
#' hierarchical clustering dendrograms
#'
#' \code{\link{mergeCloseModules}} for merging of close modules.
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
recutBlockwiseTrees = function(datExpr,
                      goodSamples, goodGenes,
                      blocks,
                      TOMFiles,
                      dendrograms,
                      corType = "pearson",
                      networkType = "unsigned",
                      deepSplit = 2,
                      detectCutHeight = 0.995, minModuleSize = min(20, ncol(datExpr)/2 ),
                      maxCoreScatter = NULL, minGap = NULL,
                      maxAbsCoreScatter = NULL, minAbsGap = NULL,
                      minSplitHeight = NULL, minAbsSplitHeight = NULL,

                      useBranchEigennodeDissim = FALSE,
                      minBranchEigennodeDissim = mergeCutHeight,

                      pamStage = TRUE, pamRespectsDendro = TRUE,
                      # minKMEtoJoin =0.7,
                      minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
                      minKMEtoStay = 0.3,
                      reassignThreshold = 1e-6,
                      mergeCutHeight = 0.15, impute = TRUE,
                      trapErrors = FALSE, numericLabels = FALSE,
                      verbose = 0, indent = 0, ...)
{
  spaces = indentSpaces(indent)

  #if (verbose>0)
  #   printFlush(paste(spaces, "Calculating module eigengenes block-wise from all genes"))
  cutreeLabels = list()
  intCorType = pmatch(corType, .corTypes)
  if (is.na(intCorType))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  dimEx = dim(datExpr)
  if (length(dimEx)!=2) stop("datExpr has incorrect dimensions.")

  nGenes = dimEx[2]
  nSamples = dimEx[1]
  allLabels = rep(0, nGenes)
  AllMEs = NULL
  allLabelIndex = NULL

  if (length(blocks)!=nGenes)
    stop("Input error: the length of 'geneRank' does not equal the number of genes in given 'datExpr'.")

  # Check data for genes and samples that have too many missing values

  nGGenes = sum(goodGenes)
  nGSamples = sum(goodSamples)
  gsg = list(goodSamples = goodSamples, goodGenes = goodGenes,
             allOK =(sum(!goodSamples) + sum(!goodGenes) == 0))

  if (!gsg$allOK) datExpr = datExpr[goodSamples, goodGenes]

  gBlocks = blocks[gsg$goodGenes]
  blockLevels = as.numeric(levels(factor(gBlocks)))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  datExpr.scaled.imputed = t(impute.knn(t(scale(datExpr)))$data)
  if (any(is.na(datExpr)))
     datExpr.scaled.imputed = t(impute.knn(t(scale(datExpr)))$data)

  corFnc = .corFnc[intCorType]
  corOptions = list(use = 'p')

  signed = networkType %in% c("signed", "signed hybrid")
  # Set up advanced tree cut methods

  otherArgs = list(...)

  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = list("branchEigengeneDissim")
    externalSplitOptions = list(list( corFnc = corFnc, corOptions = corOptions,
                                      signed = signed))
    externalSplitFncNeedsDistance = FALSE
    nExternalBranchSplitFnc = 1
    minExternalSplit = minBranchEigennodeDissim
  } else {
    branchSplitFnc = list()
    externalSplitOptions = list(list())
    externalSplitFncNeedsDistance = logical(0)
    nExternalBranchSplitFnc = 0
    minExternalSplit = numeric(0)
  }

  if ("useBranchSplit" %in% names(otherArgs))
  {
    if (otherArgs$useBranchSplit)
    {
      nExternalBranchSplitFnc = nExternalBranchSplitFnc + 1
      branchSplitFnc[[nExternalBranchSplitFnc]] = "branchSplit"
      externalSplitOptions[[nExternalBranchSplitFnc]] = list(discardProp = 0.08, minCentralProp = 0.75,
                       nConsideredPCs = 3, signed = signed, getDetails = FALSE)
      minExternalSplit[ nExternalBranchSplitFnc] = otherArgs$minBranchSplit
      externalSplitFncNeedsDistance[ nExternalBranchSplitFnc] = FALSE
    }
  }

  # Initialize various variables

  blockNo = 1
  maxUsedLabel = 0
  while (blockNo <= nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))

    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    selExpr = as.matrix(datExpr[, block])
    nBlockGenes = length(block)

    TOM = NULL;	# Will be loaded below; this gets rid of warnings form Rcheck.

    xx = try(load(TOMFiles[blockNo]), silent = TRUE)
    if (class(xx)=='try-error') {
      printFlush(paste("************\n File name", TOMFiles[blockNo],
                       "appears invalid: the load function returned the following error:\n     ",
                       xx))
      stop()
    }
    if (xx!='TOM')
      stop(paste("The file", TOMFiles[blockNo], "does not contain the appopriate variable."))

    if (class(TOM)!="dist")
      stop(paste("The file", TOMFiles[blockNo], "does not contain the appopriate distance structure."))

    dissTom = as.matrix(1-TOM)

    if (verbose>2) printFlush(paste(spaces, "....detecting modules.."))

    datExpr.scaled.imputed.block = datExpr.scaled.imputed[, block]
    if (nExternalBranchSplitFnc > 0) for (extBSFnc in 1:nExternalBranchSplitFnc)
        externalSplitOptions[[extBSFnc]]$expr = datExpr.scaled.imputed.block

    blockLabels = try(cutreeDynamic(dendro = dendrograms[[blockNo]],
                           deepSplit = deepSplit,
                           cutHeight = detectCutHeight, minClusterSize = minModuleSize,
                           method ="hybrid",
                           maxCoreScatter = maxCoreScatter, minGap = minGap,
                           maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                           minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,

                           externalBranchSplitFnc = branchSplitFnc,
                           minExternalSplit = minExternalSplit,
                           externalSplitOptions = externalSplitOptions,
                           externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                           assumeSimpleExternalSpecification = FALSE,

                           pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                           distM = dissTom,
                           verbose = verbose-3, indent = indent + 2), silent = TRUE)
    collectGarbage()
    cutreeLabels[[blockNo]] = blockLabels

    if (class(blockLabels)=='try-error') {
      if (verbose>0) {
        printFlush(paste(spaces, "*** cutreeDynamic returned the following error:\n",
                         spaces, blockLabels, spaces,
                         "Stopping the module detection here."))
      } else
        warning(paste("blockwiseModules: cutreeDynamic returned the following error:\n",
                      "      ", blockLabels, "---> Continuing with next block. "))
      next
    }
    if (sum(blockLabels>0)==0) {
      if (verbose>1) {
          printFlush(paste(spaces, "No modules detected in block", blockNo))
      }
      blockNo = blockNo + 1
      next
    }

    blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel
    maxUsedLabel = max(blockLabels)

    if (verbose>2) printFlush(paste(spaces, "....calculating module eigengenes.."))
    MEs = try(moduleEigengenes(selExpr[, blockLabels!=0], blockLabels[blockLabels!=0], impute = impute,
                           # subHubs = TRUE, trapErrors = FALSE,
                           verbose = verbose - 3, indent = indent + 2), silent = TRUE)
    if (class(MEs)=='try-error') {
      if (trapErrors){
        if (verbose>0) {
          printFlush(paste(spaces, "*** moduleEigengenes failed with the following message:"))
          printFlush(paste(spaces, "       ", MEs))
          printFlush(paste(spaces, "    ---> Stopping module detection here."))
        } else
          warning(paste("blockwiseModules: moduleEigengenes failed with the following message:",
                        "\n     ", MEs, "---> Continuing with next block. "))
        next
      } else stop(MEs)
    }

    #propMEs = as.data.frame(MEs$eigengenes[, names(MEs$eigengenes)!="ME0"])

    propMEs = MEs$eigengenes
    blockLabelIndex = as.numeric(substring(names(propMEs), 3))

    deleteModules = NULL
    changedModules = NULL

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff.

    if (verbose>2) printFlush(paste(spaces, "....checking modules for statistical meaningfulness.."))
    for (mod in 1:ncol(propMEs)) {
      modGenes = (blockLabels==blockLabelIndex[mod])
      corEval = parse(text = paste(corFnc, "(selExpr[, modGenes], propMEs[, mod]",
                                   prepComma(.corOptions[intCorType]), ")"))
      KME = as.vector(eval(corEval))
      if (intNetworkType==1) KME = abs(KME)
      if (sum(KME>minCoreKME) < minCoreKMESize) {
        blockLabels[modGenes] = 0
        deleteModules = c(deleteModules, mod)
        if (verbose>3)
          printFlush(paste0(spaces, "    ..deleting module ", mod, ": of ", sum(modGenes),
                     " total genes in the module\n       only ",  sum(KME>minCoreKME),
                     " have the requisite high correlation with the eigengene."))
      } else if (sum(KME<minKMEtoStay)>0) {
        # Remove genes whose KME is too low:
        if (verbose > 2)
          printFlush(paste(spaces, "    ..removing", sum(KME<minKMEtoStay),
                           "genes from module", mod, "because their KME is too low."))
        blockLabels[modGenes][KME < minKMEtoStay] = 0
        if (sum(blockLabels[modGenes]>0) < minModuleSize) {
          deleteModules = c(deleteModules, mod)
          blockLabels[modGenes] = 0
          if (verbose>3)
            printFlush(paste0(spaces,
                             "    ..deleting module ",blockLabelIndex[mod],
                     ": not enough genes in the module after removal of low KME genes."))
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod])
        }
      }
    }

    # Remove modules that are to be removed

    if (!is.null(deleteModules)) {
       propMEs = propMEs[, -deleteModules, drop = FALSE]
       modGenes = is.finite(match(blockLabels,
                                  blockLabelIndex[deleteModules]))
       blockLabels[modGenes] = 0
       modAllGenes = is.finite(match(allLabels,
                                     blockLabelIndex[deleteModules]))
       allLabels[modAllGenes] = 0
       blockLabelIndex = blockLabelIndex[-deleteModules]
    }

    # Check if any modules are left

    if (sum(blockLabels>0)==0) {
      if (verbose>1) {
        printFlush(paste(spaces, "No significant modules detected in block",
                         blockNo))
      }
      blockNo = blockNo + 1
      next
    }

    # Update allMEs and allLabels

    if (is.null(AllMEs)){
      AllMEs = propMEs
    } else {
      AllMEs = cbind(AllMEs, propMEs)
    }
    allLabelIndex = c(allLabelIndex, blockLabelIndex)

    assigned = block[blockLabels!=0]
    allLabels[assigned] = blockLabels[blockLabels!=0]

    rm(dissTom)
    collectGarbage()

    blockNo = blockNo + 1
  }

  # Check whether any of the already assigned genes should be re-assigned

  deleteModules = NULL
  goodLabels = allLabels[gsg$goodGenes]
  if (sum(goodLabels!=0) > 0) {
     propLabels = goodLabels[goodLabels!=0]
     assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels!=0]
     corEval = parse(text = paste(corFnc, "(datExpr[, goodLabels!=0], AllMEs",
                                  prepComma(.corOptions[intCorType]), ")"))
     KME = eval(corEval)
     if (intNetworkType == 1) {
         KME = abs(KME)
     }
     nMods = ncol(AllMEs)
     for (mod in 1:nMods) {
       modGenes = c(1:length(propLabels))[propLabels==allLabelIndex[mod]]
       KMEmodule = KME[modGenes, mod]
       KMEbest = apply(KME[modGenes, , drop = FALSE], 1, max)
       candidates = (KMEmodule < KMEbest)
       candidates[!is.finite(candidates)] = FALSE
       if (sum(candidates) > 0) {
         pModule = corPvalueFisher(KMEmodule[candidates], nSamples)
         whichBest = apply(KME[modGenes[candidates], , drop = FALSE], 1, which.max)
         pBest = corPvalueFisher(KMEbest[candidates], nSamples)
         reassign = ifelse(is.finite(pBest/pModule), (pBest/pModule < reassignThreshold), FALSE)
         if (sum(reassign)>0) {
           if (verbose > 2) {
             printFlush(paste(spaces, " ..reassigning", sum(reassign),
                              "genes from module", mod,
                              "to modules with higher KME."))
           }
           allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
           changedModules = union(changedModules, whichBest[reassign])
           if (sum(modGenes)-sum(reassign) < minModuleSize) {
             deleteModules = c(deleteModules, mod)
           } else {
             changedModules = union(changedModules, mod)
           }
         }
       }
     }
  }

  # Remove modules that are to be removed

  if (!is.null(deleteModules)) {
     AllMEs = AllMEs[, -deleteModules, drop = FALSE]
     genes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
     allLabels[genes] = 0
     allLabelIndex = allLabelIndex[-deleteModules]
     goodLabels = allLabels[gsg$goodGenes]
  }

  if (verbose>1) {
      printFlush(paste(spaces, "..merging modules that are too close.."))
  }
  if (numericLabels) {
    colors = allLabels
  } else {
    colors = labels2colors(allLabels)
  }
  mergedAllColors = colors
  MEsOK = TRUE
  mergedMods = try(mergeCloseModules(datExpr, colors[gsg$goodGenes], cutHeight = mergeCutHeight,
                                 relabel = TRUE, # trapErrors = FALSE,
                                 impute = impute,
                                 verbose = verbose-2, indent = indent + 2), silent = TRUE)
  if (class(mergedMods)=='try-error') {
    warning(paste("blockwiseModules: mergeCloseModules failed with the following error message:\n    ",
                  mergedMods, "\n--> returning unmerged colors.\n"))
    MEs = try(moduleEigengenes(datExpr, colors[gsg$goodGenes], # subHubs = TRUE, trapErrors = FALSE,
                               impute = impute,
                               verbose = verbose-3, indent = indent+3), silent = TRUE)
    if (class(MEs) == 'try-error') {
      if (!trapErrors) stop(MEs)
      if (verbose>0) {
        printFlush(paste(spaces, "*** moduleEigengenes failed with the following error message:"))
        printFlush(paste(spaces, "     ", MEs))
        printFlush(paste(spaces, "*** returning no module eigengenes.\n"))
      } else {
        warning("blockwiseModules: moduleEigengenes failed with the following error message:\n    ",
                      MEs, "\n--> returning no module eigengenes.\n")
      }
      allSampleMEs = NULL
      MEsOK = FALSE
    } else {
      if (sum(!MEs$validMEs)>0) {
        colors[gsg$goodGenes] = MEs$validColors
        MEs = MEs$eigengenes[, MEs$validMEs]
      } else {
          MEs = MEs$eigengenes
      }
      allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples,
                                          ncol = ncol(MEs)))
      allSampleMEs[gsg$goodSamples, ] = MEs[,]
      names(allSampleMEs) = names(MEs)
    }
  } else {
    mergedAllColors[gsg$goodGenes] = mergedMods$colors
    allSampleMEs = as.data.frame(matrix(NA, nrow = nSamples,
                                        ncol = ncol(mergedMods$newMEs)))
    allSampleMEs[gsg$goodSamples, ] = mergedMods$newMEs[,]
    names(allSampleMEs) = names(mergedMods$newMEs)
  }

  list(colors = mergedAllColors,
       unmergedColors = colors,
       cutreeLabels = cutreeLabels,
       MEs = allSampleMEs,
       #goodSamples = gsg$goodSamples,
       #goodGenes = gsg$goodGenes,
       #dendrograms = dendrograms,
       #TOMFiles = TOMFiles,
       #blockGenes = blockGenes,
       MEsOK = MEsOK)
}

#==========================================================================================================
#
# blockwiseIndividualTOMs
#
#==========================================================================================================

# This function calculates and saves blockwise topological overlaps for a given multi expression data. The
# argument blocks can be given to specify blocks, or the blocks can be omitted and will be calculated if
# necessary.

# Note on naming of output files: %s will translate into set number, %N into set name (if given in
# multiExpr), %b into block number.

.substituteTags = function(format, tags, replacements) {
  nTags = length(tags)
  if (length(replacements)!= nTags) {
      stop("Length of tags and replacements must be the same.")
  }

  for (t in 1:nTags) {
    format = gsub(as.character(tags[t]),
                  as.character(replacements[t]), format, fixed = TRUE)
    }

  format
}

.processFileName = function(format, setNumber, setNames, blockNumber) {
  # The following is a workaround around empty (NULL) setNames.
  # Replaces the name with the setNumber.
  if (is.null(setNames)) setNames = rep(setNumber, setNumber)

  .substituteTags(format, c("%s", "%N", "%b"),
                  c(setNumber, setNames[setNumber], blockNumber))
}

#' Calculation of block-wise topological overlaps
#'
#' Calculates topological overlaps in the given (expression) data. If the
#' number of variables (columns) in the input data is too large, the data is
#' first split using pre-clustering, then topological overlaps are calculated
#' in each block.
#'
#' The function starts by optionally filtering out samples that have too many
#' missing entries and genes that have either too many missing entries or zero
#' variance in at least one set. Genes that are filtered out are excluded from
#' the TOM calculations.
#'
#' If \code{blocks} is not given and the number of genes exceeds
#' \code{maxBlockSize}, genes are pre-clustered into blocks using the function
#' \code{\link{consensusProjectiveKMeans}}; otherwise all genes are treated in
#' a single block.
#'
#' For each block of genes, the network is constructed and (if requested)
#' topological overlap is calculated in each set. The topological overlaps can
#' be saved to disk as RData files, or returned directly within the return
#' value (see below). Note that the matrices can be big and returning them
#' within the return value can quickly exhaust the system's memory. In
#' particular, if the block-wise calculation is necessary, it is nearly certain
#' that returning all matrices via the return value will be impossible.
#'
#' @inheritParams adjacency
#' @inheritParams blockwiseModules
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param checkPower logical: should basic sanity check be performed on the
#' supplied \code{power}? If you would like to experiment with unusual powers,
#' set the argument to \code{FALSE} and proceed with caution.
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param individualTOMFileNames character string giving the file names to save
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @return A list with the following components:
#'
#' \item{actualTOMFileNames}{Only returned if input \code{saveTOMs} is
#' \code{TRUE}. A matrix of character strings giving the file names in which
#' each block TOM is saved. Rows correspond to data sets and columns to
#' blocks.}
#'
#' \item{TOMSimilarities}{Only returned if input \code{saveTOMs} is
#' \code{FALSE}. A list in which each component corresponds to one block. Each
#' component is a matrix of dimensions (N times (number of sets)), where N is
#' the length of a distance structure corresponding to the block. That is, if
#' the block contains n genes, N=n*(n-1)/2. Each column of the matrix contains
#' the topological overlap of variables in the corresponding set ( and the
#' corresponding block), arranged as a distance structure. Do note however that
#' the topological overlap is a similarity (not a distance). }
#'
#' \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector
#' of length equal number of genes giving the block label for each gene. Note
#' that block labels are not necessarilly sorted in the order in which the
#' blocks were processed (since we do not require this for the input
#' \code{blocks}). See \code{blockOrder} below. }
#'
#' \item{blockGenes}{a list with one component for each block of genes. Each
#' component is a vector giving the indices (relative to the input
#' \code{multiExpr}) of genes in the corresponding block. }
#'
#' \item{goodSamplesAndGenes}{if input \code{checkMissingData} is \code{TRUE},
#' the output of the function \code{\link{goodSamplesGenesMS}}.  A list with
#' components \code{goodGenes} (logical vector indicating which genes passed
#' the missing data filters), \code{goodSamples} (a list of logical vectors
#' indicating which samples passed the missing data filters in each set), and
#' \code{allOK} (a logical indicating whether all genes and all samples passed
#' the filters). See \code{\link{goodSamplesGenesMS}} for more details. If
#' \code{checkMissingData} is \code{FALSE}, \code{goodSamplesAndGenes} contains
#' a list of the same type but indicating that all genes and all samples passed
#' the missing data filters.}
#'
#' The following components are present mostly to streamline the interaction of
#' this function with \code{\link{blockwiseConsensusModules}}.
#'
#' \item{nGGenes}{ Number of genes that passed missing data filters (if input
#' \code{checkMissingData} is \code{TRUE}), or the number of all genes (if
#' \code{checkMissingData} is \code{FALSE}).}
#'
#' \item{gBlocks}{ the vector \code{blocks} (above), restricted to good genes
#' only. }
#'
#' \item{nThreads}{ number of threads used to calculate correlation and TOM
#' matrices. }
#'
#' \item{saveTOMs}{ logical: were calculated matrices saved in files
#' (\code{TRUE}) or returned in the return value (\code{FALSE})?}
#'
#' \item{intNetworkType, intCorType}{integer codes for network and correlation
#' type. }
#'
#' \item{nSets}{number of sets in input data.}
#'
#' \item{setNames}{the \code{names} attribute of input \code{multiExpr}.}
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{blockwiseConsensusModules}}
#' @references For a general discussion of the weighted network formalism, see
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' The blockwise approach is briefly described in the article describing this
#' package,
#'
#' Langfelder P, Horvath S (2008) "WGCNA: an R package for weighted correlation
#' network analysis".  BMC Bioinformatics 2008, 9:559
#' @keywords misc
blockwiseIndividualTOMs = function(multiExpr,
                            # Data checking options
                            checkMissingData = TRUE,
                            # Blocking options
                            blocks = NULL,
                            maxBlockSize = 5000,
                            blockSizePenaltyPower = 5,
                            nPreclusteringCenters = NULL,
                            randomSeed = 12345,
                            # Network construction arguments: correlation options
                            corType = "pearson",
                            maxPOutliers = 1,
                            quickCor = 0,
                            pearsonFallback = "individual",
                            cosineCorrelation = FALSE,
                            # Adjacency function options
                            power = 6,
                            networkType = "unsigned",
                            checkPower = TRUE,
                            replaceMissingAdjacencies = FALSE,
                            # Topological overlap options
                            TOMType = "unsigned",
                            TOMDenom = "min",
                            # Save individual TOMs? If not, they will be returned in the session.
                            saveTOMs = TRUE,
                            individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
                            # General options
                            nThreads = 0,
                            verbose = 2, indent = 0) {
  spaces = indentSpaces(indent)

  dataSize = checkSets(multiExpr, checkStructure = TRUE)
  if (dataSize$structureOK) {
    nSets = dataSize$nSets
    nGenes = dataSize$nGenes
    multiFormat = TRUE
  } else {
    multiExpr = multiSet(multiExpr)
    nSets = dataSize$nSets
    nGenes = dataSize$nGenes
    multiFormat = FALSE
  }


  if (length(power)!=1)
  {
    if (length(power)!=nSets)
      stop("Invalid arguments: Length of 'power' must equal number of sets given in 'multiExpr'.")
  } else {
    power = rep(power, nSets)
  }

  seedSaved = FALSE
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       seedSaved = TRUE
       savedSeed = .Random.seed
    }
    set.seed(randomSeed)
  }

  if (maxBlockSize >= floor(sqrt(2^31)) )
    stop("'maxBlockSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (!is.null(blocks) && (length(blocks)!=nGenes))
    stop("Input error: length of 'blocks' must equal number of genes in 'multiExpr'.")

  if (verbose>0)
     printFlush(paste(spaces, "Calculating topological overlaps block-wise from all genes"))

  intCorType = pmatch(corType, .corTypes)
  if (is.na(intCorType))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  intTOMType = pmatch(TOMType, .TOMTypes)
  if (is.na(intTOMType))
    stop(paste("Invalid 'TOMType'. Recognized values are", paste(.TOMTypes, collapse = ", ")))

  TOMDenomC = pmatch(TOMDenom, .TOMDenoms)-1
  if (is.na(TOMDenomC))
    stop(paste("Invalid 'TOMDenom'. Recognized values are", paste(.TOMDenoms, collapse = ", ")))

  if ( checkPower & ((sum(power<1)>0) | (sum(power>50)>0) ) ) stop("power must be between 1 and 50.")

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.")
  if (quickCor < 0) stop("quickCor must be positive.")

  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(fallback))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))

  if (nThreads < 0) stop("nThreads must be positive.")
  if (is.null(nThreads) || (nThreads==0)) nThreads = .useNThreads()

  nSamples = dataSize$nSamples

  # Check data for genes and samples that have too many missing values

  if (checkMissingData)
  {
    gsg = goodSamplesGenesMS(multiExpr, verbose = verbose - 1, indent = indent + 1)
    if (!gsg$allOK)
    {
      for (set in 1:nSets)
        multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]
    }
  } else {
    gsg = list(goodGenes = rep(TRUE, nGenes), goodSamples = list())
    for (set in 1:nSets)
      gsg$goodSamples[[set]] = rep(TRUE, nSamples[set])
    gsg$allOK = TRUE
  }

  nGGenes = sum(gsg$goodGenes)
  nGSamples = rep(0, nSets)
  for (set in 1:nSets) nGSamples[set] = sum(gsg$goodSamples[[set]])

  if (is.null(blocks))
  {
    if (nGGenes > maxBlockSize)
    {
      if (verbose>1) printFlush(paste(spaces, "....pre-clustering genes to determine blocks.."))
      clustering = consensusProjectiveKMeans(multiExpr, preferredSize = maxBlockSize,
                                         sizePenaltyPower = blockSizePenaltyPower, checkData = FALSE,
                                         nCenters = nPreclusteringCenters,
                                         verbose = verbose-2, indent = indent + 1)
      gBlocks = .orderLabelsBySize(clustering$clusters)
    } else
      gBlocks = rep(1, nGGenes)
    blocks = rep(NA, nGenes)
    blocks[gsg$goodGenes] = gBlocks
  } else {
    gBlocks = blocks[gsg$goodGenes]
  }

  blockLevels = as.numeric(levels(factor(gBlocks)))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  # check file names for uniqueness

  actualFileNames = NULL
  if (saveTOMs)
  {
    actualFileNames = matrix("", nSets, nBlocks)
    for (set in 1:nSets) for (b in 1:nBlocks)
      actualFileNames[set, b] = .processFileName(individualTOMFileNames, set, names(multiExpr), b)

    rownames(actualFileNames) = paste0("Set.", c(1:nSets))
    colnames(actualFileNames) = paste0("Block.", c(1:nBlocks))
    if (length(unique(as.vector(actualFileNames))) < nSets * nBlocks)
    {
      printFlush("Error: File names for (some) set/block combinations are not unique:")
      print(actualFileNames)
      stop("File names must be unique.")
    }
  }

  # Initialize various variables

  blockGenes = list()
  blockNo = 1
  collectGarbage()
  setTomDS = list()
  # Here's where the analysis starts
  for (blockNo in 1:nBlocks)
  {
    if (verbose>1 && nBlocks > 1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select the block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    nBlockGenes = length(block)
    blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]]
    errorOccurred = FALSE

    # Set up file names or memory space to hold the set TOMs
    if (!saveTOMs)
    {
      setTomDS[[blockNo]] = array(0, dim = c(nBlockGenes*(nBlockGenes-1)/2, nSets))
    }

    # For each set: calculate and save TOM

    for (set in 1:nSets)
    {
      if (verbose>2) printFlush(paste(spaces, "....Working on set", set))
      selExpr = as.matrix(multiExpr[[set]]$data[, block])

      # Calculate TOM by calling a custom C function:
      callVerb = max(0, verbose - 1); callInd = indent + 2
      CcorType = intCorType - 1
      CnetworkType = intNetworkType - 1
      CTOMType = intTOMType - 1
      # tempExpr = as.double(as.matrix(selExpr))
      warn = 0L

      tom = .Call("tomSimilarity_call", selExpr,
          as.integer(CcorType), as.integer(CnetworkType), as.double(power[set]), as.integer(CTOMType),
          as.integer(TOMDenomC),
          as.double(maxPOutliers),
          as.double(quickCor),
          as.integer(fallback),
          as.integer(cosineCorrelation),
          as.integer(replaceMissingAdjacencies),
          warn, as.integer(nThreads),
          as.integer(callVerb), as.integer(callInd), PACKAGE = "WGCNA")

      # FIXME: warn if necessary

      tomDS = as.dist(tom)
      # dim(tom) = c(nBlockGenes, nBlockGenes)
      rm(tom)

      # Save the calculated TOM either to disk in chunks or to memory.
      if (saveTOMs)
      {
        save(tomDS, file = actualFileNames[set, blockNo])
      } else {
        setTomDS[[blockNo]] [, set] = tomDS[]
      }
    }
    rm(tomDS); collectGarbage()
  }

  if (!multiFormat)
  {
    gsg$goodSamples = gsg$goodSamples[[1]]
  }

  # Re-set random number generator if necessary
  if (seedSaved) .Random.seed <<- savedSeed

  list(actualTOMFileNames = actualFileNames,
       TOMSimilarities = if(!saveTOMs) setTomDS else NULL,
       blocks = blocks,
       blockGenes = blockGenes,
       goodSamplesAndGenes = gsg,
       nGGenes = nGGenes,
       gBlocks = gBlocks,
       nThreads = nThreads,
       saveTOMs = saveTOMs,
       intNetworkType = intNetworkType,
       intCorType = intCorType,
       nSets = nSets,
       setNames = names(multiExpr)
       )
}

#==========================================================================================================
#
# lowerTri2matrix
#
#==========================================================================================================



#' Reconstruct a symmetric matrix from a distance (lower-triangular)
#' representation
#'
#' Assuming the input vector contains a vectorized form of the distance
#' representation of a symmetric matrix, this function creates the
#' corresponding matrix. This is useful when re-forming symmetric matrices that
#' have been vectorized to save storage space.
#'
#' The function assumes that \code{x} contains the vectorized form of the
#' distance representation of a symmetric matrix. In particular, \code{x} must
#' have a length that can be expressed as n*(n-1)/2, with n an integer. The
#' result of the function is then an n times n matrix.
#'
#' @param x a numeric vector
#' @param diag value to be put on the diagonal. Recycled if necessary.
#' @return A symmetric matrix whose lower triangle is given by \code{x}.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#'   # Create a symmetric matrix
#'   m = matrix(c(1:16), 4,4)
#'   mat = (m + t(m))
#'   diag(mat) = 0
#'
#'   # Print the matrix
#'   mat
#'
#'   # Take the lower triangle and vectorize it (in two ways)
#'   x1 = mat[lower.tri(mat)]
#'   x2 = as.vector(as.dist(mat))
#'
#'   all.equal(x1, x2) # The vectors are equal
#'
#'   # Turn the vectors back into matrices
#'   new.mat = lowerTri2matrix(x1, diag = 0)
#'
#'   # Did we get back the same matrix?
#'
#'   all.equal(mat, new.mat)
#'
#'
lowerTri2matrix = function(x, diag = 1) {
  if (class(x)=="dist") {
    mat = as.matrix(x)
  } else {
    n = length(x)
    n1 = (1 + sqrt(1 + 8*n))/2
    if (floor(n1)!=n1) {
        stop("Input length does not translate into matrix")
    }
    mat = matrix(0, n1, n1)
    mat[lower.tri(mat)] = x
    mat = mat + t(mat)
  }
  diag(mat) = diag
  mat
}

#==========================================================================================================
#
# blockwiseConsensusModules
#
#==========================================================================================================


.checkComponents = function(object, names) {
  objNames = names(object)
  inObj = names %in% objNames
  if (!all(inObj))
    stop(".checkComponents: object is missing the following components:\n",
         paste(names[!inObj], collapse = ", "))
}

# Function to calculate consensus modules and eigengenes from all genes.



#' Find consensus modules across several datasets.
#'
#' Perform network construction and consensus module detection across several
#' datasets.
#'
#' The function starts by optionally filtering out samples that have too many
#' missing entries and genes that have either too many missing entries or zero
#' variance in at least one set. Genes that are filtered out are left
#' unassigned by the module detection. Returned eigengenes will contain
#' \code{NA} in entries corresponding to filtered-out samples.
#'
#' If \code{blocks} is not given and the number of genes exceeds
#' \code{maxBlockSize}, genes are pre-clustered into blocks using the function
#' \code{\link{consensusProjectiveKMeans}}; otherwise all genes are treated in
#' a single block.
#'
#' For each block of genes, the network is constructed and (if requested)
#' topological overlap is calculated in each set. To minimize memory usage,
#' calculated topological overlaps are optionally saved to disk in chunks until
#' they are needed again for the calculation of the consensus network
#' topological overlap.
#'
#' Before calculation of the consensus Topological Overlap, individual TOMs are
#' optionally calibrated. Calibration methods include single quantile scaling
#' and full quantile normalization.
#'
#' Single quantile scaling raises individual TOM in sets 2,3,... to a power
#' such that the quantiles given by \code{calibrationQuantile} agree with the
#' quantile in set 1. Since the high TOMs are usually the most important for
#' module identification, the value of \code{calibrationQuantile} is close to
#' (but not equal) 1. To speed up quantile calculation, the quantiles can be
#' determined on a randomly-chosen component subset of the TOM matrices.
#'
#' Full quantile normalization, implemented in
#' \code{\link[preprocessCore]{normalize.quantiles}}, adjusts the TOM matrices
#' such that all quantiles equal each other (and equal to the quantiles of the
#' component-wise average of the individual TOM matrices).
#'
#' Note that network calibration is performed separately in each block, i.e.,
#' the normalizing transformation may differ between blocks. This is necessary
#' to avoid manipulating a full TOM in memory.
#'
#' The consensus TOM is calculated as the component-wise
#' \code{consensusQuantile} quantile of the individual (set) TOMs; that is, for
#' each gene pair (TOM entry), the \code{consensusQuantile} quantile across all
#' input sets. Alternatively, one can also use (weighted) component-wise mean
#' across all imput data sets. If requested, the consensus topological overlaps
#' are saved to disk for later use.
#'
#' Genes are then clustered using average linkage hierarchical clustering and
#' modules are identified in the resulting dendrogram by the Dynamic Hybrid
#' tree cut. Found modules are trimmed of genes whose consensus module
#' membership kME (that is, correlation with module eigengene) is less than
#' \code{minKMEtoStay}. Modules in which fewer than \code{minCoreKMESize} genes
#' have consensus KME higher than \code{minCoreKME} are disbanded, i.e., their
#' constituent genes are pronounced unassigned.
#'
#' After all blocks have been processed, the function checks whether there are
#' genes whose KME in the module they assigned is lower than KME to another
#' module. If p-values of the higher correlations are smaller than those of the
#' native module by the factor \code{reassignThresholdPS} (in every set), the
#' gene is re-assigned to the closer module.
#'
#' In the last step, modules whose eigengenes are highly correlated are merged.
#' This is achieved by clustering module eigengenes using the dissimilarity
#' given by one minus their correlation, cutting the dendrogram at the height
#' \code{mergeCutHeight} and merging all modules on each branch. The process is
#' iterated until no modules are merged. See \code{\link{mergeCloseModules}}
#' for more details on module merging.
#'
#' The argument \code{quick} specifies the precision of handling of missing
#' data in the correlation calculations. Zero will cause all calculations to be
#' executed precisely, which may be significantly slower than calculations
#' without missing data. Progressively higher values will speed up the
#' calculations but introduce progressively larger errors. Without missing
#' data, all column means and variances can be pre-calculated before the
#' covariances are calculated. When missing data are present, exact
#' calculations require the column means and variances to be calculated for
#' each covariance. The approximate calculation uses the pre-calculated mean
#' and variance and simply ignores missing data in the covariance calculation.
#' If the number of missing data is high, the pre-calculated means and
#' variances may be very different from the actual ones, thus potentially
#' introducing large errors.  The \code{quick} value times the number of rows
#' specifies the maximum difference in the number of missing entries for mean
#' and variance calculations on the one hand and covariance on the other hand
#' that will be tolerated before a recalculation is triggered. The hope is that
#' if only a few missing data are treated approximately, the error introduced
#' will be small but the potential speedup can be significant.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param checkMissingData logical: should data be checked for excessive
#' numbers of missing entries in genes and samples, and for genes with zero
#' variance? See details.
#' @param blocks optional specification of blocks in which hierarchical
#' clustering and module detection should be performed. If given, must be a
#' numeric vector with one entry per gene of \code{multiExpr} giving the number
#' of the block to which the corresponding gene belongs.
#' @param maxBlockSize integer giving maximum block size for module detection.
#' Ignored if \code{blocks} above is non-NULL. Otherwise, if the number of
#' genes in \code{datExpr} exceeds \code{maxBlockSize}, genes will be
#' pre-clustered into blocks whose size should not exceed \code{maxBlockSize}.
#' @param blockSizePenaltyPower number specifying how strongly blocks should be
#' penalized for exceeding the maximum size. Set to a lrge number or \code{Inf}
#' if not exceeding maximum block size is very important.
#' @param nPreclusteringCenters number of centers to be used in the
#' preclustering. Defaults to smaller of \code{nGenes/20} and
#' \code{100*nGenes/maxBlockSize}, where \code{nGenes} is the nunber of genes
#' (variables) in \code{multiExpr}.
#' @param randomSeed integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit. If \code{NULL} is given, the function will not save and
#' restore the seed.
#' @param individualTOMInfo Optional data for TOM matrices in individual data
#' sets. This object is returned by the function
#' \code{\link{blockwiseIndividualTOMs}}. If not given, appropriate topological
#' overlaps will be calculated using the network contruction options below.
#' @param useIndivTOMSubset If \code{individualTOMInfo} is given, this argument
#' allows to only select a subset of the individual set networks contained in
#' \code{individualTOMInfo}. It should be a numeric vector giving the indices
#' of the individual sets to be used. Note that this argument is NOT applied to
#' \code{multiExpr}.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pariwise.complete.obs} option.
#' @param maxPOutliers only used for \code{corType=="bicor"}. Specifies the
#' maximum percentile of data that can be considered outliers on either side of
#' the median separately. For each side of the median, if higher percentile
#' than \code{maxPOutliers} is considered an outlier by the weight function
#' based on \code{9*mad(x)}, the width of the weight function is increased such
#' that the percentile of outliers on that side of the median equals
#' \code{maxPOutliers}. Using \code{maxPOutliers=1} will effectively disable
#' all weight function broadening; using \code{maxPOutliers=0} will give
#' results that are quite similar (but not equal to) Pearson correlation.
#' @param quickCor real number between 0 and 1 that controls the handling of
#' missing data in the calculation of correlations. See details.
#' @param pearsonFallback Specifies whether the bicor calculation, if used,
#' should revert to Pearson when median absolute deviation (mad) is zero.
#' Recongnized values are (abbreviations of) \code{"none", "individual",
#' "all"}. If set to \code{"none"}, zero mad will result in \code{NA} for the
#' corresponding correlation. If set to \code{"individual"}, Pearson
#' calculation will be used only for columns that have zero mad. If set to
#' \code{"all"}, the presence of a single zero mad will cause the whole
#' variable to be treated in Pearson correlation manner (as if the
#' corresponding \code{robust} option was set to \code{FALSE}). Has no effect
#' for Pearson correlation. See \code{\link{bicor}}.
#' @param cosineCorrelation logical: should the cosine version of the
#' correlation calculation be used? The cosine calculation differs from the
#' standard one in that it does not subtract the mean.
#' @param power soft-thresholding power for network construction.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param checkPower logical: should basic sanity check be performed on the
#' supplied \code{power}? If you would like to experiment with unusual powers,
#' set the argument to \code{FALSE} and proceed with caution.
#' @param replaceMissingAdjacencies logical: should missing values in the
#' calculation of adjacency be replaced by 0?
#' @param TOMType one of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If
#' \code{"none"}, adjacency will be used for clustering. If \code{"unsigned"},
#' the standard TOM will be used (more generally, TOM function will receive the
#' adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
#' correlations between neighbors.
#' @param TOMDenom a character string specifying the TOM variant to be used.
#' Recognized values are \code{"min"} giving the standard TOM described in
#' Zhang and Horvath (2005), and \code{"mean"} in which the \code{min} function
#' in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
#' better results but at this time should be considered experimental.
#' @param saveIndividualTOMs logical: should individual TOMs be saved to disk
#' for later use?.
#' @param individualTOMFileNames character string giving the file names to save
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param networkCalibration network calibration method. One of "single
#' quantile", "full quantile", "none" (or a unique abbreviation of one of
#' them).
#' @param calibrationQuantile if \code{networkCalibration} is \code{"single
#' quantile"}, topological overlaps (or adjacencies if TOMs are not computed)
#' will be scaled such that their \code{calibrationQuantile} quantiles will
#' agree.
#' @param sampleForCalibration if \code{TRUE}, calibration quantiles will be
#' determined from a sample of network similarities. Note that using all data
#' can double the memory footprint of the function and the function may fail.
#' @param sampleForCalibrationFactor determines the number of samples for
#' calibration: the number is \code{1/calibrationQuantile *
#' sampleForCalibrationFactor}. Should be set well above 1 to ensure accuracy
#' of the sampled quantile.
#' @param getNetworkCalibrationSamples logical: should samples used for TOM
#' calibration be saved for future analysis? This option is only available when
#' \code{sampleForCalibration} is \code{TRUE}.
#' @param consensusQuantile quantile at which consensus is to be defined. See
#' details.
#' @param useMean logical: should the consensus be determined from a (possibly
#' weighted) mean across the data sets rather than a quantile?
#' @param setWeights Optional vector (one component per input set) of weights
#' to be used for weighted mean consensus. Only used when \code{useMean} above
#' is \code{TRUE}.
#' @param saveConsensusTOMs logical: should the consensus topological overlap
#' matrices for each block be saved and returned?.
#' @param consensusTOMFileNames character string containing the file namefiles
#' containing the consensus topological overlaps. The tag \code{\%b} will be
#' replaced by the block number. If the resulting file names are non-unique
#' (for example, because the user gives a file name without a \code{\%b} tag),
#' an error will be generated. These files are standard R data files and can be
#' loaded using the \code{\link{load}} function.
#' @param useDiskCache should calculated network similarities in individual
#' sets be temporarilly saved to disk? Saving to disk is somewhat slower than
#' keeping all data in memory, but for large blocks and/or many sets the memory
#' footprint may be too big.
#' @param chunkSize network similarities are saved in smaller chunks of size
#' \code{chunkSize}.
#' @param cacheBase character string containing the desired name for the cache
#' files. The actual file names will consists of \code{cacheBase} and a suffix
#' to make the file names unique.
#' @param cacheDir character string containing the desired path for the cache
#' files.
#' @param consensusTOMInfo optional list summarizing consensus TOM, output of
#' \code{\link{consensusTOM}}. It contains information about pre-calculated
#' consensus TOM. Supplying this argument replaces TOM calculation, so none of
#' the individual or consensus TOM calculation arguments are taken into
#' account.
#' @param deepSplit integer value between 0 and 4. Provides a simplified
#' control over how sensitive module detection should be to module splitting,
#' with 0 least and 4 most sensitive. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param detectCutHeight dendrogram cut height for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minModuleSize minimum module size for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param checkMinModuleSize logical: should sanity checks be performed on
#' \code{minModuleSize}?
#' @param maxCoreScatter maximum scatter of the core for a branch to be a
#' cluster, given as the fraction of \code{cutHeight} relative to the 5th
#' percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minGap minimum cluster gap given as the fraction of the difference
#' between \code{cutHeight} and the 5th percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxAbsCoreScatter maximum scatter of the core for a branch to be a
#' cluster given as absolute heights. If given, overrides
#' \code{maxCoreScatter}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for
#' more details.
#' @param minAbsGap minimum cluster gap given as absolute height difference. If
#' given, overrides \code{minGap}. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minSplitHeight Minimum split height given as the fraction of the
#' difference between \code{cutHeight} and the 5th percentile of joining
#' heights. Branches merging below this height will automatically be merged.
#' Defaults to zero but is used only if \code{minAbsSplitHeight} below is
#' \code{NULL}.
#' @param minAbsSplitHeight Minimum split height given as an absolute height.
#' Branches merging below this height will automatically be merged. If not
#' given (default), will be determined from \code{minSplitHeight} above.
#' @param useBranchEigennodeDissim Logical: should branch eigennode (eigengene)
#' dissimilarity be considered when merging branches in Dynamic Tree Cut?
#' @param minBranchEigennodeDissim Minimum consensus branch eigennode
#' (eigengene) dissimilarity for branches to be considerd separate. The branch
#' eigennode dissimilarity in individual sets is simly 1-correlation of the
#' eigennodes; the consensus is defined as quantile with probability
#' \code{consensusQuantile}.
#' @param stabilityLabels Optional matrix of cluster labels that are to be used
#' for calculating branch dissimilarity based on split stability. The number of
#' rows must equal the number of genes in \code{multiExpr}; the number of
#' columns (clusterings) is arbitrary. See
#' \code{\link{branchSplitFromStabilityLabels}} for details.
#' @param minStabilityDissim Minimum stability dissimilarity criterion for two
#' branches to be considered separate. Should be a number between 0
#' (essentially no dissimilarity required) and 1 (perfect dissimilarity or
#' distinguishability based on \code{stabilityLabels}). See
#' \code{\link{branchSplitFromStabilityLabels}} for details.
#' @param pamStage logical.  If TRUE, the second (PAM-like) stage of module
#' detection will be performed.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param pamRespectsDendro Logical, only used when \code{pamStage} is
#' \code{TRUE}. If \code{TRUE}, the PAM stage will respect the dendrogram in
#' the sense an object can be PAM-assigned only to clusters that lie below it
#' on the branch that the object is merged into.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param reassignThresholdPS per-set p-value ratio threshold for reassigning
#' genes between modules.  See Details.
#' @param trimmingConsensusQuantile a number between 0 and 1 specifying the
#' consensus quantile used for kME calculation that determines module trimming
#' according to the arguments below.
#' @param minCoreKME a number between 0 and 1. If a detected module does not
#' have at least \code{minModuleKMESize} genes with eigengene connectivity at
#' least \code{minCoreKME}, the module is disbanded (its genes are unlabeled
#' and returned to the pool of genes waiting for mofule detection).
#' @param minCoreKMESize see \code{minCoreKME} above.
#' @param minKMEtoStay genes whose eigengene connectivity to their module
#' eigengene is lower than \code{minKMEtoStay} are removed from the module.
#' @param impute logical: should imputation be used for module eigengene
#' calculation? See \code{\link{moduleEigengenes}} for more details.
#' @param trapErrors logical: should errors in calculations be trapped?
#' @param equalizeQuantilesForModuleMerging Logical: equalize quantiles of the
#' module eigengene networks before module merging? If \code{TRUE}, the
#' quantiles of the eigengene correlation matrices (interpreted as a single
#' vectors of non-redundant components) will be equalized across the input data
#' sets. Note that although this seems like a reasonable option, it should be
#' considered experimental and not necessarily recommended.
#' @param quantileSummaryForModuleMerging One of \code{"mean"} or
#' \code{"median"}.  If quantile equalization of the module eigengene networks
#' is performed, the resulting "normal" quantiles will be given by this
#' function of the corresponding quantiles across the input data sets.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param mergeConsensusQuantile consensus quantile for module merging. See
#' \code{mergeCloseModules} for details.
#' @param numericLabels logical: should the returned modules be labeled by
#' colors (\code{FALSE}), or by numbers (\code{TRUE})?
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @param \dots Other arguments to be passed.
#' \code{reproduceBranchEigennodeQuantileError} that instructs the function to
#' reproduce a bug in branch eigennode dissimilarity calculations for purposes
#' if reproducing old reults.
#' @return A list with the following components:
#'
#' \item{colors}{ module assignment of all input genes. A vector containing
#' either character strings with module colors (if input \code{numericLabels}
#' was unset) or numeric module labels (if \code{numericLabels} was set to
#' \code{TRUE}). The color "grey" and the numeric label 0 are reserved for
#' unassigned genes.  }
#'
#' \item{unmergedColors }{ module colors or numeric labels before the module
#' merging step. }
#'
#' \item{multiMEs}{ module eigengenes corresponding to the modules returned in
#' \code{colors}, in multi-set format. A vector of lists, one per set,
#' containing eigengenes, proportion of variance explained and other
#' information. See \code{\link{multiSetMEs}} for a detailed description. }
#'
#' \item{goodSamples}{ a list, with one component per input set. Each component
#' is a logical vector with one entry per sample from the corresponding set.
#' The entry indicates whether the sample in the set passed basic quality
#' control criteria. }
#'
#' \item{goodGenes}{a logical vector with one entry per input gene indicating
#' whether the gene passed basic quality control criteria in all sets.}
#'
#' \item{dendrograms}{a list with one component for each block of genes. Each
#' component is the hierarchical clustering dendrogram obtained by clustering
#' the consensus gene dissimilarity in the corresponding block. }
#'
#' \item{TOMFiles}{ if \code{saveConsensusTOMs==TRUE}, a vector of character
#' strings, one string per block, giving the file names of files (relative to
#' current directory) in which blockwise topological overlaps were saved. }
#'
#' \item{blockGenes}{a list with one component for each block of genes. Each
#' component is a vector giving the indices (relative to the input
#' \code{multiExpr}) of genes in the corresponding block. }
#'
#' \item{blocks}{if input \code{blocks} was given, its copy; otherwise a vector
#' of length equal number of genes giving the block label for each gene. Note
#' that block labels are not necessarilly sorted in the order in which the
#' blocks were processed (since we do not require this for the input
#' \code{blocks}). See \code{blockOrder} below. }
#'
#' \item{blockOrder}{ a vector giving the order in which blocks were processed
#' and in which \code{blockGenes} above is returned. For example,
#' \code{blockOrder[1]} contains the label of the first-processed block. }
#'
#' \item{originCount}{if the input \code{consensusQuantile==0}, this vector
#' will contain counts of how many times each set contributed the consensus
#' gene similarity value. If the counts are highly unbalanced, the consensus
#' may be biased. }
#'
#' \item{networkCalibrationSamples}{if the input
#' \code{getNetworkCalibrationSamples} is \code{TRUE}, this component is a list
#' with one component per block. Each component is again a list with two
#' components: \code{sampleIndex} contains indices of the distance structure in
#' which TOM is stored that were sampled, and \code{TOMSamples} is a matrix
#' whose rows correspond to TOM samples and columns to individual set. Hence,
#' \code{networkCalibrationSamples[[blockNo]]$TOMSamples[index, setNo]}
#' contains the TOM entry that corresponds to element
#' \code{networkCalibrationSamples[[blockNo]]$sampleIndex[index]} of the TOM
#' distance structure in block \code{blockNo} and set \code{setNo}. (For
#' details on the distance structure, see \code{\link{dist}}.)}
#' @note If the input datasets have large numbers of genes, consider carefully
#' the \code{maxBlockSize} as it significantly affects the memory footprint
#' (and whether the function will fail with a memory allocation error). From a
#' theoretical point of view it is advantageous to use blocks as large as
#' possible; on the other hand, using smaller blocks is substantially faster
#' and often the only way to work with large numbers of genes. As a rough
#' guide, it is unlikely a standard desktop computer with 4GB memory or less
#' will be able to work with blocks larger than 7000 genes.
#'
#' %Topological overlap calculations can be speeded up substantially (several
#' 10-fold times on multi-core %systems) if R is compiled with a dedicated BLAS
#' (Basic Linear Algebra Subroutines) %library such as ATLAS or GotoBLAS and
#' the package is compiled on your target system (which is always the %case for
#' Unix, Unix-like and Mac systems, but is normally not the case on Windows
#' systems).
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{goodSamplesGenesMS}} for basic quality control and filtering
#'
#' \code{\link{adjacency}}, \code{\link{TOMsimilarity}} for network
#' construction
#'
#' \code{\link{hclust}} for hierarchical clustering
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in
#' hierarchical clustering dendrograms
#'
#' \code{\link{mergeCloseModules}} for merging of close modules.
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
blockwiseConsensusModules <- function(multiExpr,
         # Data checking options
         checkMissingData = TRUE,
         # Blocking options
         blocks = NULL,
         maxBlockSize = 5000,
         blockSizePenaltyPower = 5,
         nPreclusteringCenters = NULL,
         randomSeed = 12345,
         # individual TOM information
         individualTOMInfo = NULL,
         useIndivTOMSubset = NULL,
         # Network construction arguments: correlation options
         corType = "pearson",
         maxPOutliers = 1,
         quickCor = 0,
         pearsonFallback = "individual",
         cosineCorrelation = FALSE,
         # Adjacency function options
         power = 6,
         networkType = "unsigned",
         checkPower = TRUE,
         replaceMissingAdjacencies = FALSE,
         # Topological overlap options
         TOMType = "unsigned",
         TOMDenom = "min",
         # Save individual TOMs?
         saveIndividualTOMs = TRUE,
         individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
         # Consensus calculation options: network calibration
         networkCalibration = c("single quantile", "full quantile", "none"),
         # Simple quantile calibration options
         calibrationQuantile = 0.95,
         sampleForCalibration = TRUE,
         sampleForCalibrationFactor = 1000,
         getNetworkCalibrationSamples = FALSE,
         # Consensus definition
         consensusQuantile = 0,
         useMean = FALSE,
         setWeights = NULL,
         # Saving the consensus TOM
         saveConsensusTOMs = FALSE,
         consensusTOMFileNames = "consensusTOM-block.%b.RData",
         # Internal handling of TOMs
         useDiskCache = TRUE, chunkSize = NULL,
         cacheBase = ".blockConsModsCache",
         cacheDir = ".",
         # Alternative consensus TOM input from a previous calculation
         consensusTOMInfo = NULL,
         # Basic tree cut options
         deepSplit = 2,
         detectCutHeight = 0.995, minModuleSize = 20,
         checkMinModuleSize = TRUE,
         # Advanced tree cut opyions
         maxCoreScatter = NULL, minGap = NULL,
         maxAbsCoreScatter = NULL, minAbsGap = NULL,
         minSplitHeight = NULL, minAbsSplitHeight = NULL,
         useBranchEigennodeDissim = FALSE,
         minBranchEigennodeDissim = mergeCutHeight,
         stabilityLabels = NULL,
         minStabilityDissim = NULL,
         pamStage = TRUE,  pamRespectsDendro = TRUE,
         # Gene joining and removal from a module, and module "significance" criteria
         reassignThresholdPS = 1e-4,
         trimmingConsensusQuantile = consensusQuantile,
         # minKMEtoJoin =0.7,
         minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
         minKMEtoStay = 0.2,
         # Module eigengene calculation options
         impute = TRUE,
         trapErrors = FALSE,
         # Module merging options
         equalizeQuantilesForModuleMerging = FALSE,
         quantileSummaryForModuleMerging = "mean",
         mergeCutHeight = 0.15,
         mergeConsensusQuantile = consensusQuantile,
         # Output options
         numericLabels = FALSE,
         # General options
         nThreads = 0,
         verbose = 2, indent = 0,
         ...)
{
  spaces = indentSpaces(indent)

  dataSize = checkSets(multiExpr)
  nSets = dataSize$nSets
  nGenes = dataSize$nGenes

  if (length(power)!=1)
  {
    if (length(power)!=nSets)
      stop("Invalid arguments: Length of 'power' must equal number of sets given in 'multiExpr'.")
  } else {
    power = rep(power, nSets)
  }

  seedSaved = FALSE
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       seedSaved = TRUE
       savedSeed = .Random.seed
    }
    set.seed(randomSeed)
  }

  if ( (consensusQuantile < 0) | (consensusQuantile > 1) )
    stop("'consensusQuantile' must be between 0 and 1.")

  if (checkMinModuleSize & (minModuleSize > nGenes/2))
  {
    minModuleSize = nGenes/2
    warning("blockwiseConsensusModules: minModuleSize appeared too large and was lowered to",
            minModuleSize,
            ". If you insist on keeping the original minModuleSize, set checkMinModuleSize = FALSE.")
  }

  if (verbose>0)
     printFlush(paste(spaces, "Calculating consensus modules and module eigengenes",
                      "block-wise from all genes"))

  # prepare scaled and imputed multiExpr.
  multiExpr.scaled = apply(multiExpr, scale)
  hasMissing = unlist(multiSet2list(apply(multiExpr, function(x) { any(is.na(x)) })))
  # Impute those that have missing data
  multiExpr.scaled.imputed = multiSet.mapply(function(x, doImpute){
      if (doImpute) t(impute.knn(t(x))$data) else x },
                                   multiExpr.scaled, hasMissing)
  branchSplitFnc = NULL
  minBranchDissimilarities = numeric(0)
  externalSplitFncNeedsDistance = logical(0)
  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = "multiSet.branchEigengeneDissim"
    minBranchDissimilarities = minBranchEigennodeDissim
    externalSplitFncNeedsDistance = FALSE
  }

  if (!is.null(stabilityLabels))
  {
    branchSplitFnc = c(branchSplitFnc, "branchSplitFromStabilityLabels")
    minBranchDissimilarities = c(minBranchDissimilarities, minStabilityDissim)
    externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, FALSE)
  }

  # Basic checks on consensusTOMInfo

  if (!is.null(consensusTOMInfo))
  {

    .checkComponents(consensusTOMInfo, c("saveConsensusTOMs", "individualTOMInfo", "goodSamplesAndGenes"))

    if (length(consensusTOMInfo$individualTOMInfo$blocks)!=nGenes)
      stop("Inconsistent number of genes in 'consensusTOMInfo$individualTOMInfo$blocks'.")

    if (!is.null(consensusTOMInfo$consensusQuantile) &&
              (consensusQuantile!=consensusTOMInfo$consensusQuantile) )
       warning(immediate. = TRUE,
              "blockwiseConsensusModules: given (possibly default) 'consensusQuantile' is different\n",
              "from the value recorded in 'consensusTOMInfo'. This is normally undesirable and may\n",
              "indicate a mistake in the function call.")
  }

  # Handle "other arguments"

  args = list(...)
  if (is.null(args$reproduceBranchEigennodeQuantileError))
  {
     reproduceBranchEigennodeQuantileError = FALSE
  } else reproduceBranchEigennodeQuantileError = args$reproduceBranchEigennodeQuantileError

  # If topological overlaps weren't calculated yet, calculate them.

  removeIndividualTOMsOnExit = FALSE
  nBlocks.0 = length(unique(blocks))

  if (is.null(individualTOMInfo))
  {
    if (is.null(consensusTOMInfo))
    {
      individualTOMInfo = blockwiseIndividualTOMs(multiExpr = multiExpr,
                           checkMissingData = checkMissingData,
                           blocks = blocks,
                           maxBlockSize = maxBlockSize,
                           blockSizePenaltyPower = blockSizePenaltyPower,
                           nPreclusteringCenters = nPreclusteringCenters,
                           randomSeed = NULL,
                           corType = corType,
                           maxPOutliers = maxPOutliers,
                           quickCor = quickCor,
                           pearsonFallback = pearsonFallback,
                           cosineCorrelation = cosineCorrelation,
                           power = power,
                           networkType = networkType,
                           replaceMissingAdjacencies= replaceMissingAdjacencies,
                           TOMType = TOMType,
                           TOMDenom = TOMDenom,
                           saveTOMs = useDiskCache | nBlocks.0>1,
                           individualTOMFileNames = individualTOMFileNames,
                           nThreads = nThreads,
                           verbose = verbose, indent = indent)
      removeIndividualTOMsOnExit = TRUE
    } else
      individualTOMInfo = consensusTOMInfo$individualTOMInfo
  }

  if (is.null(useIndivTOMSubset))
  {
    if (individualTOMInfo$nSets != nSets)
      stop(paste("Number of sets in individualTOMInfo and in multiExpr do not agree.\n",
                 "  To use a subset of individualTOMInfo, set useIndivTOMSubset appropriately."))

    useIndivTOMSubset = c(1:nSets)
  }

  if (length(useIndivTOMSubset)!=nSets)
    stop("Length of 'useIndivTOMSubset' must equal the number of sets in 'multiExpr'")

  if (length(unique(useIndivTOMSubset))!=nSets)
    stop("Entries of 'useIndivTOMSubset' must be unique")

  if (any(useIndivTOMSubset<1) | any(useIndivTOMSubset>individualTOMInfo$nSets))
    stop("All entries of 'useIndivTOMSubset' must be between 1 and the number of sets in individualTOMInfo")

  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  intNetworkType = individualTOMInfo$intNetworkType
  intCorType = individualTOMInfo$intCorType

  corFnc = match.fun(.corFnc[intCorType])
  corOptions = list(use = 'p')

  fallback = pmatch(pearsonFallback, .pearsonFallbacks)

  nSamples = dataSize$nSamples

  allLabels = rep(0, nGenes)
  allLabelIndex = NULL

  # Restrict data to goodSamples and goodGenes

  gsg = individualTOMInfo$goodSamplesAndGenes

  # Restrict gsg to used sets

  gsg$goodSamples = gsg$goodSamples[useIndivTOMSubset]
  if (!gsg$allOK)
    multiExpr = multiExpr[gsg$goodSamples, gsg$goodGenes]

  nGGenes = sum(gsg$goodGenes)
  nGSamples = rep(0, nSets)
  for (set in 1:nSets) nGSamples[set] = sum(gsg$goodSamples[[ set ]])

  blocks = individualTOMInfo$blocks
  gBlocks = individualTOMInfo$gBlocks

  blockLevels = sort(unique(gBlocks))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/nSets)

  reassignThreshold = reassignThresholdPS^nSets

  consMEs = vector(mode = "list", length = nSets)
  dendros = list()

  maxUsedLabel = 0
  collectGarbage()
  # Here's where the analysis starts

  removeConsensusTOMOnExit = FALSE
  if (is.null(consensusTOMInfo) && (nBlocks==1 || saveConsensusTOMs || getNetworkCalibrationSamples))
  {
    consensusTOMInfo = consensusTOM(
         individualTOMInfo = individualTOMInfo,
         useIndivTOMSubset = useIndivTOMSubset,

         networkCalibration = networkCalibration,
         saveCalibratedIndividualTOMs = FALSE,

         calibrationQuantile = calibrationQuantile,
         sampleForCalibration = sampleForCalibration,
         sampleForCalibrationFactor = sampleForCalibrationFactor,
         getNetworkCalibrationSamples = getNetworkCalibrationSamples,

         consensusQuantile = consensusQuantile,
         useMean = useMean,
         setWeights = setWeights,

         # Return options
         saveConsensusTOMs = saveConsensusTOMs,
         consensusTOMFileNames = consensusTOMFileNames,
         returnTOMs = nBlocks==1,

         # Internal handling of TOMs
         useDiskCache = useDiskCache,
         chunkSize = chunkSize,
         cacheBase = cacheBase,
         cacheDir = cacheDir,
         verbose = verbose, indent = indent)
     removeConsensusTOMOnExit = !saveConsensusTOMs
  }

  blockwiseConsensusCalculation = is.null(consensusTOMInfo)

  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    nBlockGenes = length(block)

    selExpr = multiExpr[ , block]
    errorOccurred = FALSE

    if (blockwiseConsensusCalculation)
    {
       # This code is only reached if input saveConsensusTOMs is FALSE and there are at least 2 blocks.
       consensusTOMInfo = consensusTOM(
         individualTOMInfo = individualTOMInfo,
         useIndivTOMSubset = useIndivTOMSubset,

         useBlocks = blockNo,

         networkCalibration = networkCalibration,
         saveCalibratedIndividualTOMs = FALSE,

         calibrationQuantile = calibrationQuantile,
         sampleForCalibration = sampleForCalibration,
         sampleForCalibrationFactor = sampleForCalibrationFactor,
         getNetworkCalibrationSamples = FALSE,

         consensusQuantile = consensusQuantile,
         useMean = useMean,
         setWeights = setWeights,

         saveConsensusTOMs = FALSE,
         returnTOMs = TRUE,

         useDiskCache = useDiskCache,
         chunkSize = chunkSize,
         cacheBase = cacheBase,
         cacheDir = cacheDir)

       consTomDS = consensusTOMInfo$consensusTOM[[1]]
       # Remove the consensus TOM from the structure.
       consensusTOMInfo$consensusTOM[[1]] = NULL
       consensusTOMInfo$consensusTOM = NULL
    } else {
       if (consensusTOMInfo$saveConsensusTOMs)
       {
         consTomDS = .loadObject(file = consensusTOMInfo$TOMFiles[blockNo],
                                 size = nBlockGenes * (nBlockGenes-1)/2)
       } else
         consTomDS = consensusTOMInfo$consensusTOM[[blockNo]]
    }

    # Temporary "cast" so fastcluster::hclust doesn't complain about non-integer size.

    attr(consTomDS, "Size") = as.integer(attr(consTomDS, "Size"))

    consTomDS = 1-consTomDS

    collectGarbage()

    if (verbose>2) printFlush(paste(spaces, "....clustering and detecting modules.."))
    errorOccured = FALSE
    dendros[[blockNo]] = fastcluster::hclust(consTomDS, method = "average")
    if (verbose > 8)
    {
      if (interactive())
        plot(dendros[[blockNo]], labels = FALSE, main = paste("Block", blockNo))
    }

    externalSplitOptions = list()
    e.index = 1
    if (useBranchEigennodeDissim)
    {
      externalSplitOptions[[e.index]] = list(multiExpr = multiExpr.scaled.imputed[, block],
                                      corFnc = corFnc, corOptions = corOptions,
                                      consensusQuantile = consensusQuantile,
                                      signed = networkType %in% c("signed", "signed hybrid"),
                                      reproduceQuantileError = reproduceBranchEigennodeQuantileError)
      e.index = e.index +1
    }
    if (!is.null(stabilityLabels))
    {
      externalSplitOptions[[e.index]] = list(stabilityLabels = stabilityLabels)
      e.index = e.index + 1
    }
    collectGarbage()
    #blockLabels = try(cutreeDynamic(dendro = dendros[[blockNo]],
    blockLabels = cutreeDynamic(dendro = dendros[[blockNo]],
                    distM = as.matrix(consTomDS),
                    deepSplit = deepSplit,
                    cutHeight = detectCutHeight, minClusterSize = minModuleSize,
                    method ="hybrid",
                    maxCoreScatter = maxCoreScatter, minGap = minGap,
                    maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                    minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,

                    externalBranchSplitFnc = branchSplitFnc,
                    minExternalSplit = minBranchDissimilarities,
                    externalSplitOptions = externalSplitOptions,
                    externalSplitFncNeedsDistance = externalSplitFncNeedsDistance,
                    assumeSimpleExternalSpecification = FALSE,

                    pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                    verbose = verbose-3, indent = indent + 2)
                    #verbose = verbose-3, indent = indent + 2), silent = TRUE)
    if (verbose > 8)
    {
      print(table(blockLabels))
      if (interactive())
        plotDendroAndColors(dendros[[blockNo]], labels2colors(blockLabels), dendroLabels = FALSE,
           main = paste("Block", blockNo))
    }
    if (class(blockLabels)=='try-error')
    {
      (if (verbose>0) printFlush else warning)
           (paste(spaces, "blockwiseConsensusModules: cutreeDynamic failed:\n    ", spaces,
                  blockLabels, "\n", spaces, "    Error occured in block", blockNo, "\n",
                  spaces, "   Continuing with next block. "))
      next
    } else {
      blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel
      maxUsedLabel = max(blockLabels)
    }
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
          printFlush(paste(spaces, "No modules detected in block", blockNo,
                           "--> continuing with next block."))
      }
      next
    }

    # Calculate eigengenes for this batch

    if (verbose>2) printFlush(paste(spaces, "....calculating eigengenes.."))
    blockAssigned = c(1:nBlockGenes)[blockLabels!=0]
    blockLabelIndex = as.numeric(levels(as.factor(blockLabels[blockAssigned])))
    blockConsMEs = try(multiSetMEs(selExpr, universalColors = blockLabels,
                                   excludeGrey = TRUE, grey = 0, impute = impute,
                                   # trapErrors = TRUE, returnValidOnly = TRUE,
                                   verbose = verbose-4, indent = indent + 3), silent = TRUE)
    if (class(blockConsMEs)=='try-error')
    {
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** multiSetMEs failed with the message:"))
        printFlush(paste(spaces, "     ", blockConsMEs))
        printFlush(paste(spaces, "*** --> Ending module detection here"))
      } else warning(paste("blocwiseConsensusModules: multiSetMEs failed with the message: \n",
               "      ", blockConsMEs, "\n--> continuing with next block."))
      next
    }

    deleteModules = NULL
    changedModules = NULL

    # find genes whose closest module eigengene has cor higher than minKMEtoJoin and assign them
    # Removed - should not change blocks before clustering them
    #unassGenes = c(c(1:nGGenes)[-block][allLabels[-block]==0], block[blockLabels==0])
    #if (length(unassGenes) > 0)
    #{
      #blockKME = array(0, dim = c(length(unassGenes), ncol(blockConsMEs[[1]]$data), nSets))
      #corEval = parse(text = paste(.corFnc[intCorType],
                         #"( multiExpr[[set]]$data[, unassGenes], blockConsMEs[[set]]$data,",
                         #.corOptions[intCorType], ")"))
      #for (set in 1:nSets) blockKME[, , set] = eval(corEval)
      #if (intNetworkType==1) blockKME = abs(blockKME)
      #consKME = as.matrix(apply(blockKME, c(1,2), min))
      #consKMEmax = apply(consKME, 1, max)
      #closestModule = blockLabelIndex[apply(consKME, 1, which.max)]
      #assign = (consKMEmax >= minKMEtoJoin )
      #if (sum(assign>0))
      #{
        #allLabels[unassGenes[assign]] = closestModule[assign]
        #changedModules = union(changedModules, closestModule[assign])
      #}
      #rm(blockKME, consKME, consKMEmax)
    #}

    collectGarbage()

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
    # the required minimum consensus KME

    if (verbose>2)
      printFlush(paste(spaces, "....checking consensus modules for statistical meaningfulness.."))

    for (mod in 1:ncol(blockConsMEs[[1]]$data))
    {
      modGenes = (blockLabels==blockLabelIndex[mod])
      KME = matrix(0, nrow = sum(modGenes), ncol = nSets)
      corEval = parse(text = paste(.corFnc[intCorType],
                       "( selExpr[[set]]$data[, modGenes], blockConsMEs[[set]]$data[, mod]",
                      prepComma(.corOptions[intCorType]), ")"))
      for (set in 1:nSets) KME[, set] = as.vector(eval(corEval))
      if (intNetworkType==1) KME = abs(KME)
      consKME = apply(KME, 1, quantile, probs = trimmingConsensusQuantile, names = FALSE, na.rm = TRUE)
      if (sum(consKME>minCoreKME) < minCoreKMESize)
      {
        blockLabels[modGenes] = 0
        deleteModules = union(deleteModules, mod)
        if (verbose>3)
          printFlush(paste0(spaces, "    ..deleting module ",blockLabelIndex[mod],
                           ": of ", sum(modGenes),
                     " total genes in the module only ",  sum(consKME>minCoreKME),
                     " have the requisite high correlation with the eigengene in all sets."))
      } else if (sum(consKME<minKMEtoStay)>0)
      {
        if (verbose > 3)
          printFlush(paste(spaces, "    ..removing", sum(consKME<minKMEtoStay),
                           "genes from module", blockLabelIndex[mod], "because their KME is too low."))
        blockLabels[modGenes][consKME < minKMEtoStay] = 0
        if (sum(blockLabels[modGenes]>0) < minModuleSize)
        {
          deleteModules = union(deleteModules, mod)
          blockLabels[modGenes] = 0
          if (verbose>3)
            printFlush(paste0(spaces, "    ..deleting module ",blockLabelIndex[mod],
                     ": not enough genes in the module after removal of low KME genes."))
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod])
        }
      }
    }

    # Remove marked modules

    if (!is.null(deleteModules))
    {
       for (set in 1:nSets) blockConsMEs[[set]]$data = blockConsMEs[[set]]$data[, -deleteModules]
       modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
       blockLabels[modGenes] = 0
       modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
       allLabels[modAllGenes] = 0
       blockLabelIndex = blockLabelIndex[-deleteModules]
    }

    # Check whether there's anything left
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
        printFlush(paste(spaces, "  ..No significant modules detected in block", blockNo))
        printFlush(paste(spaces, "  ..continuing with next block."))
      }
      next
    }

    # Update module eigengenes

    for (set in 1:nSets)
      if (is.null(dim(blockConsMEs[[set]]$data)))
        dim(blockConsMEs[[set]]$data) = c(length(blockConsMEs[[set]]$data), 1)

    if (is.null(consMEs[[1]]))
    {
       for (set in 1:nSets) consMEs[[set]] = list(data = blockConsMEs[[set]]$data)
    } else for (set in 1:nSets)
       consMEs[[set]]$data = cbind(consMEs[[set]]$data, blockConsMEs[[set]]$data)

    # Update allLabels

    allLabelIndex = c(allLabelIndex, blockLabelIndex)
    allLabels[gsg$goodGenes][block[blockAssigned]] = blockLabels[blockAssigned]

    collectGarbage()

  }

  # Check whether any of the already assigned genes (in this or previous blocks) should be re-assigned

  if (verbose>2)
    printFlush(paste(spaces, "....checking for genes that should be reassigned.."))

  deleteModules = NULL
  goodLabels = allLabels[gsg$goodGenes]
  if (sum(goodLabels!=0) > 0)
  {
     propLabels = goodLabels[goodLabels!=0]
     assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels!=0]
     corEval = parse(text = paste(.corFnc[intCorType],
                                  "(multiExpr[[set]]$data[, goodLabels!=0], consMEs[[set]]$data",
                                  prepComma(.corOptions[intCorType]), ")"))
     nMods = ncol(consMEs[[1]]$data)
     lpValues = array(0, dim = c(length(propLabels), nMods, nSets))
     sumSign = array(0, dim = c(length(propLabels), nMods))
     if (verbose>3)
       printFlush(paste(spaces, "......module membership p-values.."))
     for (set in 1:nSets)
     {
       KME = eval(corEval)
       if (intNetworkType == 1) KME = abs(KME)
       lpValues[,,set] = -2*log(corPvalueFisher(KME, nGSamples[set], twoSided = FALSE))
       sumSign = sumSign + sign(KME)
     }
     if (verbose>3)
       printFlush(paste(spaces, "......module membership scores.."))
     scoreAll = as.matrix(apply(lpValues, c(1,2), sum)) * (nSets + sumSign)/(2*nSets)
     scoreAll[!is.finite(scoreAll)] = 0.001 # This low should be enough
     bestScore = apply(scoreAll, 1, max)
     if (intNetworkType==1) sumSign = abs(sumSign)
     if (verbose>3)
     {
       cat(paste(spaces, "......individual modules.."))
       pind = initProgInd()
     }
     for (mod in 1:nMods)
     {
       modGenes = c(1:length(propLabels))[propLabels==allLabelIndex[mod]]
       scoreMod = scoreAll[modGenes, mod]
       candidates = (bestScore[modGenes] > scoreMod)
       candidates[!is.finite(candidates)] = FALSE
       if (sum(candidates) > 0)
       {
         pModule = pchisq(scoreMod[candidates], nSets, log.p = TRUE)
         whichBest = apply(scoreAll[modGenes[candidates], ], 1, which.max)
         pBest = pchisq(bestScore[modGenes[candidates]], nSets, log.p = TRUE)
         reassign =  ifelse(is.finite(pBest - pModule),
                            ( (pBest - pModule) < log(reassignThreshold) ),
                            FALSE)
         if (sum(reassign)>0)
         {
           allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
           changedModules = union(changedModules, whichBest[reassign])
           if (sum(modGenes)-sum(reassign) < minModuleSize)
           {
             deleteModules = union(deleteModules, mod)
           } else
             changedModules = union(changedModules, mod)
         }
       }
       if (verbose > 3) pind = updateProgInd(mod/nMods, pind)
     }
     rm(lpValues, sumSign, scoreAll)
     if (verbose > 3) printFlush("")
  }

  # Remove marked modules

  if (!is.null(deleteModules))
  {
     # for (set in 1:nSets) consMEs[[set]]$data = consMEs[[set]]$data[, -deleteModules]
     modGenes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
     allLabels[modGenes] = 0
     # allLabelIndex = allLabelIndex[-deleteModules]
  }

  if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."))
  #print(table(allLabels))
  #print(is.numeric(allLabels))
  if (numericLabels) {
    colors = allLabels
  } else {
    colors = labels2colors(allLabels)
  }
  mergedColors = colors
  mergedMods = try(mergeCloseModules(multiExpr, colors[gsg$goodGenes],
                                     equalizeQuantiles = equalizeQuantilesForModuleMerging,
                                     quantileSummary = quantileSummaryForModuleMerging,
                                     consensusQuantile = mergeConsensusQuantile,
                                     cutHeight = mergeCutHeight,
                                     relabel = TRUE, impute = impute,
                                     verbose = verbose-2, indent = indent + 2), silent = TRUE)
  if (class(mergedMods)=='try-error')
  {
    if (verbose>0)
    {
      printFlush(paste(spaces, 'blockwiseConsensusModules: mergeCloseModule failed with this message:\n',
            spaces, '    ', mergedMods, spaces,
            '---> returning unmerged consensus modules'))
    } else warning(paste('blockwiseConsensusModules: mergeCloseModule failed with this message:\n     ',
                          mergedMods, '---> returning unmerged consensus modules'))
    MEs = try(multiSetMEs(multiExpr, universalColors = colors[gsg$goodGenes]
                          # trapErrors = TRUE, returnValidOnly = TRUE
                          ), silent = TRUE)
    if (class(MEs)=='try-error')
    {
      warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
            MEs, '---> returning empty module eigengenes'))
      allSampleMEs = NULL
    } else {
      if (!MEs[[1]]$allOK) mergedColors[gsg$goodGenes] = MEs[[1]]$validColors
      allSampleMEs = vector(mode = "list", length = nSets)
      for (set in 1:nSets)
      {
        allSampleMEs[[set]] =
           list(data = as.data.frame(matrix(NA, nrow = nGSamples[set], ncol = ncol(MEs[[set]]$data))))
        allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = MEs[[set]]$data[,]
        names(allSampleMEs[[set]]$data) = names(MEs[[set]]$data)
      }
    }
  } else {
    mergedColors[gsg$goodGenes] = mergedMods$colors
    allSampleMEs = vector(mode = "list", length = nSets)
    for (set in 1:nSets)
    {
      allSampleMEs[[set]] =
         list(data = as.data.frame(matrix(NA, nrow = nGSamples[set],
                                          ncol = ncol(mergedMods$newMEs[[1]]$data))))
      allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = mergedMods$newMEs[[set]]$data[,]
      names(allSampleMEs[[set]]$data) = names(mergedMods$newMEs[[set]]$data)
    }
  }

  if (seedSaved) .Random.seed <<- savedSeed

  if (removeConsensusTOMOnExit)
  {
    .checkAndDelete(consensusTOMInfo$TOMFiles)
    consensusTOMInfo$TOMFiles = NULL
  }

  if (removeIndividualTOMsOnExit)
  {
    .checkAndDelete(individualTOMInfo$actualTOMFileNames)
    individualTOMInfo$actualTOMFileNames = NULL
  }

  # Under no circumstances return consensus TOM or individual TOM similarities within the returned list.
  consensusTOMInfo$consensusTOM = NULL
  individualTOMInfo$TOMSimilarities = NULL

  list(colors = mergedColors,
       unmergedColors = colors,
       multiMEs = allSampleMEs,
       goodSamples = gsg$goodSamples,
       goodGenes = gsg$goodGenes,
       dendrograms = dendros,
       TOMFiles = consensusTOMInfo$TOMFiles,
       blockGenes = individualTOMInfo$blockGenes,
       blocks = blocks,
       originCount = consensusTOMInfo$originCount,
       networkCalibrationSamples = consensusTOMInfo$networkCalibrationSamples,
       individualTOMInfo = individualTOMInfo,
       consensusTOMInfo = if (saveConsensusTOMs) consensusTOMInfo else NULL,
       consensusQuantile = consensusQuantile
      )
}


#==========================================================================================================
#
# recutConsensusTrees
#
#==========================================================================================================



#' Repeat blockwise consensus module detection from pre-calculated data
#'
#' Given consensus networks constructed for example using
#' \code{\link{blockwiseConsensusModules}}, this function (re-)detects modules
#' in them by branch cutting of the corresponding dendrograms. If repeated
#' branch cuts of the same gene network dendrograms are desired, this function
#' can save substantial time by re-using already calculated networks and
#' dendrograms.
#'
#'
#' For details on blockwise consensus module detection, see
#' \code{\link{blockwiseConsensusModules}}. This function implements the module
#' detection subset of the functionality of
#' \code{\link{blockwiseConsensusModules}}; network construction and clustering
#' must be performed in advance. The primary use of this function is to
#' experiment with module detection settings without having to re-execute long
#' network and clustering calculations whose results are not affected by the
#' cutting parameters.
#'
#' This function takes as input the networks and dendrograms that are produced
#' by \code{\link{blockwiseConsensusModules}}.  Working block by block, modules
#' are identified in the dendrograms by the Dynamic Hybrid tree cut.  Found
#' modules are trimmed of genes whose consensus module membership kME (that is,
#' correlation with module eigengene) is less than \code{minKMEtoStay}. Modules
#' in which fewer than \code{minCoreKMESize} genes have consensus KME higher
#' than \code{minCoreKME} are disbanded, i.e., their constituent genes are
#' pronounced unassigned.
#'
#' After all blocks have been processed, the function checks whether there are
#' genes whose KME in the module they assigned is lower than KME to another
#' module. If p-values of the higher correlations are smaller than those of the
#' native module by the factor \code{reassignThresholdPS} (in every set), the
#' gene is re-assigned to the closer module.
#'
#' In the last step, modules whose eigengenes are highly correlated are merged.
#' This is achieved by clustering module eigengenes using the dissimilarity
#' given by one minus their correlation, cutting the dendrogram at the height
#' \code{mergeCutHeight} and merging all modules on each branch. The process is
#' iterated until no modules are merged. See \code{\link{mergeCloseModules}}
#' for more details on module merging.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param goodSamples a list with one component per set. Each component is a
#' logical vector specifying which samples are considered "good" for the
#' analysis. See \code{\link{goodSamplesGenesMS}}.
#' @param goodGenes a logical vector with length equal number of genes in
#' \code{multiExpr} that specifies which genes are considered "good" for the
#' analysis. See \code{\link{goodSamplesGenesMS}}.
#' @param blocks specification of blocks in which hierarchical clustering and
#' module detection should be performed. A numeric vector with one entry per
#' gene of \code{multiExpr} giving the number of the block to which the
#' corresponding gene belongs.
#' @param TOMFiles a vector of character strings specifying file names in which
#' the block-wise topological overlaps are saved.
#' @param dendrograms a list of length equal the number of blocks, in which
#' each component is a hierarchical clustering dendrograms of the genes that
#' belong to the block.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pariwise.complete.obs} option.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}. Note that while no new networks are computed in
#' this function, this parameter affects the interpretation of correlations in
#' this function.
#' @param deepSplit integer value between 0 and 4. Provides a simplified
#' control over how sensitive module detection should be to module splitting,
#' with 0 least and 4 most sensitive. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param detectCutHeight dendrogram cut height for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minModuleSize minimum module size for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param checkMinModuleSize logical: should sanity checks be performed on
#' \code{minModuleSize}?
#' @param maxCoreScatter maximum scatter of the core for a branch to be a
#' cluster, given as the fraction of \code{cutHeight} relative to the 5th
#' percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minGap minimum cluster gap given as the fraction of the difference
#' between \code{cutHeight} and the 5th percentile of joining heights. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param maxAbsCoreScatter maximum scatter of the core for a branch to be a
#' cluster given as absolute heights. If given, overrides
#' \code{maxCoreScatter}. See \code{\link[dynamicTreeCut]{cutreeDynamic}} for
#' more details.
#' @param minAbsGap minimum cluster gap given as absolute height difference. If
#' given, overrides \code{minGap}. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minSplitHeight Minimum split height given as the fraction of the
#' difference between \code{cutHeight} and the 5th percentile of joining
#' heights. Branches merging below this height will automatically be merged.
#' Defaults to zero but is used only if \code{minAbsSplitHeight} below is
#' \code{NULL}.
#' @param minAbsSplitHeight Minimum split height given as an absolute height.
#' Branches merging below this height will automatically be merged. If not
#' given (default), will be determined from \code{minSplitHeight} above.
#' @param useBranchEigennodeDissim Logical: should branch eigennode (eigengene)
#' dissimilarity be considered when merging branches in Dynamic Tree Cut?
#' @param minBranchEigennodeDissim Minimum consensus branch eigennode
#' (eigengene) dissimilarity for branches to be considerd separate. The branch
#' eigennode dissimilarity in individual sets is simly 1-correlation of the
#' eigennodes; the consensus is defined as quantile with probability
#' \code{consensusQuantile}.
#' @param pamStage logical.  If TRUE, the second (PAM-like) stage of module
#' detection will be performed.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param pamRespectsDendro Logical, only used when \code{pamStage} is
#' \code{TRUE}. If \code{TRUE}, the PAM stage will respect the dendrogram in
#' the sense an object can be PAM-assigned only to clusters that lie below it
#' on the branch that the object is merged into.  See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param trimmingConsensusQuantile a number between 0 and 1 specifying the
#' consensus quantile used for kME calculation that determines module trimming
#' according to the arguments below.
#' @param minCoreKME a number between 0 and 1. If a detected module does not
#' have at least \code{minModuleKMESize} genes with eigengene connectivity at
#' least \code{minCoreKME}, the module is disbanded (its genes are unlabeled
#' and returned to the pool of genes waiting for mofule detection).
#' @param minCoreKMESize see \code{minCoreKME} above.
#' @param minKMEtoStay genes whose eigengene connectivity to their module
#' eigengene is lower than \code{minKMEtoStay} are removed from the module.
#' @param reassignThresholdPS per-set p-value ratio threshold for reassigning
#' genes between modules.  See Details.
#' @param mergeCutHeight dendrogram cut height for module merging.
#' @param mergeConsensusQuantile consensus quantile for module merging. See
#' \code{mergeCloseModules} for details.
#' @param impute logical: should imputation be used for module eigengene
#' calculation? See \code{\link{moduleEigengenes}} for more details.
#' @param trapErrors logical: should errors in calculations be trapped?
#' @param numericLabels logical: should the returned modules be labeled by
#' colors (\code{FALSE}), or by numbers (\code{TRUE})?
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{colors}{ module assignment of all input genes. A vector containing
#' either character strings with module colors (if input \code{numericLabels}
#' was unset) or numeric module labels (if \code{numericLabels} was set to
#' \code{TRUE}). The color "grey" and the numeric label 0 are reserved for
#' unassigned genes.  }
#'
#' \item{unmergedColors }{ module colors or numeric labels before the module
#' merging step. }
#'
#' \item{multiMEs}{ module eigengenes corresponding to the modules returned in
#' \code{colors}, in multi-set format. A vector of lists, one per set,
#' containing eigengenes, proportion of variance explained and other
#' information. See \code{\link{multiSetMEs}} for a detailed description. }
#' @note Basic sanity checks are performed on given arguments, but it is left
#' to the user's responsibility to provide valid input.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{blockwiseConsensusModules}} for the full blockwise modules
#' calculation. Parts of its output are natural input for this function.
#'
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for adaptive branch cutting in
#' hierarchical clustering dendrograms
#'
#' \code{\link{mergeCloseModules}} for merging of close modules.
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
recutConsensusTrees = function(multiExpr,
                            goodSamples, goodGenes,
                            blocks,
                            TOMFiles,
                            dendrograms,
                            corType = "pearson",
                            networkType = "unsigned",
                            deepSplit = 2,
                            detectCutHeight = 0.995, minModuleSize = 20,
                            checkMinModuleSize = TRUE,
                            maxCoreScatter = NULL, minGap = NULL,
                            maxAbsCoreScatter = NULL, minAbsGap = NULL,
                            minSplitHeight = NULL, minAbsSplitHeight = NULL,

                            useBranchEigennodeDissim = FALSE,
                            minBranchEigennodeDissim = mergeCutHeight,

                            pamStage = TRUE,  pamRespectsDendro = TRUE,
                            # minKMEtoJoin =0.7,
                            trimmingConsensusQuantile = 0,
                            minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
                            minKMEtoStay = 0.2,
                            reassignThresholdPS = 1e-4,
                            mergeCutHeight = 0.15,
                            mergeConsensusQuantile = trimmingConsensusQuantile,
                            impute = TRUE,
                            trapErrors = FALSE,
                            numericLabels = FALSE,
                            verbose = 2, indent = 0)
{
  spaces = indentSpaces(indent)

  dataSize = checkSets(multiExpr)
  nSets = dataSize$nSets
  nGenes = dataSize$nGenes
  nSamples = dataSize$nSamples

  if (length(blocks)!=nGenes)
    stop("Input error: length of 'blocks' must equal number of genes in 'multiExpr'.")

  #if (verbose>0)
  #   printFlush(paste(spaces, "Calculating consensus modules and module eigengenes",
  #                    "block-wise from all genes"))

  # If we're merging branches by correlation within cutreeHybrid, prepare scaled and imputed multiExpr.

  if (useBranchEigennodeDissim)
  {
    multiExpr.scaled = apply(multiExpr, scale)
    hasMissing = unlist(multiSet2list(apply(multiExpr, function(x) { any(is.na(x)) })))
    # Impute those that have missing data
    multiExpr.scaled.imputed = multiSet.mapply(function(x, doImpute)
                           { if (doImpute) t(impute.knn(t(x))$data) else x },
                                     multiExpr.scaled, hasMissing)
    branchSplitFnc = "multiSet.branchEigengeneDissim"
  }


  intCorType = pmatch(corType, .corTypes)
  if (is.na(intCorType))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  allLabels = rep(0, nGenes)
  allLabelIndex = NULL

  corFnc = match.fun(.corFnc[intCorType])
  corOptions = list(use = 'p')

  # Get rid of bad genes and bad samples

  gsg = list(goodGenes = goodGenes, goodSamples = goodSamples, allOK = TRUE)

  gsg$allOK = (sum(!gsg$goodGenes)==0)
  nGGenes = sum(gsg$goodGenes)
  nGSamples = rep(0, nSets)
  for (set in 1:nSets)
  {
    nGSamples[set] = sum(gsg$goodSamples[[set]])
    gsg$allOK  = gsg$allOK && (sum(!gsg$goodSamples[[set]])==0)
  }

  if (!gsg$allOK)
    for (set in 1:nSets)
      multiExpr[[set]]$data = multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes]

  gBlocks = blocks[gsg$goodGenes]

  blockLevels = as.numeric(levels(factor(gBlocks)))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  reassignThreshold = reassignThresholdPS^nSets

  # prepare scaled and imputed multiExpr.
  multiExpr.scaled = apply(multiExpr, scale)
  hasMissing = unlist(multiSet2list(apply(multiExpr, function(x) { any(is.na(x)) })))
  # Impute those that have missing data
  multiExpr.scaled.imputed = multiSet.mapply(function(x, doImpute)
                         { if (doImpute) t(impute.knn(t(x))$data) else x },
                                   multiExpr.scaled, hasMissing)
  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = "multiSet.branchEigengeneDissim"
  } else
    branchSplitFnc = NULL

  # Initialize various variables

  consMEs = vector(mode = "list", length = nSets)

  blockNo = 1
  maxUsedLabel = 0
  collectGarbage()
  # Here's where the analysis starts

  while (blockNo <= nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select most connected genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    nBlockGenes = length(block)
    selExpr = vector(mode = "list", length = nSets)
    for (set in 1:nSets)
      selExpr[[set]] = list(data = multiExpr[[set]]$data[, block])

    # This is how TOMs are saved:
    #if (saveTOMs)
    #{
    #   TOMFiles[blockNo] = paste(saveTOMFileBase, "-block.", blockNo, ".RData", sep="")
    #   save(consTomDS, file = TOMFiles[blockNo])
    #}
    #consTomDS = 1-consTomDS
    #collectGarbage()

    xx = try(load(TOMFiles[blockNo]), silent = TRUE)
    if (class(xx)=='try-error')
    {
      printFlush(paste("************\n File name", TOMFiles[blockNo],
                       "appears invalid: the load function returned the following error:\n     ",
                       xx))
      stop()
    }
    if (xx!='consTomDS')
      stop(paste("The file", TOMFiles[blockNo], "does not contain the appopriate variable."))

    if (class(consTomDS)!="dist")
      stop(paste("The file", TOMFiles[blockNo], "does not contain the appopriate distance structure."))

    consTomDS = 1-consTomDS
    collectGarbage()

    if (verbose>2) printFlush(paste(spaces, "....clustering and detecting modules.."))
    errorOccured = FALSE

    blockLabels = try(cutreeDynamic(dendro = dendrograms[[blockNo]],
                    deepSplit = deepSplit,
                    cutHeight = detectCutHeight, minClusterSize = minModuleSize,
                    method ="hybrid",
                    maxCoreScatter = maxCoreScatter, minGap = minGap,
                    maxAbsCoreScatter = maxAbsCoreScatter, minAbsGap = minAbsGap,
                    minSplitHeight = minSplitHeight, minAbsSplitHeight = minAbsSplitHeight,

                    externalBranchSplitFnc = if (useBranchEigennodeDissim)
                                                branchSplitFnc else NULL,
                    minExternalSplit = minBranchEigennodeDissim,
                    externalSplitOptions = list(multiExpr = multiExpr.scaled.imputed[ , block],
                                                corFnc = corFnc, corOptions = corOptions,
                                                consensusQuantile = mergeConsensusQuantile,
                                                signed = networkType %in% c("signed", "signed hybrid")),
                    externalSplitFncNeedsDistance = FALSE,

                    pamStage = pamStage, pamRespectsDendro = pamRespectsDendro,
                    distM = as.matrix(consTomDS),
                    verbose = verbose-3, indent = indent + 2), silent = TRUE)
    if (class(blockLabels)=='try-error')
    {
      (if (verbose>0) printFlush else warning)
           (paste(spaces, "blockwiseConsensusModules: cutreeDynamic failed:\n    ",
                  blockLabels, "\nError occured in block", blockNo, "\nContinuing with next block."))
      next
    } else {
      blockLabels[blockLabels>0] = blockLabels[blockLabels>0] + maxUsedLabel
      maxUsedLabel = max(blockLabels)
    }
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
          printFlush(paste(spaces, "No modules detected in block", blockNo))
          printFlush(paste(spaces, "  Continuing with next block."))
      }
      next
    }

    # Calculate eigengenes for this batch

    if (verbose>2) printFlush(paste(spaces, "....calculating eigengenes.."))
    blockAssigned = c(1:nBlockGenes)[blockLabels!=0]
    blockLabelIndex = as.numeric(levels(as.factor(blockLabels[blockAssigned])))
    blockConsMEs = try(multiSetMEs(selExpr, universalColors = blockLabels,
                                   excludeGrey = TRUE, grey = 0, impute = impute,
                                   # trapErrors = TRUE, returnValidOnly = TRUE,
                                   verbose = verbose-4, indent = indent + 3), silent = TRUE)
    if (class(blockConsMEs)=='try-error')
    {
      if (verbose>0)
      {
        printFlush(paste(spaces, "*** multiSetMEs failed with the message:"))
        printFlush(paste(spaces, "     ", blockConsMEs))
        printFlush(paste(spaces, "*** --> Ending module detection here"))
      } else warning(paste("blocwiseConsensusModules: multiSetMEs failed with the message: \n",
               "      ", blockConsMEs, "\n--> Continuing with next block."))
      next
    }

    deleteModules = NULL
    changedModules = NULL

    # find genes whose closest module eigengene has cor higher than minKMEtoJoin and assign them
    # Removed - should not change blocks before clustering them
    #unassGenes = c(c(1:nGGenes)[-block][allLabels[-block]==0], block[blockLabels==0])
    #if (length(unassGenes) > 0)
    #{
      #blockKME = array(0, dim = c(length(unassGenes), ncol(blockConsMEs[[1]]$data), nSets))
      #corEval = parse(text = paste(.corFnc[intCorType],
                         #"( multiExpr[[set]]$data[, unassGenes], blockConsMEs[[set]]$data,",
                         #.corOptions[intCorType], ")"))
      #for (set in 1:nSets) blockKME[, , set] = eval(corEval)
      #if (intNetworkType==1) blockKME = abs(blockKME)
      #consKME = as.matrix(apply(blockKME, c(1,2), min))
      #consKMEmax = apply(consKME, 1, max)
      #closestModule = blockLabelIndex[apply(consKME, 1, which.max)]
      #assign = (consKMEmax >= minKMEtoJoin )
      #if (sum(assign>0))
      #{
        #allLabels[unassGenes[assign]] = closestModule[assign]
        #changedModules = union(changedModules, closestModule[assign])
      #}
      #rm(blockKME, consKME, consKMEmax)
    #}

    collectGarbage()

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
    # the required minimum consensus KME

    if (verbose>2)
      printFlush(paste(spaces, "....checking consensus modules for statistical meaningfulness.."))

    for (mod in 1:ncol(blockConsMEs[[1]]$data))
    {
      modGenes = (blockLabels==blockLabelIndex[mod])
      KME = matrix(0, nrow = sum(modGenes), ncol = nSets)
      corEval = parse(text = paste(.corFnc[intCorType],
                       "( selExpr[[set]]$data[, modGenes], blockConsMEs[[set]]$data[, mod]",
                      prepComma(.corOptions[intCorType]), ")"))
      for (set in 1:nSets) KME[, set] = as.vector(eval(corEval))
      if (intNetworkType==1) KME = abs(KME)
      consKME = apply(KME, 1, quantile, probs = trimmingConsensusQuantile, na.rm = TRUE)
      if (sum(consKME>minCoreKME) < minCoreKMESize)
      {
        blockLabels[modGenes] = 0
        deleteModules = union(deleteModules, mod)
        if (verbose>3)
          printFlush(paste0(spaces, "    ..deleting module ",blockLabelIndex[mod],
                           ": of ", sum(modGenes),
                     " total genes in the module only ",  sum(consKME>minCoreKME),
                     " have the requisite high correlation with the eigengene in all sets."))
      } else if (sum(consKME<minKMEtoStay)>0)
      {
        if (verbose > 3)
          printFlush(paste(spaces, "    ..removing", sum(consKME<minKMEtoStay),
                           "genes from module", blockLabelIndex[mod], "because their KME is too low."))
        blockLabels[modGenes][consKME < minKMEtoStay] = 0
        if (sum(blockLabels[modGenes]>0) < minModuleSize)
        {
          deleteModules = union(deleteModules, mod)
          blockLabels[modGenes] = 0
          if (verbose>3)
            printFlush(paste0(spaces, "    ..deleting module ",blockLabelIndex[mod],
                     ": not enough genes in the module after removal of low KME genes."))
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod])
        }
      }
    }

    # Remove marked modules

    if (!is.null(deleteModules))
    {
       for (set in 1:nSets) blockConsMEs[[set]]$data = blockConsMEs[[set]]$data[, -deleteModules]
       modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
       blockLabels[modGenes] = 0
       modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
       allLabels[modAllGenes] = 0
       blockLabelIndex = blockLabelIndex[-deleteModules]
    }

    # Check whether there's anything left
    if (sum(blockLabels>0)==0)
    {
      if (verbose>1)
      {
        printFlush(paste(spaces, " No significant modules detected in block", blockNo))
        printFlush(paste(spaces, " Continuing with next block."))
      }
      next
    }

    # Update module eigengenes

    for (set in 1:nSets)
      if (is.null(dim(blockConsMEs[[set]]$data)))
        dim(blockConsMEs[[set]]$data) = c(length(blockConsMEs[[set]]$data), 1)

    if (is.null(consMEs[[1]]))
    {
       for (set in 1:nSets) consMEs[[set]] = list(data = blockConsMEs[[set]]$data)
    } else for (set in 1:nSets)
       consMEs[[set]]$data = cbind(consMEs[[set]]$data, blockConsMEs[[set]]$data)

    # Update allLabels

    allLabelIndex = c(allLabelIndex, blockLabelIndex)
    allLabels[block[blockAssigned]] = blockLabels[blockAssigned]

    collectGarbage()

    blockNo = blockNo + 1

  }

  # Check whether any of the already assigned genes (in this or previous blocks) should be re-assigned

  if (verbose>2)
    printFlush(paste(spaces, "....checking for genes that should be reassigned.."))

  deleteModules = NULL
  goodLabels = allLabels[gsg$goodGenes]
  if (sum(goodLabels!=0) > 0)
  {
     propLabels = goodLabels[goodLabels!=0]
     assGenes = (c(1:nGenes)[gsg$goodGenes])[goodLabels!=0]
     corEval = parse(text = paste(.corFnc[intCorType],
                                  "(multiExpr[[set]]$data[, goodLabels!=0], consMEs[[set]]$data",
                                  prepComma(.corOptions[intCorType]), ")"))
     nMods = ncol(consMEs[[1]]$data)
     lpValues = array(0, dim = c(length(propLabels), nMods, nSets))
     sumSign = array(0, dim = c(length(propLabels), nMods))
     if (verbose>3)
       printFlush(paste(spaces, "......module membership p-values.."))
     for (set in 1:nSets)
     {
       KME = eval(corEval)
       if (intNetworkType == 1) KME = abs(KME)
       lpValues[,,set] = -2*log(corPvalueFisher(KME, nGSamples[set], twoSided = FALSE))
       sumSign = sumSign + sign(KME)
     }
     if (verbose>3)
       printFlush(paste(spaces, "......module membership scores.."))
     scoreAll = as.matrix(apply(lpValues, c(1,2), sum)) * (nSets + sumSign)/(2*nSets)
     scoreAll[!is.finite(scoreAll)] = 0.001 # This low should be enough
     bestScore = apply(scoreAll, 1, max)
     if (intNetworkType==1) sumSign = abs(sumSign)
     if (verbose>3)
     {
       cat(paste(spaces, "......individual modules.."))
       pind = initProgInd()
     }
     for (mod in 1:nMods)
     {
       modGenes = c(1:length(propLabels))[propLabels==allLabelIndex[mod]]
       scoreMod = scoreAll[modGenes, mod]
       candidates = (bestScore[modGenes] > scoreMod)
       candidates[!is.finite(candidates)] = FALSE
       if (sum(candidates) > 0)
       {
         pModule = pchisq(scoreMod[candidates], nSets, log.p = TRUE)
         whichBest = apply(scoreAll[modGenes[candidates], ], 1, which.max)
         pBest = pchisq(bestScore[modGenes[candidates]], nSets, log.p = TRUE)
         reassign =  ifelse(is.finite(pBest - pModule),
                            ( (pBest - pModule) < log(reassignThreshold) ),
                            FALSE)
         if (sum(reassign)>0)
         {
           allLabels[assGenes[modGenes[candidates][reassign]]] = whichBest[reassign]
           changedModules = union(changedModules, whichBest[reassign])
           if (sum(modGenes)-sum(reassign) < minModuleSize)
           {
             deleteModules = union(deleteModules, mod)
           } else
               changedModules = union(changedModules, mod)
         }
       }
       if (verbose > 3) pind = updateProgInd(mod/nMods, pind)
     }
     rm(lpValues, sumSign, scoreAll)
     if (verbose > 3) printFlush("")
  }

  # Remove marked modules

  if (!is.null(deleteModules))
  {
     # for (set in 1:nSets) consMEs[[set]]$data = consMEs[[set]]$data[, -deleteModules]
     modGenes = is.finite(match(allLabels, allLabelIndex[deleteModules]))
     allLabels[modGenes] = 0
     # allLabelIndex = allLabelIndex[-deleteModules]
  }

  if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."))
  #print(table(allLabels))
  #print(is.numeric(allLabels))
  if (numericLabels) {
    colors = allLabels
  } else {
    colors = labels2colors(allLabels)
  }
  mergedColors = colors
  mergedMods = try(mergeCloseModules(multiExpr, colors[gsg$goodGenes],
                                     consensusQuantile = mergeConsensusQuantile,
                                     cutHeight = mergeCutHeight,
                                     relabel = TRUE, impute = impute,
                                     verbose = verbose-2, indent = indent + 2), silent = TRUE)
  if (class(mergedMods)=='try-error')
  {
    if (verbose>0)
    {
      printFlush(paste(spaces, 'blockwiseConsensusModules: mergeCloseModule failed with this message:\n',
            spaces, '    ', mergedMods, spaces,
            '---> returning unmerged consensus modules'))
    } else warning(paste('blockwiseConsensusModules: mergeCloseModule failed with this message:\n     ',
                          mergedMods, '---> returning unmerged consensus modules'))
    MEs = try(multiSetMEs(multiExpr, universalColors = colors[gsg$goodGenes]
                          # trapErrors = TRUE, returnValidOnly = TRUE
                          ), silent = TRUE)
    if (class(MEs)=='try-error')
    {
      warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
            MEs, '---> returning empty module eigengenes'))
      allSampleMEs = NULL
    } else {
      if (!MEs[[1]]$allOK) mergedColors[gsg$goodGenes] = MEs[[1]]$validColors
      allSampleMEs = vector(mode = "list", length = nSets)
      for (set in 1:nSets)
      {
        allSampleMEs[[set]] =
           list(data = as.data.frame(matrix(NA, nrow = nGSamples[set], ncol = ncol(MEs[[set]]$data))))
        allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = MEs[[set]]$data[,]
        names(allSampleMEs[[set]]$data) = names(MEs[[set]]$data)
      }
    }
  } else {
    mergedColors[gsg$goodGenes] = mergedMods$colors
    allSampleMEs = vector(mode = "list", length = nSets)
    for (set in 1:nSets)
    {
      allSampleMEs[[set]] =
         list(data = as.data.frame(matrix(NA, nrow = nGSamples[set],
                                          ncol = ncol(mergedMods$newMEs[[1]]$data))))
      allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = mergedMods$newMEs[[set]]$data[,]
      names(allSampleMEs[[set]]$data) = names(mergedMods$newMEs[[set]]$data)
    }
  }

  list(colors = mergedColors,
       unmergedColors = colors,
       multiMEs = allSampleMEs
  #     goodSamples = gsg$goodSamples,
  #     goodGenes = gsg$goodGenes,
  #     dendrograms = dendros,
  #     blockGenes = blockGenes,
  #     originCount = originCount,
  #     TOMFiles = TOMFiles
      )
}



#======================================================================================================
#
# preliminary partitioning
#
#======================================================================================================



#' Projective K-means (pre-)clustering of expression data
#'
#' Implementation of a variant of K-means clustering for expression data.
#'
#'
#' The principal aim of this function within WGCNA is to pre-cluster a large
#' number of genes into smaller blocks that can be handled using standard WGCNA
#' techniques.
#'
#' This function implements a variant of K-means clustering that is suitable
#' for co-expression analysis. Cluster centers are defined by the first
#' principal component, and distances by correlation (more precisely,
#' 1-correlation). The distance between a gene and a cluster is multiplied by a
#' factor of \eqn{max(clusterSize/preferredSize,
#' 1)^{sizePenaltyPower}}{\code{max(clusterSize/preferredSize,
#' 1)^sizePenaltyPower}}, thus penalizing clusters whose size exceeds
#' \code{preferredSize}. The function starts with randomly generated cluster
#' assignment (hence the need to set the random seed for repeatability) and
#' executes interations of calculating new centers and reassigning genes to
#' nearest center until the clustering becomes stable. Before returning, nearby
#' clusters are iteratively combined if their combined size is below
#' \code{preferredSize}.
#'
#' The standard principal component calculation via the function \code{svd}
#' fails from time to time (likely a convergence problem of the underlying
#' lapack functions). Such errors are trapped and the principal component is
#' approximated by a weighted average of expression profiles in the cluster. If
#' \code{verbose} is set above 2, an informational message is printed whenever
#' this approximation is used.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. NAs are allowed, but not too many.
#' @param preferredSize preferred maximum size of clusters.
#' @param nCenters number of initial clusters. Empirical evidence suggests that
#' more centers will give a better preclustering; the default is an attempt to
#' arrive at a reasonable number.
#' @param sizePenaltyPower parameter specifying how severe is the penalty for
#' clusters that exceed \code{preferredSize}.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param randomSeed integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit.
#' @param checkData logical: should data be checked for genes with zero
#' variance and genes and samples with excessive numbers of missing samples?
#' Bad samples are ignored; returned cluster assignment for bad genes will be
#' \code{NA}.
#' @param imputeMissing logical: should missing values in \code{datExpr} be
#' imputed before the calculations start? The early imputation makes the code
#' run faster but may produce slightly different results if re-running older
#' calculations.
#' @param maxIterations maximum iterations to be attempted.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{clusters}{ a numerical
#' vector with one component per input gene, giving the cluster number in which
#' the gene is assigned. } \item{centers}{ cluster centers, that is their first
#' principal components. }
#' @author Peter Langfelder
#' @keywords cluster
projectiveKMeans = function (
      datExpr,
      preferredSize = 5000,
      nCenters = as.integer(min(ncol(datExpr)/20, preferredSize^2/ncol(datExpr))),
      sizePenaltyPower = 4,
      networkType = "unsigned",
      randomSeed = 54321,
      checkData = TRUE,
      imputeMissing = TRUE,
      maxIterations = 1000,
      verbose = 0, indent = 0 )
{

  centerGeneDist = function(centers, oldDst = NULL, changed = c(1:nCenters),
                          blockSize = 50000, verbose = 0, spaces = "")
  {
    if (is.null(oldDst))
    {
      oldDst = array(0, c(nCenters, nGenes))
      changed = c(1:nCenters)
    }
    dstAll = oldDst
    nChanged = length(changed)

    nBlocks = ceiling(ncol(datExpr)/blockSize)
    blocks = allocateJobs(ncol(datExpr), nBlocks)

    if (verbose > 5)
      pind = initProgInd(paste0(spaces, "   ..centerGeneDist: "))

    for (b in 1:nBlocks)
    {
      if (intNetworkType==1)
      {
         dst = 1-abs(cor(centers[, changed], datExpr[, blocks[[b]] ]))
      } else {
         dst = 1-cor(centers[, changed], datExpr[, blocks[[b]] ])
      }
      dstAll[changed, blocks[[b]]] = dst
      if (verbose > 5) pind = updateProgInd(b/nBlocks, pind)
    }
    dstAll
  }

  memberCenterDist = function(dst, membership, sizes = NULL, changed = c(1:nCenters), oldMCDist = NULL)
  {
    if (is.null(oldMCDist))
    {
       changed = c(1:nCenters)
       oldMCDist = rep(0, nGenes)
    }
    centerDist = oldMCDist
    if (!is.null(changed))
    {
      if (is.null(sizes)) sizes = table(membership)
      if (length(sizes)!=nCenters)
      {
        sizes2 = rep(0, nCenters)
        sizes2[as.numeric(names(sizes))] = sizes
        sizes = sizes2
      }
      if (is.finite(sizePenaltyPower))
      {
        sizeCorrections = (sizes/preferredSize)^sizePenaltyPower
        sizeCorrections[sizeCorrections < 1] = 1
      } else {
        sizeCorrections = rep(1, length(sizes))
        sizeCorrections[sizes > preferredSize] = Inf
      }
      for (cen in changed) if (sizes[cen]!=0)
      {
        if (is.finite(sizeCorrections[cen]))
        {
          centerDist[membership==cen] = dst[cen, membership==cen] * sizeCorrections[cen]
        } else
          centerDist[membership==cen] = 10 + dst[cen, membership==cen]
      }
    }
    centerDist
  }

  spaces = indentSpaces(indent)
  if (verbose > 0)
    printFlush(paste(spaces, "Projective K-means:"))

  datExpr = scale(as.matrix(datExpr))

  if (any(is.na(datExpr)))
  {
    if (imputeMissing)
    {
      printFlush(paste0(spaces, "projectiveKMeans: imputing missing data in 'datExpr'.\n",
                  "To reproduce older results, use 'imputeMissing = FALSE'. "))
      datExpr = t(impute.knn(t(datExpr))$data)
    } else {
      printFlush(paste0(spaces, "projectiveKMeans: there are missing data in 'datExpr'.\n",
          "SVD will not work; will use a weighted mean approximation."))
    }
  }

  #if (preferredSize >= floor(sqrt(2^31)) )
  #  stop("'preferredSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (exists(".Random.seed"))
  {
     seedSaved = TRUE
     savedSeed = .Random.seed
  } else seedSaved = FALSE
  set.seed(randomSeed)

  nAllSamples = nrow(datExpr)
  nAllGenes = ncol(datExpr)

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  if (verbose > 1)
    printFlush(paste(spaces, "..using", nCenters, "centers."))

  # Check data for genes and samples that have too many missing values

  if (checkData)
  {
    if (verbose > 0)
      printFlush(paste(spaces, "..checking data for excessive number of missing values.."))
    gsg = goodSamplesGenes(datExpr, verbose = verbose -1, indent = indent + 1)
    if (!gsg$allOK) datExpr = datExpr[gsg$goodSamples, gsg$goodGenes]
  }
  nGenes = ncol(datExpr)
  nSamples = nrow(datExpr)

  datExpr[is.na(datExpr)] = 0

  centers = matrix(0, nSamples, nCenters)

  randGeneIndex = sample(nGenes, size = nGenes)
  temp = rep(c(1:nCenters), times = ceiling(nGenes/nCenters))
  membership = temp[randGeneIndex]

  if (verbose > 0)
    printFlush(paste(spaces, "..k-means clustering.."))

  changed = c(1:nCenters); dst = NULL; centerDist = NULL
  iteration = 0
  while (!is.null(changed) && (iteration < maxIterations))
  {
    iteration = iteration + 1
    if (verbose > 1) printFlush(paste(spaces, " ..iteration", iteration))
    clusterSizes = table(membership)
    if (verbose > 5) pind = initProgInd(paste(spaces, " ....calculating centers: "))
    for (cen in sort(changed))
    {
        centers[, cen] = .alignedFirstPC(datExpr[, membership==cen], verbose = verbose-2,
                                         indent = indent+2)
        if (verbose > 5) pind = updateProgInd(cen/nCenters, pind)
    }
    if (verbose > 5) { pind = updateProgInd(1, pind); printFlush("")}

    if (verbose > 5) printFlush(paste(spaces, " ....calculating center to gene distances"))
    dst = centerGeneDist(centers, dst, changed, verbose = verbose, spaces = spaces)
    centerDist = memberCenterDist(dst, membership, clusterSizes, changed, centerDist)

    nearestDist = rep(0, nGenes)
    nearest = rep(0, nGenes)
    if (verbose > 5) printFlush(paste(spaces, " ....finding nearest center for each gene"))
    minRes = .C("minWhichMin", as.double(dst), as.integer(nCenters), as.integer(nGenes),
                as.double(nearestDist), as.double(nearest), PACKAGE = "WGCNA")
    nearestDist = minRes[[4]]
    nearest = minRes[[5]]+1
    rm(minRes); collectGarbage()

    if (sum(centerDist>nearestDist)>0)
    {
      proposedMemb = nearest
      accepted = FALSE
      while (!accepted && (sum(proposedMemb!=membership)>0))
      {
        if (verbose > 2)
          cat(paste(spaces, "   ..proposing to move", sum(proposedMemb!=membership), "genes"))
        moved = c(1:nGenes)[proposedMemb!=membership]
        newCentDist = memberCenterDist(dst, proposedMemb)
        gotWorse = newCentDist[moved] > centerDist[moved]
        if (sum(!is.finite(gotWorse))>0)
           warning("Have NAs in gotWorse.")
        if (sum(gotWorse)==0)
        {
          accepted = TRUE
          if (verbose > 2) printFlush(paste("..move accepted."))
        } else {
          if (verbose > 2) printFlush(paste("..some genes got worse. Trying again."))
          ord = order(centerDist[moved[gotWorse]] - newCentDist[moved[gotWorse]])
          n = ceiling(length(ord)*3/5)
          proposedMemb[moved[gotWorse][ord[c(1:n)]]] = membership[moved[gotWorse][ord[c(1:n)]]]
        }
      }
      if (accepted)
      {
        propSizes = table(proposedMemb)
        keep = as.numeric(names(propSizes))
        centers = centers[, keep]
        dst = dst[keep, ]
        changedAll = union(membership[moved], proposedMemb[moved])
        changedKeep = changedAll[is.finite(match(changedAll, keep))]
        changed = rank(changedKeep);    # This is a way to make say 1,3,4 into 1,2,3
        membership = as.numeric(as.factor(proposedMemb))
        if ( (verbose > 1) && (sum(keep) < nCenters))
          printFlush(paste(spaces, "  ..dropping", nCenters - sum(keep),
                           "centers because ther clusters are empty."))
        nCenters = length(keep)
      } else {
        changed = NULL
        if (verbose > 2) printFlush(paste("Could not find a suitable move to improve the clustering."))
      }
    } else {
      changed = NULL
      if (verbose > 2)
        printFlush(paste("Clustering is stable: all genes are closest to their assigned center."))
    }
    if (verbose > 5)
    {
      printFlush("Sizes of biggest preliminary clusters:")
      order = order(-as.numeric(clusterSizes))
      print(as.numeric(clusterSizes)[order[1:min(100, length(order))]])
    }
  }

  if (verbose > 2 & verbose <6)
  {
    printFlush("Sizes of preliminary clusters:")
    print(clusterSizes)
  }



  # merge nearby clusters if their sizes allow merging
  if (verbose > 0) printFlush(paste(spaces, "..merging smaller clusters..."))
  small = (clusterSizes < preferredSize)
  done = FALSE
  while (!done & (sum(small)>1))
  {
    smallIndex = c(1:nCenters)[small]
    nSmall = sum(small)
    if (intNetworkType==1)
    {
       clustDist = 1-abs(cor(centers[, small]))
    } else {
       clustDist = 1-cor(centers[, small])
    }

    diag(clustDist) = 10
    distOrder = order(as.vector(clustDist))[seq(from=2, to = nSmall * (nSmall-1), by=2)]
    i = 1; canMerge = FALSE
    while (i <= length(distOrder) && (!canMerge))
    {
      whichJ = smallIndex[as.integer( (distOrder[i]-1)/nSmall + 1)]
      whichI = smallIndex[distOrder[i] - (whichJ-1) * nSmall]
      canMerge = sum(clusterSizes[c(whichI, whichJ)]) < preferredSize
      i = i + 1
    }
    if (canMerge)
    {
      membership[membership==whichJ] = whichI
      clusterSizes[whichI] = clusterSizes[whichI] + clusterSizes[whichJ]
      centers[, whichI] = .alignedFirstPC(datExpr[, membership==whichI], verbose = verbose-2,
                                          indent = indent+2)
      nCenters = nCenters -1
      if (verbose > 3)
        printFlush(paste(spaces, "  ..merged clusters", whichI, "and", whichJ,
                   "whose combined size is", clusterSizes[whichI]))
      membership[membership>whichJ] = membership[membership>whichJ] - 1
      centers = centers[,-whichJ]
      clusterSizes = clusterSizes[-whichJ]
      small = (clusterSizes < preferredSize)
    } else done = TRUE
  }

  if (checkData)
  {
    membershipAll = rep(NA, nAllGenes)
    membershipAll[gsg$goodGenes] = membership
  } else
    membershipAll = membership

  if (seedSaved) .Random.seed <<- savedSeed

  if (verbose > 2)
  {
    printFlush("Sorted sizes of final clusters:")
    print(sort(as.numeric(table(membership))))
  }

  return( list(clusters = membershipAll, centers = centers) )
}

#======================================================================================================
#
# Consensus preliminary partitioning
#
#======================================================================================================



#' Consensus projective K-means (pre-)clustering of expression data
#'
#' Implementation of a consensus variant of K-means clustering for expression
#' data across multiple data sets.
#'
#'
#' The principal aim of this function within WGCNA is to pre-cluster a large
#' number of genes into smaller blocks that can be handled using standard WGCNA
#' techniques.
#'
#' This function implements a variant of K-means clustering that is suitable
#' for co-expression analysis. Cluster centers are defined by the first
#' principal component, and distances by correlation. Consensus distance across
#' several sets is defined as the maximum of the corresponding distances in
#' individual sets; however, if \code{useMean} is set, the mean distance will
#' be used instead of the maximum.  The distance between a gene and a center of
#' a cluster is multiplied by a factor of \eqn{max(clusterSize/preferredSize,
#' 1)^{sizePenaltyPower}}{\code{max(clusterSize/preferredSize,
#' 1)^sizePenaltyPower}}, thus penalizing clusters whose size exceeds
#' \code{preferredSize}. The function starts with randomly generated cluster
#' assignment (hence the need to set the random seed for repeatability) and
#' executes interations of calculating new centers and reassigning genes to
#' nearest (in the consensus sense) center until the clustering becomes stable.
#' Before returning, nearby clusters are iteratively combined if their combined
#' size is below \code{preferredSize}.
#'
#' Consensus distance defined as maximum of distances in all sets is consistent
#' with the approach taken in \code{\link{blockwiseConsensusModules}}, but the
#' procedure may not converge. Hence it is advisable to use the mean as
#' consensus in cases where there are multiple data sets (4 or more, say)
#' and/or if the input data sets are very different.
#'
#' The standard principal component calculation via the function \code{svd}
#' fails from time to time (likely a convergence problem of the underlying
#' lapack functions). Such errors are trapped and the principal component is
#' approximated by a weighted average of expression profiles in the cluster. If
#' \code{verbose} is set above 2, an informational message is printed whenever
#' this approximation is used.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param preferredSize preferred maximum size of clusters.
#' @param nCenters number of initial clusters. Empirical evidence suggests that
#' more centers will give a better preclustering; the default is
#' \code{as.integer(min(nGenes/20, 100*nGenes/preferredSize))} and is an
#' attempt to arrive at a reasonable number given the resources available.
#' @param sizePenaltyPower parameter specifying how severe is the penalty for
#' clusters that exceed \code{preferredSize}.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param randomSeed integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit.
#' @param checkData logical: should data be checked for genes with zero
#' variance and genes and samples with excessive numbers of missing samples?
#' Bad samples are ignored; returned cluster assignment for bad genes will be
#' \code{NA}.
#' @param imputeMissing logical: should missing values in \code{datExpr} be
#' imputed before the calculations start? If the missing data are not imputed,
#' they will be replaced by 0 which can be problematic.
#' @param useMean logical: should mean distance across sets be used instead of
#' maximum? See details.
#' @param maxIterations maximum iterations to be attempted.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{clusters}{ a numerical
#' vector with one component per input gene, giving the cluster number in which
#' the gene is assigned. }
#'
#' \item{centers}{ a vector of lists, one list per set. Each list contains a
#' component \code{data} that contains a matrix whose columns are the cluster
#' centers in the corresponding set. }
#'
#' \item{unmergedClusters}{ a numerical vector with one component per input
#' gene, giving the cluster number in which the gene was assigned before the
#' final merging step. }
#'
#' \item{unmergedCenters}{ a vector of lists, one list per set. Each list
#' contains a component \code{data} that contains a matrix whose columns are
#' the cluster centers before merging in the corresponding set. }
#' @author Peter Langfelder
#' @seealso \code{\link{projectiveKMeans}}
#' @keywords cluster
consensusProjectiveKMeans = function (
      multiExpr,
      preferredSize = 5000,
      nCenters = NULL,
      sizePenaltyPower = 4,
      networkType = "unsigned",
      randomSeed = 54321,
      checkData = TRUE,
      imputeMissing = TRUE,
      useMean = (length(multiExpr) > 3),
      maxIterations = 1000,
      verbose = 0, indent = 0 ) {

  centerGeneDist = function(centers, oldDst = NULL, changed = c(1:nCenters))
  {
    if (is.null(oldDst))
    {
      oldDst = array(0, c(nCenters, nGenes))
      changed = c(1:nCenters)
    }
    dstAll = oldDst
    nChanged = length(changed)
    if (nChanged!=0)
    {
      dstX = array(0, c(nSets, nChanged, nGenes))
      for (set in 1:nSets)
        if (intNetworkType==1)
        {
           dstX[set, , ] = -1+abs(cor(centers[[set]]$data[, changed], multiExpr[[set]]$data))
        } else {
           dstX[set, , ] = -1+cor(centers[[set]]$data[, changed], multiExpr[[set]]$data)
        }
      dst = array(0, c(nChanged, nGenes))
      if (useMean)
      {
        minRes = .C("mean", as.double(dstX), as.integer(nSets), as.integer(nGenes * nChanged),
                    as.double(dst), PACKAGE = "WGCNA")
      } else {
        which = array(0, c(nChanged, nGenes))
        minRes = .C("minWhichMin", as.double(dstX), as.integer(nSets), as.integer(nGenes * nChanged),
                    as.double(dst), as.double(which), PACKAGE = "WGCNA")
      }
      dstAll[changed, ] = -minRes[[4]]
    }
    dstAll
  }

  memberCenterDist = function(dst, membership, sizes = NULL, changed = c(1:nCenters), oldMCDist = NULL)
  {
    if (is.null(oldMCDist))
    {
       changed = c(1:nCenters)
       oldMCDist = rep(0, nGenes)
    }
    centerDist = oldMCDist
    if (!is.null(changed))
    {
      if (is.null(sizes)) sizes = table(membership)
      if (length(sizes)!=nCenters)
      {
        sizes2 = rep(0, nCenters)
        sizes2[as.numeric(names(sizes))] = sizes
        sizes = sizes2
      }
      if (is.finite(sizePenaltyPower))
      {
        sizeCorrections = (sizes/preferredSize)^sizePenaltyPower
        sizeCorrections[sizeCorrections < 1] = 1
      } else {
        sizeCorrections = rep(1, length(sizes))
        sizeCorrections[sizes > preferredSize] = Inf
      }
      for (cen in changed) if (sizes[cen]!=0)
      {
        if (is.finite(sizeCorrections[cen]))
        {
          centerDist[membership==cen] = dst[cen, membership==cen] * sizeCorrections[cen]
        } else
          centerDist[membership==cen] = 10 + dst[cen, membership==cen]
      }
    }
    centerDist
 }

  spaces = indentSpaces(indent)

  if (verbose > 0)
    printFlush(paste(spaces, "Consensus projective K-means:"))

  allSize = checkSets(multiExpr)
  nSamples = allSize$nSamples
  nGenes = allSize$nGenes
  nSets = allSize$nSets

  if (preferredSize >= floor(sqrt(2^31)) )
    stop("'preferredSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (exists(".Random.seed"))
  {
     seedSaved = TRUE
     savedSeed = .Random.seed
  } else seedSaved = FALSE
  set.seed(randomSeed)

  if (is.null(nCenters)) nCenters = as.integer(min(nGenes/20, 100 * nGenes/preferredSize))

  if (verbose > 1)
    printFlush(paste(spaces, "..using", nCenters, "centers."))

  intNetworkType = charmatch(networkType, .networkTypes)
  if (is.na(intNetworkType))
    stop(paste("Unrecognized networkType argument.",
         "Recognized values are (unique abbreviations of)", paste(.networkTypes, collapse = ", ")))

  for (set in 1:nSets)
     multiExpr[[set]]$data = scale(as.matrix(multiExpr[[set]]$data))

  # Check data for genes and samples that have too many missing values
  if (checkData)
  {
    if (verbose > 0)
      printFlush(paste(spaces, "..checking data for excessive number of missing values.."))
    gsg = goodSamplesGenesMS(multiExpr, verbose = verbose - 1, indent = indent + 1)
    for (set in 1:nSets)
    {
      if (!gsg$allOK)
        multiExpr[[set]]$data = scale(multiExpr[[set]]$data[gsg$goodSamples[[set]], gsg$goodGenes])
    }
  }

  anyNA = apply(multiExpr, function(x) any(is.na(x)), mdaSimplify = TRUE)
  if (imputeMissing && any(anyNA))
  {
    if (verbose > 0)
      printFlush(paste(spaces, "..imputing missing data.."))
    multiExpr[anyNA] = apply(multiExpr[anyNA], function(x) t(impute.knn(t(x))$data), mdaVerbose = verbose>1)
  } else if (any(anyNA))
  {
    printFlush(paste(spaces, "Found missing data. These will be replaced by zeros;\n",
                     spaces, "  for a better replacement, use 'imputeMissing=TRUE'."))
    for (set in 1:nSets)
      multiExpr[[set]]$data[is.na(multiExpr[[set]]$data)] = 0
  }

  setSize = checkSets(multiExpr)
  nGenes = setSize$nGenes
  nSamples = setSize$nSamples

  centers = vector(mode="list", length = nSets)
  for (set in 1:nSets)
    centers[[set]] = list(data = matrix(0, nSamples[set], nCenters))

  randGeneIndex = sample(nGenes, size = nGenes)
  temp = rep(c(1:nCenters), times = ceiling(nGenes/nCenters))
  membership = temp[randGeneIndex]

  if (verbose > 0)
    printFlush(paste(spaces, "..consensus k-means clustering.."))

  changed = c(1:nCenters); dst = NULL
  iteration = 0
  centerDist = NULL
  while (!is.null(changed) && (iteration < maxIterations)) {
    iteration = iteration + 1
    if (verbose > 1) printFlush(paste(spaces, " ..iteration", iteration))
    clusterSizes = table(membership)
    for (set in 1:nSets) for (cen in changed)
        centers[[set]]$data[, cen] = .alignedFirstPC(multiExpr[[set]]$data[, membership==cen],
                                                     verbose = verbose-2, indent = indent+2)
    dst = centerGeneDist(centers, dst, changed)
    centerDist = memberCenterDist(dst, membership, clusterSizes, changed, centerDist)

    nearestDist = rep(0, nGenes)
    nearest = rep(0, nGenes)
    minRes = .C("minWhichMin", as.double(dst), as.integer(nCenters), as.integer(nGenes),
                as.double(nearestDist), as.double(nearest), PACKAGE = "WGCNA")
    nearestDist = minRes[[4]]
    nearest = minRes[[5]]+1
    changed = NULL
    rm(minRes); collectGarbage()
    if (sum(centerDist>nearestDist)>0) {
      proposedMemb = nearest
      accepted = FALSE
      while (!accepted && (sum(proposedMemb!=membership)>0)) {
        if (verbose > 2)
          cat(paste(spaces, "   ..proposing to move", sum(proposedMemb!=membership), "genes"))
        moved = c(1:nGenes)[proposedMemb!=membership]
        newCentDist = memberCenterDist(dst, proposedMemb)
        gotWorse = newCentDist[moved] > centerDist[moved]
        if (sum(gotWorse)==0){
          accepted = TRUE
          if (verbose > 2) printFlush(paste("..move accepted."))
        } else {
          if (verbose > 2) printFlush(paste("..some genes got worse. Trying again."))
          ord = order(centerDist[moved[gotWorse]] - newCentDist[moved[gotWorse]])
          n = ceiling(length(ord)*3/5)
          proposedMemb[moved[gotWorse][ord[c(1:n)]]] = membership[moved[gotWorse][ord[c(1:n)]]]
        }
      }
      if (accepted) {
        propSizes = table(proposedMemb)
        keep = as.numeric(names(propSizes))
        for (set in 1:nSets)
          centers[[set]]$data = centers[[set]]$data[, keep]
        dst = dst[keep, ]
        changedAll = union(membership[moved], proposedMemb[moved])
        changedKeep = changedAll[is.finite(match(changedAll, keep))]
        changed = rank(changedKeep);	# This is a way to make say 1,3,4 into 1,2,3
        membership = as.numeric(as.factor(proposedMemb))
        if ( (verbose > 1) && (sum(keep) < nCenters))
          printFlush(paste(spaces, "  ..dropping", nCenters - sum(keep),
                           "centers because ther clusters are empty."))
        nCenters = length(keep)
      } else {
        changed = NULL
        if (verbose > 2) printFlush(paste("Could not find a suitable move to improve the clustering."))
      }
    } else {
      changed = NULL
      if (verbose > 2)
        printFlush(paste("Clustering is stable: all genes are closest to their assigned center."))
    }
  }

  consCenterDist = function(centers, select) {
    if (is.logical(select)) {
      nC = sum(select)
    } else {
      nC = length(select)
    }
    distX = array(0, dim=c(nSets, nC, nC))
    for (set in 1:nSets)
    {
      if (intNetworkType==1)
      {
         distX[set, , ] = -1+abs(cor(centers[[set]]$data[, select]))
      } else {
         distX[set, , ] = -1+cor(centers[[set]]$data[, select])
      }
    }
    dst = matrix(0, nC, nC)
    which = matrix(0, nC, nC)
    minRes = .C("minWhichMin", as.double(distX), as.integer(nSets), as.integer(nC*nC),
                as.double(dst), as.double(which), PACKAGE = "WGCNA")
    dst[,] = -minRes[[4]]
    dst
  }

  unmergedMembership = membership
  unmergedCenters = centers
  # merge nearby clusters if their sizes allow merging
  if (verbose > 0) printFlush(paste(spaces, "..merging smaller clusters..."))
  clusterSizes = as.vector(table(membership))
  small = (clusterSizes < preferredSize)
  done = FALSE
  while (!done & (sum(small)>1)) {
    smallIndex = c(1:nCenters)[small]
    nSmall = sum(small)
    clustDist = consCenterDist(centers, smallIndex)

    diag(clustDist) = 10
    distOrder = order(as.vector(clustDist))[seq(from=2, to = nSmall * (nSmall-1), by=2)]
    i = 1; canMerge = FALSE
    while (i <= length(distOrder) && (!canMerge)) {
      whichJ = smallIndex[as.integer( (distOrder[i]-1)/nSmall + 1)]
      whichI = smallIndex[distOrder[i] - (whichJ-1) * nSmall]
      canMerge = sum(clusterSizes[c(whichI, whichJ)]) < preferredSize
      i = i + 1
    }
    if (canMerge) {
      membership[membership==whichJ] = whichI
      clusterSizes[whichI] = sum(clusterSizes[c(whichI, whichJ)])
      #if (verbose > 1)
      #  printFlush(paste(spaces, "  ..merged clusters", whichI, "and", whichJ,
      #             "whose combined size is", clusterSizes[whichI]))
      for (set in 1:nSets)
        centers[[set]]$data[, whichI] = .alignedFirstPC(multiExpr[[set]]$data[, membership==whichI],
                                                        verbose = verbose-2, indent = indent+2)
      for (set in 1:nSets)
        centers[[set]]$data = centers[[set]]$data[,-whichJ]
      membership[membership>whichJ] = membership[membership>whichJ] - 1
      nCenters = nCenters -1
      clusterSizes = as.vector(table(membership))
      small = (clusterSizes < preferredSize)
      if (verbose > 3)
        printFlush(paste(spaces, "  ..merged clusters", whichI, "and", whichJ,
                   "whose combined size is", clusterSizes[whichI]))
    } else done = TRUE
  }

  if (checkData){
    membershipAll = rep(NA, allSize$nGenes)
    membershipAll[gsg$goodGenes] = membership
  } else
    membershipAll = membership

  if (seedSaved) .Random.seed <<- savedSeed

  return(list(clusters = membershipAll, centers = centers, unmergedClusters = unmergedMembership,
               unmergedCenters = unmergedCenters))
}
