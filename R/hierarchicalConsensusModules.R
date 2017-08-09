# Hierarchical consensus modules

#' Hierarchical consensus network construction and module identification
#'
#' Hierarchical consensus network construction and module identification across
#' multiple data sets.
#'
#' This function calculates a consensus network with a flexible, possibly
#' hierarchical consensus specification, identifies (consensus) modules in the
#' network, and calculates their eigengenes. "Blockwise" calculation is
#' available for large data sets for which a full network (TOM or adjacency
#' matrix) would not fit into avilable RAM.
#'
#' The input can be either several numerical data sets (expression etc) in the
#' argument \code{multiExpr} together with all necessary network construction
#' options, or a pre-calculated network, typically the result of a call to
#' \code{\link{hierarchicalConsensusTOM}}.
#'
#' Steps in the network construction include the following: (1) optional
#' filtering of variables (genes) and observations (samples) that contain too
#' many missing values or have zero variance; (2) optional pre-clustering to
#' split data into blocks of manageable size; (3) calculation of adjacencies
#' and optionally of TOMs in each individual data set; (4) calculation of
#' consensus network from the individual networks; (5) hierarchical clustering
#' and module identification; (6) trimming of modules by removing genes with
#' low correlation with the eigengene of the module; and (7) merging of modules
#' whose eigengenes are strongly correlated.
#'
#' Steps 1-4 (up to and including the calculation of consensus network from the
#' individual networks) are handled by the function
#' \code{\link{hierarchicalConsensusTOM}}.
#'
#' Variables (genes) are clustered using average-linkage hierarchical
#' clustering and modules are identified in the resulting dendrogram by the
#' Dynamic Hybrid tree cut. Found modules are trimmed of genes whose consensus
#' module membership kME (that is, correlation with module eigengene) is less
#' than \code{minKMEtoStay}. Modules in which fewer than \code{minCoreKMESize}
#' genes have consensus KME higher than \code{minCoreKME} are disbanded, i.e.,
#' their constituent genes are pronounced unassigned.
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
#' @param multiExpr Expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param multiExpr.imputed If \code{multiExpr} contain missing data, this
#' argument can be used to supply the expression data with missing data
#' imputed. If not given, the \code{\link[impute]{impute.knn}} function will be
#' used to impute the missing data.
#' @param checkMissingData Logical: should data be checked for excessive
#' numbers of missing entries in genes and samples, and for genes with zero
#' variance? See details.
#' @param blocks Optional specification of blocks in which hierarchical
#' clustering and module detection should be performed. If given, must be a
#' numeric vector with one entry per gene of \code{multiExpr} giving the number
#' of the block to which the corresponding gene belongs.
#' @param maxBlockSize Integer giving maximum block size for module detection.
#' Ignored if \code{blocks} above is non-NULL. Otherwise, if the number of
#' genes in \code{datExpr} exceeds \code{maxBlockSize}, genes will be
#' pre-clustered into blocks whose size should not exceed \code{maxBlockSize}.
#' @param blockSizePenaltyPower Number specifying how strongly blocks should be
#' penalized for exceeding the maximum size. Set to a lrge number or \code{Inf}
#' if not exceeding maximum block size is very important.
#' @param nPreclusteringCenters Number of centers to be used in the
#' preclustering. Defaults to smaller of \code{nGenes/20} and
#' \code{100*nGenes/maxBlockSize}, where \code{nGenes} is the nunber of genes
#' (variables) in \code{multiExpr}.
#' @param randomSeed Integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit. If \code{NULL} is given, the function will not save and
#' restore the seed.
#' @param networkOptions A single list of class \code{\link{NetworkOptions}}
#' giving options for network calculation for all of the networks, or a
#' \code{\link{multiSet}} structure containing one such list for each input
#' data set.
#' @param saveIndividualTOMs Logical: should individual TOMs be saved to disk
#' (\code{TRUE}) or retuned directly in the return value (\code{FALSE})?
#' @param individualTOMFileNames Character string giving the file names to save
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param keepIndividualTOMs Logical: should individual TOMs be retained after
#' the calculation is finished?
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param saveConsensusTOM Logical: should the consensus TOM be saved to disk?
#' @param consensusTOMFilePattern Character string giving the file names to
#' save consensus TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param keepConsensusTOM Logical: should consensus TOM be retained after the
#' calculation ends? Depending on \code{saveConsensusTOM}, the retained TOM is
#' either saved to disk or returned within the return value.
#' @param useDiskCache Logical: should disk cache be used for consensus
#' calculations? The disk cache can be used to store chunks of calibrated data
#' that are small enough to fit one chunk from each set into memory (blocks may
#' be small enough to fit one block of one set into memory, but not small
#' enough to fit one block from all sets in a consensus calculation into memory
#' at the same time). Using disk cache is slower but lessens the memory
#' footprint of the calculation. As a general guide, if individual data are
#' split into blocks, we recommend setting this argument to \code{TRUE}. If
#' this argument is \code{NULL}, the function will decide whether to use disk
#' cache based on the number of sets and block sizes.
#' @param chunkSize Integer giving the chunk size. If left \code{NULL}, a
#' suitable size will be chosen automatically.
#' @param cacheDir Directory in which to save cache files. The files are
#' deleted on normal exit but persist if the function terminates abnormally.
#' @param cacheBase Base for the file names of cache files.
#' @param consensusTOMInfo If the consensus TOM has been pre-calculated using
#' function \code{\link{hierarchicalConsensusTOM}}, this argument can be used
#' to supply it. If given, the consensus TOM calculation options above are
#' ignored.
#' @param deepSplit Numeric value between 0 and 4. Provides a simplified
#' control over how sensitive module detection should be to module splitting,
#' with 0 least and 4 most sensitive. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param detectCutHeight Dendrogram cut height for module detection. See
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for more details.
#' @param minModuleSize Minimum module size for module detection. See
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
#' @param stabilityCriterion One of \code{c("Individual fraction", "Common
#' fraction")}, indicating which method for assessing stability similarity of
#' two branches should be used. We recommend \code{"Individual fraction"} which
#' appears to perform better; the \code{"Common fraction"} method is provided
#' for backward compatibility since it was the (only) method available prior to
#' WGCNA version 1.60.
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
#' @param calibrateMergingSimilarities Logical: should module eigengene
#' similarities be calibrataed before calculating the consensus? Although
#' calibration is in principle desirable, the calibration methods currently
#' available assume large data and do not work very well on eigengene
#' similarities.
#' @param mergeCutHeight Dendrogram cut height for module merging.
#' @param collectGarbage Logical: should garbage be collected after some of the
#' memory-intensive steps?
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @param \dots Other arguments. Currently ignored.
#' @return List with the following components: \item{labels}{A numeric vector
#' with one component per variable (gene), giving the module label of each
#' variable (gene).  Label 0 is reserved for unassigned variables; module
#' labels are sequential and smaller numbers are used for larger modules.}
#'
#' \item{unmergedLabels}{A numeric vector with one component per variable
#' (gene), giving the unmerged module label of each variable (gene), i.e.,
#' module labels before the call to module merging.}
#'
#' \item{colors}{A character vector with one component per variable (gene),
#' giving the module colors. The labels are mapped to colors using
#' \code{\link{labels2colors}}.}
#'
#' \item{unmergedColors}{A character vector with one component per variable
#' (gene), giving the unmerged module colors.}
#'
#' \item{multiMEs}{Module eigengenes corresponding to the modules returned in
#' \code{colors}, in multi-set format. A vector of lists, one per set,
#' containing eigengenes, proportion of variance explained and other
#' information. See \code{\link{multiSetMEs}} for a detailed description.}
#'
#' \item{dendrograms}{A list with one component for each block of genes. Each
#' component is the hierarchical clustering dendrogram obtained by clustering
#' the consensus gene dissimilarity in the corresponding block. }
#'
#' \item{consensusTOMInfo}{A list detailing various aspects of the consensus
#' TOM. See \code{\link{hierarchicalConsensusTOM}} for details.}
#'
#' \item{blockInfo}{A list with information about blocks as well as the
#' vriables and observations (genes and samples) retained after filtering out
#' those with zero variance and too many missing values.}
#'
#' \item{moduleIdentificationArguments}{A list with the module identification
#' arguments supplied to this function. Contains \code{deepSplit},
#' \code{detectCutHeight}, \code{minModuleSize}, \code{maxCoreScatter},
#' \code{minGap}, \code{maxAbsCoreScatter}, \code{minAbsGap},
#' \code{minSplitHeight}, \code{useBranchEigennodeDissim},
#' \code{minBranchEigennodeDissim}, \code{minStabilityDissim}, \code{pamStage},
#' \code{pamRespectsDendro}, \code{minCoreKME}, \code{minCoreKMESize},
#' \code{minKMEtoStay}, \code{calibrateMergingSimilarities}, and
#' \code{mergeCutHeight}.}
#' @note If the input datasets have large numbers of genes, consider carefully
#' the \code{maxBlockSize} as it significantly affects the memory footprint
#' (and whether the function will fail with a memory allocation error). From a
#' theoretical point of view it is advantageous to use blocks as large as
#' possible; on the other hand, using smaller blocks is substantially faster
#' and often the only way to work with large numbers of genes. As a rough
#' guide, when 4GB of memory are available, blocks should be no larger than
#' 8,000 genes; with 8GB one can handle some 13,000 genes; with 16GB around
#' 20,000; and with 32GB around 30,000. Depending on the operating system and
#' its setup, these numbers may vary substantially.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{hierarchicalConsensusTOM}} for calculation of hierarchical
#' consensus networks (adjacency and TOM), and a more detailed description of
#' the calculation;
#'
#' \code{\link[fastcluster]{hclust}} and
#' \code{\link[dynamicTreeCut]{cutreeHybrid}} for hierarchical clustering and
#' the Dynamic Tree Cut branch cutting method;
#'
#' \code{\link{mergeCloseModules}} for module merging;
#'
#' \code{\link{blockwiseModules}} for an analogous analysis on a single data
#' set.
#' @references
#'
#' Non-hierarchical consensus networks are described in Langfelder P, Horvath S
#' (2007), Eigengene networks for studying the relationships between
#' co-expression modules. BMC Systems Biology 2007, 1:54.
#'
#' More in-depth discussion of selected topics can be found at
#' http://www.peterlangfelder.com/ , and an FAQ at
#' https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/faq.html
#' .
#' @keywords misc
#' @export
hierarchicalConsensusModules = function(multiExpr,

                                        # Optional: multiExpr wigth imputed missing data
                                        multiExpr.imputed = NULL,

                                        # Data checking options

                                        checkMissingData = TRUE,

                                        # Blocking options

                                        blocks = NULL,
                                        maxBlockSize = 5000,
                                        blockSizePenaltyPower = 5,
                                        nPreclusteringCenters = NULL,
                                        randomSeed = 12345,

                                        # ...or information needed to construct individual networks

                                        # Network construction options. This can be a single object of class NetworkOptions, or a multiSet
                                        # structure of NetworkOptions objects, one per element of multiExpr.

                                        networkOptions,

                                        # Save individual TOMs?

                                        saveIndividualTOMs = TRUE,
                                        individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",
                                        keepIndividualTOMs = FALSE,

                                        # Consensus calculation options

                                        consensusTree = NULL,  # if not given, the one in consensusTOMInfo will be used.

                                        # Return options
                                        saveConsensusTOM = TRUE,
                                        consensusTOMFilePattern = "consensusTOM-%a-Block%b.RData",

                                        # Keep the consensus? Note: I will not have an option to keep intermediate results here.
                                        keepConsensusTOM = saveConsensusTOM,

                                        # Internal handling of TOMs

                                        useDiskCache = NULL, chunkSize = NULL,
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
                                        stabilityCriterion = c("Individual fraction", "Common fraction"),
                                        minStabilityDissim = NULL,

                                        pamStage = TRUE,  pamRespectsDendro = TRUE,

                                        # Gene joining and removal from a module, and module "significance" criteria
                                        # reassignThresholdPS = 1e-4, ## For now do not do gene reassignment - have to think more about how
                                        # to do it.

                                        minCoreKME = 0.5, minCoreKMESize = minModuleSize/3,
                                        minKMEtoStay = 0.2,

                                        # Module eigengene calculation options

                                        impute = TRUE,
                                        trapErrors = FALSE,

                                        # Module merging options

                                        calibrateMergingSimilarities = FALSE,
                                        mergeCutHeight = 0.15,

                                        # General options
                                        collectGarbage = TRUE,
                                        verbose = 2, indent = 0,
                                        ...)
{
  spaces = indentSpaces(indent)

  dataSize = checkSets(multiExpr)
  nSets = dataSize$nSets
  nGenes = dataSize$nGenes
  # nSamples = dataSize$nSamples

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
      savedSeed = .Random.seed
      on.exit(.Random.seed <<-savedSeed, add = FALSE)
    }
    set.seed(randomSeed)
  }

  if (checkMinModuleSize & (minModuleSize > nGenes/5))
  {
    minModuleSize = nGenes/5
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
  if (is.null(multiExpr.imputed)) {
    multiExpr.imputed = multiSet.mapply(function(x, doImpute)
    { if (doImpute) t(impute.knn(t(x))$data) else x },
    multiExpr.scaled, hasMissing)
  } else {
    size.imp = checkSets(multiExpr.imputed)
    if (!isTRUE(all.equal(size.imp, dataSize)))
      stop("If given, 'multiExpr.imputed' must have the same dimensions in each set as 'multiExpr'.")
  }
  branchSplitFnc = NULL
  minBranchDissimilarities = numeric(0)
  externalSplitFncNeedsDistance = logical(0)
  if (useBranchEigennodeDissim)
  {
    branchSplitFnc = "hierarchicalBranchEigengeneDissim"
    minBranchDissimilarities = minBranchEigennodeDissim
    externalSplitFncNeedsDistance = FALSE
  }

  if (!is.null(stabilityLabels))
  {
    stabilityCriterion = match.arg(stabilityCriterion)
    branchSplitFnc = c(branchSplitFnc,
                       if (stabilityCriterion=="Individual fraction")
                         "branchSplitFromStabilityLabels.individualFraction" else "branchSplitFromStabilityLabels")
    minBranchDissimilarities = c(minBranchDissimilarities, minStabilityDissim)
    externalSplitFncNeedsDistance = c(externalSplitFncNeedsDistance, FALSE)
  }

  otherArgs = list(...)

  getDetails = FALSE
  if ("getDetails" %in% names(otherArgs)) getDetails = otherArgs$getDetails


  if (is.null(consensusTOMInfo))
  {
    if (is.null(consensusTree))
      stop("Either 'consensusTOMInfo' or 'consensusTree' must be given.")

    consensusTOMInfo = hierarchicalConsensusTOM(
      multiExpr = multiExpr,
      checkMissingData = checkMissingData,
      blocks = blocks,
      maxBlockSize = maxBlockSize,
      blockSizePenaltyPower = blockSizePenaltyPower,
      nPreclusteringCenters = nPreclusteringCenters,
      randomSeed = NULL,

      networkOptions = networkOptions,

      keepIndividualTOMs = keepIndividualTOMs,
      individualTOMFileNames = individualTOMFileNames,

      consensusTree = consensusTree,

      saveCalibratedIndividualTOMs = FALSE,
      getCalibrationSamples = FALSE,

      # Return options
      saveConsensusTOM = saveConsensusTOM,
      consensusTOMFilePattern = consensusTOMFilePattern,

      keepIntermediateResults = FALSE,

      # Internal handling of TOMs
      useDiskCache = useDiskCache,
      chunkSize = chunkSize,
      cacheBase = cacheBase,
      cacheDir = cacheDir,
      collectGarbage = collectGarbage,
      verbose = verbose, indent = indent)

    removeConsensusTOMOnExit = !keepConsensusTOM
  } else {
    # Basic checks on consensusTOMInfo
    .checkComponents(consensusTOMInfo, c("individualTOMInfo", "consensusData", "consensusTree"))

    if (length(consensusTOMInfo$individualTOMInfo$blockInfo$blocks)!=nGenes)
      stop("Inconsistent number of genes in 'consensusTOMInfo$individualTOMInfo$blockInfo$blocks'.")

    if (!is.null(consensusTree) && !isTRUE(all.equal(consensusTree, consensusTOMInfo$consensusTree)))
      warning(immediate. = TRUE,
              "hierarchicalConsensusModules: given 'consensusTree' is different\n",
              "from the 'consensusTree' component in 'consensusTOMInfo'. \n",
              "This is normally undesirable and may\n",
              "indicate a mistake in the function call.")

    if (is.null(consensusTree)) consensusTree = consensusTOMInfo$consensusTree

    removeConsensusTOMOnExit = FALSE
    networkOptions = consensusTOMInfo$individualTOMInfo$networkOptions
  }

  allLabels = rep(0, nGenes)
  allLabelIndex = NULL

  # Restrict data to goodSamples and goodGenes

  gsg = consensusTOMInfo$individualTOMInfo$blockInfo$goodSamplesAndGenes

  if (!gsg$allOK)
    multiExpr = subset(multiExpr, gsg$goodSamples, gsg$goodGenes)

  nGGenes = sum(gsg$goodGenes)
  nGSamples = sapply(gsg$goodSamples, sum)

  blocks = consensusTOMInfo$individualTOMInfo$blockInfo$blocks
  gBlocks = consensusTOMInfo$individualTOMInfo$blockInfo$gBlocks

  blockLevels = sort(unique(gBlocks))
  blockSizes = table(gBlocks)
  nBlocks = length(blockLevels)

  # reassignThreshold = reassignThresholdPS^nSets

  consMEs = vector(mode = "list", length = nSets)
  dendros = list()

  cutreeLabels = list()
  maxUsedLabel = 0
  # Here's where the analysis starts

  for (blockNo in 1:nBlocks)
  {
    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    nBlockGenes = length(block)

    selExpr = subset(multiExpr, , block)
    errorOccurred = FALSE
    consTomDS = BD.getData(consensusTOMInfo$consensusData, blockNo)
    # Temporary "cast" so fastcluster::hclust doesn't complain about non-integer size.
    # attr(consTomDS, "Size") = as.integer(attr(consTomDS, "Size")); ## This should not be needed now.

    consTomDS = 1-consTomDS

    if (collectGarbage) gc()

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
      externalSplitOptions[[e.index]] = list(multiExpr = subset(multiExpr.imputed,, block),
                                             networkOptions = networkOptions,
                                             consensusTree = consensusTree)
      e.index = e.index +1
    }
    if (!is.null(stabilityLabels))
    {
      externalSplitOptions[[e.index]] = list(stabilityLabels = stabilityLabels)
      e.index = e.index + 1
    }

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
    if (getDetails)
    {
      cutreeLabels[[blockNo]] = blockLabels
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
    blockLabelIndex = sort(unique(blockLabels[blockAssigned]))
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

    if (collectGarbage) gc()

    # Check modules: make sure that of the genes present in the module, at least a minimum number
    # have a correlation with the eigengene higher than a given cutoff, and that all member genes have
    # the required minimum consensus KME

    if (verbose>2)
      printFlush(paste(spaces, "....checking consensus modules for statistical meaningfulness.."))

    KME = multiSet.mapply(function(expr, me, netOpt)
    {
      # printFlush("=============================================================")
      # print(netOpt$corOptions)
      kme = do.call(netOpt$corFnc, c(list(x = expr, y = me), netOpt$corOptions))
      if (!grepl("signed", netOpt$networkType)) kme = abs(kme)
      kme
    }, selExpr, blockConsMEs, networkOptions, returnList = TRUE)
    consKME = simpleHierarchicalConsensusCalculation(KME, consensusTree)

    for (mod in 1:ncol(blockConsMEs[[1]]$data)) {
      modGenes = (blockLabels==blockLabelIndex[mod])
      consKME1 = consKME[modGenes, mod]
      if (sum(consKME1>minCoreKME) < minCoreKMESize) {
        blockLabels[modGenes] = 0
        deleteModules = union(deleteModules, mod)
        if (verbose>3)
          printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod],
                           ": of ", sum(modGenes),
                           " total genes in the module only ",  sum(consKME1>minCoreKME),
                           " have the requisite high correlation with the eigengene in all sets.", sep=""))
      } else if (sum(consKME1<minKMEtoStay)>0) {
        if (verbose > 3)
          printFlush(paste(spaces, "    ..removing", sum(consKME1<minKMEtoStay),
                           "genes from module", blockLabelIndex[mod], "because their KME is too low."))
        blockLabels[modGenes][consKME1 < minKMEtoStay] = 0
        if (sum(blockLabels[modGenes]>0) < minModuleSize) {
          deleteModules = union(deleteModules, mod)
          blockLabels[modGenes] = 0
          if (verbose>3)
            printFlush(paste(spaces, "    ..deleting module ",blockLabelIndex[mod],
                             ": not enough genes in the module after removal of low KME genes.", sep=""))
        } else {
          changedModules = union(changedModules, blockLabelIndex[mod])
        }
      }
    }

    # Remove marked modules

    if (!is.null(deleteModules)) {
      for (set in 1:nSets) blockConsMEs[[set]]$data = blockConsMEs[[set]]$data[, -deleteModules]
      modGenes = is.finite(match(blockLabels, blockLabelIndex[deleteModules]))
      blockLabels[modGenes] = 0
      modAllGenes = is.finite(match(allLabels, blockLabelIndex[deleteModules]))
      allLabels[modAllGenes] = 0
      blockLabelIndex = blockLabelIndex[-deleteModules]
    }

    # Check whether there's anything left
    if (sum(blockLabels>0)==0) {
      if (verbose>1) {
        printFlush(paste(spaces, "  ..No significant modules detected in block", blockNo))
        printFlush(paste(spaces, "  ..continuing with next block."))
      }
      next
    }

    # Update module eigengenes

    for (set in 1:nSets)
      if (is.null(dim(blockConsMEs[[set]]$data)))
        dim(blockConsMEs[[set]]$data) = c(length(blockConsMEs[[set]]$data), 1)

    if (is.null(consMEs[[1]])) {
      for (set in 1:nSets) consMEs[[set]] = list(data = blockConsMEs[[set]]$data)
    } else for (set in 1:nSets)
      consMEs[[set]]$data = cbind(consMEs[[set]]$data, blockConsMEs[[set]]$data)

    # Update allLabels

    allLabelIndex = c(allLabelIndex, blockLabelIndex)
    allLabels[gsg$goodGenes][block[blockAssigned]] = blockLabels[blockAssigned]

    collectGarbage()

  }

  if (verbose>1) printFlush(paste(spaces, "..merging consensus modules that are too close.."))

  #print(table(allLabels))
  #print(is.numeric(allLabels))

  mergedLabels = rep(NA, nGenes)

  mergedMods = try(hierarchicalMergeCloseModules(multiExpr, allLabels[gsg$goodGenes],
                                                 networkOptions = networkOptions, consensusTree = consensusTree,
                                                 calibrateMESimilarities = calibrateMergingSimilarities,
                                                 cutHeight = mergeCutHeight,
                                                 relabel = TRUE,
                                                 verbose = verbose-2, indent = indent + 2), silent = TRUE)
  if (class(mergedMods)=='try-error') {
    if (verbose>0) {
      printFlush(paste(spaces, 'blockwiseConsensusModules: mergeCloseModule failed with this message:\n',
                       spaces, '    ', mergedMods, spaces,
                       '---> returning unmerged consensus modules'))
    } else warning(paste('blockwiseConsensusModules: mergeCloseModule failed with this message:\n     ',
                         mergedMods, '---> returning unmerged consensus modules'))
    MEs = try(multiSetMEs(multiExpr, universalColors = allLabels[gsg$goodGenes]
                          # trapErrors = TRUE, returnValidOnly = TRUE
    ), silent = TRUE)
    if (class(MEs)=='try-error') {
      warning(paste('blockwiseConsensusModules: ME calculation failed with this message:\n     ',
                    MEs, '---> returning empty module eigengenes'))
      allSampleMEs = NULL
    } else {
      mergedLabels[gsg$goodGenes] = allLabels[gsg$goodGenes]
      allSampleMEs = vector(mode = "list", length = nSets)
      for (set in 1:nSets) {
        allSampleMEs[[set]] =
          list(data = as.data.frame(matrix(NA, nrow = nGSamples[set], ncol = ncol(MEs[[set]]$data))))
        allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = MEs[[set]]$data[,]
        names(allSampleMEs[[set]]$data) = names(MEs[[set]]$data)
      }
    }
  } else {
    mergedLabels[gsg$goodGenes] = mergedMods$labels
    allSampleMEs = vector(mode = "list", length = nSets)
    for (set in 1:nSets) {
      allSampleMEs[[set]] =
        list(data = as.data.frame(matrix(NA, nrow = nGSamples[set],
                                         ncol = ncol(mergedMods$newMEs[[1]]$data))))
      allSampleMEs[[set]]$data[gsg$goodSamples[[set]], ] = mergedMods$newMEs[[set]]$data[,]
      names(allSampleMEs[[set]]$data) = names(mergedMods$newMEs[[set]]$data)
    }
  }

  names(allSampleMEs) = names(multiExpr)

  if (removeConsensusTOMOnExit) {
    BD.checkAndDeleteFiles(consensusTOMInfo$consensusData)
    consensusTOMInfo$consensusData = NULL
  }

  list(labels = mergedLabels,
       unmergedLabels = allLabels,
       colors = labels2colors(mergedLabels),
       unmergedColors = labels2colors(allLabels),
       multiMEs = allSampleMEs,
       dendrograms = dendros,
       consensusTOMInfo = consensusTOMInfo,
       blockInfo = consensusTOMInfo$individualTOMInfo$blockInfo,
       moduleIdentificationArguments = list(
         deepSplit = deepSplit,
         detectCutHeight = detectCutHeight,
         minModuleSize = minModuleSize,
         maxCoreScatter = maxCoreScatter,
         minGap = minGap,
         maxAbsCoreScatter = maxAbsCoreScatter,
         minAbsGap = minAbsGap,
         minSplitHeight = minAbsSplitHeight,
         useBranchEigennodeDissim = useBranchEigennodeDissim,
         minBranchEigennodeDissim = minBranchEigennodeDissim,
         minStabilityDissim = minStabilityDissim,
         pamStage = pamStage,
         pamRespectsDendro = pamRespectsDendro,
         minCoreKME = minCoreKME,
         minCoreKMESize = minCoreKMESize,
         minKMEtoStay = minKMEtoStay,
         calibrateMergingSimilarities = calibrateMergingSimilarities,
         mergeCutHeight = mergeCutHeight),
       details = if(getDetails) list(cutreeLabels = cutreeLabels) else NULL
  )
}
