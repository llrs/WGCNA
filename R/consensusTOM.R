# New multilevel specification of consensus: a hierarchical list.

# The consensus calculation needs to be general enough so it can be used for module merging (consensus
# eigengene network) as well as for calculation of consensus kMEs, and possibly for other stuff as well.

# Consensus specification for a single operation:
#   - inputs: a list specifying the input of the consensus. This should be general enough to handle
#   blockwise data but also not specific to adjacencies.
#     The inut should be a multiSet structure, with each
#   - calibrationMethod: currently one of "none", "simple quantile", "full quantile"
#   - consensusQuantile
#   - saveCalibratedIndividualTOMs (= FALSE)
#   - calibratedIndividualTOMFilePattern (= "calibratedIndividualTOM-Set%s-Block%b.RData")

# The consensus calculation itself does not need information about correlation type etc; the consensus
# calculation is very general.

# For network analysis applications of the consensus: will also have to keep information about how the
# individual network adjacencies (TOMs) are to be created or were created - correlation type, network type,
# TOM type, correlation options etc.


# So we will keep 2 different types of information around:
# 1. network construction options. A simple list giving the necessary network construction options.
# 2. ConsensusOptions: This will also be a list giving the options but not holding the actual consensus
# inputs or outputs.

# outputs of network construction and consensus construction should be held in separate lists.


# Output of individualTOMs:
#   - adjacency data, a list of blockwiseData instances
#   - block information
#   - network construction options, separately for each adjacency (network construction could be different)

# For consensusTOM:
#   Inputs:
#      - adjacency data: either from individual TOMs or from other consensus TOMs
#          . note that constructing a complicated consensus manually (i.e., using consensus TOMs
#            as inputs to higher-level consensus) is of limited use
#            since consensus modules would need the full consensus tree anyway.
#      - consensus tree
#   Not really needed: block information
#   outputs:
#      - consensus adjacency data
#      - copy of consensus options
#      - other (diagnostic etc) information

# For consensus modules
#   Inputs:
#      - undelying expression data
#      - optional consensus TOM
#      - block information
#      - network options for each expression data set
#      - consensus tree


.checkPower = function(power) {
    if (any(!is.finite(power)) | (sum(power<1)>0) | (sum(power>50)>0) )
        stop("power must be between 1 and 50.")
}


#==========================================================================================================
#
# Defining a single consensus operation
#
#==========================================================================================================



#' Create a new consensus tree
#'
#' This function creates a new consensus tree, a class for representing
#' "recipes" for hierarchical consensus calculations.
#'
#' Consensus trees specify a "recipe" for the calculation of hierarchical
#' consensus in \code{\link{hierarchicalConsensusCalculation}} and other
#' functions.
#'
#' @aliases newConsensusTree ConsensusTree
#' @param consensusOptions An object of class \code{ConsensusOptions}, usually
#' obtained by calling \code{\link{newConsensusOptions}}.
#' @param inputs A vector (or list) of inputs. Each component can be either a
#' character string giving a names of a data set, or another
#' \code{ConsensusTree} object.
#' @param analysisName Optional specification of a name for this consensus
#' analysis. While this has no effect on the actual consensus calculation, some
#' functions use this character string to make certain file names unique.
#' @return A list with class set to \code{"ConsensusTree"} with these
#' components: \code{consensusOptions}A copy of the input
#' \code{consensusOptions}. \code{inputs}A copy of the input \code{inputs}.
#' \code{analysisName}A copy of the input \code{analysisName}.
#' @author Peter Langfelder
#' @seealso \code{\link{hierarchicalConsensusCalculation}} for hierarchical
#' consensus calculation for which a \code{ConsensusTree} object specifies the
#' recipe
#' @keywords misc
newConsensusTree = function(consensusOptions = newConsensusOptions(),
                            inputs,
                            analysisName = NULL) {
  if (!inherits(consensusOptions, "ConsensusOptions"))
    stop("'consensusOptions' must be of class ConsensusOptions.")
  out = list(consensusOptions = consensusOptions,
       inputs = inputs,
       analysisName = analysisName)
  class(out) = c("ConsensusTree", class(out))
  out
}



#' Create a list holding consensus calculation options.
#'
#' This function creates a list of class \code{ConsensusOptions} that holds
#' options for consensus calculations. This list holds options for a
#' single-level analysis.
#'
#'
#' @aliases newConsensusOptions ConsensusOptions
#' @param calibration Calibration method. One of \code{"full quantile", "single
#' quantile", "none"} (or a unique abbreviation of one of them).
#' @param calibrationQuantile if \code{calibration} is \code{"single
#' quantile"}, input data to a consensus calculation will be scaled such that
#' their \code{calibrationQuantile} quantiles will agree.
#' @param sampleForCalibration if \code{TRUE}, calibration quantiles will be
#' determined from a sample of network similarities. Note that using all data
#' can double the memory footprint of the function and the function may fail.
#' @param sampleForCalibrationFactor Determines the number of samples for
#' calibration: the number is \code{1/calibrationQuantile *
#' sampleForCalibrationFactor}. Should be set well above 1 to ensure accuracy
#' of the sampled quantile.
#' @param consensusQuantile Quantile at which consensus is to be defined. See
#' details.
#' @param useMean Logical: should the consensus be calculated using (weighted)
#' mean rather than a quantile?
#' @param setWeights Optional specification of weights when \code{useMean} is
#' \code{TRUE}.
#' @param analysisName Optional character string naming the consensus analysis.
#' Useful for identifying partial consensus calculation in hierarchical
#' consensus analysis.
#' @return A list of type \code{ConsensusOptions} that holds copies of the
#' input arguments.
#' @author Peter Langfelder
#' @keywords misc
newConsensusOptions = function(
      calibration = c("full quantile", "single quantile", "none"),

      # Simple quantile scaling options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE,
      sampleForCalibrationFactor = 1000,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,
      # Name to prevent clashing of files
      analysisName = "")
{
  calibration = match.arg(calibration)
  if (any(!is.finite(setWeights))) stop("Entries of 'setWeights' must all be finite.")
  if ( (consensusQuantile < 0) | (consensusQuantile > 1) )
      stop("'consensusQuantile' must be between 0 and 1.")
  out = list( calibration = calibration,
              calibrationQuantile = calibrationQuantile,
              sampleForCalibration = sampleForCalibration,
              sampleForCalibrationFactor = sampleForCalibrationFactor,
              consensusQuantile = consensusQuantile,
              useMean = useMean,
              setWeights = setWeights,
              analysisName = analysisName)
  class(out) = c("ConsensusOptions", class(out))
  out
}



#' Creates a list of correlation options.
#'
#' Convenience function to create a re-usable list of correlation options.
#'
#'
#' @aliases newCorrelationOptions CorrelationOptions
#' @param corType Character specifying the type of correlation function.
#' Currently supported options are \code{"pearson", "bicor"}.
#' @param maxPOutliers Maximum proportion of outliers for biweight
#' mid-correlation. See \code{\link{bicor}}.
#' @param quickCor Real number between 0 and 1 that controls the handling of
#' missing data in the calculation of correlations. See \code{\link{bicor}}.
#' @param pearsonFallback Specifies whether the bicor calculation should revert
#' to Pearson when median absolute deviation (mad) is zero. Recongnized values
#' are (abbreviations of) \code{"none", "individual", "all"}. If set to
#' \code{"none"}, zero mad will result in \code{NA} for the corresponding
#' correlation. If set to \code{"individual"}, Pearson calculation will be used
#' only for columns that have zero mad. If set to \code{"all"}, the presence of
#' a single zero mad will cause the whole variable to be treated in Pearson
#' correlation manner (as if the corresponding \code{robust} option was set to
#' \code{FALSE}).
#' @param cosineCorrelation Logical: calculate cosine biweight midcorrelation?
#' Cosine bicorrelation is similar to standard bicorrelation but the median
#' subtraction is not performed.
#' @param nThreads A non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param corFnc Correlation function to be called in R code. Should
#' correspoind to the value of \code{corType} above.
#' @param corOptions A list of options to be supplied to the correlation
#' function (in addition to appropriate arguments \code{x} and \code{y}).
#' @return A list containing a copy of the input arguments. The output has
#' class \code{CorrelationOptions}.
#' @author Peter Langfelder
#' @keywords misc
newCorrelationOptions = function(
      corType = c("pearson", "bicor"),
      maxPOutliers = 0.05,
      quickCor = 0,
      pearsonFallback = "individual",
      cosineCorrelation = FALSE,
      nThreads = 0,
      corFnc = if (corType=="bicor") "bicor" else "cor",
      corOptions = c(
        list(use = 'p',
             cosine = cosineCorrelation,
             quick = quickCor,
             nThreads = nThreads),
        if (corType=="bicor")
           list(maxPOutliers = maxPOutliers,
                pearsonFallback = pearsonFallback) else NULL))
{
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.")
  if (quickCor < 0) stop("quickCor must be positive.")
  if (nThreads < 0) stop("nThreads must be positive.")
  corType = match.arg(corType)
  if ( (maxPOutliers < 0) | (maxPOutliers > 1)) stop("maxPOutliers must be between 0 and 1.")
  if (quickCor < 0) stop("quickCor must be positive.")
  fallback = pmatch(pearsonFallback, .pearsonFallbacks)
  if (is.na(fallback))
      stop(paste("Unrecognized 'pearsonFallback'. Recognized values are (unique abbreviations of)\n",
           paste(.pearsonFallbacks, collapse = ", ")))
  out = list(
    corType = corType,
    maxPOutliers = maxPOutliers,
    quickCor = quickCor,
    pearsonFallback = pearsonFallback,
    pearsonFallback.code = match(pearsonFallback, .pearsonFallbacks),
    cosineCorrelation = cosineCorrelation,
    corFnc = corFnc,
    corOptions = corOptions,
    corType.code = match(corType, .corTypes))
  class(out) = c("CorrelationOptions", class(out))
  out
}



#' Create a list of network construction arguments (options).
#'
#' This function creates a reusable list of network calculation
#' arguments/options.
#'
#'
#' @aliases newNetworkOptions NetworkOptions
#' @param correlationOptions A list of correlation options. See
#' \code{\link{newCorrelationOptions}}.
#' @param replaceMissingAdjacencies Logical: should missing adjacencies be
#' replaced by zero?
#' @param power Soft-thresholding power for network construction.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param checkPower Logicel: should the power be checked for sanity?
#' @param TOMType One of \code{"none"}, \code{"unsigned"}, \code{"signed"}. If
#' \code{"none"}, adjacency will be used for clustering. If \code{"unsigned"},
#' the standard TOM will be used (more generally, TOM function will receive the
#' adjacency as input). If \code{"signed"}, TOM will keep track of the sign of
#' correlations between neighbors.
#' @param TOMDenom Character string specifying the TOM variant to be used.
#' Recognized values are \code{"min"} giving the standard TOM described in
#' Zhang and Horvath (2005), and \code{"mean"} in which the \code{min} function
#' in the denominator is replaced by \code{mean}. The \code{"mean"} may produce
#' better results but at this time should be considered experimental.
#' @return A list of class \code{NetworkOptions}.
#' @author Peter Langfelder
#' @seealso code\link{newCorrelationOptions}
#' @keywords misc
#' @export
newNetworkOptions = function(
    correlationOptions = newCorrelationOptions(),

    # Adjacency options
    replaceMissingAdjacencies = TRUE,
    power = 6,
    networkType = c("signed hybrid", "signed", "unsigned"),
    checkPower = TRUE,

    # Topological overlap options
    TOMType = c("signed", "unsigned", "none"),
    TOMDenom = c("mean", "min")) {
  if (checkPower) .checkPower(power)
  networkType = match.arg(networkType)
  TOMType = match.arg(TOMType)
  TOMDenom = match.arg(TOMDenom)
  out = c(correlationOptions,
      list(replaceMissingAdjacencies = replaceMissingAdjacencies,
           power = power,
           networkType = networkType,
           TOMType = TOMType,
           TOMDenom = TOMDenom,
           networkType.code = match(networkType, .networkTypes),
           TOMType.code = match(TOMType, .TOMTypes),
           TOMDenom.code = match(TOMDenom, .TOMDenoms)))
  class(out) = c("NetworkOptions", class(correlationOptions))
  out
}

#====================================================================================================
#
# cor, network, and consensus calculations
#
#====================================================================================================

.corCalculation = function(x, y = NULL, correlationOptions) {
  if (!inherits(correlationOptions, "CorrelationOptions"))
    stop(".corCalculation: 'correlationOptions' does not have correct type.")
  do.call(match.fun(correlationOptions$corFnc), c(list(x = x, y= y), correlationOptions$corOptions))
}


# network calculation. Returns the resulting topological overlap or
.networkCalculation = function(data, networkOptions,
                verbose = 1, indent = 0) {
  if (!inherits(networkOptions, "NetworkOptions"))
    stop(".networkCalculation: 'networkOptions' does not have correct type.")

   callVerb = max(0, verbose - 1); callInd = indent + 2
   CcorType = networkOptions$corType.code - 1
   CnetworkType = networkOptions$networkType.code - 1
   CTOMType = networkOptions$TOMType.code -1
   CTOMDenom = networkOptions$TOMDenom.code -1

   warn = as.integer(0)
   # For now return the full matrix; eventually we may return just the dissimilarity.
   # To make full use of the lower triagle space saving we'd have to modify the underlying C code
   # dramatically, otherwise it will still need to allocate the full matrix for the matrix multiplication.
   .Call("tomSimilarity_call", data,
          as.integer(CcorType), as.integer(CnetworkType),
          as.double(networkOptions$power), as.integer(CTOMType),
          as.integer(CTOMDenom),
          as.double(networkOptions$maxPOutliers),
          as.double(networkOptions$quickCor),
          as.integer(networkOptions$pearsonFallback.code),
          as.integer(networkOptions$cosineCorrelation),
          as.integer(networkOptions$replaceMissingAdjacencies),
          warn, as.integer(min(1, networkOptions$nThreads)),
          as.integer(callVerb), as.integer(callInd), PACKAGE = "WGCNA")
}

# the following is contained in consensusOptions:
#  out = list( calibration = calibration,
#              calibrationQuantile = calibrationQuantile,
#              sampleForCalibration = sampleForCalibration,
#              sampleForCalibrationFactor = sampleForCalibrationFactor,
#              consensusQuantile = consensusQuantile,
#              useMean = useMean,
#              setWeights = setWeights)




#' Calculation of a (single) consenus with optional data calibration.
#'
#' This function calculates a single consensus from given individual data,
#' optionally first calibrating the individual data to make them comparable.
#'
#' Consensus is defined as the element-wise (also known as "parallel") quantile
#' of the individual data at probability given by the \code{consensusQuantile}
#' element of \code{consensusOptions}. Depending on the value of component
#' \code{calibration} of \code{consensusOptions}, the individual data are first
#' calibrated. For \code{consensusOptions$calibration="full quantile"}, the
#' individual data are quantile normalized using
#' \code{\link{normalize.quantiles}}. For
#' \code{consensusOptions$calibration="single quantile"}, the individual data
#' are raised to a power such that the quantiles at probability
#' \code{consensusOptions$calibrationQuantile} are the same.  For
#' \code{consensusOptions$calibration="none"}, the individual data are not
#' calibrated.
#'
#' @param individualData Individual data from which the consensus is to be
#' calculated. It can be either a list or a \code{\link{multiSet}} structure.
#' Each element in \code{individulData} can in turn either be a numeric obeject
#' (vector, matrix or array) or a \code{\link{BlockwiseData}} structure.
#' @param consensusOptions A list of class \code{ConsensusOptions} that
#' contains options for the consensus calculation. A suitable list can be
#' obtained by calling function \code{\link{newConsensusOptions}}.
#' @param useBlocks When \code{individualData} contains
#' \code{\link{BlockwiseData}}, this argument can be an integer vector with
#' indices of blocks for which the calculation should be performed.
#' @param randomSeed If non-\code{NULL}, the function will save the current
#' state of the random generator, set the given seed, and restore the random
#' seed to its original state upon exit. If \code{NULL}, the seed is not set
#' nor is it restored on exit.
#' @param saveCalibratedIndividualData Logical: should calibrated individual
#' data be saved?
#' @param calibratedIndividualDataFilePattern Pattern from which file names for
#' saving calibrated individual data are determined. The conversions \code{\%a},
#' \code{\%s} and \code{\%b} will be replaced by analysis name, set number and
#' block number, respectively.
#' @param saveConsensusData Logical: should final consensus be saved
#' (\code{TRUE}) or returned in the return value (\code{FALSE})?
#' @param consensusDataFileNames Pattern from which file names for saving the
#' final consensus are determined. The conversions \code{\%a} and \code{\%b} will
#' be replaced by analysis name and block number, respectively.
#' @param getCalibrationSamples When calibration method in the
#' \code{consensusOptions} component of \code{ConsensusTree} is \code{"single
#' quantile"}, this logical argument determines whether the calibration samples
#' should be retuned within the return value.
#' @param useDiskCache Logical: should disk cache be used for consensus
#' calculations? The disk cache can be used to sture chunks of calibrated data
#' that are small enough to fit one chunk from each set into memory (blocks may
#' be small enough to fit one block of one set into memory, but not small enogh
#' to fit one block from all sets in a consensus calculation into memory at the
#' same time). Using disk cache is slower but lessens the memry footprint of
#' the calculation. As a general guide, if individual data are split into
#' blocks, we recommend setting this argument to \code{TRUE}. If this argument
#' is \code{NULL}, the function will decide whether to use disk cache based on
#' the number of sets and block sizes.
#' @param chunkSize Integer giving the chunk size. If left \code{NULL}, a
#' suitable size will be chosen automatically.
#' @param cacheDir Directory in which to save cache files. The files are
#' deleted on normal exit but persist if the function terminates abnormally.
#' @param cacheBase Base for the file names of cache files.
#' @param collectGarbage Logical: should garbage collection be forced after
#' each major calculation?
#' @param verbose Integer level of verbosity of diagnostic messages. Zero means
#' silent, higher values make the output progressively more and more verbose.
#' @param indent Indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{consensusData}{A
#' \code{\link{BlockwiseData}} list containing the consensus.}
#' \item{nSets}{Number of input data sets.}
#' \item{saveCalibratedIndividualData}{Copy of the input
#' \code{saveCalibratedIndividualData}.} \item{calibratedIndividualData}{If
#' input \code{saveCalibratedIndividualData} is \code{TRUE}, a list in which
#' each component is a \code{\link{BlockwiseData}} structure containing the
#' calibrated individual data for the corresponding input individual data set.}
#' \item{calibrationSamples}{If \code{consensusOptions$calibration} is
#' \code{"single quantile"} and \code{getCalibrationSamples} is \code{TRUE}, a
#' list in which eahc component contains the calibration samples for the
#' corresponding input individual data set.} \code{originCount}When
#' \code{consensusOptions$consensusQuantile} is 0, this vector of length
#' \code{nSets} contains the counts of the minima from each set.
#' @author Peter Langfelder
#' @seealso \code{\link[preprocessCore]{normalize.quantiles}} for quantile
#' normalization.
#' @references Consensus network analysis was originally described in
#' Langfelder P, Horvath S. Eigengene networks for studying the relationships
#' between co-expression modules. BMC Systems Biology 2007, 1:54
#' http://www.biomedcentral.com/1752-0509/1/54
#' @keywords misc
consensusCalculation = function(
      # a list or multiSet structure of either numeric vectors (possibly arrays) or blockwiseAdj objects
      individualData,
      consensusOptions,

      useBlocks = NULL,
      randomSeed = NULL,
      saveCalibratedIndividualData = FALSE,
      calibratedIndividualDataFilePattern = "calibratedIndividualData-%a-Set%s-Block%b.RData",

      # Return options: the data can be either saved or returned but not both.
      saveConsensusData = TRUE,
      consensusDataFileNames = "consensusData-%a-Block%b.RData",
      getCalibrationSamples= FALSE,

      # Internal handling of data
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      # Behaviour


#' Iterative garbage collection.
#'
#' Performs garbage collection until free memory idicators show no change.
#'
#'
#' @return None.
#' @author Steve Horvath
#' @keywords utilities
      collectGarbage = FALSE,
      verbose = 1, indent = 0)
{
  nSets = length(individualData)

  if (! is.multiSet(individualData))
     individualData = list2multiSet(individualData)

  setNames = names(individualData)
  if (is.null(setNames)) setNames = rep("", nSets)

  blockwise = inherits(individualData[[1]]$data, "BlockwiseData")

  if (!blockwise)
  {
    blockDimnames = .mtd.checkDimConsistencyAndGetDimnames(individualData)
    blockLengths = length(individualData[[1]]$data)
    blockAttributes = attributes(individualData[[1]]$data)
    metaData = list()
  } else {
    blockLengths = BD.blockLengths(individualData[[1]]$data)
    blockAttributes = individualData[[1]]$data$attributes
    metaData = BD.getMetaData(individualData[[1]]$data, blocks = 1)
  }
  nBlocks = length(blockLengths)

  spaces = indentSpaces(indent)

  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed)
    }
    set.seed(randomSeed)
  }

  setWeights = consensusOptions$setWeights
  if (is.null(setWeights)) setWeights = rep(1, nSets)

  if (length(setWeights)!=nSets)
    stop("Length of 'setWeights' must equal the number of sets.")

  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))
  if (is.null(useDiskCache)) useDiskCache = .useDiskCache(individualData, chunkSize = chunkSize)

  # Initialize various variables

  if (getCalibrationSamples)
  {
    if (!consensusOptions$sampleForCalibration)
      stop(paste("Incompatible input options: calibrationSamples can only be returned",
                 "if sampleForCalibration is TRUE."))
    calibrationSamples = list()
  }

  blockLevels = 1:nBlocks
  if (is.null(useBlocks)) useBlocks = blockLevels
  useBlockIndex = match(useBlocks, blockLevels)

  if (!all(useBlocks %in% blockLevels))
    stop("All entries of 'useBlocks' must be valid block levels.")

  if (any(duplicated(useBlocks)))
    stop("Entries of 'useBlocks' must be unique.")

  nUseBlocks = length(useBlocks)
  if (nUseBlocks==0)
    stop("'useBlocks' cannot be non-NULL and empty at the same time.")

  consensus.out = list()

  consensusFiles = rep("", nUseBlocks)
  originCount = rep(0, nSets)

  calibratedIndividualDataFileNames = NULL
  if (saveCalibratedIndividualData)
  {
    calibratedIndividualDataFileNames = matrix("", nSets, nBlocks)
    for (set in 1:nSets) for (b in 1:nBlocks)
      calibratedIndividualDataFileNames[set, b] =
                .processFileName(calibratedIndividualDataFilePattern, setNumber = set, setNames = setNames,
                           blockNumber = b, analysisName = consensusOptions$analysisName)
  }
  if (collectGarbage) gc()

  calibratedIndividualData.saved = list()
  consensusData = NULL
  dataFiles = character(nUseBlocks)

  # Here's where the analysis starts
  for (blockIndex in 1:nUseBlocks)
  {
    block = useBlockIndex[blockIndex]

    if (verbose>1) printFlush(paste0(spaces, "..Working on block ", block, "."))

    scaleQuant = rep(1, nSets)
    scalePowers = rep(1, nSets)

    useDiskCache1 = useDiskCache && nSets > 1;  ### No need to use disk cache when there is only 1 set.
    # Set up file names or memory space to hold the set Data
    if (useDiskCache1)
    {
      nChunks = ceiling(blockLengths[block]/chunkSize)
      chunkFileNames = array("", dim = c(nChunks, nSets))
      on.exit(.checkAndDelete(chunkFileNames), add = TRUE)
    } else nChunks = 1

    if (nChunks==1) useDiskCache1 = FALSE
    if (!useDiskCache1)
    {
      # Note: setTomDS will contained the scaled set Data matrices.
      calibratedData = array(0, dim = c(blockLengths[block], nSets))
    }

    # sample entry indices from the distance structure for Data scaling, if requested

    if (consensusOptions$calibration=="single quantile" &&
           consensusOptions$sampleForCalibration)
    {
      qx = min(consensusOptions$calibrationQuantile, 1-consensusOptions$calibrationQuantile)
      nScGenes = min(consensusOptions$sampleForCalibrationFactor * 1/qx, blockLengths[block])
      scaleSample = sample(blockLengths[block], nScGenes)
      if (getCalibrationSamples)
        calibrationSamples[[blockIndex]] = list(sampleIndex = scaleSample,
                                            TOMSamples = matrix(NA, nScGenes, nSets))
    }
    if (consensusOptions$calibration %in% c("single quantile", "none"))
    {
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste0(spaces, "....Working on set ", set, " (", setNames[set], ")"))

        # We need to drop dimensions here but we may need them later. Keep that in mind.
        tomDS = as.numeric(.getData(individualData[[set]]$data, block, simplify = TRUE))

        if (consensusOptions$calibration=="single quantile")
        {
          # Scale Data so that calibrationQuantile agree in each set
          if (consensusOptions$sampleForCalibration)
          {
            if (consensusOptions$getCalibrationSamples)
            {
              calibrationSamples[[blockIndex]]$dataSamples[, set] = tomDS[scaleSample]
              scaleQuant[set] = quantile(calibrationSamples[[blockIndex]]$dataSamples[, set],
                                         probs = consensusOptions$calibrationQuantile, type = 8)
            } else {
              scaleQuant[set] = quantile(tomDS[scaleSample], probs = consensusOptions$calibrationQuantile,
                                         type = 8)
            }
          } else
            scaleQuant[set] = quantile(x = tomDS, probs = consensusOptions$calibrationQuantile, type = 8)
          if (set>1)
          {
             scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
             tomDS = tomDS^scalePowers[set]
          }
          if (saveCalibratedIndividualData)
             calibratedIndividualData.saved[[set]] =
                  addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(tomDS, blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, block])
        }

        # Save the calculated Data either to disk in chunks or to memory.

        if (useDiskCache1)
        {
          if (verbose > 3) printFlush(paste(spaces, "......saving Data similarity to disk cache.."))
          sc = .saveChunks(tomDS, chunkSize, cacheBase, cacheDir = cacheDir)
          chunkFileNames[, set] = sc$files
          chunkLengths = sc$chunkLengths
        } else {
          calibratedData[, set] = tomDS
        }
      }
      if (collectGarbage) gc()
    } else if (consensusOptions$calibration=="full quantile")
    {
      # Step 1: load each data set, get order, split Data into chunks according to order, and save.
      if (verbose>1) printFlush(paste0(spaces, "..working on quantile normalization"))
      if (useDiskCache1)
      {
        orderFiles = rep("", nSets)
        on.exit(.checkAndDelete(orderFiles),add = TRUE)
      }
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste0(spaces, "....Working on set ", set, " (", setNames[set], ")"))
        tomDS = as.numeric(.getData(individualData[[set]]$data, block, simplify = TRUE))

        if (useDiskCache1)
        {
          # Order Data (this may take a long time...)
          if (verbose > 3) printFlush(paste0(spaces, "......ordering Data"))
          time = system.time({order1 = .qorder(tomDS)})
          if (verbose > 1) { printFlush("Time to order Data:"); print(time); }
          # save the order
          orderFiles[set] = tempfile(pattern = paste0(".orderForSet", set), tmpdir = cacheDir)
          if (verbose > 3) printFlush(paste0(spaces, "......saving order and ordered Data"))
          save(order1, file = orderFiles[set])
          # Save ordered tomDS into chunks
          tomDS.ordered = tomDS[order1]
          sc = .saveChunks(tomDS.ordered, chunkSize, cacheBase, cacheDir = cacheDir)
          chunkFileNames[, set] = sc$files
          chunkLengths = sc$chunkLengths
        } else {
           calibratedData[, set] = tomDS
        }
      }
      if (useDiskCache1)
      {
        # Step 2: Load chunks one by one and quantile normalize
        if (verbose > 2) printFlush(paste0(spaces, "....quantile normalizing chunks"))
        for (c in 1:nChunks)
        {
          if (verbose > 3) printFlush(paste0(spaces, "......QN for chunk ", c, " of ", nChunks))
          chunkData = matrix(NA, chunkLengths[c], nSets)
          for (set in 1:nSets)
            chunkData[, set] = .loadObject(chunkFileNames[c, set])

          time = system.time({ chunk.norm = normalize.quantiles(chunkData, copy = FALSE);})
          if (verbose > 1) { printFlush("Time to QN chunk:"); print(time); }
          # Save quantile normalized chunks
          for (set in 1:nSets)
          {
            temp = chunk.norm[, set]
            save(temp, file = chunkFileNames[c, set])
          }
        }
        if (verbose > 2) printFlush(paste0(spaces, "....putting together full QN'ed Data"))
        # Put together full Data
        for (set in 1:nSets)
        {
           load(orderFiles[set])
           start = 1
           for (c in 1:nChunks)
           {
             end = start + chunkLengths[c] - 1
             tomDS[order1[start:end]] = .loadObject(chunkFileNames[c, set], size = chunkLengths[c])
             start = start + chunkLengths[c]
           }
           if (saveCalibratedIndividualData)
              calibratedIndividualData.saved[[set]] = addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(tomDS, blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, blockIndex])
           .saveChunks(tomDS, chunkSize, fileNames = chunkFileNames[, set])
           unlink(orderFiles[set])
        }
      } else {
        # If disk cache is not being used, simply call normalize.quantiles on the full set.
        if (nSets > 1) calibratedData = normalize.quantiles(calibratedData)
        if (saveCalibratedIndividualData) for (set in 1:nSets)
        {
           calibratedIndividualData.saved[[set]] = addBlockToBlockwiseData(
                        calibratedIndividualData.saved[[set]],
                        .setAttrFromList(calibratedData[, set], blockAttributes[[blockIndex]]),
                        external = TRUE,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = calibratedIndividualDataFileNames[set, blockIndex])

        }
      }
    } else stop("Unrecognized value of 'calibration' in consensusOptions: ", consensusOptions$calibration)

    # Calculate consensus
    if (verbose > 2)
      printFlush(paste(spaces, "....Calculating consensus"))

    # create an empty consTomDS distance structure.
    consTomDS = numeric(blockLengths[block])

    if (useDiskCache1)
    {
      start = 1
      for (chunk in 1:nChunks)
      {
        if (verbose > 3) printFlush(paste(spaces, "......working on chunk", chunk))
        end = start + chunkLengths[chunk] - 1
        setChunks = array(0, dim = c(chunkLengths[chunk], nSets))
        for (set in 1:nSets)
        {
          load(file = chunkFileNames[chunk, set])
          setChunks[, set] = temp
          file.remove(chunkFileNames[chunk, set])
        }
        setWeightMat = matrix(setWeights, chunkLengths[chunk], nSets, byrow = TRUE)
        tmp = .consensusCalculation.base(setChunks, useMean = consensusOptions$useMean,
                                         setWeightMat = setWeightMat,
                                         consensusQuantile = consensusOptions$consensusQuantile)
        consTomDS[start:end] = tmp$consensus
        if (!is.null(tmp$originCount))
        {
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
        }
        start = end + 1
      }
    } else {
      setWeightMat = matrix(setWeights, blockLengths[block], nSets, byrow = TRUE)
      tmp = .consensusCalculation.base(calibratedData, useMean = consensusOptions$useMean,
                                       setWeightMat = setWeightMat,
                                       consensusQuantile = consensusOptions$consensusQuantile)
      consTomDS[] = tmp$consensus
      if (!is.null(tmp$originCount))
      {
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
      }
    }

    # Save the consensus Data if requested
    if (saveConsensusData)
    {
       if (!grepl("%b", consensusDataFileNames))
         stop(paste("File name for consensus data must contain the tag %b somewhere in the file name -\n",
                    "   - this tag will be replaced by the block number. "))
       dataFiles[blockIndex] = .substituteTags(consensusDataFileNames, c("%b", "%a"),
                                              c(block, consensusOptions$analysisName[1]))
    }
    consensusData = addBlockToBlockwiseData(
                        consensusData,
                        .setAttrFromList(consTomDS, blockAttributes[[blockIndex]]),
                        external = saveConsensusData,
                        recordAttributes = TRUE,
                        metaData = metaData,
                        blockFile = if (saveConsensusData) dataFiles[blockIndex] else NULL)
    if (collectGarbage) gc()
  }

  list(
       #blockwiseData
       consensusData = consensusData,
       # int
       nSets = nSets,
       # Logical
       saveCalibratedIndividualData = saveCalibratedIndividualData,
       # List of blockwise data of length nSets
       calibratedIndividualData = calibratedIndividualData.saved,
       # List with one component per block
       calibrationSamples = if (getCalibrationSamples) calibrationSamples else NULL,
       # Numeric vector with nSets components
       originCount = originCount,

       consensusOptions = consensusOptions

      )
}



#==========================================================================================================
#
# Hierarchical consensus calculation
#
#==========================================================================================================

#  hierarchical consensus tree: a list with the following components:
#    inputs: either an atomic character vector whose entries match names of individualData, or a list in
#      which each component can either be a single character string giving a name in individualDara, or
#      another hierarchical consensus tree.
#    consensusOptions: a list of class ConsensusOptions
#    analysisName: optional, analysis name used for naming files when saving to disk.

# Here individualData is a list or multiSet in which every component is either a blockwiseData instance or
# a numeric object (matrix, vector etc). Function consensusCalculation handles both.



#' Hierarchical consensus calculation
#'
#' Hierarchical consensus calculation with optional data calibration.
#'
#' This function calculates consensus in a hierarchical manner, using a
#' separate (and possibly different) set of consensus options at each step. The
#' "recipe" for the consensus calculation is supplied in the argument
#' \code{consensusTree}.
#'
#' The argument \code{consensusTree} should have the following components: (1)
#' \code{inputs} must be either a character vector whose components match
#' \code{names(inputData)}, or consensus trees in the own right. (2)
#' \code{consensusOptions} must be a list of class \code{"ConsensusOptions"}
#' that specifies options for calculating the consensus. A suitable set of
#' options can be obtained by calling \code{\link{newConsensusOptions}}. (3)
#' Optionally, the component \code{analysisName} can be a single character
#' string giving the name for the analysis. When intermediate results are
#' returned, they are returned in a list whose names will be set from
#' \code{analysisName} components, if they exist.
#'
#' The actual consensus calculation at each level of the consensus tree is
#' carried out in function \code{\link{consensusCalculation}}. The consensus
#' options for each individual consensus calculation are independent from one
#' another, i.e., the consensus options for different steps can be different.
#'
#' @param individualData Individual data from which the consensus is to be
#' calculated. It can be either a list or a \code{\link{multiSet}} structure.
#' Each element in \code{individulData} can in turn either be a numeric object
#' (vector, matrix or array) or a \code{\link{BlockwiseData}} structure.
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param level Integer which the user should leave at 1.  This serves to keep
#' default set names unique.
#' @param useBlocks When \code{individualData} contains
#' \code{\link{BlockwiseData}}, this argument can be an integer vector with
#' indices of blocks for which the calculation should be performed.
#' @param randomSeed If non-\code{NULL}, the function will save the current
#' state of the random generator, set the given seed, and restore the random
#' seed to its original state upon exit. If \code{NULL}, the seed is not set
#' nor is it restored on exit.
#' @param saveCalibratedIndividualData Logical: should calibrated individual
#' data be saved?
#' @param calibratedIndividualDataFilePattern Pattern from which file names for
#' saving calibrated individual data are determined. The conversions \code{\%a},
#' \code{\%s} and \code{\%b} will be replaced by analysis name, set number and
#' block number, respectively.
#' @param saveConsensusData Logical: should final consensus be saved
#' (\code{TRUE}) or returned in the return value (\code{FALSE})?
#' @param consensusDataFileNames Pattern from which file names for saving the
#' final consensus are determined. The conversions \code{\%a} and \code{\%b} will
#' be replaced by analysis name and block number, respectively.
#' @param getCalibrationSamples When calibration method in the
#' \code{consensusOptions} component of \code{ConsensusTree} is \code{"single
#' quantile"}, this logical argument determines whether the calibration samples
#' should be returned within the return value.
#' @param keepIntermediateResults Logical: should results of intermediate
#' consensus calculations (if any) be kept? These are always returned as
#' \code{BlockwiseData} whose data are saved to disk.
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
#' @param collectGarbage Logical: should garbage collection be forced after
#' each major calculation?
#' @param verbose Integer level of verbosity of diagnostic messages. Zero means
#' silent, higher values make the output progressively more and more verbose.
#' @param indent Indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list containing the output of the top level call to
#' \code{\link{consensusCalculation}}; if \code{keepIntermediateResults} is
#' \code{TRUE}, component \code{inputs} contains a (possibly recursive) list of
#' the results of intermediate consensus calculations. Names of the
#' \code{inputs} list are taken from the corresponding \code{analysisName}
#' components if they exist, otherwise from names of the corresponding
#' \code{inputs} components of the supplied \code{consensusTree}.  See example
#' below for an example of a relatively simple consensus tree.
#' @author Peter Langfelder
#' @seealso \code{\link{newConsensusOptions}} for obtaining a suitable list of
#' consensus options
#'
#' \code{\link{consensusCalculation}} for the actual calculation of a consensus
#' that underpins this function.
#' @keywords misc
#' @examples
#' # We generate 3 simple matrices
#' set.seed(5)
#' data = replicate(3, matrix(rnorm(10*100), 10, 100))
#' names(data) = c("Set1", "Set2", "Set3")
#' # Put together a consensus tree. In this example the final consensus uses
#' # as input set 1 and a consensus of sets 2 and 3.
#'
#' # First define the consensus of sets 2 and 3:
#' consTree.23 = newConsensusTree(
#'            inputs = c("Set2", "Set3"),
#'            consensusOptions = newConsensusOptions(calibration = "none",
#'                                consensusQuantile = 0.25),
#'            analysisName = "Consensus of sets 1 and 2")
#'
#' # Now define the final consensus
#' consTree.final = newConsensusTree(
#'    inputs = list("Set1", consTree.23),
#'    consensusOptions = newConsensusOptions(calibration = "full quantile",
#'                                consensusQuantile = 0),
#'    analysisName = "Final consensus")
#' # FIXME
#' \dontrun{
#' consensus = hierarchicalConsensusCalculation(
#'   individualData = data,
#'   consensusTree = consTree.final,
#'   saveConsensusData = FALSE,
#'   keepIntermediateResults = FALSE)
#'
#' names(consensus)
#' }
hierarchicalConsensusCalculation = function(
   individualData,

   consensusTree,

   level = 1,
   useBlocks = NULL,
   randomSeed = NULL,
   saveCalibratedIndividualData = FALSE,
   calibratedIndividualDataFilePattern = "calibratedIndividualData-%a-Set%s-Block%b.RData",

   # Return options: the data can be either saved or returned but not both.
   saveConsensusData = TRUE,
   consensusDataFileNames = "consensusData-%a-Block%b.RData",
   getCalibrationSamples= FALSE,

   # Return the intermediate results as well?
   keepIntermediateResults = FALSE,

   # Internal handling of data
   useDiskCache = NULL, chunkSize = NULL,
   cacheDir = ".",
   cacheBase = ".blockConsModsCache",

   # Behaviour
   collectGarbage = FALSE,
   verbose = 1, indent = 0)
{
  spaces = indentSpaces(indent)
  individualNames = names(individualData)
  if (is.null(individualNames))
    stop("'individualData' must be a named list.")

  if (!is.multiSet(individualData))
    individualData = list2multiSet(individualData)

  if (!"inputs" %in% names(consensusTree))
    stop("'consensusTree' must contain component 'inputs'.")

  if (!"consensusOptions" %in% names(consensusTree))
    stop("'consensusTree' must contain component 'consensusOptions'.")

  if (!is.null(randomSeed)) {
    if (exists(".Random.seed")) {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed)
    }
    set.seed(randomSeed)
  }

  # Set names for consensusTree$inputs so that the names are informative.

  if (is.null(names(consensusTree$inputs))) {
    names(consensusTree$inputs) = paste0("Level.", level, ".Input.",
                                                  1:length(consensusTree$inputs))
    validInputNames = FALSE
  } else
    validInputNames = TRUE

  isChar = base::sapply(consensusTree$inputs, is.character)
  names(consensusTree$inputs)[isChar] = consensusTree$inputs[isChar]
  if (!is.null(consensusTree$analysisName))
    consensusTree$consensusOptions$analysisName = consensusTree$analysisName

  # Recursive step if necessary
  if (verbose > 0)
    printFlush(
      paste0(
        spaces,
        "------------------------------------------------------------------\n",
        spaces, "   Working on ", consensusTree$consensusOptions$analysisName,
        "\n", spaces,
        "------------------------------------------------------------------"))
  names(consensusTree$inputs) = make.unique(make.names(names(consensusTree$inputs)))
  inputs0 = multiSet.mapply(function(inp1, name) {
     if (is.character(inp1)) {
        if (!inp1 %in% names(individualData))
          stop("Element '", inp1, "' is not among names of 'individualData'.")
        inp1
     } else {
        if ("analysisName" %in% names(inp1))
          name1 = inp1$analysisName
        else name1 = name
        inp1$consensusOptions$analysisName = name1
        hierarchicalConsensusCalculation(individualData, inp1,
                   useBlocks = useBlocks,
                   level = level + 1,
                   randomSeed = NULL,
                   saveCalibratedIndividualData = saveCalibratedIndividualData,
                   calibratedIndividualDataFilePattern =calibratedIndividualDataFilePattern,
                   saveConsensusData = saveConsensusData,
                   consensusDataFileNames = consensusDataFileNames,
                   getCalibrationSamples = getCalibrationSamples,
                   keepIntermediateResults = keepIntermediateResults,
                   useDiskCache = useDiskCache,
                   chunkSize = chunkSize,
                   cacheDir = cacheDir,
                   cacheBase = cacheBase,
                   collectGarbage = collectGarbage,
                   verbose = verbose -2, indent = indent + 2)
     }
  }, consensusTree$inputs, names(consensusTree$inputs))

  names(inputs0) = names(consensusTree$inputs)

  inputData = apply(inputs0, function(inp1) {
    if (is.character(inp1)) {
       individualData[[inp1]]$data
    } else
       inp1$consensusData
  })

  inputIsIntermediate = !base::sapply(consensusTree$inputs, is.character)

  # Need to check that all inputData have the same format.
  # In particular, some could be plain numeric data and
  # some could be BlockwiseData.

  nInputs1 = length(inputData)
  isBlockwise = apply(inputData, inherits, "BlockwiseData", mdaSimplify = TRUE)
  if (any(!isBlockwise)) for (i in which(!isBlockwise))
     inputData[[i]]$data = newBlockwiseData(list(inputData[[i]]$data), external = FALSE)

  names(inputData) = names(consensusTree$inputs)

  # Calculate the consensus

  if (verbose > 0)
    printFlush(paste0(spaces, "..Final consensus calculation.."))
  consensus = consensusCalculation(
      individualData = inputData,
      consensusOptions = consensusTree$consensusOptions,
      randomSeed = NULL,
      saveCalibratedIndividualData = saveCalibratedIndividualData,
      calibratedIndividualDataFilePattern =calibratedIndividualDataFilePattern,
      saveConsensusData = saveConsensusData,
      consensusDataFileNames = consensusDataFileNames,
      getCalibrationSamples = getCalibrationSamples,
      useDiskCache = useDiskCache,
      chunkSize = chunkSize,
      cacheDir = cacheDir,
      cacheBase = cacheBase,
      collectGarbage = collectGarbage,
      verbose = verbose-1, indent = indent+1)

  if (saveConsensusData && !keepIntermediateResults && any(inputIsIntermediate))
    apply(inputData[inputIsIntermediate], BD.checkAndDeleteFiles)

  out = c(consensus, ifelse(keepIntermediateResults,
                            list(inputs = inputs0), NULL))

  out
}


#==========================================================================================================
#
# Simple hierarchical consensus calculation from numeric data, with minimum checking and no calibration.
#
#==========================================================================================================

# Simpler version of consensus calculation, suitable for small data where calibration is not
# necessary.



#' Simple calculation of a single consenus
#'
#' This function calculates a single consensus from given individual data.
#'
#' Consensus is defined as the element-wise (also known as "parallel") quantile
#' of of the individual data at probability given by the
#' \code{consensusQuantile} element of \code{consensusOptions}.
#'
#' @param individualData Individual data from which the consensus is to be
#' calculated. It can be either a list or a \code{\link{multiSet}} structure
#' in which each element is a numeric vector or array.
#' @param consensusOptions A list of class \code{ConsensusOptions} that
#' contains options for the consensus calculation. A suitable list can be
#' obtained by calling function \code{\link{newConsensusOptions}}.
#' @param verbose Integer level of verbosity of diagnostic messages. Zero means
#' silent, higher values make the output progressively more and more verbose.
#' @param indent Indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A numeric vector or array of the same dimensions as each element of
#' \code{individualData}
#' @author Peter Langfelder
#' @seealso \code{\link{consensusCalculation}} for consensus calculation that
#' can work with \code{\link{BlockwiseData}} and can calibrate data before
#' calculating consensus.
#' @references Consensus network analysis was originally described in
#' Langfelder P, Horvath S. Eigengene networks for studying the relationships
#' between co-expression modules. BMC Systems Biology 2007, 1:54
#' http://www.biomedcentral.com/1752-0509/1/54
#' @keywords misc
simpleConsensusCalculation = function(
   # multiSet or  list of numeric vectors
   individualData,
   consensusOptions,
   verbose = 1, indent = 0)
{
  nSets = length(individualData)

  if (is.multiSet(individualData))
    individualData = multiSet2list(individualData)

  if (consensusOptions$useMean) {
    setWeights = consensusOptions$setWeights
    if (is.null(setWeights)) setWeights = rep(1, nSets)
    if (length(setWeights)!=nSets)
      stop("Length of 'setWeights' must equal the number of sets.")
  } else setWeightMat = NULL

  .consensusCalculation.base.FromList(individualData, useMean = consensusOptions$useMean,
                                   setWeights = setWeights,
                                   consensusQuantile = consensusOptions$consensusQuantile)$consensus
}


# Simple hierarchical consensus



#' Simple hierarchical consensus calculation
#'
#' Hierarchical consensus calculation without calibration.
#'
#' This function calculates consensus in a hierarchical manner, using a
#' separate (and possibly different) set of consensus options at each step. The
#' "recipe" for the consensus calculation is supplied in the argument
#' \code{consensusTree}.
#'
#' The argument \code{consensusTree} should have the following components: (1)
#' \code{inputs} must be either a character vector whose components match
#' \code{names(inputData)}, or consensus trees in the own right. (2)
#' \code{consensusOptions} must be a list of class \code{"ConsensusOptions"}
#' that specifies options for calculating the consensus. A suitable set of
#' options can be obtained by calling \code{\link{newConsensusOptions}}. (3)
#' Optionally, the component \code{analysisName} can be a single character
#' string giving the name for the analysis. When intermediate results are
#' returned, they are returned in a list whose names will be set from
#' \code{analysisName} components, if they exist.
#'
#' Unlike the similar function \code{\link{hierarchicalConsensusCalculation}},
#' this function ignores the calibration settings in the
#' \code{consensusOptions} component of \code{consensusTree}; no calibration of
#' input data is performed.
#'
#' The actual consensus calculation at each level of the consensus tree is
#' carried out in function \code{\link{simpleConsensusCalculation}}. The
#' consensus options for each individual consensus calculation are independent
#' from one another, i.e., the consensus options for different steps can be
#' different.
#'
#' @param individualData Individual data from which the consensus is to be
#' calculated. It can be either a list or a \code{\link{multiSet}} structure.
#' Each element in \code{individulData} should be a numeric object (vector,
#' matrix or array).
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param level Integer which the user should leave at 1.  This serves to keep
#' default set names unique.
#' @return A list with a single component \code{consensus}, containing the
#' consensus data of the same dimensions as the individual entries in the input
#' \code{individualData}. This perhaps somewhat cumbersome convention is used
#' to make the output compatible with that of
#' \code{\link{hierarchicalConsensusCalculation}}.
#' @author Peter Langfelder
#' @seealso \code{\link{simpleConsensusCalculation}} for a "single-level"
#' consensus calculation
#'
#' \code{\link{hierarchicalConsensusCalculation}} for hierarchical consensus
#' calculation with calibration
#' @keywords misc
simpleHierarchicalConsensusCalculation = function(
   # multiSet or  list of numeric vectors
   individualData,
   consensusTree,
   level = 1)

{
  individualNames = names(individualData)
  if (is.null(individualNames))
    stop("'individualData' must be named.")

  if (is.null(names(consensusTree$inputs)))
    names(consensusTree$inputs) = paste0("Level.", level, ".Input.",
                                                  1:length(consensusTree$inputs))

  if (is.multiSet(individualData))
    individualData = multiSet2list(individualData)

  isChar = sapply(consensusTree$inputs, is.character)
  names(consensusTree$inputs)[isChar] = consensusTree$inputs[isChar]

  # Recursive step if necessary

  names(consensusTree$inputs) = make.unique(make.names(names(consensusTree$inputs)))
  inputData = mapply(function(inp1, name)
   {
     if (is.character(inp1))
     {
        if (!inp1 %in% names(individualData))
          stop("Element '", inp1, "' is not among names of 'individualData'.")
        individualData[[inp1]]
     } else {
        if ("analysisName" %in% names(inp1)) name1 = inp1$analysisName else name1 = name
        inp1$consensusOptions$analysisName = name1
        simpleHierarchicalConsensusCalculation(individualData, inp1,
                   level = level + 1)
     }
  }, consensusTree$inputs, names(consensusTree$inputs), SIMPLIFY = FALSE)

  # Calculate the consensus

  simpleConsensusCalculation(
      individualData = inputData,
      consensusOptions = consensusTree$consensusOptions)
}



#==========================================================================================================
#
# Utility functions for handling possibly disk-backed blockwise data.
#
#==========================================================================================================

.getAttributesOrEmptyList = function(object)
{
  att = attributes(object)
  if (is.null(att)) list() else att
}



#' Create, merge and expand BlockwiseData objects
#'
#' These functions create, merge and expand BlockwiseData objects for holding
#' in-memory or disk-backed blockwise data. Blockwise here means that the data
#' is too large to be loaded or processed in one piece and is therefore split
#' into blocks that can be handled one by one in a divide-and-conquer manner.
#'
#' Several functions in this package use the concept of blockwise, or
#' "divide-and-conquer", analysis. The BlockwiseData class is meant to hold the
#' blockwise data, or all necessary information about blockwise data that is
#' saved in disk files.
#'
#' The data can be stored in disk files (one file per block) or in-memory. In
#' memory storage is provided so that same code can be used for both smaller
#' (single-block) data where disk storage could slow down operations as well as
#' larger data sets where disk storage and block by block analysis are
#' necessary.
#'
#' @aliases newBlockwiseData BlockwiseData mergeBlockwiseData
#' addBlockToBlockwiseData
#' @param data A list in which each component carries the data of a single
#' block.
#' @param external Logical: should the data be disk-backed (\code{TRUE}) or
#' in-memory (\code{FALSE})?
#' @param fileNames When \code{external} is \code{TRUE}, this argument must be
#' a character vector of the same length as \code{data}, giving the file names
#' for the data to be saved to, or where the data is already located.
#' @param doSave Logical: should data be saved? If this is \code{FALSE}, it is
#' the user's responsibility to ensure the files supplied in \code{fileNames}
#' already exist and contain the expected data.
#' @param recordAttributes Logical: should \code{attributes} of the given data
#' be recorded within the object?
#' @param metaData A list giving any additional meta-data for \code{data} that
#' should be attached to the object.
#' @param bwData An existing \code{BlockwiseData} object.
#' @param blockData A vector, matrix or array carrying the data of a single
#' block.
#' @param blockFile File name where data contained in \code{blockData} should
#' be saved.
#' @param ... One or more objects of class \code{BlockwiseData}.
#' @return All three functions return a list with the class set to
#' \code{"BlockwiseData"}, containing the following components:
#' \item{external}{Copy of the input argument \code{external}} \item{data}{If
#' \code{external} is \code{TRUE}, an empty list, otherwise a copy of the input
#' \code{data}.} \item{fileNames}{Copy of the input argument \code{fileNames}.}
#' \item{lengths}{A vector of lengths (results of \code{\link{length}}) of
#' elements of \code{data}.} \item{attributes}{If input \code{recordAttributes}
#' is \code{TRUE}, a list with one component per block (component of
#' \code{data}); each component is in turn a list of attributes of that
#' component of \code{data}.} \item{metaData}{A copy of the input
#' \code{metaData}.}
#' @section Warning: The definition of \code{BlockwiseData} should be
#' considered experimental and may change in the future.
#' @author Peter Langfelder
#' @seealso Other functions on \code{BlockwiseData}:
#'
#' \code{\link{BD.getData}} for retrieving data
#'
#' \code{\link{BD.actualFileNames}} for retrieving file names of files
#' containing data
#'
#' \code{\link{BD.nBlocks}} for retrieving the number of blocks
#'
#' \code{\link{BD.blockLengths}} for retrieving block lengths
#'
#' \code{\link{BD.getMetaData}} for retrieving metadata
#'
#' \code{\link{BD.checkAndDeleteFiles}} for deleting files of an unneeded
#' object.
#' @keywords misc
newBlockwiseData = function(data, external = FALSE, fileNames = NULL,
                             doSave = external,
                             recordAttributes = TRUE,
                             metaData = list())
{
  if (length(external)==0)
     stop("'external' must be logical of length 1.")

  if (!is.null(dim(data)) || !is.list(data))
      stop("'data' must be a list without dimensions.")

  if (recordAttributes)
  {
    attributes = lapply(data, .getAttributesOrEmptyList)
  } else
    attributes = NULL

  nBlocks = length(data)

  if (length(metaData) > 0)
  {
     if (length(metaData)!=nBlocks)
         stop("If 'metaData' are given, it must be a list with one component per component of 'data'.")
  } else {
    metaData = .listRep(list(), nBlocks)
  }

  lengths = sapply(data, length)

  if (doSave && !external)
    warning("newBlockwiseData: Cannot save when 'external' is not TRUE. Data will not be written to disk.")

  if (external)
  {
    if (is.null(fileNames)) stop("When 'external' is TRUE, 'fileNames' must be given.")
  } else
    fileNames = NULL

  out = list(external = external, data = if (external) list() else data, fileNames = fileNames,
               lengths = lengths, attributes = attributes, metaData = metaData)
  if (doSave && external)
  {
    if (nBlocks!=length(fileNames)) stop("Length of 'data' and 'fileNames' must be the same.")
    mapply(function(object, f) save(object, file = f), data, fileNames)
  }

  class(out) = "BlockwiseData"
  out
}


mergeBlockwiseData = function(...)
{
  args = list(...)
  args = args[ sapply(args, length) > 0]
  if (!all(sapply(args, inherits, "BlockwiseData")))
    stop("All arguments must be of class 'BlockwiseData'.")

  external1 = .checkLogicalConsistency(args, "external")
  .checkListNamesConsistency(lapply(args, getElement, "attributes"), "attributes")
  .checkListNamesConsistency(lapply(args, getElement, "metaData"), "metaData")

  out = list(external = external1, data = do.call(c, lapply(args, .getElement, "data")),
             fileNames = do.call(c, lapply(args, .getElement, "fileNames")),
             lengths = do.call(c, lapply(args, .getElement, "lengths")),
             attributes = do.call(c, lapply(args, .getElement, "attributes")),
             metaData = do.call(c, lapply(args, .getElement, "metaData")))
  class(out) = "BlockwiseData"
  out
}


# Under normal circumstance arguments external, dist and diag should not be set by the calling fnc, but this
# function can also be used to start a new instance of blockwise data.

addBlockToBlockwiseData = function(bwData,
               blockData,
               external = bwData$external,
               blockFile = NULL,
               doSave = external,
               recordAttributes = !is.null(bwData$attributes),
               metaData = NULL)
{
  badj1 = newBlockwiseData(external = external,
                           data = if (is.null(blockData)) NULL else list(blockData),
                           fileNames = blockFile,
                           recordAttributes = recordAttributes,
                           metaData = list(metaData),
                           doSave = doSave)
  mergeBlockwiseData(bwData, badj1)
}

BD.actualFileNames = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.")
  if (bwData$external) bwData$fileNames else character(0)
}

BD.nBlocks = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.")
  length(bwData$lengths)
}


BD.blockLengths = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise structure.")
  bwData$lengths
}

BD.getMetaData = function(bwData, blocks = NULL, simplify = TRUE)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise structure.")
  if (is.null(blocks)) blocks = 1:BD.nBlocks(bwData)
  if ( (length(blocks)==0) | any(!is.finite(blocks)))
    stop("'block' must be present and finite.")

  if (any(blocks<1) | (blocks > BD.nBlocks(bwData)))
    stop("All entries in 'block' must be between 1 and ", BD.nBlocks(bwData))

  out = bwData$metaData[blocks]
  if (length(blocks)==1 && simplify)
    out= out[[1]]
  out
}



#' Various basic operations on \code{BlockwiseData} objects.
#'
#' These functions implement basic operations on \code{\link{BlockwiseData}}
#' objects. Blockwise here means that the data is too large to be loaded or
#' processed in one piece and is therefore split into blocks that can be
#' handled one by one in a divide-and-conquer manner.
#'
#' Several functions in this package use the concept of blockwise, or
#' "divide-and-conquer", analysis. The BlockwiseData class is meant to hold the
#' blockwise data, or all necessary information about blockwise data that is
#' saved in disk files.
#'
#' @aliases BD.actualFileNames BD.nBlocks BD.blockLengths BD.getMetaData
#' BD.getData BD.checkAndDeleteFiles
#' @param bwData A \code{BlockwiseData} object.
#' @param blocks Optional vector of integers specifying the blocks on which to
#' execute the operation.
#' @param simplify Logical: if the \code{blocks} argument above is of length 1,
#' should the returned list be simplified by removing the redundant outer
#' \code{list} structure?
#' @return \item{BD.actualFileNames}{returns a vector of character strings
#' giving the file names in which the files are saved, or \code{NULL} if the
#' data are held in-memory.}
#'
#' \item{BD.nBlocks}{returns the number of blocks in the input object.}
#'
#' \item{BD.blockLengths}{returns the block lengths (results of applying
#' \code{\link{length}} to the data in each block).}
#'
#' \item{BD.getMetaData}{returns a list with one component per block. Each
#' component is in turn a list containing the stored meta-data for the
#' corresponding block. If \code{blocks} is of length 1 and \code{simplify} is
#' \code{TRUE}, the outer (redundant) \code{list} is removed.}
#'
#' \item{BD.getData}{returns a list with one component per block. Each
#' component is in turn a list containing the stored data for the corresponding
#' block. If \code{blocks} is of length 1 and \code{simplify} is \code{TRUE},
#' the outer (redundant) \code{list} is removed.}
#'
#' \item{BD.checkAndDeleteFiles}{deletes the files referenced in the input
#' \code{bwData} if they exist.}
#' @section Warning: The definition of \code{BlockwiseData} and the functions
#' here should be considered experimental and may change in the future.
#' @author Peter Langfelder
#' @seealso Definition of and other functions on \code{\link{BlockwiseData}}:
#'
#' \code{\link{newBlockwiseData}} for creating new \code{BlockwiseData}
#' objects
#'
#' \code{\link{mergeBlockwiseData}} for merging blockwise data structure
#'
#' \code{\link{addBlockToBlockwiseData}} for adding a new block to existing
#' blockwise data
#' @keywords misc
BD.getData = function(bwData, blocks = NULL, simplify = TRUE)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure.")

  if (is.null(blocks)) blocks = 1:BD.nBlocks(bwData)
  if ( (length(blocks)==0) | any(!is.finite(blocks)))
    stop("'block' must be present and finite.")

  if (any(blocks<1) | (blocks > BD.nBlocks(bwData)))
    stop("All entries in 'block' must be between 1 and ", BD.nBlocks(bwData))

  if (bwData$external)
  {
    lengths = BD.blockLengths(bwData)
    out = mapply(.loadObject, bwData$fileNames[blocks], name = 'object', size = lengths[blocks],
                 SIMPLIFY = FALSE)
  } else
    out = bwData$data[blocks]
  if (length(blocks)==1 && simplify)
    out= out[[1]]
  out
}

BD.checkAndDeleteFiles = function(bwData)
{
  if (!inherits(bwData, "BlockwiseData")) stop("'bwData' is not a blockwise data structure")
  if (bwData$external)
    .checkAndDelete(bwData$fileNames)
}

.getData = function(x, ...)
{
  if (inherits(x, "BlockwiseData")) return(BD.getData(x, ...))
  x
}

.setAttr = function(object, name, value)
{
  attr(object, name) = value
  object
}

.setAttrFromList = function(object, valueList)
{
  if (length(valueList) > 0) for (i in 1:length(valueList))
      attr(object, names(valueList)[i]) = valueList[[i]]
  object
}

# A version of getElement that returns NULL if name does not name a valid object
.getElement = function(lst, name)
{
  if (name %in% names(lst)) lst[[name, exact = TRUE]] else NULL
}

.checkLogicalConsistency = function(objects, name)
{
  vals = sapply(objects, getElement, name)
  if (!all(vals) && !all(!vals))
    stop("All arguments must have the same value of '", name, "'.")
  vals[1]
}

.checkListNamesConsistency = function(lst, name)
{
  names = lapply(lst, names)
  if (!all(sapply(names, function(x) isTRUE(all.equal(x, names[[1]])))))
    stop("Not all names agree in ", name)
}

#==============================================================================================
#
# general utility functions
#
#==============================================================================================

# Try to guess whether disk cache should be used
# this should work for both multiSet as well as for simple lists of arrays.

.useDiskCache = function(multiExpr, blocks = NULL, chunkSize = NULL)
{
  mtd = is.multiSet(multiExpr)

  nSets = length(multiExpr)
  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))

  if (length(blocks) == 0) {
    blockLengths = if (mtd) checkSets(multiExpr)$nGenes else ncol(multiExpr[[1]])
  } else
    blockLengths = as.numeric(table(blocks))

  max(blockLengths) > chunkSize
}

.dimensions = function(x)
{
   if (is.null(dim(x))) return(length(x))
   return(dim(x))
}

.shiftList = function(c0, lst)
{
  if (length(lst)==0) return(list(c0))
  ll = length(lst)
  out = list()
  out[[1]] = c0
  for (i in 1:ll)
    out[[i+1]] = lst[[i]]
  out
}

.checkListDimConsistencyAndGetDimnames = function(dataList)
{
  nPars = length(dataList)
  dn = NULL
  for (p in 1:nPars)
  {
      if (mode(dataList[[p]])!="numeric")
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(dataList[[p]])
       } else {
          if (!isTRUE(all.equal(.dimensions(dataList[[p]]), dim)))
             stop("Argument dimensions are not consistent.")
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."))
       if (is.null(dn)) dn = dimnames(dataList[[p]])
   }
   dn
}

.mtd.checkDimConsistencyAndGetDimnames = function(mtd)
{
  nPars = length(mtd)
  dn = NULL
  for (p in 1:nPars)
  {
      if (mode(mtd[[p]]$data)!="numeric")
          stop(paste("Argument number", p, " is not numeric."))
       if (p==1) {
          dim = .dimensions(mtd[[p]]$data)
       } else {
          if (!isTRUE(all.equal(.dimensions(mtd[[p]]$data), dim)))
             stop("Argument dimensions are not consistent.")
       }
       if (prod(dim)==0) stop(paste("Argument has zero dimension."))
       if (is.null(dn)) dn = dimnames(mtd[[p]]$data)
   }
   dn
}



#==============================================================================================
#
# general utility functions for working with disk cache
#
#==============================================================================================

.saveChunks = function(data, chunkSize, fileBase, cacheDir, fileNames = NULL)
{
   ld = length(data)
   nChunks = ceiling(ld/chunkSize)
   if (is.null(fileNames))
   {
     if (length(fileBase)!=1)
       stop("Internal error: length of 'fileBase' must be 1.")

     fileNames = rep("", nChunks)
     x = 1
     for (c in 1:nChunks)
     {
       fileNames[c] = tempfile(pattern = fileBase, tmpdir = cacheDir)
       # This effectively reserves the file name
       save(x, file = fileNames[c])
     }
   } else {
     if (length(fileNames)!=nChunks)
       stop("Internal error: length of 'fileNames' must equal the number of chunks.")
   }

   chunkLengths = rep(0, nChunks)
   start = 1
   for (c in 1:nChunks)
   {
     end = min(start + chunkSize-1, ld)
     chunkLengths[c] = end - start + 1
     temp = data[start:end]
     save(temp, file = fileNames[c])
     start = end + 1
   }
   rm(temp)
   collectGarbage()
   list(files = fileNames, chunkLengths = chunkLengths)
}

.loadAsList = function(file) {
  env = new.env()
  load(file = file, envir = env)
  as.list(env)
}


.loadObject = function(file, name = NULL, size = NULL) {
  x = .loadAsList(file)
  if (!is.null(name) && (names(x)[1]!=name)) {
    stop("File ", file, " does not contain object '", name, "'.")
  }

  obj = x[[1]]
  if (!is.null(size) && (length(obj)!=size))
    stop("Object '", name, "' does not have the correct length.")
  obj
}

.vector2dist = function(x) {
  n = length(x)
  n1 = (1 + sqrt(1 + 8*n))/2
  if (floor(n1)!=n1) stop("Input length not consistent with a distance structure.")
  attributes(x) = list(Size = as.integer(n1), Diag = FALSE, Upper = FALSE)
  class(x) = "dist"
  x
}

.emptyDist = function(nObjects, fill = 0) {
  n = (nObjects * (nObjects-1))/2
  .vector2dist(rep(fill, n))
}

.checkAndDelete = function(files) {
  if (length(files)>0) lapply(as.character(files), function(file) if (file.exists(file)) file.remove(file))
  NULL
}

.qorder = function(data)
{
  data = as.numeric(data)
  .Call("qorder", data, PACKAGE = "WGCNA")
}

# Actual consensus calculation distilled into one function. data is assumed to have sets in columns
# and samples/observations/whatever in rows. setWeightMat should be a matrix of dimensions (nSets, 1)
# and be normalized to sum=1.

.consensusCalculation.base = function(data, useMean, setWeightMat, consensusQuantile) {
  nSets = ncol(data)
  if (nSets==1) {
    out.list = list(consensus = data)
    if (consensusQuantile==0) out.list$originCount = c(`1`=nrow(data))
  } else {
    if (useMean) {
       if (any(is.na(data))) {
         finiteMat = 1-is.na(data)
         data[is.na(data)] = 0
         out = data %*% setWeightMat / finiteMat%*%setWeightMat
       } else {
         out = data %*% setWeightMat
       }
       out.list = list(consensus = out)
    } else if (consensusQuantile == 0)  {
        #min = rep(0, nrow(data))
        #which = rep(0, nrow(data))
        #whichmin = .C("minWhichMin_row", as.double(data),
        #              as.integer(nrow(data)), as.integer(ncol(data)),
        #              as.double(min), as.double(which), PACKAGE = "WGCNA")
        #min = whichmin[[4]]
        #which = whichmin[[5]] + 1
        #rm(whichmin)
        whichmin = .Call("minWhich_call", data, 1L, PACKAGE = "WGCNA")
        out.list = list(consensus = whichmin$min, originCount = table(as.integer(whichmin$which)))
    } else {
       out.list = list(consensus = rowQuantileC(data, p = consensusQuantile))
    }
  }
  out.list
}

.consensusCalculation.base.FromList = function(dataList, useMean, setWeights, consensusQuantile) {
  nSets = length(dataList)
  if (nSets==1) {
    out.list = list(consensus = dataList[[1]])
    if (consensusQuantile == 0) out.list$originCount = c(`1`=length(dataList[[1]]))
  } else {
    if (useMean) {
       out.list = list(consensus = pmean(dataList, setWeights))
    } else if (consensusQuantile == 0) {
        whichmin = pminWhich.fromList(dataList)
        min = whichmin$min
        which = whichmin$which
        out.list = list(consensus = min, originCount = table(as.integer(which)))
    } else {
       out.list = list(consensus = pquantile.fromList(dataList, prob = consensusQuantile))
    }
  }
  out.list
}

#==============================================================================================
#
# utility functions for working with multiple blockwise adjacencies.
#
#==============================================================================================



#' Create a list holding information about dividing data into blocks
#'
#' This function creates a list storing information about dividing data into
#' blocks, as well as about possibly excluding genes or samples with excessive
#' numbers of missing data.
#'
#'
#' @aliases newBlockInformation BlockInformation
#' @param blocks A vector giving block labels. It is assumed to be a numeric
#' vector with block labels consecutive integers starting at 1.
#' @param goodSamplesAndGenes A list returned by \code{\link{goodSamplesGenes}}
#' or \code{\link{goodSamplesGenesMS}}.
#' @return A list with \code{class} attribute set to \code{BlockInformation},
#' with the following componens: \item{blocks}{A copy of the input
#' \code{blocks}.} \item{blockGenes}{A list with one component per block,
#' giving the indices of elements in \code{block} whose value is the same.}
#' \item{goodSamplesAndGenes}{A copy of input \code{goodSamplesAndGenes}.}
#' \item{nGGenes}{Number of `good' genes in \code{goodSamplesAndGenes}.}
#' \item{gBlocks}{The input \code{blocks} restricted to `good' genes in
#' \code{goodSamplesAndGenes}.}
#' @author Peter Langfelder
#' @seealso \code{\link{goodSamplesGenes}}, \code{\link{goodSamplesGenesMS}}.
#' @keywords misc
newBlockInformation = function(
    blocks,
    goodSamplesAndGenes) {
  blockGenes = tapply(1:length(blocks), blocks, identity)
  names(blockGenes) = sort(unique(blocks))
  nGGenes = sum(goodSamplesAndGenes$goodGenes)
  gBlocks = blocks[goodSamplesAndGenes$goodGenes]
  out = list(blocks = blocks,
             blockGenes = blockGenes,
             goodSamplesAndGenes = goodSamplesAndGenes,
             nGGenes = nGGenes,
             gBlocks = gBlocks)
  class(out) = c("BlockInformation", class(out))
  out
}

#=======================================================================================================
#
# individualTOMs
# This is essentially a re-badged blockwiseIndividualTOMs with a different input and ouptut format.
#
#=======================================================================================================

# The following is contained in networkOptions:

    #corType = corType,
    #maxPOutliers = maxPOutliers,
    #quickCor = quickCor,
    #pearsonFallback = pearsonFallback,
    #cosineCorrelation = cosineCorrelation,
    #corFnc = corFnc,
    #corOptions = corOptions,
    #corType.code = match(corType, .corTypes),

    # Adjacency options
    #replaceMissingAdjacencies = TRUE,
    #power = 6,
    #networkType = c("signed hybrid", "signed", "unsigned"),

    # Topological overlap options

    #TOMType = c("signed", "unsigned"),
    #TOMDenom = c("mean", "min"))




#' Calculate individual correlation network matrices
#'
#' This function calculates correlation network matrices (adjacencies or
#' topological overlaps), after optionally first pre-clustering input data into
#' blocks.
#'
#' The function starts by optionally filtering out samples that have too many
#' missing entries and genes that have either too many missing entries or zero
#' variance in at least one set. Genes that are filtered out are excluded from
#' the network calculations.
#'
#' If \code{blocks} is not given and the number of genes (columns) in
#' \code{multiExpr} exceeds \code{maxBlockSize}, genes are pre-clustered into
#' blocks using the function \code{\link{consensusProjectiveKMeans}}; otherwise
#' all genes are treated in a single block. Any missing data in
#' \code{multiExpr} will be imputed; if imputed data are already available,
#' they can be supplied separately.
#'
#' For each block of genes, the network adjacency is constructed and (if
#' requested) topological overlap is calculated in each set. The topological
#' overlaps can be saved to disk as RData files, or returned directly within
#' the return value (see below). Note that the matrices can be big and
#' returning them within the return value can quickly exhaust the system's
#' memory. In particular, if the block-wise calculation is necessary, it is
#' usually impossible to return all matrices in the return value.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param multiExpr.imputed Optional version of \code{multiExpr} with missing
#' data imputed. If not given and \code{multiExpr} contains missing data, they
#' will be imputed using the function \code{\link[impute]{impute.knn}}.
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
#' @param networkOptions A single list of class \code{\link{NetworkOptions}}
#' giving options for network calculation for all of the networks, or a
#' \code{\link{multiSet}} structure containing one such list for each input
#' data set.
#' @param saveTOMs logical: should individual TOMs be saved to disk
#' (\code{TRUE}) or retuned directly in the return value (\code{FALSE})?
#' @param individualTOMFileNames character string giving the file names to save
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param collectGarbage Logical: should garbage collection be called after
#' each block calculation? This can be useful when the data are large, but
#' could unnecessarily slow down calculation with small data.
#' @param verbose Integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent Indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{blockwiseAdjacencies}{A
#' \code{\link{multiSet}} structure containing (possibly blockwise) network
#' matrices for each input data set. The network matrices are stored as
#' \code{\link{BlockwiseData}} objects.} \item{setNames}{A copy of
#' \code{names(multiExpr)}.} \item{nSets}{Number of sets in \code{multiExpr}}
#' \item{blockInfo}{A list of class \code{\link{BlockInformation}}, giving
#' information about blocks and gene and sample filtering.}
#' \item{networkOptions}{The input \code{networkOptions}, returned as a
#' \code{\link{multiSet}} structure with one entry per input data set.}
#' @author Peter Langfelder
#' @seealso Input arguments and output components of this function use
#' \code{\link{multiSet}}, \code{\link{NetworkOptions}},
#' \code{\link{BlockwiseData}}, and \code{\link{BlockInformation}}.
#'
#' Underlying functions of interest include
#' \code{\link{consensusProjectiveKMeans}},
#' \code{\link{TOMsimilarityFromExpr}}.
#' @keywords misc
individualTOMs = function(
   multiExpr,

   multiExpr.imputed = NULL,  ## Optional, useful for pre-clustering if preclustering is needed.

   # Data checking options

   checkMissingData = TRUE,

   # Blocking options

   blocks = NULL,
   maxBlockSize = 5000,
   blockSizePenaltyPower = 5,
   nPreclusteringCenters = NULL,
   randomSeed = 12345,

   # Network construction options. This can be a single object of class NetworkOptions, or a multiSet
   # structure of NetworkOptions objects, one per element of multiExpr.

   networkOptions,

   # Save individual TOMs? This is equivalent to using external = TRUE in blockwiseData

   saveTOMs = TRUE,
   individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

   # Behaviour options
   collectGarbage = TRUE,
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

  if (inherits(networkOptions, "NetworkOptions")) {
    networkOptions = list2multiSet(.listRep(networkOptions, nSets))
  } else {
    if (!all(apply(networkOptions, inherits, "NetworkOptions", mdaSimplify = TRUE)))
       stop("'networkOptions' must be of class 'NetworkOptions' or a multiSet structure\n",
            "   of objects of the class.\n",
            "   See newNetworkOptions for creating valid network options.")
  }

  if (is.null(names(multiExpr)))
    names(multiExpr) = paste0("Set", 1:nSets)

  if (!is.null(randomSeed)) {
    if (exists(".Random.seed")) {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed)
    }
    set.seed(randomSeed)
  }

  #if (maxBlockSize >= floor(sqrt(2^31)) )
  #  stop("'maxBlockSize must be less than ", floor(sqrt(2^31)), ". Please decrease it and try again.")

  if (!is.null(blocks) && (length(blocks)!=nGenes))
    stop("Input error: length of 'blocks' must equal number of genes in 'multiExpr'.")

  if (verbose>0)
     printFlush(paste(spaces, "Calculating topological overlaps block-wise from all genes"))

  nSamples = dataSize$nSamples

  # Check data for genes and samples that have too many missing values

  # Check that multiExpr has valid (mtd)column names. If column names are missing, generate them.

  colIDs = colnames(multiExpr)
  if (is.null(colIDs)) colIDs = c(1:dataSize$nGenes)

  if (checkMissingData) {
    gsg = goodSamplesGenesMS(multiExpr, verbose = verbose - 1, indent = indent + 1)
    if (!gsg$allOK) {
      multiExpr = subset(multiExpr, gsg$goodSamples, gsg$goodGenes)
      if (!is.null(multiExpr.imputed))
        multiExpr.imputed = subset(multiExpr.imputed, gsg$goodSamples, gsg$goodGenes)
      colIDs = colIDs[gsg$goodGenes]
    }
  } else {
    gsg = list(goodGenes = rep(TRUE, nGenes), goodSamples = lapply(nSamples, function(n) rep(TRUE, n)))
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
      clustering = consensusProjectiveKMeans(
                       if (!is.null(multiExpr.imputed)) multiExpr.imputed else multiExpr,
                       preferredSize = maxBlockSize,
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

  if (any(blockSizes > sqrt(2^31)-1))
    printFlush(paste0(spaces,
            "Found block(s) with size(s) larger than limit of 'int' indexing.\n",
            spaces, " Support for such large blocks is experimental; please report\n",
            spaces, " any issues to Peter.Langfelder@gmail.com."))

  # check file names for uniqueness

  actualFileNames = matrix("", nSets, nBlocks)
  if (saveTOMs)
  {
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
  setTomDS = replicate(nSets, list())
  # Here's where the analysis starts
  for (blockNo in 1:nBlocks)
  {
    if (verbose>1 && nBlocks > 1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select the block genes
    block = c(1:nGGenes)[gBlocks==blockLevels[blockNo]]
    #nBlockGenes = length(block)
    #blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]]
    #errorOccurred = FALSE

    # For each set: calculate and save TOM

    for (set in 1:nSets)
    {
      if (verbose>2) printFlush(paste(spaces, "....Working on set", set))
      selExpr = as.matrix(multiExpr[[set]]$data[, block])

      tomDS = as.dist(.networkCalculation(selExpr, networkOptions[[set]]$data,
                      verbose = verbose -2, indent = indent+2))

      setTomDS[[set]]$data = addBlockToBlockwiseData(if (blockNo==1) NULL else setTomDS[[set]]$data,
                            external = saveTOMs,
                            blockData = tomDS,
                            blockFile = actualFileNames[set, blockNo],
                            recordAttributes = TRUE,
                            metaData = list(IDs = colIDs[block]))
    }
    if (collectGarbage) { rm(tomDS); gc(); }
  }

  names(setTomDS) = names(multiExpr)

  if (!multiFormat)
  {
    gsg$goodSamples = gsg$goodSamples[[1]]
    setTomDS = setTomDS[[1]]$data
  }

  blockInformation = newBlockInformation(blocks, gsg)

  list(blockwiseAdjacencies = setTomDS,
       setNames = names(multiExpr),
       nSets = length(multiExpr),
       blockInfo = blockInformation,
       networkOptions = networkOptions
       )
}


#=====================================================================================================
#
# hierarchical consensus TOM
#
#=====================================================================================================



#' Calculation of hierarchical consensus topological overlap matrix
#'
#' This function calculates consensus topological overlap in a hierarchical
#' manner.
#'
#' This function is essentially a wrapper for
#' \code{\link{hierarchicalConsensusCalculation}}, with a few additional
#' operations specific to calculations of topological overlaps.
#'
#' @param multiExpr Expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
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
#' @param keepIndividualTOMs Logical: should individual TOMs be retained after
#' the calculation is finished?
#' @param individualTOMFileNames Character string giving the file names to save
#' individual TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param individualTOMInfo A list, typically returned by
#' \code{\link{individualTOMs}}, containing information about the topological
#' overlap matrices in the individual data sets in \code{multiExpr}. See the
#' output of \code{\link{individualTOMs}} for details on the content of the
#' list.
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param useBlocks Optional vector giving the blocks that should be used for
#' the calcualtions. If \code{NULL}, all all blocks will be used.
#' @param saveCalibratedIndividualTOMs Logical: should the calibrated
#' individual TOMs be saved?
#' @param calibratedIndividualTOMFilePattern Specification of file names in
#' which calibrated individual TOMs should be saved.
#' @param saveConsensusTOM Logical: should the consensus TOM be saved to disk?
#' @param consensusTOMFilePattern Character string giving the file names to
#' save consensus TOMs into. The following tags should be used to make the file
#' names unique for each set and block: \code{\%s} will be replaced by the set
#' number; \code{\%N} will be replaced by the set name (taken from
#' \code{names(multiExpr)}) if it exists, otherwise by set number; \code{\%b}
#' will be replaced by the block number. If the file names turn out to be
#' non-unique, an error will be generated.
#' @param getCalibrationSamples Logical: should the sampled values used for
#' network calibration be returned?
#' @param keepIntermediateResults Logical: should intermediate consensus TOMs
#' be saved as well?
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
#' @param chunkSize network similarities are saved in smaller chunks of size
#' \code{chunkSize}. If \code{NULL}, an appropriate chunk size will be
#' determined from an estimate of available memory. Note that if the chunk size
#' is greater than the memory required for storing intemediate results, disk
#' cache use will automatically be disabled.
#' @param cacheDir character string containing the directory into which cache
#' files should be written. The user should make sure that the filesystem has
#' enough free space to hold the cache files which can get quite large.
#' @param cacheBase character string containing the desired name for the cache
#' files. The actual file names will consists of \code{cacheBase} and a suffix
#' to make the file names unique.
#' @param collectGarbage Logical: should garbage be collected after
#' memory-intensive operations?
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list that contains the output of
#' \code{\link{hierarchicalConsensusCalculation}} and two extra components:
#' \item{individualTOMInfo}{A copy of the input \code{individualTOMInfo} if it
#' was non-\code{NULL}, or the result of \code{\link{individualTOMs}}. }
#' \item{consensusTree}{A copy of the input \code{consensusTree}.}
#' @author Peter Langfelder
#' @seealso \code{\link{hierarchicalConsensusCalculation}} for the actual
#' hierarchical consensus calculation
#'
#' \code{\link{individualTOMs}} for the calculation of individual TOMs in a
#' format suitable for consensus calculation.
#' @keywords misc
hierarchicalConsensusTOM = function(
      # Supply either ...
      # ... information needed to calculate individual TOMs

      multiExpr,

      # Data checking options
      checkMissingData = TRUE,

      # Blocking options
      blocks = NULL,
      maxBlockSize = 20000,
      blockSizePenaltyPower = 5,
      nPreclusteringCenters = NULL,
      randomSeed = 12345,

      # Network construction options. This can be a single object of class NetworkOptions, or a multiSet
      # structure of NetworkOptions objects, one per element of multiExpr.

      networkOptions,

      # Save individual TOMs?

      keepIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

      # ... or information about individual (more precisely, input) TOMs

      individualTOMInfo = NULL,

      # Consensus calculation options
      consensusTree,

      useBlocks = NULL,

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",

      # Return options
      saveConsensusTOM = TRUE,
      consensusTOMFilePattern = "consensusTOM-%a-Block%b.RData",
      getCalibrationSamples = FALSE,

      # Return the intermediate results as well?
      keepIntermediateResults = saveConsensusTOM,

      # Internal handling of TOMs
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      # Behavior
      collectGarbage = TRUE,
      verbose = 1,
      indent = 0)
{
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed)
    }
    set.seed(randomSeed)
  }

  localIndividualTOMCalculation = is.null(individualTOMInfo)
  if (is.null(individualTOMInfo))
  {
    if (missing(multiExpr)) stop("Either 'individualTOMInfo' or 'multiExpr' must be given.")
    time = system.time({individualTOMInfo = individualTOMs(multiExpr = multiExpr,
                         checkMissingData = checkMissingData,

                         blocks = blocks,
                         maxBlockSize = maxBlockSize,
                         blockSizePenaltyPower = blockSizePenaltyPower,
                         nPreclusteringCenters = nPreclusteringCenters,
                         randomSeed = NULL,

                         networkOptions = networkOptions,

                         saveTOMs = useDiskCache | keepIndividualTOMs,
                         individualTOMFileNames = individualTOMFileNames,

                         collectGarbage = collectGarbage,
                         verbose = verbose, indent = indent);})
    if (verbose > 1) { printFlush("Timimg for individual TOMs:"); print(time); }
  }
  consensus = hierarchicalConsensusCalculation(individualTOMInfo$blockwiseAdjacencies,
           consensusTree,
           level = 1,
           useBlocks = useBlocks,
           randomSeed = NULL,
           saveCalibratedIndividualData = saveCalibratedIndividualTOMs,
           calibratedIndividualDataFilePattern = calibratedIndividualTOMFilePattern,
           saveConsensusData = saveConsensusTOM,
           consensusDataFileNames = consensusTOMFilePattern,
           getCalibrationSamples= getCalibrationSamples,
           keepIntermediateResults = keepIntermediateResults,
           useDiskCache = useDiskCache,
           chunkSize = chunkSize,
           cacheDir = cacheDir,
           cacheBase = cacheBase,
           collectGarbage = collectGarbage,
           verbose = verbose,
           indent = indent)

  if (localIndividualTOMCalculation)
  {
     if (!keepIndividualTOMs)
     {
        # individual TOMs are contained in the individualTOMInfo list; remove them.
        apply(individualTOMInfo$blockwiseAdjacencies, BD.checkAndDeleteFiles)
        individualTOMInfo$blockwiseAdjacencies = NULL
     }
  }

  c( consensus,
     list(individualTOMInfo = individualTOMInfo,
          consensusTree = consensusTree)
   )
}


#========================================================================================================
#
# Merge consensusTOMInfo lists
#
#========================================================================================================
#
# Caution: at present the function does not check that the inputs are compatible and compatible with the
# supplied blocks.

.mergeConsensusTOMInformationLists = function(blockInformation, consensusTOMInfoList)
{

  blocks = blockInformation$blocks

  blockLevels = unique(blocks)
  nBlocks = length(blockLevels)

  out = consensusTOMInfoList[[1]]

  # Merging consensus information
  out$consensusData = do.call(mergeBlockwiseData, lapply(consensusTOMInfoList, getElement, "consensusData"))
  if (out$saveCalibratedIndividualData)
  {
    out$calibratedIndividualData = lapply(1:out$nSets, function(set)
           do.call(mergeBlockwiseData,
                lapply(consensusTOMInfoList, function(ct) ct$calibratedIndividualData[[set]])))
  }

  if (!is.null(out$calibrationSamples))
  {
    out$calibrationSamples = do.call(c, lapply(consensusTOMInfoList, getElement, "calibrationSamples"))
  }
  out$originCount = rowSums(do.call(cbind, lapply(consensusTOMInfoList, getElement, "originCount")))

  # Merging information in individualTOMInfo

  out$individualTOMInfo$blockwiseAdjacencies = lapply(1:out$nSets, function(set)
      list(data = do.call(mergeBlockwiseData, lapply(consensusTOMInfoList, function(ct)
                     ct$individualTOMInfo$blockwiseAdjacencies[[set]]$data))))

  out$individualTOMInfo$blockInfo = blockInformation

  out
}

#========================================================================================================
#
# consensusTOM: old, single-layer consensus.
#
#========================================================================================================



#' Consensus network (topological overlap).
#'
#' Calculation of a consensus network (topological overlap).
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
#' @param nPreclusteringCenters number of centers for pre-clustering. Larger
#' numbers typically results in better but slower pre-clustering. The default
#' is \code{as.integer(min(nGenes/20, 100*nGenes/preferredSize))} and is an
#' attempt to arrive at a reasonable number given the resources available.
#' @param randomSeed integer to be used as seed for the random number generator
#' before the function starts. If a current seed exists, it is saved and
#' restored upon exit. If \code{NULL} is given, the function will not save and
#' restore the seed.
#' @param corType character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pariwise.complete.obs} option.
#' @param maxPOutliers only used for \code{corType == "bicor"}. Specifies the
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
#' @param replaceMissingAdjacencies logical: should missing values in the
#' calculation of adjacency be replaced by 0?
#' @param power soft-thresholding power for network construction.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param checkPower logical: should basic sanity check be performed on the
#' supplied \code{power}? If you would like to experiment with unusual powers,
#' set the argument to \code{FALSE} and proceed with caution.
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
#' for later use?
#' @param individualTOMFileNames Load a previously saved TOM?
#' @param individualTOMInfo Optional data for TOM matrices in individual data
#' sets. This object is returned by the function
#' \code{\link{blockwiseIndividualTOMs}}. If not given, appropriate topological
#' overlaps will be calculated using the network contruction options below.
#' @param useIndivTOMSubset If \code{individualTOMInfo} is given, this argument
#' allows to only select a subset of the individual set networks contained in
#' \code{individualTOMInfo}. It should be a numeric vector giving the indices
#' of the individual sets to be used. Note that this argument is NOT applied to
#' \code{multiExpr}.
#' @param useBlocks optional specification of blocks that should be used for
#' the calcualtions. The default is to use all blocks.
#' @param networkCalibration network calibration method. One of "single
#' quantile", "full quantile", "none" (or a unique abbreviation of one of
#' them).
#' @param saveCalibratedIndividualTOMs logical: should the calibrated
#' individual TOMs be saved?
#' @param calibratedIndividualTOMFilePattern pattern of file names for saving
#' calibrated individual TOMs.
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
#' @param getNetworkCalibrationSamples logical: should the sampled values used
#' for network calibration be returned?
#' @param consensusQuantile quantile at which consensus is to be defined. See
#' details.
#' @param useMean logical: should the consensus be determined from a (possibly
#' weighted) mean across the data sets rather than a quantile?
#' @param setWeights Optional vector (one component per input set) of weights
#' to be used for weighted mean consensus. Only used when \code{useMean} above
#' is \code{TRUE}.
#' @param saveConsensusTOMs logical: should the consensus topological overlap
#' matrices for each block be saved and returned?
#' @param consensusTOMFileNames Load a previously saved TOM?
#' @param consensusTOMFilePattern character string containing the file
#' namefiles containing the consensus topological overlaps. The tag \code{\%b}
#' will be replaced by the block number. If the resulting file names are
#' non-unique (for example, because the user gives a file name without a
#' \code{\%b} tag), an error will be generated. These files are standard R data
#' files and can be loaded using the \code{\link{load}} function.
#' @param returnTOMs logical: should calculated consensus TOM(s) be returned?
#' @param useDiskCache should calculated network similarities in individual
#' sets be temporarilly saved to disk? Saving to disk is somewhat slower than
#' keeping all data in memory, but for large blocks and/or many sets the memory
#' footprint may be too big. See \code{chunkSize} below for additional
#' information.
#' @param useDiskCache should calculated network similarities in individual
#' sets be temporarilly saved to disk? Saving to disk is somewhat slower than
#' keeping all data in memory, but for large blocks and/or many sets the memory
#' footprint may be too big. If not given (the default), the function will
#' determine the need of caching based on the size of the data. See
#' \code{chunkSize} below for additional information.
#' @param chunkSize network similarities are saved in smaller chunks of size
#' \code{chunkSize}. If \code{NULL}, an appropriate chunk size will be
#' determined from an estimate of available memory. Note that if the chunk size
#' is greater than the memory required for storing intemediate results, disk
#' cache use will automatically be disabled.
#' @param cacheDir character string containing the directory into which cache
#' files should be written. The user should make sure that the filesystem has
#' enough free space to hold the cache files which can get quite large.
#' @param cacheBase character string containing the desired name for the cache
#' files. The actual file names will consists of \code{cacheBase} and a suffix
#' to make the file names unique.
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
#' @return List with the following components:
#'
#' \item{consensusTOM}{only present if input \code{returnTOMs} is \code{TRUE}.
#' A list containing consensus TOM for each block, stored as a distance
#' structure.}
#'
#' \item{TOMFiles}{only present if input \code{saveConsensusTOMs} is
#' \code{TRUE}. A vector of file names, one for each block, in which the TOM
#' for the corresponding block is stored. TOM is saved as a distance structure
#' to save space.}
#'
#' \item{saveConsensusTOMs}{a copy of the inputsaveConsensusTOMs.}
#'
#' \item{individualTOMInfo}{information about individual set TOMs. A copy of
#' the input \code{individualTOMInfo} if given; otherwise the result of calling
#' \code{blockwiseIndividualTOMs}. See \code{blockwiseIndividualTOMs} for
#' details.}
#'
#' Further components are retained for debugging and/or convenience.
#'
#' \item{useIndivTOMSubset}{a copy of the input \code{useIndivTOMSubset}.}
#'
#' \item{goodSamplesAndGenes}{a list containing information about which samples
#' and genes are "good" in the sense that they do not contain more than a
#' certain fraction of missing data and (for genes) have non-zero variance. See
#' \code{\link{goodSamplesGenesMS}} for details.}
#'
#' \item{nGGenes}{number of "good" genes in \code{goodSamplesGenes} above. }
#'
#' \item{nSets}{number of input sets.}
#'
#' \item{saveCalibratedIndividualTOMs}{a copy of the input
#' \code{saveCalibratedIndividualTOMs}.}
#'
#' \item{calibratedIndividualTOMFileNames}{if input
#' \code{saveCalibratedIndividualTOMs} is \code{TRUE}, this component will
#' contain the file names of calibrated individual networks. The file names are
#' arranged in a character matrix with each row corresponding to one input set
#' and each column to one block.}
#'
#' \item{networkCalibrationSamples}{if input
#' \code{getNetworkCalibrationSamples} is \code{TRUE}, a list with one
#' component per block. Each component is in turn a list with two components:
#' \code{sampleIndex} is a vector contain the indices of the TOM samples (the
#' indices refer to a flattened distance structure), and \code{TOMSamples} is a
#' matrix of TOM samples with each row corresponding to a sample in
#' \code{sampleIndex}, and each column to one input set.}
#'
#' \item{consensusQuantile}{a copy of the input \code{consensusQuantile}.}
#'
#' \item{originCount}{a vector with one component per input set.  When
#' \code{consensusQuantile} equals zero, \code{originCount} contains the number
#' of entries in the consensus TOM that come from each set (i.e., the number of
#' times the TOM in the set was the minimum). When \code{consensusQuantile} is
#' not zero or the "mean" consensus is used, this vector contains zeroes.}
#' @author Peter Langfelder
#' @seealso \code{\link{blockwiseIndividualTOMs}} for calculation of
#' topological overlaps across multiple sets.
#' @references WGCNA methodology has been described in
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17 PMID: 16646834
#'
#' The original reference for the WGCNA package is
#'
#' Langfelder P, Horvath S (2008) WGCNA: an R package for weighted correlation
#' network analysis. BMC Bioinformatics 2008, 9:559 PMID: 19114008
#'
#' For consensus modules, see
#'
#' Langfelder P, Horvath S (2007) "Eigengene networks for studying the
#' relationships between co-expression modules", BMC Systems Biology 2007, 1:54
#'
#' This function uses quantile normalization described, for example, in
#'
#' Bolstad BM1, Irizarry RA, Astrand M, Speed TP (2003) "A comparison of
#' normalization methods for high density oligonucleotide array data based on
#' variance and bias", Bioinformatics. 2003 Jan 22;19(2):1
#' @keywords misc
consensusTOM = function(
      # Supply either ...
      # ... information needed to calculate individual TOMs

      multiExpr,

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
      replaceMissingAdjacencies = FALSE,

      # Adjacency function options

      power = 6,
      networkType = "unsigned",
      checkPower = TRUE,

      # Topological overlap options

      TOMType = "unsigned",
      TOMDenom = "min",

      # Save individual TOMs?

      saveIndividualTOMs = TRUE,
      individualTOMFileNames = "individualTOM-Set%s-Block%b.RData",

      # ... or individual TOM information

      individualTOMInfo = NULL,
      useIndivTOMSubset = NULL,

   ##### Consensus calculation options

      useBlocks = NULL,

      networkCalibration = c("single quantile", "full quantile", "none"),

      # Save calibrated TOMs?
      saveCalibratedIndividualTOMs = FALSE,
      calibratedIndividualTOMFilePattern = "calibratedIndividualTOM-Set%s-Block%b.RData",

      # Simple quantile scaling options
      calibrationQuantile = 0.95,
      sampleForCalibration = TRUE, sampleForCalibrationFactor = 1000,
      getNetworkCalibrationSamples = FALSE,

      # Consensus definition
      consensusQuantile = 0,
      useMean = FALSE,
      setWeights = NULL,

      # Return options
      saveConsensusTOMs = TRUE,
      consensusTOMFilePattern = "consensusTOM-Block%b.RData",
      returnTOMs = FALSE,

      # Internal handling of TOMs
      useDiskCache = NULL, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      nThreads = 1,

      # Diagnostic messages
      verbose = 1,
      indent = 0)
{
  spaces = indentSpaces(indent)
  networkCalibration = match.arg(networkCalibration)

  seedSaved = FALSE
  if (!is.null(randomSeed))
  {
    if (exists(".Random.seed"))
    {
       savedSeed = .Random.seed
       on.exit(.Random.seed <<-savedSeed)
    }
    set.seed(randomSeed)
  }

  if (is.null(useDiskCache)) useDiskCache = .useDiskCache(multiExpr, blocks, chunkSize)

  if (any(!is.finite(setWeights))) stop("Entries of 'setWeights' must all be finite.")

  localIndividualTOMCalculation = is.null(individualTOMInfo)
  if (is.null(individualTOMInfo))
  {
    if (missing(multiExpr)) stop("Either 'individualTOMInfo' or 'multiExpr' must be given.")

    dataSize = checkSets(multiExpr)
    nSets.all = dataSize$nSets
    nGenes = dataSize$nGenes

    if (length(power)!=1)
    {
      if (length(power)!=nSets.all)
        stop("Invalid arguments: Length of 'power' must equal number of sets given in 'multiExpr'.")
    } else {
      power = rep(power, nSets.all)
    }

    if ( (consensusQuantile < 0) | (consensusQuantile > 1) )
      stop("'consensusQuantile' must be between 0 and 1.")

    time = system.time({individualTOMInfo = blockwiseIndividualTOMs(multiExpr = multiExpr,
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
                         replaceMissingAdjacencies = replaceMissingAdjacencies,
                         power = power,
                         networkType = networkType,
                         TOMType = TOMType,
                         TOMDenom = TOMDenom,
                         saveTOMs = useDiskCache | saveIndividualTOMs,
                         individualTOMFileNames = individualTOMFileNames,
                         nThreads = nThreads,
                         verbose = verbose, indent = indent);})
    if (verbose > 1) { printFlush("Timimg for individual TOMs:"); print(time); }

    if (!saveIndividualTOMs & useDiskCache)
       on.exit(.checkAndDelete(individualTOMInfo$actualTOMFileNames), add = TRUE)

  } else {
    nSets.all = if (individualTOMInfo$saveTOMs)
             nrow(individualTOMInfo$actualTOMFileNames) else ncol(individualTOMInfo$TOMSimilarities[[1]])
    nGenes = length(individualTOMInfo$blocks)
  }
  nGoodGenes = length(individualTOMInfo$gBlocks)

  if (is.null(setWeights)) setWeights = rep(1, nSets.all)
  if (length(setWeights)!=nSets.all)
    stop("Length of 'setWeights' must equal the number of sets.")

  setWeightMat = as.matrix(setWeights/sum(setWeights))

  if (is.null(useIndivTOMSubset))
  {
    if (individualTOMInfo$nSets != nSets.all)
      stop(paste("Number of sets in individualTOMInfo and in multiExpr do not agree.\n",
                 "  To use a subset of individualTOMInfo, set useIndivTOMSubset appropriately."))

    useIndivTOMSubset = c(1:nSets.all)
  }

  nSets = length(useIndivTOMSubset)


  if (length(unique(useIndivTOMSubset))!=nSets)
    stop("Entries of 'useIndivTOMSubset' must be unique")

  if (any(useIndivTOMSubset<1) | any(useIndivTOMSubset>individualTOMInfo$nSets))
    stop("All entries of 'useIndivTOMSubset' must be between 1 and the number of sets in individualTOMInfo")

  # if ( (minKMEtoJoin >1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  gsg = individualTOMInfo$goodSamplesAndGenes

  # Restrict gsg to used sets

  gsg$goodSamples = gsg$goodSamples[useIndivTOMSubset]

  if (is.null(chunkSize)) chunkSize = as.integer(.largestBlockSize/(2*nSets))

  # Initialize various variables

  if (getNetworkCalibrationSamples)
  {
    if (!sampleForCalibration)
      stop(paste("Incompatible input options: networkCalibrationSamples can only be returned",
                 "if sampleForCalibration is TRUE."))
    networkCalibrationSamples = list()
  }

  blockLevels = sort(unique(individualTOMInfo$gBlocks))
  nBlocks = length(blockLevels)

  if (is.null(useBlocks)) useBlocks = blockLevels

  useBlockIndex = match(useBlocks, blockLevels)

  if (!all(useBlocks %in% blockLevels))
    stop("All entries of 'useBlocks' must be valid block levels.")

  if (any(duplicated(useBlocks)))
    stop("Entries of 'useBlocks' must be unique.")

  nUseBlocks = length(useBlocks)
  if (nUseBlocks==0)
    stop("'useBlocks' cannot be non-NULL and empty at the same time.")

  consensusTOM.out = list()

  TOMFiles = rep("", nUseBlocks)
  originCount = rep(0, nSets)

  calibratedIndividualTOMFileNames = NULL
  if (saveCalibratedIndividualTOMs)
  {
    calibratedIndividualTOMFileNames = matrix("", nSets, nBlocks)
    for (set in 1:nSets) for (b in 1:nBlocks)
      calibratedIndividualTOMFileNames[set, b] = .processFileName(calibratedIndividualTOMFilePattern,
                                 setNumber = set, blockNumber = b)

  }
  gc()

  # Here's where the analysis starts

  for (blockIndex in 1:nUseBlocks)
  {
    blockNo = useBlockIndex[blockIndex]

    if (verbose>1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select block genes
    block = c(1:nGoodGenes)[individualTOMInfo$gBlocks==blockLevels[blockNo]]
    nBlockGenes = length(block)
    # blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks==blockLevels[blockNo]]
    scaleQuant = rep(1, nSets)
    scalePowers = rep(1, nSets)

    # Set up file names or memory space to hold the set TOMs
    if (useDiskCache)
    {
      nChunks = ceiling(nBlockGenes * (nBlockGenes-1)/2/chunkSize)
      chunkFileNames = array("", dim = c(nChunks, nSets))
      on.exit(.checkAndDelete(chunkFileNames), add = TRUE)
    } else nChunks = 1

    if (nChunks==1) useDiskCache = FALSE
    if (!useDiskCache)
    {
      # Note: setTomDS will contained the scaled set TOM matrices.
      setTomDS = array(0, dim = c(nBlockGenes*(nBlockGenes-1)/2, nSets))
    }

    # create an empty consTomDS distance structure.

    consTomDS = .emptyDist(nBlockGenes)

    # sample entry indices from the distance structure for TOM scaling, if requested

    if (networkCalibration=="single quantile" && sampleForCalibration)
    {
      qx = min(calibrationQuantile, 1-calibrationQuantile)
      nScGenes = min(sampleForCalibrationFactor * 1/qx, length(consTomDS))
      nTOMEntries = length(consTomDS)
      scaleSample = sample(nTOMEntries, nScGenes)
      if (getNetworkCalibrationSamples)
        networkCalibrationSamples[[blockIndex]] = list(sampleIndex = scaleSample,
                                            TOMSamples = matrix(NA, nScGenes, nSets))
    }
    if (networkCalibration %in% c("single quantile", "none"))
    {
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste(spaces, "....Working on set", useIndivTOMSubset[set]))
        if (individualTOMInfo$saveTOMs)
        {
           tomDS = .loadObject(individualTOMInfo$ actualTOMFileNames[useIndivTOMSubset[set], blockNo],
                               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2)
        } else {
          tomDS = consTomDS
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]] [, useIndivTOMSubset[set]]
        }

        if (networkCalibration=="single quantile")
        {
          # Scale TOMs so that calibrationQuantile agree in each set
          if (sampleForCalibration)
          {
            if (getNetworkCalibrationSamples)
            {
              networkCalibrationSamples[[blockIndex]]$TOMSamples[, set] = tomDS[scaleSample]
              scaleQuant[set] = quantile(networkCalibrationSamples[[blockIndex]]$TOMSamples[, set],
                                         probs = calibrationQuantile, type = 8)
            } else {
              scaleQuant[set] = quantile(tomDS[scaleSample], probs = calibrationQuantile, type = 8)
            }
          } else
            scaleQuant[set] = quantile(x = tomDS, probs = calibrationQuantile, type = 8)
          if (set>1)
          {
             scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
             tomDS = tomDS^scalePowers[set]
          }
          if (saveCalibratedIndividualTOMs)
             save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
        }

        # Save the calculated TOM either to disk in chunks or to memory.

        if (useDiskCache)
        {
          if (verbose > 3) printFlush(paste(spaces, "......saving TOM similarity to disk cache.."))
          sc = .saveChunks(tomDS, chunkSize, cacheBase, cacheDir = cacheDir)
          chunkFileNames[, set] = sc$files
          chunkLengths = sc$chunkLengths
        } else {
          setTomDS[, set] = tomDS[]
        }
        rm(tomDS)
        gc()
      }
    } else if (networkCalibration=="full quantile")
    {
      # Step 1: load each TOM, get order, split TOM into chunks according to order, and save.
      if (verbose>1) printFlush(paste0(spaces, "..working on quantile normalization"))
      if (useDiskCache)
      {
        orderFiles = rep("", nSets)
        on.exit(.checkAndDelete(orderFiles),add = TRUE)
      }
      for (set in 1:nSets)
      {
        if (verbose>2) printFlush(paste(spaces, "....Working on set", useIndivTOMSubset[set]))
        if (individualTOMInfo$saveTOMs)
        {
           tomDS = .loadObject(individualTOMInfo$ actualTOMFileNames[useIndivTOMSubset[set], blockNo],
                               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2)
        } else {
          tomDS = consTomDS
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]] [, useIndivTOMSubset[set]]
        }
        if (useDiskCache)
        {
          # Order TOM (this may take a long time...)
          if (verbose > 3) printFlush(paste0(spaces, "......ordering TOM"))
          time = system.time({order1 = .qorder(tomDS)})
          if (verbose > 3) { printFlush("Time to order TOM:"); print(time); }
          # save the order
          orderFiles[set] = tempfile(pattern = paste0(".orderForSet", set), tmpdir = cacheDir)
          if (verbose > 3) printFlush(paste0(spaces, "......saving order and ordered TOM"))
          save(order1, file = orderFiles[set])
          # Save ordered tomDS into chunks
          tomDS.ordered = tomDS[order1]
          sc = .saveChunks(tomDS.ordered, chunkSize, cacheBase, cacheDir = cacheDir)
          chunkFileNames[, set] = sc$files
          chunkLengths = sc$chunkLengths
        } else {
          setTomDS[, set] = tomDS[]
        }
      }
      if (useDiskCache)
      {
        # Step 2: Load chunks one by one and quantile normalize
        if (verbose > 2) printFlush(paste0(spaces, "....quantile normalizing chunks"))
        for (c in 1:nChunks)
        {
          if (verbose > 3) printFlush(paste0(spaces, "......QN for chunk ", c, " of ", nChunks))
          chunkData = matrix(NA, chunkLengths[c], nSets)
          for (set in 1:nSets)
            chunkData[, set] = .loadObject(chunkFileNames[c, set])

          time = system.time({ chunk.norm = normalize.quantiles(chunkData, copy = FALSE);})
          if (verbose > 1) { printFlush("Time to QN chunk:"); print(time); }
          # Save quantile normalized chunks
          for (set in 1:nSets)
          {
            temp = chunk.norm[, set]
            save(temp, file = chunkFileNames[c, set])
          }
        }

        if (verbose > 2) printFlush(paste0(spaces, "....putting together full QN'ed TOMs"))
        # Put together full TOMs
        for (set in 1:nSets)
        {
           load(orderFiles[set])
           start = 1
           for (c in 1:nChunks)
           {
             end = start + chunkLengths[c] - 1
             tomDS[order1[start:end]] = .loadObject(chunkFileNames[c, set], size = chunkLengths[c])
             start = start + chunkLengths[c]
           }
           if (saveCalibratedIndividualTOMs)
              save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
           .saveChunks(tomDS, chunkSize, fileNames = chunkFileNames[, set])
           unlink(orderFiles[set])
        }
      } else {
        # If disk cache is not being used, simply call normalize.quantiles on the full set.
        setTomDS = normalize.quantiles(setTomDS)
        if (saveCalibratedIndividualTOMs) for (set in 1:nSets)
        {
           tomDS = .vector2dist(setTomDS[, set])
           save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
        }
      }
    } else stop("Unrecognized value of 'networkCalibration': ", networkCalibration)

    # Calculate consensus network
    if (verbose > 2)
      printFlush(paste(spaces, "....Calculating consensus network"))
    if (useDiskCache)
    {
      start = 1
      for (chunk in 1:nChunks)
      {
        if (verbose > 3) printFlush(paste(spaces, "......working on chunk", chunk))
        end = start + chunkLengths[chunk] - 1
        setChunks = array(0, dim = c(chunkLengths[chunk], nSets))
        for (set in 1:nSets)
        {
          load(file = chunkFileNames[chunk, set])
          setChunks[, set] = temp
          file.remove(chunkFileNames[chunk, set])
        }
        if (useMean | consensusQuantile > 0)
        {
          consTomDS[start:end] = .consensusCalculation.base(
                     setChunks, useMean = useMean, setWeightMat = setWeightMat,
                     consensusQuantile = consensusQuantile)$consensus
        } else {
          tmp = .consensusCalculation.base(setChunks, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)
          consTomDS[start:end] = tmp$consensus
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
          rm(tmp)
        }
        start = end + 1
      }
    } else {
      if (useMean | consensusQuantile > 0)
      {
         consTomDS[] = .consensusCalculation.base(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus
      } else {
          tmp = .consensusCalculation.base(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)
          consTomDS[] = tmp$consensus
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
          rm(tmp)
      }
    }

    # Save the consensus TOM if requested

    if (saveConsensusTOMs)
    {
       TOMFiles[blockIndex] = .substituteTags(consensusTOMFilePattern, "%b", blockNo)
       if (TOMFiles[blockIndex]==consensusTOMFilePattern)
         stop(paste("File name for consensus TOM must contain the tag %b somewhere in the file name -\n",
                    "   - this tag will be replaced by the block number. "))
       save(consTomDS, file = TOMFiles[blockIndex])
    }

    if (returnTOMs) consensusTOM.out[[blockIndex]] = consTomDS

    gc()
  }

  if (!saveConsensusTOMs) TOMFiles = NULL
  if (!returnTOMs)  consensusTOM.out = NULL

  if (localIndividualTOMCalculation)
  {
     if (!individualTOMInfo$saveTOMs)
     {
        # individual TOMs are contained in the individualTOMInfo list; remove them.
        individualTOMInfo$TOMSimilarities = NULL
     }
  }

  list(consensusTOM = consensusTOM.out,
       TOMFiles = TOMFiles,
       saveConsensusTOMs = saveConsensusTOMs,

       individualTOMInfo = individualTOMInfo,
       useIndivTOMSubset = useIndivTOMSubset,
       goodSamplesAndGenes = gsg,
       nGGenes = nGoodGenes,
       nSets = nSets,

       saveCalibratedIndividualTOMs = saveCalibratedIndividualTOMs,
       calibratedIndividualTOMFileNames = calibratedIndividualTOMFileNames,
       networkCalibrationSamples = if (getNetworkCalibrationSamples) networkCalibrationSamples else NULL,

       consensusQuantile = consensusQuantile,
       originCount = originCount

      )
}
