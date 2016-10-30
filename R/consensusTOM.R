.saveChunks <- function(data, chunkSize, fileBase, cacheDir, fileNames = NULL) {
   ld = length(data)
   nChunks = ceiling(ld/chunkSize)
   if (is.null(fileNames)) {
     if (length(fileBase) != 1) {
       stop("Internal error: length of 'fileBase' must be 1.")
     }

     fileNames = rep("", nChunks)
     x = 1
     for (c in 1:nChunks) {
       fileNames[c] = tempfile(pattern = fileBase, tmpdir = cacheDir)
       # This effectively reserves the file name
       save(x, file = fileNames[c])
     }
   } else {
     if (length(fileNames) != nChunks) {
       stop("Internal error: length of 'fileNames' must equal the number of chunks.")
     }
   }

   chunkLengths <- rep(0, nChunks)
   start <- 1
   for (c in 1:nChunks) {
     end <- min(start + chunkSize-1, ld)
     chunkLengths[c] <- end - start + 1
     temp = data[start:end]
     save(temp, file = fileNames[c])
     start <- end + 1
   }
   rm(temp)
   collectGarbage()
   list(files = fileNames, chunkLengths = chunkLengths)
}

.loadObject <- function(file, name = NULL, size = NULL) {
  x <- load(file)
  if (!is.null(name) && (x != name)) {
    stop("File ", file, " does not contain object '", name, "'.")
  }

  obj <- get(x)
  if (!is.null(size) && (length(obj) != size)) {
    stop("Object '", name, "' does not have the correct length.")
  }
  obj
}

.qorder <- function(data) {
  data <- as.numeric(data)
  .Call("qorder", data)
}

# Actual consensus calculation distilled into one function. setTomMat is assumed to have sets in columns
# and gene pairs in rows. setWeightMat should be a matrix of dimensions (nSets, 1) and be normalized to sum=1.

.consensusCalculation <- function(setTomMat, useMean, setWeightMat,
                                 consensusQuantile) {
  if (useMean) {
     if (anyNA(setTomMat)) {
       finiteMat <- 1 - is.na(setTomMat)
       setTomMat[is.na(setTomMat)] = 0
       out <- setTomMat %*% setWeightMat / finiteMat %*% setWeightMat
     } else {
       out <- setTomMat %*% setWeightMat
     }
     out.list <- list(consensus = out)
  } else if (consensusQuantile == 0) {
      min <- rep(0, nrow(setTomMat))
      which <-  rep(0, nrow(setTomMat))
      whichmin<- .C("minWhichMin_row", as.double(setTomMat),
                    as.integer(nrow(setTomMat)), as.integer(ncol(setTomMat)),
                    as.double(min), as.double(which), PACKAGE = "WGCNA")
      min <- whichmin[[4]]
      which <- whichmin[[5]] + 1
      rm(whichmin)
      out.list <- list(consensus = min, originCount = table(as.integer(which)))
  } else {
     out.list <- list(consensus = rowQuantileC(setTomMat, p = consensusQuantile))
  }
  out.list
}

.vector2dist <- function(x) {
  n <- length(x)
  n1 <- (1 + sqrt(1 + 8*n))/2
  if (floor(n1) != n1) {
      stop("Input length not consistent with a distance structure.")
  }
  attributes(x) <- list(Size = as.integer(n1), Diag = FALSE, Upper = FALSE)
  class(x) <- "dist"
  x
}

.emptyDist <- function(nObjects, fill = 0) {
  n <- (nObjects * (nObjects-1))/2
  .vector2dist(rep(fill, n))
}

.checkAndDelete <- function(files) {
  if (length(files) > 0) {
      lapply(as.character(files), function(file) {
          if (file.exists(file)) {
              file.remove(file)
          }
      })
  }
    NULL
}



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
#' for later use?
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
#' @param returnTOMs logical: should calculated consensus TOM(s) be returned?
#' @param useDiskCache should calculated network similarities in individual
#' sets be temporarilly saved to disk? Saving to disk is somewhat slower than
#' keeping all data in memory, but for large blocks and/or many sets the memory
#' footprint may be too big. See \code{chunkSize} below for additional
#' information.
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
#' @param individualTOMFileNames Load a previously saved TOM?
#' @param consensusTOMFileNames Load a previously saved TOM?
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
consensusTOM <- function(
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
      consensusTOMFileNames = "consensusTOM-Block%b.RData",
      returnTOMs = FALSE,

      # Internal handling of TOMs
      useDiskCache = TRUE, chunkSize = NULL,
      cacheDir = ".",
      cacheBase = ".blockConsModsCache",

      nThreads = 1,

      # Diagnostic messages
      verbose = 1,
      indent = 0) {
  spaces = indentSpaces(indent)
  networkCalibration = match.arg(networkCalibration)

  seedSaved = FALSE
  if (!is.null(randomSeed)) {
    if (exists(".Random.seed")) {
       seedSaved = TRUE
       savedSeed = .Random.seed
    }
    set.seed(randomSeed)
  }

  if (any(!is.finite(setWeights))) {
      stop("Entries of 'setWeights' must all be finite.")
  }

  localIndividualTOMCalculation = is.null(individualTOMInfo)
  if (is.null(individualTOMInfo)){
    if (missing(multiExpr)) {
        stop("Either 'individualTOMInfo' or 'multiExpr' must be given.")
    }

    dataSize = checkSets(multiExpr)
    nSets.all = dataSize$nSets
    nGenes = dataSize$nGenes

    if (length(power) != 1) {
      if (length(power) != nSets.all)
        stop("Invalid arguments: Length of 'power' must equal number ",
             "of sets given in 'multiExpr'.")
    } else {
      power = rep(power, nSets.all)
    }

    if ( (consensusQuantile < 0) | (consensusQuantile > 1)) {
      stop("'consensusQuantile' must be between 0 and 1.")
    }
    time = system.time({individualTOMInfo = blockwiseIndividualTOMs(
        multiExpr = multiExpr,
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
        verbose = verbose, indent = indent)
    })
    if (verbose > 1) {
        printFlush("Timimg for individual TOMs:")
        print(time)
    }

    if (!saveIndividualTOMs & useDiskCache) {
       on.exit(.checkAndDelete(individualTOMInfo$actualTOMFileNames),
               add = TRUE)
    }
  } else {
    nSets.all = if (individualTOMInfo$saveTOMs) {
        nrow(individualTOMInfo$actualTOMFileNames)
    } else {
        ncol(individualTOMInfo$TOMSimilarities[[1]])
    }
    nGenes = length(individualTOMInfo$blocks)
  }
  nGoodGenes = length(individualTOMInfo$gBlocks)

  if (is.null(setWeights)) {
      setWeights = rep(1, nSets.all)
  }
  if (length(setWeights) != nSets.all)
    stop("Length of 'setWeights' must equal the number of sets.")

  setWeightMat = as.matrix(setWeights/sum(setWeights))

  if (is.null(useIndivTOMSubset)) {
    if (individualTOMInfo$nSets != nSets.all) {
      stop("Number of sets in individualTOMInfo and in ",
           "multiExpr do not agree.\n To use a subset of individualTOMInfo, ",
           "set useIndivTOMSubset appropriately.")
    }
    useIndivTOMSubset = c(1:nSets.all)
  }

  nSets = length(useIndivTOMSubset)


  if (length(unique(useIndivTOMSubset)) != nSets) {
      stop("Entries of 'useIndivTOMSubset' must be unique")
  }

  if (any(useIndivTOMSubset < 1) | any(
      useIndivTOMSubset > individualTOMInfo$nSets)) {
    stop("All entries of 'useIndivTOMSubset' must be between 1 and the ",
         "number of sets in individualTOMInfo")
  }
  # if ( (minKMEtoJoin  > 1) | (minKMEtoJoin  <0) ) stop("minKMEtoJoin  must be between 0 and 1.")

  gsg = individualTOMInfo$goodSamplesAndGenes

  # Restrict gsg to used sets

  gsg$goodSamples = gsg$goodSamples[useIndivTOMSubset]

  if (is.null(chunkSize)) {
      chunkSize = as.integer(.largestBlockSize/(2*nSets))
  }

  # Initialize various variables

  if (getNetworkCalibrationSamples) {
    if (!sampleForCalibration) {
      stop("Incompatible input options: networkCalibrationSamples can only be",
           " returned if sampleForCalibration is TRUE.")
    }
    networkCalibrationSamples = list()
  }

  blockLevels = sort(unique(individualTOMInfo$gBlocks))
  nBlocks = length(blockLevels)

  if (is.null(useBlocks)) {
      useBlocks = blockLevels
  }

  useBlockIndex = match(useBlocks, blockLevels)

  if (!all(useBlocks %in% blockLevels)) {
    stop("All entries of 'useBlocks' must be valid block levels.")
  }

  if (any(duplicated(useBlocks))) {
    stop("Entries of 'useBlocks' must be unique.")
  }

  nUseBlocks = length(useBlocks)
  if (nUseBlocks == 0) {
    stop("'useBlocks' cannot be non-NULL and empty at the same time.")
  }

  consensusTOM.out = list()

  TOMFiles = rep("", nUseBlocks)
  originCount = rep(0, nSets)

  calibratedIndividualTOMFileNames = NULL
  if (saveCalibratedIndividualTOMs) {
    calibratedIndividualTOMFileNames = matrix("", nSets, nBlocks)
    for (set in 1:nSets) {
        for (b in 1:nBlocks) {
            calibratedIndividualTOMFileNames[set, b] = .processFileName(
                calibratedIndividualTOMFilePattern, set,
                individualTOMInfo$setNames, b)
        }
    }
  }
  collectGarbage()

  # Here's where the analysis starts

  for (blockIndex in 1:nUseBlocks) {
    blockNo = useBlockIndex[blockIndex]

    if (verbose > 1) printFlush(paste(spaces, "..Working on block", blockNo, "."))
    # Select block genes
    block = c(1:nGoodGenes)[individualTOMInfo$gBlocks == blockLevels[blockNo]]
    nBlockGenes = length(block)
    # blockGenes[[blockNo]] = c(1:nGenes)[gsg$goodGenes][gBlocks == blockLevels[blockNo]]
    scaleQuant = rep(1, nSets)
    scalePowers = rep(1, nSets)

    # Set up file names or memory space to hold the set TOMs
    if (useDiskCache) {
      nChunks = ceiling(nBlockGenes * (nBlockGenes-1)/2/chunkSize)
      chunkFileNames = array("", dim = c(nChunks, nSets))
      on.exit(.checkAndDelete(chunkFileNames), add = TRUE)
    } else {
        nChunks = 1
    }

    if (nChunks == 1) {
        useDiskCache = FALSE
    }
    if (!useDiskCache) {
      # Note: setTomDS will contained the scaled set TOM matrices.
      setTomDS = array(0, dim = c(nBlockGenes*(nBlockGenes-1)/2, nSets))
    }

    # create an empty consTomDS distance structure.

    consTomDS = .emptyDist(nBlockGenes)

    # sample entry indices from the distance structure for TOM scaling, if requested

    if (networkCalibration == "single quantile" && sampleForCalibration) {
      qx = min(calibrationQuantile, 1-calibrationQuantile)
      nScGenes = min(sampleForCalibrationFactor * 1/qx, length(consTomDS))
      nTOMEntries = length(consTomDS)
      scaleSample = sample(nTOMEntries, nScGenes)
      if (getNetworkCalibrationSamples) {
        networkCalibrationSamples[[blockIndex]] = list(
            sampleIndex = scaleSample,
            TOMSamples = matrix(NA, nScGenes, nSets)
        )
      }
    }
    if (networkCalibration %in% c("single quantile", "none")) {
      for (set in 1:nSets){
        if (verbose > 2) {
            printFlush(paste(spaces, "....Working on set",
                             useIndivTOMSubset[set]))
        }
        if (individualTOMInfo$saveTOMs) {
           tomDS = .loadObject(
               individualTOMInfo$actualTOMFileNames[useIndivTOMSubset[set], blockNo],
               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2)
        } else {
          tomDS = consTomDS
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]][, useIndivTOMSubset[set]]
        }

        if (networkCalibration == "single quantile") {
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
          if (set > 1){
             scalePowers[set] = log(scaleQuant[1])/log(scaleQuant[set])
             tomDS = tomDS^scalePowers[set]
          }
          if (saveCalibratedIndividualTOMs) {
             save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
          }
        }

        # Save the calculated TOM either to disk in chunks or to memory.

        if (useDiskCache) {
          if (verbose > 3) {
              printFlush(paste(spaces,
                               "......saving TOM similarity to disk cache.."))
          }
          sc = .saveChunks(tomDS, chunkSize, cacheBase, cacheDir = cacheDir)
          chunkFileNames[, set] = sc$files
          chunkLengths = sc$chunkLengths
        } else {
          setTomDS[, set] = tomDS[]
        }
        rm(tomDS)
        collectGarbage()
      }
    } else if (networkCalibration == "full quantile") {
      # Step 1: load each TOM, get order, split TOM into chunks according to order, and save.
      if (verbose > 1) {
          printFlush(paste0(spaces, "..working on quantile normalization"))
      }
      if (useDiskCache) {
        orderFiles = rep("", nSets)
        on.exit(.checkAndDelete(orderFiles),add = TRUE)
      }
      for (set in 1:nSets) {
        if (verbose > 2) {
            printFlush(paste(spaces, "....Working on set",
                             useIndivTOMSubset[set]))
        }
          if (individualTOMInfo$saveTOMs) {
           tomDS = .loadObject(individualTOMInfo$ actualTOMFileNames[useIndivTOMSubset[set], blockNo],
                               name = "tomDS", size = nBlockGenes*(nBlockGenes-1)/2)
        } else {
          tomDS = consTomDS
          tomDS[] = individualTOMInfo$TOMSimilarities[[blockNo]] [, useIndivTOMSubset[set]]
        }
        if (useDiskCache) {
          # Order TOM (this may take a long time...)
          if (verbose > 3) {
              printFlush(paste0(spaces, "......ordering TOM"))
          }
          time = system.time({order1 = .qorder(tomDS)})
          if (verbose > 1) {
              printFlush("Time to order TOM:")
              print(time)
          }
          # save the order
          orderFiles[set] = tempfile(pattern = paste0("orderForSet", set),
                                     tmpdir = cacheDir)
          if (verbose > 3) {
              printFlush(paste0(spaces, "......saving order and ordered TOM"))
          }
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
      if (useDiskCache) {
        # Step 2: Load chunks one by one and quantile normalize
        if (verbose > 2) {
            printFlush(paste0(spaces, "....quantile normalizing chunks"))
        }
        for (c in 1:nChunks) {
          if (verbose > 3) {
              printFlush(paste0(spaces, "......QN for chunk ", c, " of ", nChunks))
          }
          chunkData = matrix(NA, chunkLengths[c], nSets)
          for (set in 1:nSets) {
            chunkData[, set] = .loadObject(chunkFileNames[c, set])
          }
          time = system.time({
              chunk.norm = normalize.quantiles(chunkData, copy = FALSE)})
          if (verbose > 1) {
              printFlush("Time to QN chunk:")
              print(time)
          }
          # Save quantile normalized chunks
          for (set in 1:nSets) {
            temp = chunk.norm[, set]
            save(temp, file = chunkFileNames[c, set])
          }
        }

        if (verbose > 2) {
            printFlush(paste0(spaces, "....putting together full QN'ed TOMs"))
        }
        # Put together full TOMs
        for (set in 1:nSets) {
           load(orderFiles[set])
           start = 1
           for (c in 1:nChunks) {
             end = start + chunkLengths[c] - 1
             tomDS[order1[start:end]] = .loadObject(chunkFileNames[c, set], size = chunkLengths[c])
             start = start + chunkLengths[c]
           }
           if (saveCalibratedIndividualTOMs) {
              save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
           }
           .saveChunks(tomDS, chunkSize, fileNames = chunkFileNames[, set])
           unlink(orderFiles[set])
        }
      } else {
        # If disk cache is not being used, simply call normalize.quantiles on the full set.
        setTomDS = normalize.quantiles(setTomDS)
        if (saveCalibratedIndividualTOMs) {
            for (set in 1:nSets) {
                tomDS = .vector2dist(setTomDS[, set])
                save(tomDS, file = calibratedIndividualTOMFileNames[set, blockNo])
            }
        }
      }
    } else {
        stop("Unrecognized value of 'networkCalibration': ",
             networkCalibration)
    }

    # Calculate consensus network
    if (verbose > 2) {
      printFlush(paste(spaces, "....Calculating consensus network"))
    }
    if (useDiskCache) {
      start = 1
      for (chunk in 1:nChunks) {
        if (verbose > 3) {
            printFlush(paste(spaces, "......working on chunk", chunk))
        }
        end = start + chunkLengths[chunk] - 1
        setChunks = array(0, dim = c(chunkLengths[chunk], nSets))
        for (set in 1:nSets) {
          load(file = chunkFileNames[chunk, set])
          setChunks[, set] = temp
          file.remove(chunkFileNames[chunk, set])
        }
        if (useMean | consensusQuantile > 0) {
          consTomDS[start:end] = .consensusCalculation(setChunks, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus
        } else {
          tmp = .consensusCalculation(setChunks, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)
          consTomDS[start:end] = tmp$consensus
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
          rm(tmp)
        }
        start = end + 1
      }
    } else {
      if (useMean | consensusQuantile > 0) {
         consTomDS[] = .consensusCalculation(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)$consensus
      } else {
          tmp = .consensusCalculation(setTomDS, useMean = useMean, setWeightMat = setWeightMat,
                                             consensusQuantile = consensusQuantile)
          consTomDS[] = tmp$consensus
          countIndex = as.numeric(names(tmp$originCount))
          originCount[countIndex] = originCount[countIndex] + tmp$originCount
          rm(tmp)
      }
    }

    # Save the consensus TOM if requested

    if (saveConsensusTOMs) {
       TOMFiles[blockIndex] = .substituteTags(consensusTOMFileNames, "%b", blockNo)
       if (TOMFiles[blockIndex] == consensusTOMFileNames) {
         stop("File name for consensus TOM must contain the tag %b somewhere ",
              "in the file name -\n   - this tag will be replaced by the ",
              "block number. ")
       }
       save(consTomDS, file = TOMFiles[blockIndex])
    }

    if (returnTOMs) {
        consensusTOM.out[[blockIndex]] = consTomDS
    }

    collectGarbage()
  }

  if (!saveConsensusTOMs) {
      TOMFiles = NULL
  }
  if (!returnTOMs) {
      consensusTOM.out = NULL
  }

  if (localIndividualTOMCalculation) {
     if (!individualTOMInfo$saveTOMs) {
        # individual TOMs are contained in the individualTOMInfo list; remove them.
        individualTOMInfo$TOMSimilarities = NULL
     }
  }


  if (seedSaved) {
      .Random.seed <<- savedSeed
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
       networkCalibrationSamples = ifelse(getNetworkCalibrationSamples,
                                          networkCalibrationSamples, NULL),
       consensusQuantile = consensusQuantile,
       originCount = originCount
      )
}
