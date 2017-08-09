# Call blockwise modules several times with sampled data, collect the results
# -02: first entry in the result will be the base modules.

# Also: add functions for determining module stability and whether modules should be merged.

#===================================================================================================
#
# sampledBlockwiseModules
#
#===================================================================================================



#' Blockwise module identification in sampled data
#'
#' This function repeatedly resamples the samples (rows) in supplied data and
#' identifies modules on the resampled data.
#'
#' For each run, samples (but not genes) are randomly sampled to obtain a
#' perturbed data set; a full network analysis and module identification is
#' carried out, and the results are returned in a list with one component per
#' run.
#'
#' For each run, the soft-thresholding power can optionally be adjusted such
#' that the mean adjacency in the re-sampled data set equals the mean adjacency
#' in the original data.
#'
#' @param datExpr Expression data. A matrix (preferred) or data frame in which
#' columns are genes and rows ar samples.
#' @param nRuns Number of network construction and module identification runs.
#' @param startRunIndex Number to be assigned to the start run. The run number
#' or index is used to make saved files unique; it has no effect on the actual
#' results of the run.
#' @param endRunIndex Number (index) of the last run. If given, \code{nRuns} is
#' ignored.
#' @param replace Logical: should samples (observations or rows in entries in
#' \code{multiExpr}) be sampled with replacement?
#' @param fraction Fraction of samples to sample for each run.
#' @param randomSeed Integer specifying the random seed. If non-NULL, the
#' random number generator state is saved before the seed is set and restored
#' at the end of the function. If \code{NULL}, the random number generator
#' state is not changed nor saved at the start, and not restored at the end.
#' @param checkSoftPower Logical: should the soft-tresholding power be adjusted
#' to approximately match the connectivity distribution of the sampled data set
#' and the full data set?
#' @param nPowerCheckSamples Number of genes to be sampled from the full data
#' set to calculate connectivity and match soft-tresholding powers.
#' @param skipUnsampledCalculation Logical: should a calculation on original
#' (not resampled) data be skipped?
#' @param corType Character string specifying the correlation to be used.
#' Allowed values are (unique abbreviations of) \code{"pearson"} and
#' \code{"bicor"}, corresponding to Pearson and bidweight midcorrelation,
#' respectively. Missing values are handled using the
#' \code{pairwise.complete.obs} option.
#' @param power Soft-thresholding power for network construction.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param saveTOMs Logical: should the networks (topological overlaps) be saved
#' for each run? Note that for large data sets (tens of thousands of nodes) the
#' TOM files are rather large.
#' @param saveTOMFileBase Character string giving the base of the file names
#' for TOMs. The actual file names will consist of a concatenation of
#' \code{saveTOMFileBase} and \code{"-run-<run number>-Block-<block
#' number>.RData"}.
#' @param \dots Other arguments to \code{\link{blockwiseModules}}.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with one component per run. Each component is a list with the
#' following components: \item{mods}{The output of the function
#' \code{\link{blockwiseModules}} applied to a resampled data set.}
#' \item{samples}{Indices of the samples selected for the resampled data step
#' for this run.} \item{powers}{Actual soft-thresholding powers used in this
#' run.}
#' @author Peter Langfelder
#' @seealso \code{\link{blockwiseModules}} for the underlying network analysis
#' and module identification
#'
#' \code{\link{sampledHierarchicalConsensusModules}} for a similar resampling
#' analysis of consensus networks.
#' @references An application of this function is described in the motivational
#' example section of
#'
#' Langfelder P, Horvath S (2012) Fast R Functions for Robust Correlations and
#' Hierarchical Clustering. Journal of Statistical Software 46(11) 1-17; PMID:
#' 23050260 PMCID: PMC3465711
#' @keywords misc
#' @export
sampledBlockwiseModules = function(
  datExpr,
  nRuns,
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE,
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  skipUnsampledCalculation = FALSE,
  corType = "pearson",
  power = 6,
  networkType = "unsigned",
  saveTOMs = FALSE,
  saveTOMFileBase = "TOM",
  ...,
  verbose = 2, indent = 0)

{

  spaces = indentSpaces(indent)

  result = list()
  runTOMFileBase = saveTOMFileBase
  nSamples = nrow(datExpr)
  nGenes = ncol(datExpr)

  corTypeI = pmatch(corType, .corTypes)
  if (is.na(corTypeI))
    stop(paste("Invalid 'corType'. Recognized values are", paste(.corTypes, collapse = ", ")))

  corFnc = .corFnc[corTypeI]

  seedSaved = FALSE
  if (!is.null(randomSeed)) {
    if (exists(".Random.seed")) {
       seedSaved = TRUE
       savedSeed = .Random.seed
    }
    set.seed(randomSeed)
  }

  if (checkSoftPower) {
    if (verbose > 0) printFlush(paste(spaces, "...calculating reference mean adjacencies.."))
    useGenes = sample(nGenes, nPowerCheckSamples, replace = FALSE)
    adj = adjacency(datExpr[, useGenes], power = power, type = networkType,
                      corFnc = corFnc)
    refAdjMeans = mean(as.dist(adj))
  }

  for (run in startRunIndex:endRunIndex) {
    set.seed(randomSeed + 2*run + 1)
    if (verbose > 0) printFlush(paste(spaces, "...working on run", run, ".."))
    if (saveTOMs)
      runTOMFileBase = paste(saveTOMFileBase, "-run-", run, sep = "")

    if (run > startRunIndex || skipUnsampledCalculation) {
      useSamples = sample(nSamples, as.integer(nSamples * fraction), replace = replace)
    } else
      useSamples = c(1:nSamples)

    if (verbose > 2) {
      printFlush(paste(spaces, "Using the following samples: "))
      print(useSamples)
    }
    samExpr = as.matrix(datExpr[useSamples, ])
    samPowers = power
    if (checkSoftPower) {
      if (verbose > 1) printFlush(paste(spaces, "  ...calculating mean adjacencies in sampled data.."))
      adj = adjacency(samExpr[, useGenes], power = power, type = networkType,
                      corFnc = corFnc)
      sampledAdjMeans = mean(as.dist(adj))
      samPowers = power * log( refAdjMeans) / log( sampledAdjMeans)
      if (!is.finite(samPowers)) samPowers = power
    }

    mods = blockwiseModules(
      datExpr = samExpr,
      randomSeed = NULL,
      power = samPowers,
      corType = corType,
      networkType = networkType,
      saveTOMs = saveTOMs,
      saveTOMFileBase = runTOMFileBase,
      ...,
      verbose = verbose-2, indent = indent+2)

    result[[run]] = list(mods = mods, samples = useSamples, powers = samPowers)
  }

  if (seedSaved) .Random.seed <<- savedSeed

  result
}

#===================================================================================================
#
# sampledHierarchicalConsensusModules
#
#===================================================================================================



#' Hierarchical consensus module identification in sampled data
#'
#' This function repeatedly resamples the samples (rows) in supplied data and
#' identifies hierarchical consensus modules on the resampled data.
#'
#' For each run, samples (but not genes) are randomly sampled to obtain a
#' perturbed data set; a full network analysis and module identification is
#' carried out, and the results are returned in a list with one component per
#' run.
#'
#' For each run, the soft-thresholding power can optionally be adjusted such
#' that the mean adjacency in the re-sampled data set equals the mean adjacency
#' in the original data.
#'
#' @param multiExpr Expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param networkOptions A single list of class \code{\link{NetworkOptions}}
#' giving options for network calculation for all of the networks, or a
#' \code{\link{multiSet}} structure containing one such list for each input
#' data set.
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param nRuns Number of network construction and module identification runs.
#' @param startRunIndex Number to be assigned to the start run. The run number
#' or index is used to make saved files unique; it has no effect on the actual
#' results of the run.
#' @param endRunIndex Number (index) of the last run. If given, \code{nRuns} is
#' ignored.
#' @param replace Logical: should samples (observations or rows in entries in
#' \code{multiExpr}) be sampled with replacement?
#' @param fraction Fraction of samples to sample for each run.
#' @param randomSeed Integer specifying the random seed. If non-NULL, the
#' random number generator state is saved before the seed is set and restored
#' at the end of the function. If \code{NULL}, the random number generator
#' state is not changed nor saved at the start, and not restored at the end.
#' @param checkSoftPower Logical: should the soft-tresholding power be adjusted
#' to approximately match the connectivity distribution of the sampled data set
#' and the full data set?
#' @param nPowerCheckSamples Number of genes to be sampled from the full data
#' set to calculate connectivity and match soft-tresholding powers.
#' @param individualTOMFileNames Pattern for file names for files holding
#' individual TOMs. The tags \code{\%r, \%a, \%b} are replaced by run number,
#' analysis name and block number, respectively. The TOM files are usually
#' temporary but can be retained, see \code{keepConsensusTOM} below.
#' @param keepConsensusTOMs Logical: should the (final) consensus TOMs of each
#' sampled calculation be retained after the run ends? Note that for large data
#' sets (tens of thousands of nodes) the TOM files are rather large.
#' @param consensusTOMFilePattern Patter of the file to be stored
#' @param skipUnsampledCalculation Logical: should a calculation on original
#' (not resampled) data be skipped?
#' @param \dots Other arguments to \code{\link{hierarchicalConsensusModules}}.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @param saveRunningResults Logical: should the cumulative results be saved
#' after each run on resampled data?
#' @param runningResultsFile File name of file in which to save running results
#' into. In case of a parallel execution (say on several nodes of a cluster),
#' one should choose a unique name for each process to avoid overwriting the
#' same file.
#' @return
#'
#' A list with one component per run. Each component is a list with the
#' following components: \item{mods}{The output of the function
#' \code{\link{hierarchicalConsensusModules}} on the resampled data.}
#' \item{samples}{Indices of the samples selected for the resampled data step
#' for this run.} \item{powers}{Actual soft-thresholding powers used in this
#' run.}
#' @author Peter Langfelder
#' @seealso \code{\link{hierarchicalConsensusModules}} for consensus networ
#' analysis and module identification
#'
#' \code{\link{sampledBlockwiseModules}} for a similar resampling analysis for
#' a single data set.
#' @keywords misc
#' @export
#' @importFrom utils object.size
sampledHierarchicalConsensusModules = function(
  multiExpr,

  networkOptions,
  consensusTree,

  nRuns,
  startRunIndex = 1,
  endRunIndex = startRunIndex + nRuns -1,
  replace = FALSE,
  fraction = if (replace) 1.0 else 0.63,
  randomSeed = 12345,
  checkSoftPower = TRUE,
  nPowerCheckSamples = 2000,
  individualTOMFileNames = paste0("individualTOM-Run.%r-Set%s-Block%b.RData"),
  keepConsensusTOMs = FALSE,
  consensusTOMFilePattern = "consensusTOM-Run.%r-%a-Block.%b.RData",
  skipUnsampledCalculation = FALSE,
  ...,
  verbose = 2, indent = 0,
  saveRunningResults = TRUE,
  runningResultsFile = "results.tmp.RData")
{

  spaces = indentSpaces(indent)

  result = list()
  exprSize = checkSets(multiExpr)
  nSets = exprSize$nSets
  nSamples = exprSize$nSamples

  if (inherits(networkOptions, "NetworkOptions"))
    networkOptions = list2multiSet(replicate(nSets, networkOptions, simplify = FALSE))

  if (!is.null(randomSeed)) {
    if (exists(".Random.seed")) {
       savedSeed = .Random.seed
       on.exit({ .Random.seed <<- savedSeed }, add = FALSE)
    }
    set.seed(randomSeed)
  }

  powers = unlist(apply(networkOptions, getElement, "power"))

  if (nPowerCheckSamples > exprSize$nGenes) nPowerCheckSamples = exprSize$nGenes

  if (checkSoftPower) {
    if (verbose > 0) printFlush(paste(spaces, "...calculating reference mean adjacencies.."))
    useGenes = sample(exprSize$nGenes, nPowerCheckSamples, replace = FALSE)
    refAdjMeans = rep(0, nSets)
    for (set in 1:nSets) {
      adj = adjacency(multiExpr[[set]]$data[, useGenes],
                      power = networkOptions[[set]]$data$power, type = networkOptions[[set]]$data$networkType,
                      corFnc = networkOptions[[set]]$data$corFnc,
                      corOptions = networkOptions[[set]]$data$corOptions)
      refAdjMeans[set] = mean(as.dist(adj))
    }
  }

  for (run in startRunIndex:endRunIndex) {
    runTOMFileBase = .substituteTags(consensusTOMFilePattern, "%r", run)
    individualTOMFiles1 = .substituteTags(individualTOMFileNames, "%r", run)
    set.seed(randomSeed + 2*run + 1)

    if (verbose > 0) printFlush(paste(spaces, "Working on run", run, ".."))

    useSamples = list()
    for (set in 1:nSets) {
      if (run > startRunIndex-skipUnsampledCalculation) {
         printFlush("This run will be on sampled data.")
         useSamples[[set]] = sample(nSamples[set],
                                    as.integer(nSamples[set] * fraction),
                                    replace = replace)
      } else
         useSamples[[set]] = c(1:nSamples[set])
    }
    samExpr = subset(multiExpr, useSamples)
    samPowers = powers
    if (checkSoftPower) {
      if (verbose > 1)
        printFlush(paste(spaces,
                         "  ...calculating mean adjacencies in sampled data.."))
      sampledAdjMeans = rep(0, nSets)
      for (set in 1:nSets) {
        adj = adjacency(samExpr[[set]]$data[, useGenes],
                        power = networkOptions[[set]]$data$power,
                        type = networkOptions[[set]]$data$networkType,
                        corFnc = networkOptions[[set]]$data$corFnc,
                        corOptions = networkOptions[[set]]$data$corOptions)
        sampledAdjMeans[set] = mean(as.dist(adj))
      }
      samPowers = powers * log( refAdjMeans) / log( sampledAdjMeans)
      samPowers[!is.finite(samPowers)] = powers[!is.finite(samPowers)]
    }

    networkOptions1 = multiSet.mapply(function(x, power) { x$power = power; x; },
                       networkOptions, samPowers)

    collectGarbage()

    mods = hierarchicalConsensusModules(
      multiExpr = samExpr,
      randomSeed = NULL,

      networkOptions = networkOptions1,
      consensusTree = consensusTree,

      saveConsensusTOM = TRUE,
      consensusTOMFilePattern = runTOMFileBase,
      saveIndividualTOMs = TRUE,
      individualTOMFileNames = individualTOMFiles1,

      keepIndividualTOMs = FALSE,
      keepConsensusTOM = keepConsensusTOMs,
      ...,
      verbose = verbose-2, indent = indent+2)

    result[[run]] = list(mods = mods, samples = useSamples, powers = samPowers)

    if (saveRunningResults) save(result, file = runningResultsFile)

    print(lapply(mods, object.size))
    print(gc())

  }

  result
}

#' Hierarchical consensus calculation of module eigengene dissimilarity
#'
#' Hierarchical consensus calculation of module eigengene dissimilarities, or
#' more generally, correlation-based dissimilarities of sets of vectors.
#'
#' This function first calculates the similarities of the ME vectors from their
#' correlations, using the appropriate options in \code{networkOptions}
#' (correlation type and options, signed or unsigned dissimilarity etc). This
#' results in a similarity matrix in each of the input data sets.
#'
#' Next, a hierarchical consensus of the similarities is calculated via a call
#' to \code{\link{hierarchicalConsensusCalculation}}, using the consensus
#' specification and options in \code{consensusTree}. In typical use,
#' \code{consensusTree} contains the same consensus specification as the
#' consensus network calculation that gave rise to the consensus modules whose
#' eigengenes are contained in \code{MEs} but this is not mandatory.
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
#' In the final step, the consensus similarity is turned into a dissimilarity
#' by subtracting it from 1.
#'
#' @param MEs A \code{\link{multiSet}} structure containing vectors (usually
#' module eigengenes) whose consensus dissimilarity is to be calculated.
#' @param networkOptions A \code{\link{multiSet}} structure containing, for
#' each input data set, a list of class \code{\link{NetworkOptions}} giving
#' options for network calculation for all of the networks.
#' @param consensusTree A list specifying the consensus calculation. See
#' details.
#' @param greyName Name of the "grey" module eigengene. Currently not used.
#' @param calibrate Logical: should the dissimilarities be calibrated using the
#' calibration method specified in \code{consensusTree}? See details.
#' @return A matrix with rows and columns corresponding to the variables
#' (modules) in MEs, containing the consensus dissimilarities.
#' @author Peter Langfelder
#' @seealso \code{\link{hierarchicalConsensusCalculation}} for the actual
#' consensus calculation.
#' @keywords misc
#' @export
hierarchicalConsensusMEDissimilarity = function(MEs, networkOptions,
                                                consensusTree,
                                                greyName = "ME0",
                                                calibrate = FALSE){
  nSets = checkSets(MEs)$nSets

  if (inherits(networkOptions, "NetworkOptions"))
    networkOptions = list2multiSet(.listRep(networkOptions, nSets))

  .hierarchicalConsensusMEDissimilarity(MEs, networkOptions, consensusTree,
                                        greyName = greyName, calibrate = calibrate)
}


.hierarchicalConsensusMEDissimilarity = function(multiMEs,
                                                 networkOptions,
                                                 consensusTree,
                                                 greyName,
                                                 calibrate)
{
  nSets = checkSets(multiMEs)$nSets
  useMEs = which(colnames(multiMEs)!=greyName)
  useNames = colnames(multiMEs)[useMEs]
  nUseMEs = length(useMEs)
  #  if (nUseMEs<2)
  #    stop("Something is wrong: there are two or more proper modules, but less than two proper",
  #         "eigengenes. Please check that the grey color label and module eigengene label",
  #         "are correct.")

  if (!is.multiSet(networkOptions, strict = FALSE))
    stop("'networkOptions' must be either a single list of class 'NetworkOptions'\n",
         "or a MultiData structure containing one such list per input set. ")

  if (length(networkOptions)!=nSets)
    stop("Number of sets in 'multiMEs' and 'networkOptions' must be the same.")

  MEDiss = multiSet.mapply(function(me, netOpt) {
    cor.me = do.call(netOpt$corFnc,
                     c(list(x = me), netOpt$corOptions))
    if (!grepl("signed", netOpt$networkType)) cor.me = abs(cor.me)
    cor.me
  }, subset(multiMEs, , useMEs), networkOptions, returnList = TRUE)

  if (calibrate) {
    cons = hierarchicalConsensusCalculation(MEDiss,
                                            consensusTree = consensusTree,
                                            level = 1,
                                            # Return options: the data can be either saved or returned but not both.
                                            saveConsensusData = FALSE,
                                            keepIntermediateResults = FALSE,
                                            # Internal handling of data
                                            useDiskCache = FALSE,
                                            # Behaviour
                                            collectGarbage = FALSE,
                                            verbose = 0, indent = 0)$consensusData
    cons = BD.getData(cons, blocks = 1)
  } else
    cons = simpleHierarchicalConsensusCalculation(MEDiss, consensusTree)
  consDiss = 1-cons
  colnames(consDiss) = rownames(consDiss) = useNames
  consDiss
}


#' Merge close (similar) hierarchical consensus modules
#'
#' Merges hierarchical consensus modules that are too close as measured by the
#' correlation of their eigengenes.
#'
#' This function merges input modules that are closely related. The
#' similarities are quantified by correlations of module eigengenes; a
#' ``consensus'' similarity is calculated using
#' \code{hierarchicalConsensusMEDissimilarity} according to the recipe in
#' \code{consensusTree}. Once the (dis-)similarities are calculated, average
#' linkage hierarchical clustering of the module eigengenes is performed, the
#' dendrogram is cut at the height \code{cutHeight} and modules on each branch
#' are merged. The process is (optionally) repeated until no more modules are
#' merged.
#'
#' If, for a particular module, the module eigengene calculation fails, a
#' hubgene approximation will be used.
#'
#' The user should be aware that if a computational error occurs and
#' \code{trapErrors==TRUE}, the returned list (see below) will not contain all
#' of the components returned upon normal execution.
#'
#' @param multiExpr Expression data in the multi-set format (see
#' \code{\link{multiSet}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param labels A vector (numeric, character or a factor) giving module colors
#' for genes. The method only makes sense when genes have the same color label
#' in all sets, hence a single vector.
#' @param MEs If module eigengenes have been calculated before, the user can
#' save some computational time by inputting them. \code{MEs} should have the
#' same format as \code{multiExpr}. If they are not given, they will be
#' calculated.
#' @param unassdColor Specifies the string that labels unassigned genes. Module
#' of this color will not enter the module eigengene clustering and will not be
#' merged with other modules.
#' @param impute Should missing values be imputed in eigengene calculation? If
#' imputation is disabled, the presence of \code{NA} entries will cause the
#' eigengene calculation to fail and eigengenes will be replaced by their
#' hubgene approximation. See \code{\link{moduleEigengenes}} for more details.
#' @param networkOptions A single list of class \code{\link{NetworkOptions}}
#' giving options for network calculation for all of the networks, or a
#' \code{\link{multiSet}} structure containing one such list for each input
#' data set.
#' @param consensusTree A list specifying the consensus calculation. See
#' \code{\link{newConsensusTree}} for details.
#' @param calibrateMESimilarities Logical: should module eigengene similarities
#' be calibrated? This setting overrides the calibration options in
#' \code{consensusTree}.
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param iterate Controls whether the merging procedure should be repeated
#' until there is no change. If FALSE, only one iteration will be executed.
#' @param relabel Controls whether, after merging, color labels should be
#' ordered by module size.
#' @param colorSeq Color labels to be used for relabeling. Defaults to the
#' standard color order used in this package if \code{colors} are not numeric,
#' and to integers starting from 1 if \code{colors} is numeric.
#' @param getNewMEs Controls whether module eigengenes of merged modules should
#' be calculated and returned.
#' @param getNewUnassdME When doing module eigengene manipulations, the
#' function does not normally calculate the eigengene of the 'module' of
#' unassigned ('grey') genes. Setting this option to \code{TRUE} will force the
#' calculation of the unassigned eigengene in the returned newMEs, but not in
#' the returned oldMEs.
#' @param trapErrors Controls whether computational errors in calculating
#' module eigengenes, their dissimilarity, and merging trees should be trapped.
#' If \code{TRUE}, errors will be trapped and the function will return the
#' input colors. If \code{FALSE}, errors will cause the function to stop.
#' @param verbose Controls verbosity of printed progress messages. 0 means
#' silent, up to (about) 5 the verbosity gradually increases.
#' @param indent A single non-negative integer controlling indentation of
#' printed messages. 0 means no indentation, each unit above that adds two
#' spaces.
#' @return If no errors occurred, a list with components \item{labels}{Labels
#' for the genes corresponding to merged modules. The function attempts to
#' mimic the mode of the input \code{labels}: if the input \code{labels} is
#' numeric, character and factor, respectively, so is the output. Note,
#' however, that if the function performs relabeling, a standard sequence of
#' labels will be used: integers starting at 1 if the input \code{labels} is
#' numeric, and a sequence of color labels otherwise (see \code{colorSeq}
#' above).}
#'
#' \item{dendro}{Hierarchical clustering dendrogram (average linkage) of the
#' eigengenes of the most recently computed tree. If \code{iterate} was set
#' TRUE, this will be the dendrogram of the merged modules, otherwise it will
#' be the dendrogram of the original modules.}
#'
#' \item{oldDendro}{Hierarchical clustering dendrogram (average linkage) of the
#' eigengenes of the original modules.}
#'
#' \item{cutHeight}{The input cutHeight.}
#'
#' \item{oldMEs}{Module eigengenes of the original modules in the sets given by
#' \code{useSets}.}
#'
#' \item{newMEs}{Module eigengenes of the merged modules in the sets given by
#' \code{useSets}.}
#'
#' \item{allOK}{A logical set to \code{TRUE}.}
#'
#' If an error occurred and \code{trapErrors==TRUE}, the list only contains
#' these components:
#'
#' \item{colors}{A copy of the input colors.}
#'
#' \item{allOK}{a logical set to \code{FALSE}.}
#' @author Peter Langfelder
#' @seealso \code{\link{multiSetMEs}} for calculation of (consensus) module
#' eigengenes across multiple data sets
#'
#' \code{\link{newConsensusTree}} for information about consensus trees
#'
#' \code{\link{hierarchicalConsensusMEDissimilarity}} for calculation of
#' hierarchical consensus eigengene dissimilarity.
#' @keywords misc
#' @export
hierarchicalMergeCloseModules = function(
  # input data
  multiExpr, labels,

  # Optional starting eigengenes
  MEs = NULL,

  unassdColor = if (is.numeric(labels)) 0 else "grey",
  # If missing data are present, impute them?
  impute = TRUE,


  # Options for eigengene network construction
  networkOptions,

  # Options for constructing the consensus
  consensusTree,
  calibrateMESimilarities = FALSE,

  # Merging options
  cutHeight = 0.2,
  iterate = TRUE,

  # Output options
  relabel = FALSE,
  colorSeq = NULL,
  getNewMEs = TRUE,
  getNewUnassdME = TRUE,

  # Options controlling behaviour of the function
  trapErrors = FALSE,
  verbose = 1, indent = 0)
{

  MEsInSingleFrame = FALSE
  spaces = indentSpaces(indent)

  #numCols = is.numeric(labels)
  #facCols = is.factor(labels)
  #charCols = is.character(labels)

  origColors = labels

  labels = labels[, drop = TRUE]

  greyName = paste(moduleColor.getMEprefix(), unassdColor, sep="")

  if (verbose>0) printFlush(paste(spaces,
                                  "mergeCloseModules: Merging modules whose distance is less than", cutHeight))

  if (verbose>3) printFlush(paste(spaces,
                                  "  .. will use unassigned ME label", greyName))

  setsize = checkSets(multiExpr)
  nSets = setsize$nSets

  if (!is.null(MEs))
  {
    checkMEs = checkSets(MEs, checkStructure = TRUE)
    if (checkMEs$structureOK)
    {
      if (nSets!=checkMEs$nSets)
        stop("Input error: numbers of sets in multiExpr and MEs differ.")
      for (set in 1:nSets)
      {
        if (checkMEs$nSamples[set]!=setsize$nSamples[set])
          stop(paste("Number of samples in MEs is incompatible with subset length for set", set))
      }
    } else {
      if (MEsInSingleFrame)
      {
        MEs = fixDataStructure(MEs)
        checkMEs = checkSets(MEs)
      } else {
        stop("MEs do not have the appropriate structure (same as multiExpr). ")
      }
    }
  }

  if (inherits(networkOptions, "NetworkOptions"))
    networkOptions = list2multiSet(.listRep(networkOptions, nSets))

  if (setsize$nGenes!=length(labels))
    stop("Number of genes in multiExpr is different from the length of original labels. They must equal.")

  done = FALSE; iteration = 1

  MergedColors = labels
  #ok = try(
  #{
  while (!done)
  {
    if (is.null(MEs))
    {
      MEs = multiSetMEs(multiExpr, colors = NULL, universalColors = labels,
                        impute = impute,
                        subHubs = TRUE, trapErrors = FALSE, excludeGrey = TRUE,
                        grey = unassdColor,
                        verbose = verbose-1, indent = indent+1)
      #MEs = consensusOrderMEs(MEs, useAbs = useAbs, greyLast = FALSE)
      #collectGarbage()
    } else if (nlevels(as.factor(labels))!=checkMEs$nGenes)
    {
      if ((iteration==1) & (verbose>0)) printFlush(paste(spaces, "  Number of given module labels",
                                                         "does not match number of given MEs => recalculating the MEs."))
      MEs = multiSetMEs(multiExpr, colors = NULL, universalColors = labels,
                        impute = impute,
                        subHubs = TRUE, trapErrors = FALSE, excludeGrey = TRUE,
                        grey = unassdColor,
                        verbose = verbose-1, indent = indent+1)
      #MEs = consensusOrderMEs(MEs, useAbs = useAbs, greyLast = FALSE)
      #collectGarbage()
    }
    if (iteration==1) oldMEs = MEs

    # Check labels for number of distinct labels that are not grey

    colLevs = as.character(levels(as.factor(labels)))
    if  ( length(colLevs[colLevs!=as.character(unassdColor)])<2 )
    {
      printFlush(paste(spaces,
                       "mergeCloseModules: less than two proper modules."))
      printFlush(paste(spaces, " ..color levels are",
                       paste(colLevs, collapse = ", ")))
      printFlush(paste(spaces, " ..there is nothing to merge."))
      MergedNewColors = labels
      MergedColors = labels
      nOldMods = 1; nNewMods = 1
      oldTree = NULL; Tree = NULL
      break
    }

    # Cluster the found module eigengenes and merge ones that are too close according to the specified
    # quantile.

    nOldMods = nlevels(as.factor(labels))

    ConsDiss = .hierarchicalConsensusMEDissimilarity(MEs, networkOptions = networkOptions,
                                                     consensusTree = consensusTree,
                                                     greyName = greyName,
                                                     calibrate = calibrateMESimilarities)

    Tree = fastcluster::hclust(as.dist(ConsDiss), method = "average")
    if (iteration==1) oldTree = Tree
    TreeBranches = as.factor(moduleNumber(dendro = Tree, cutHeight = cutHeight, minSize = 1))
    UniqueBranches = levels(TreeBranches)
    nBranches = nlevels(TreeBranches)
    NumberOnBranch = table(TreeBranches)
    MergedColors = labels

    # Merge modules on the same branch

    for (branch in 1:nBranches) if (NumberOnBranch[branch]>1) {
      ModulesOnThisBranch = names(TreeBranches)[TreeBranches==UniqueBranches[branch]]
      ColorsOnThisBranch = substring(ModulesOnThisBranch, 3)
      if (is.numeric(origColors)) ColorsOnThisBranch = as.numeric(ColorsOnThisBranch)
      if (verbose>3)
        printFlush(paste(spaces, "  Merging original labels",
                         paste(ColorsOnThisBranch, collapse=", ")))
      for (color in 2:length(ColorsOnThisBranch))
        MergedColors[MergedColors==ColorsOnThisBranch[color]] = ColorsOnThisBranch[1]
    }

    MergedColors = MergedColors[, drop = TRUE]

    nNewMods = nlevels(as.factor(MergedColors))

    if (nNewMods<nOldMods & iterate) {
      labels = MergedColors
      MEs = NULL
    } else {
      done = TRUE
    }
    iteration = iteration+1
  }
  if (relabel) {
    RawModuleColors = levels(as.factor(MergedColors))
    # relabel the merged labels to the usual order based on the number of genes in each module
    if (is.null(colorSeq)) {
      if (is.numeric(origColors)) {
        colorSeq = c(1:length(table(origColors)))
      } else {
        nNewColors = length(RawModuleColors)
        colorSeq = labels2colors(c(1:nNewColors))
      }
    }

    # nGenesInModule = rep(0, nNewMods)
    # for (mod in 1:nNewMods) nGenesInModule[mod] = sum(MergedColors==RawModuleColors[mod])
    nGenesInModule = table(MergedColors)

    SortedRawModuleColors = RawModuleColors[order(-nGenesInModule)]

    # Change the color names to the standard sequence, but leave grey grey
    # (that's why rank in general does not equal color)
    MergedNewColors = MergedColors
    if (is.factor(MergedNewColors)) MergedNewColors = as.character(MergedNewColors)
    if (verbose>3) printFlush(paste(spaces, "   Changing original labels:"))
    rank = 0
    for (color in 1:length(SortedRawModuleColors)) if (SortedRawModuleColors[color]!=unassdColor)
    {
      rank = rank + 1
      if (verbose>3) printFlush(paste(spaces, "      ", SortedRawModuleColors[color],
                                      "to ", colorSeq[rank]))
      MergedNewColors[MergedColors==SortedRawModuleColors[color]] = colorSeq[rank]
    }
    if (is.factor(MergedColors)) MergedNewColors = as.factor(MergedNewColors)
  } else {
    MergedNewColors = MergedColors
  }
  MergedNewColors = MergedNewColors[, drop = TRUE]

  if (getNewMEs) {
    if (nNewMods<nOldMods | relabel | getNewUnassdME) {
      if (verbose>0) printFlush(paste(spaces, "  Calculating new MEs..."))
      NewMEs = multiSetMEs(multiExpr, colors = NULL, universalColors = MergedNewColors,
                           impute = impute, subHubs = TRUE, trapErrors = FALSE,
                           excludeGrey = !getNewUnassdME, grey = unassdColor,
                           verbose = verbose-1, indent = indent+1)
      newMEs = orderMEsByHierarchicalConsensus(NewMEs, networkOptions,
                                               consensusTree,
                                               greyName = greyName, calibrate = calibrateMESimilarities)

      ConsDiss = .hierarchicalConsensusMEDissimilarity(newMEs, networkOptions = networkOptions,
                                                       consensusTree = consensusTree,
                                                       greyName = greyName,
                                                       calibrate = calibrateMESimilarities)
      if (length(ConsDiss) > 1) {
        Tree = fastcluster::hclust(as.dist(ConsDiss), method = "average")
      } else Tree = NULL
    } else {
      newMEs = MEs
    }
  } else {
    newMEs = NULL
  }

  list(labels = MergedNewColors, dendro = Tree, oldDendro = oldTree, cutHeight = cutHeight,
       oldMEs = oldMEs, newMEs = newMEs, allOK = TRUE)
}


#' Order module eigengenes by their hierarchical consensus similarity
#'
#' This function calculates a hiearchical consensus similarity of the input
#' eigengenes, clusters the eigengenes according to the similarity and returns
#' the input module eigengenes ordered by the order of resulting dendrogram.
#'
#'
#' @param MEs Module eigengenes, or more generally, vectors, to be ordered, in
#' a \code{\link{multiSet}} format: A vector of lists, one per set. Each set
#' must contain a component \code{data} that contains the module eigenegens or
#' general vectors, with rows corresponding to samples and columns to genes or
#' probes.
#' @param networkOptions A single list of class \code{\link{NetworkOptions}}
#' giving options for network calculation for all of the networks, or a
#' \code{\link{multiSet}} structure containing one such list for each input
#' data set.
#' @param consensusTree A list specifying the consensus calculation. See
#' \code{\link{newConsensusTree}} for details.
#' @param greyName Specifies the column name of eigengene of the "module" that
#' contains unassigned genes. This eigengene (column) will be excluded from the
#' clustering and will be put last in the order.
#' @param calibrate Logical: should module eigengene similarities be
#' calibrated? This setting overrides the calibration options in
#' \code{consensusTree}.
#' @return A \code{\link{multiSet}} structure of the same format as the input
#' \code{MEs}, with columns ordered by the calculated dendrogram.
#' @author Peter Langfelder
#' @seealso \code{\link{hierarchicalConsensusMEDissimilarity}} for calculating
#' the consensus ME dissimilarity
#' @keywords misc
#' @export
orderMEsByHierarchicalConsensus = function(MEs, networkOptions, consensusTree,
                                           greyName = "ME0",
                                           calibrate = FALSE)
{
  Diss = .hierarchicalConsensusMEDissimilarity(MEs, networkOptions, consensusTree,
                                               greyName = greyName, calibrate = calibrate)
  order = .clustOrder(Diss, greyLast = TRUE, greyName = greyName)
  subset(MEs, , order)
}
