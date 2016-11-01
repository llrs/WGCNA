#======================================================================================================
#
# Aligned svd: svd plus aligning the result along the average expression of the input data.
#
#======================================================================================================

# CAUTION: this function assumes normalized x and no non-missing values.

.alignedFirstPC = function(x, power = 6, verbose = 2, indent = 0)
{
  x = as.matrix(x)
  #printFlush(paste(".alignedFirstPC: dim(x) = ", paste(dim(x), collapse = ", ")))
  pc = try( svd(x, nu = 1, nv = 0)$u[,1] , silent = TRUE)
  if (class(pc)=='try-error')
  {
    #file = "alignedFirstPC-svdError-inputData-x.RData"
    #save(x, file = file)
    #stop(paste("Error in .alignedFirstPC: svd failed with following message: \n    ",
    #                 pc, "\n. Saving the offending data into file", file))
    if (verbose > 0)
    {
      spaces = indentSpaces(indent)
      printFlush(paste(spaces, ".alignedFirstPC: FYI: svd failed, using a weighted mean instead.\n",
                       spaces, "  ...svd reported:", pc))
    }
    pc = rowMeans(x, na.rm = TRUE)
    weight = matrix(abs(cor(x, pc, use = 'p'))^power, nrow(x), ncol(x), byrow = TRUE)
    pc = scale(rowMeans(x * weight, na.rm = TRUE))
  } else {
    weight = abs(cor(x, pc))^power
    meanX = rowWeightedMeans(x, weight, na.rm = TRUE)
    cov1 = cov(pc, meanX)
    if (!is.finite(cov1)) cov1 = 0
    if (cov1 < 0) pc = -pc
  }
  pc
}


#==========================================================================================================
#
# Branch eigengene split (dissimilarity) calculation
#
#==========================================================================================================

# Assumes correct input: multiExpr is scaled to mean 0 and variance 1, branch1 and branch2 are numeric
# indices that have no overlap.



#' Branch dissimilarity based on eigennodes (eigengenes).
#'
#' Calculation of branch dissimilarity based on eigennodes (eigengenes) in
#' single set and multi-data situations. This function is used as a plugin for
#' the dynamicTreeCut package and the user should not call this function
#' directly. This function is experimental and subject to change.
#'
#'
#' @aliases branchEigengeneDissim multiSet.branchEigengeneDissim
#' @param expr Expression data.
#' @param branch1 Branch 1.
#' @param branch2 Branch 2.
#' @param corFnc Correlation function.
#' @param corOptions Other arguments to the correlation function.
#' @param signed Should the network be considered signed?
#' @param \dots Other arguments for compatibility; currently unused.
#' @return A single number or a list containing details of the calculation.
#' @author Peter Langfelder
#' @keywords misc
branchEigengeneDissim = function(expr, branch1, branch2,
                            corFnc = cor, corOptions = list(use = 'p'),
                            signed = TRUE, ...)
{
  expr.branch1 = expr[ , branch1]
  expr.branch2 = expr[ , branch2]

  corOptions$x = .alignedFirstPC(expr.branch1, verbose = 0)
  corOptions$y = .alignedFirstPC(expr.branch2, verbose = 0)

  corFnc = match.fun(corFnc)
  cor0 = as.numeric(do.call(corFnc, corOptions))
  if (length(cor0) != 1) stop ("Internal error in branchEigengeneDissim: cor has length ", length(cor0))
  if (signed)  1-cor0 else 1-abs(cor0)
}

multiSet.branchEigengeneDissim = function(multiExpr, branch1, branch2,
                            corFnc = cor, corOptions = list(use = 'p'),
                            consensusQuantile = 0, signed = TRUE, reproduceQuantileError = FALSE, ...)
{
  setSplits.list = multiSet.apply(multiExpr, branchEigengeneDissim,
                        branch1 = branch1, branch2 = branch2,
                        corFnc = corFnc, corOptions = corOptions,
                        signed = signed)
  setSplits = unlist(multiSet2list(setSplits.list))
  quantile(setSplits, prob = if (reproduceQuantileError) consensusQuantile else 1-consensusQuantile,
           na.rm = TRUE, names = FALSE)
}


#==========================================================================================================
#
# Branch split calculation
#
#==========================================================================================================
# Caution: highly experimental!

# Define a function that's supposed to decide whether two given groups of expression data are part of a
# single module, or truly two independent modules.
# assumes no missing values for now.
# assumes all data is scaled to mean zero and equal variance.

# return: criterion is zero or near zero if it looks like a single module, and is near 1 if looks like two
# modules.

# Careful: in the interest of speedy execution, the function doesn't check arguments for validity. For
# example, it assumes that expr is already scaled to the same mean and variance, branch1 and branch2 are valid
# indices, nConsideredPCs does not exceed any of the dimensions of expr etc.

.histogramsWithCommonBreaks = function(data, groups, discardProp = 0.08)
{

  if (is.list(data))
  {
    lengths = sapply(data, length)
    data = data[lengths>0]
    lengths = sapply(data, length)
    nGroups = length(lengths)
    groups = rep( c(1:nGroups), lengths)
    data = unlist(data)
  }

  if (discardProp > 0)
  {
    # get rid of outliers on either side - those are most likely not interesting.
    # The code is somewhat involved because I want to get rid of outliers that are defined with respect to
    # the combined data, but no more than a certain proportion of either of the groups.

    sizes = table(groups)
    nAll = length(data)
    order = order(data)
    ordGrp = groups[order]
    cs = rep(0, nAll)
    nGroups = length(sizes)
    for (g in 1:nGroups)
      cs[ordGrp==g] = ((1:sizes[g])-0.5)/sizes[g]

    firstKeep = min(which(cs > discardProp))
    first = data[order[firstKeep]]

    # Analogous upper quantile
    lastKeep = max(which(cs < 1-discardProp))
    last = data[order[lastKeep]]

    keep = ( (data >= first) & (data <= last) )
    data = data[keep]
    groups = groups[keep]
  } else {
    last = max(data, na.rm = TRUE)
    first = min(data, na.rm = TRUE)
  }

  # Find the smaller of the two groups and define histogram bin size from the number of elements in that
  # group; the aim is to prevent the group getting splintered into too many bins of the histogram.

  sizes = table(groups)
  smallerInd = which.min(sizes)
  smallerSize = sizes[smallerInd]
  nBins = ceiling(5 + ifelse(smallerSize > 25, sqrt(smallerSize)-4, 1 ))

  smaller = data[groups==smallerInd]
  binSize = (max(smaller) - min(smaller))/nBins

  nAllBins = ceiling((last-first)/binSize)
  breaks = first + c(0:nAllBins) * (last - first)/nAllBins

  tapply(data, groups, hist, breaks = breaks, plot = FALSE)
}





#' Branch split.
#'
#' Calculation of branch split based on expression data. This function is used
#' as a plugin for the dynamicTreeCut package and the user should not call this
#' function directly.
#'
#'
#' @param expr Expression data.
#' @param branch1 Branch 1,
#' @param branch2 Branch 2.
#' @param discardProp Proportion of data to be discarded as outliers.
#' @param minCentralProp Minimum central proportion
#' @param nConsideredPCs Number of principal components to consider.
#' @param signed Should the network be considered signed?
#' @param getDetails Should details of the calculation be returned?
#' @param \dots Other arguments. Present for compatibility; currently unusued.
#' @return A single number or a list containing detils of the calculation.
#' @author Peter Langfelder
#' @keywords misc
branchSplit = function(expr, branch1, branch2, discardProp = 0.05, minCentralProp = 0.75,
                       nConsideredPCs = 3, signed = FALSE, getDetails = TRUE, ...)
{
  nGenes = c(length(branch1), length(branch2))
  #combined = cbind(expr[, branch1], expr[, branch2])
  combinedScaled = cbind(expr[, branch1]/sqrt(length(branch1)), expr[, branch2]/sqrt(length(branch2)))
  groups = c(rep(1, nGenes[1]), rep(2, nGenes[2]) )

  # get the combination of PCs that best approximates the groups vector
  svd = svd(combinedScaled, nu = 0, nv = nConsideredPCs)
  v2 = svd$v * c( rep(sqrt(length(branch1)), length(branch1)), rep(sqrt(length(branch2)), length(branch2)))

  #svd = svd(combinedScaled, nu = nConsideredPCs, nv = 0)
  #v2 = cor(combinedScaled, svd$u)

  if (!signed) v2 = v2 * sign(v2[, 1])

  cor2 = predict(lm(groups~., data = as.data.frame(v2)))

  # get the histograms of the projections in both groups, but make sure the binning is the same for both.

  # get rid of outliers on either side - those are most likely not interesting.
  # The code is somewhat involved because I want to get rid of outliers that are defined with respect to
  # the combined data, but no more than a certain proportion of either of the groups.

  h = .histogramsWithCommonBreaks(cor2, groups, discardProp)

  maxAll = max(c(h[[1]]$counts, h[[2]]$counts))
  h[[1]]$counts = h[[1]]$counts/maxAll
  h[[2]]$counts = h[[2]]$counts/maxAll
  max1 = max(h[[1]]$counts)
  max2 = max(h[[2]]$counts)
  minMax = min(max1, max2)

  if (FALSE) {
    plot(h[[1]]$mids, h[[1]]$counts, type = "l")
    lines(h[[2]]$mids, h[[2]]$counts, type = "l", col = "red")
    lines(h[[2]]$mids, h[[1]]$counts + h[[2]]$counts, type = "l", col = "blue")
  }
  # Locate "central" bins: those whose scaled counts exceed a threshold.

  central = list()
  central[[1]] = h[[1]]$counts > minCentralProp * minMax
  central[[2]] = h[[2]]$counts > minCentralProp * minMax

  # Do central bins overlap?

  overlap = (min(h[[1]]$mids[central[[1]]]) <= max(h[[2]]$mids[central[[2]]])) &
               (min(h[[2]]$mids[central[[2]]]) <= max(h[[1]]$mids[central[[1]]]))

  if (overlap)
  {
    result = list(middleCounts = NULL, criterion = minCentralProp, split = -1, histograms = h)
  } else {
    # Locate the region between the two central regions and check whether the gap is deep enough.
    if (min(h[[1]]$mids[central[[1]]]) > max(h[[2]]$mids[central[[2]]]))
    {
      left = 2; right = 1
    } else {
      left = 1; right = 2
    }
    leftEdge = max(h[[left]]$mids[central[[left]]])
    rightEdge = min(h[[right]]$mids[central[[right]]])
    middle = ( (h[[left]]$mids > leftEdge) & (h[[left]]$mids < rightEdge) )
    nMiddle = sum(middle)
    if (nMiddle==0)
    {
      result = list(middleCounts = NULL, criterion = minCentralProp, split = -1, histograms = h)
    } else {
      # Reference level: 75th percentile of the central region of the smaller branch
      #refLevel1 = quantile(h[[1]]$counts [ central[[1]] ], prob = 0.75)
      #refLevel2 = quantile(h[[2]]$counts [ central[[2]] ], prob = 0.75)
      refLevel1 = mean(h[[1]]$counts [ central[[1]] ], na.rm = TRUE)
      refLevel2 = mean(h[[2]]$counts [ central[[2]] ], na.rm = TRUE)
      peakRefLevel = min(refLevel1, refLevel2)

      middleCounts = h[[left]]$counts[middle] + h[[right]]$counts[middle]
      #troughRefLevel = quantile(middleCounts, prob = 0.25)
      troughRefLevel = mean(middleCounts, na.rm = TRUE)
      meanCorrFactor = sqrt(min(nMiddle + 1, 3) / min(nMiddle, 3))
               # =sqrt(2, 3/2, 1), for nMiddle=1,2,3,..
      result = list(middleCounts = middleCounts,
                    criterion = troughRefLevel * meanCorrFactor,
                    split = (peakRefLevel - troughRefLevel * meanCorrFactor)/peakRefLevel,
                    histograms = h)
    }
  }
  if (getDetails) result else result$split
}



#==========================================================================================================
#
# Dissimilarity-based branch split
#
#==========================================================================================================

.meanInRange = function(mat, rangeMat)
{
  nc = ncol(mat)
  means = rep(0, nc)
  for (c in 1:nc)
  {
    col = mat[, c]
    means[c] = mean( col[col >=rangeMat[c, 1] & col <= rangeMat[c, 2]], na.rm= TRUE)
  }
  means
}

.sizeDependentQuantile = function(p, sizes, minNumber = 5)
{
  refSize = minNumber/p
  correctionFactor = pmin( rep(1, length(sizes)), sizes/refSize)

  pmin(rep(1, length(sizes)), p/correctionFactor)
}




#' Branch split based on dissimilarity.
#'
#' Calculation of branch split based on a dissimilarity matrix. This function
#' is used as a plugin for the dynamicTreeCut package and the user should not
#' call this function directly. This function is experimental and subject to
#' change.
#'
#'
#' @param dissimMat Dissimilarity matrix.
#' @param branch1 Branch 1.
#' @param branch2 Branch 2.
#' @param upperP Percentile of (closest) objects to be considered.
#' @param minNumberInSplit Minimum number of objects to be considered.
#' @param getDetails Should details of the calculation be returned?
#' @param \dots Other arguments for compatibility; currently unused.
#' @return A single number or a list containing details of the calculation.
#' @author Peter Langfelder
#' @keywords misc
branchSplit.dissim = function(dissimMat, branch1, branch2, upperP,
                              minNumberInSplit = 5, getDetails = FALSE, ...)
{
  lowerP = 0

  sizes = c(length(branch1), length(branch2))
  upperP = .sizeDependentQuantile(upperP, sizes, minNumber = minNumberInSplit)

  multiP = as.data.frame(rbind(rep(0, 2), upperP))

  outDissim = list(list(data = dissimMat[branch2, branch1]),
                     list(data = dissimMat[branch1, branch2]))

  quantiles = multiSet.mapply(colQuantiles, outDissim, probs = multiP, MoreArgs = list(drop = FALSE))
  averages = multiSet.mapply(.meanInRange, outDissim, quantiles)

  averageQuantiles = multiSet.mapply(quantile, averages, prob = multiP, MoreArgs = list(drop = FALSE))

  betweenQuantiles = multiSet.mapply(function(x, quantiles) { x>=quantiles[1] & x <=quantiles[2]},
                                averages, averageQuantiles)

  selectedDissim = list(list(data = dissimMat[branch1, branch1[betweenQuantiles[[1]]$data] ]),
                    list(data = dissimMat[branch2, branch1[ betweenQuantiles[[1]]$data] ]),
                    list(data = dissimMat[branch2, branch2[betweenQuantiles[[2]]$data]]),
                    list(data = dissimMat[branch1, branch2[betweenQuantiles[[2]]$data]]))

  #n1 = length(branch1)
  #m1 = sum(betweenQuantiles[[1]]$data)
  #indexMat = cbind((1:n1)[betweenQuantiles[[1]]$data], 1:m1)

  # Remove the points nearest to branch 2 from the distances in branch 1
  selectedDissim[[1]]$data[ betweenQuantiles[[1]]$data, ] = NA

  #n2 = length(branch2)
  #m2 = sum(betweenQuantiles[[2]]$data)
  #indexMat = cbind((1:n2)[betweenQuantiles[[2]]$data], 1:m2)

  selectedDissim[[3]]$data[ betweenQuantiles[[2]]$data, ] = NA

  multiP.ext = cbind(multiP, multiP[, c(2,1)])

  selectedDissimQuantiles = multiSet.mapply(colQuantiles, selectedDissim, probs = multiP.ext,
                                      MoreArgs = list( drop = FALSE, na.rm = TRUE))

  selectedAverages = multiSet.mapply(.meanInRange, selectedDissim, selectedDissimQuantiles)

  if (FALSE)
  {
    par(mfrow = c(1,2))
    verboseBoxplot(c(selectedAverages[[1]]$data, selectedAverages[[2]]$data),
                    c( rep("in", length(selectedAverages[[1]]$data)),
                       rep("out", length(selectedAverages[[2]]$data))), main = "branch 1",
                   xlab = "", ylab = "mean distance")
    verboseBoxplot(c(selectedAverages[[3]]$data, selectedAverages[[4]]$data),
                    c( rep("in", length(selectedAverages[[3]]$data)),
                       rep("out", length(selectedAverages[[4]]$data))), main = "branch 2",
                   xlab = "", ylab = "mean distance")
  }

  separation = function(x, y)
  {
    nx = length(x)
    ny = length(y)
    if (nx*ny==0) return(0)

    mx = mean(x, na.rm = TRUE)
    my = mean(y, na.rm = TRUE)

    if (!is.finite(mx) | !is.finite(my)) return(0)

    if (nx > 1) varx = var(x, na.rm = TRUE) else varx = 0
    if (ny > 0) vary = var(y, na.rm = TRUE) else vary = 0

    if (is.na(varx)) varx = 0
    if (is.na(vary)) vary = 0

    if (varx + vary == 0) {
      if (my==mx) return(0) else return(Inf)
    }
    out = abs(my-mx)/(sqrt(varx) + sqrt(vary))
    if (is.na(out)) out = 0
    out
  }

  separations = c(separation(selectedAverages[[1]]$data,
                             selectedAverages[[2]]$data),
                  separation(selectedAverages[[3]]$data,
                             selectedAverages[[4]]$data))

  out = max(separations, na.rm = TRUE)
  if (is.na(out)) {
      out = 0
  }
  if (out < 0) {
      out = 0
  }
  if (getDetails) {
    return(list(split = out,
                distances = list(within1 = selectedAverages[[1]]$data,
                                 from1to2 = selectedAverages[[2]]$data,
                                 within2 = selectedAverages[[3]]$data,
                                 from2to1 = selectedAverages[[4]]$data)))
  }
  out
}

#========================================================================================================
#
# Branch dissimilarity based on a series of alternate branch labels
#
#========================================================================================================

# this function measures the split of branch1 and branch2 based on alternate labels, typically derived from
# resampled or otherwise perturbed data (but could also be derived from an independent data set).

# Basic idea: if two branches are separate, their membership should be predictable from the alternate
# labels.

# For each module in the alternate labeling: assign it to the branch in which |module ^ branch|/|branch| is
# bigger; then calculate |module ^ branch|/|branch| for the other branch as a classification error. Add all
# classification errors for a single alternate labeling and average them over labelings.

# This method is invariant under splitting of alternate module as long as the branch to which the modules are
# assigned does not change. So in this sense the splitting settings in the resampling study shouldn't
# matter too much but to some degree they still do.

# stabilityLabels: a matrix of dimensions (nGenes) x (number of alternate labels)



#' Branch split (dissimilarity) statistic derived from labels determined from a
#' stability study
#'
#' This function evaluates how different two branches are based on a series of
#' cluster labels that are usually obtained in a stability study but can in
#' principle be arbitrary. The idea is to quantify how well membership on the
#' two tested branches can be predicted from clusters in the given stability
#' labels.
#'
#' The idea is to measure how well clusters in \code{stabilityLabels} can
#' distinguish the two given branches. For example, if a cluster C intersects
#' with branch1 but not branch2, it can distinguish branches 1 and 2 perfectly.
#' On the other hand, if there is a cluster C that contains both branch 1 and
#' branch 2, the two branches are indistinguishable (based on the test
#' clustering).
#'
#' Formally, for each cluster C in each clustering in \code{stabilityLabels},
#' its contribution to the branch similarity is min(r1, r2), where r1 =
#' |intersect(C, branch1)|/|branch1| and r2 = |intersect(C,
#' branch2)|/|branch2|. The statistics for clusters in each clustering are
#' added; the sums are then averaged across the clusterings. Since the result
#' is a similarity statistic, the final dissimilarity is defined as
#' 1-similarity. The dissimilarity ranges between 0 (branch1 and branch2 are
#' indistinguishable) and 1 (branch1 and branch2 are perfectly
#' distinguishable).
#'
#' This is a very simple statistic that does not attempt to correct for the
#' similarity that would be expected by chance.
#'
#' @param branch1 A vector of indices giving members of branch 1.
#' @param branch2 A vector of indices giving members of branch 1.
#' @param stabilityLabels A matrix of cluster labels. Each column corresponds
#' to one clustering and each row to one object (whose indices \code{branch1}
#' and \code{branch2} refer to).
#' @param ignoreLabels Label or labels that do not constitute proper clusters
#' in \code{stabilityLabels}, for example because they label unassigned
#' objects.
#' @param \dots Ignored.
#' @return Branch dissimilarity (a single number between 0 and 1).
#' @author Peter Langfelder
#' @seealso This function is utilized in \code{\link{blockwiseModules}} and
#' \code{\link{blockwiseConsensusModules}}.
#' @keywords misc
branchSplitFromStabilityLabels <- function( branch1, branch2, stabilityLabels,
                                           ignoreLabels = 0, ...) {
  nLabels = ncol(stabilityLabels)
  n1 = length(branch1)
  n2 = length(branch2)

  sim = 0
  for (l in 1:nLabels) {
      lab1 = stabilityLabels[branch1, l]
      lab2 = stabilityLabels[branch2, l]
      commonLevels = intersect(unique(lab1), unique(lab2))
      commonLevels = setdiff(commonLevels, ignoreLabels)
      if (length(commonLevels) > 0) {
          for (cl in commonLevels) {
              #printFlush(paste0("Common level ", cl, " in clustering ", l))
              r1 = sum(lab1==cl)/n1
              r2 = sum(lab2==cl)/n2
              sim = sim + min(r1, r2)
          }
      }
  }

  1-sim/nLabels
}
