# nearest centroid predictor

# 09: Add weighted measures of sample-centroid similarity.
#     Fix CVfold such that leave one out cross-validation works.

# 08: change the way sample clustering is called. Use a do.call.

# 07: modify bagging and add boosting
#     - revert the predictor to state where it can only take a single trait and a single number of features.
#     - make sure the predictor respects nFeatures=0

# -06:
#    - add new arguments assocCut.hi and assocCut.lo
#    - make nNetworkFeatures equal nFeatures
#    = Bagging and boosting is broken.

# 05: add bagging

# 04: add a self-tuning version

# version 03: try to build a sample network predictor. 

# Cluster the samples in each class separately and identify clusters. It's a bit of a question whether we
# can automate the cluster identification completely. Anyway, use the clusters as additional centroids and
# as prediction use the class of each centroid.


# version 02: add the option to use a quantile of class distances instead of distance from centroid




# Inverse distance between colunms of x and y
# assume that y has few columns and compute the distance matrix between columns of x and columns of y.

.euclideanDist.forNCP = function(x, y, use = 'p') 
{ 
  x = as.matrix(x)
  y = as.matrix(y)

  ny = ncol(y)
  diff = matrix(NA, ncol(x), ny)
  
  for (cy in 1:ny)
    diff[, cy] = apply( (x - y[, cy])^2, 2, sum, na.rm = TRUE)

  -diff
}

#=====================================================================================================
#
# Main predictor function.
#
#=====================================================================================================
# best suited to prediction of factors.

# Classification: for each level find nFeatures.eachSide that best distinguish the level from all other
# levels. Actually this doesn't make much sense since I will have to put all the distinguishing sets
# together to form the profiles, so features that may have no relationship to a level will be added if
# there are more than two levels. I can fix that using some roundabout analysis but for now forget it. 

# To do: check that the CVfold validation split makes sense, i.e. none of the bins contains all observations
# of any class.

# Could also add robust standardization

# Output a measure of similarity to class centroids
# Same thing for the gene voting predictor

# prediction for heterogenous cases: sample network in training cases get modules and eigensamples and
# similarly in the controls then use all centroids for classification between k nearest neighbor and
# nearest centroid

# CAUTION: the function standardizes each gene (unless standardization is turned off), so the sample
# networks may be different from what would be expected from the supplied data. 

# Work internally with just numeric entries. Corresponding levels of the response are saved and restored at
# the very end. 



#' Nearest centroid predictor
#' 
#' Nearest centroid predictor for binary (i.e., two-outcome) data. Implements a
#' whole host of options and improvements such as accounting for within-class
#' heterogeneity using sample networks, various ways of feature selection and
#' weighing etc.
#' 
#' 
#' Nearest centroid predictor works by forming a representative profile
#' (centroid) across features for each class from the training data, then
#' assigning each test sample to the class of the nearest representative
#' profile. The representative profile can be formed either as mean or as athe
#' first principal component ("eigensample"; this choice is governed by the
#' option \code{centroidMethod}).
#' 
#' When the number of features is large and only a small fraction is likely to
#' be associated with the outcome, feature selection can be used to restrict
#' the features that actually enter the centroid. Feature selection can be
#' based either on their association with the outcome calculated from the
#' training data using \code{assocFnc}, or on user-supplied feature
#' significance (e.g., derived from literature, argument
#' \code{featureSignificance}). In either case, features can be selected by
#' high and low association tresholds or by taking a fixed number of highest-
#' and lowest-associated features.
#' 
#' As an alternative to centroids, the predictor can also assign test samples
#' based on a given quantile of the distances from the training samples in each
#' class (argument \code{useQuantile}). This may be advantageous if the samples
#' in each class form irregular clusters. Note that setting
#' \code{useQuantile=0} (i.e., using minimum distance in each class)
#' essentially gives a nearest neighbor predictor: each test sample will be
#' assigned to the class of its nearest training neighbor.
#' 
#' If features exhibit non-trivial correlations among themselves (such as, for
#' example, in gene expression data), one can attempt to down-weigh features
#' that do not exhibit the same correlation in the test set. This is done by
#' using essentially the same predictor to predict _features_ from all other
#' features in the test data (using the training data to train the feature
#' predictor). Because test features are known, the prediction accuracy can be
#' evaluated. If a feature is predicted badly (meaning the error in the test
#' set is much larger than the error in the cross-validation prediction in
#' training data), it may mean that its quality in the training or test data is
#' low (for example, due to excessive noise or outliers). Such features can be
#' downweighed using the argument \code{weighByPrediction}. The extra factor is
#' min(1, (root mean square prediction error in test set)/(root mean square
#' cross-validation prediction error in the trainig data)^weighByPrediction),
#' that is it is never bigger than 1.
#' 
#' Unless the features' mean and variance can be ascribed clear meaning, the
#' (training) features should be scaled to mean 0 and variance 1 before the
#' centroids are formed.
#' 
#' The function implements a basic option for removal of spurious effects in
#' the training and test data, by removng a fixed number of leading principal
#' components from the features. This sometimes leads to better prediction
#' accuracy but should be used with caution.
#' 
#' If samples within each class are heterogenous, a single centroid may not
#' represent each class well. This function can deal with within-class
#' heterogeneity by clustering samples (separately in each class), then using a
#' one representative (mean, eigensample) or quantile for each cluster in each
#' class to assign test samples. Various similarity measures, specified by
#' \code{adjFnc}, can be used to construct the sample network adjacency.
#' Similarly, the user can specify a clustering function using
#' \code{clusteringFnc}. The requirements on the clustering function are
#' described in a separate section below.
#' 
#' @param x Training features (predictive variables). Each column corresponds
#' to a feature and each row to an observation.
#' @param y The response variable. Can be a single vector or a matrix with
#' arbitrary many columns. Number of rows (observations) must equal to the
#' number of rows (observations) in x.
#' @param xtest Optional test set data. A matrix of the same number of columns
#' (i.e., features) as \code{x}. If test set data are not given, only the
#' prediction on training data will be returned.
#' @param featureSignificance Optional vector of feature significance for the
#' response variable. If given, it is used for feature selection (see details).
#' Should preferably be signed, that is features can have high negative
#' significance.
#' @param assocFnc Character string specifying the association function. The
#' association function should behave roughly as \code{link{cor}} in that it
#' takes two arguments (a matrix and a vector) plus options and returns the
#' vector of associations between the columns of the matrix and the vector. The
#' associations may be signed (i.e., negative or positive).
#' @param assocOptions Character string specifying options to the association
#' function.
#' @param assocCut.hi Association (or featureSignificance) threshold for
#' including features in the predictor. Features with associtation higher than
#' \code{assocCut.hi} will be included. If not given, the threshold method will
#' not be used; instead, a fixed number of features will be included as
#' specified by \code{nFeatures.hi} and \code{nFeatures.lo}.
#' @param assocCut.lo Association (or featureSignificance) threshold for
#' including features in the predictor. Features with associtation lower than
#' \code{assocCut.lo} will be included. If not given, defaults to
#' \code{-assocCut.hi}. If \code{assocCut.hi} is \code{NULL}, the threshold
#' method will not be used; instead, a fixed number of features will be
#' included as specified by \code{nFeatures.hi} and \code{nFeatures.lo}.
#' @param nFeatures.hi Number of highest-associated features (or features with
#' highest \code{featureSignificance}) to include in the predictor. Only used
#' if \code{assocCut.hi} is \code{NULL}.
#' @param nFeatures.lo Number of lowest-associated features (or features with
#' highest \code{featureSignificance}) to include in the predictor. Only used
#' if \code{assocCut.hi} is \code{NULL}.
#' @param weighFeaturesByAssociation (Optional) power to downweigh features
#' that are less associated with the response. See details.
#' @param scaleFeatureMean Logical: should the training features be scaled to
#' mean zero? Unless there are good reasons not to scale, the features should
#' be scaled.
#' @param scaleFeatureVar Logical: should the training features be scaled to
#' unit variance? Again, unless there are good reasons not to scale, the
#' features should be scaled.
#' @param centroidMethod One of \code{"mean"} and \code{"eigensample"},
#' specifies how the centroid should be calculated. \code{"mean"} takes the
#' mean across all samples (or all samples within a sample module, if sample
#' networks are used), whereas \code{"eigensample"} calculates the first
#' principal component of the feature matrix and uses that as the centroid.
#' @param simFnc Character string giving the similarity function for measuring
#' the similarity between test samples and centroids. This function should
#' behave roughly like the function \code{\link{cor}} in that it takes two
#' arguments (\code{x}, \code{y}) and calculates the pair-wise similarities
#' between columns of \code{x} and \code{y}. For convenience, the value
#' \code{"dist"} is treated specially: the Euclidean distance between the
#' columns of \code{x} and \code{y} is calculated and its negative is returned
#' (so that smallest distance corresponds to highest similarity). Since values
#' of this function are only used for ranking centroids, its values are not
#' restricted to be positive or within certain bounds.
#' @param simOptions Character string specifying the options to the similarity
#' function.
#' @param useQuantile If non-NULL, the "nearest quantiloid" will be used
#' instead of the nearest centroid. See details.
#' @param sampleWeights Optional specification of sample weights. Useful for
#' example if one wants to explore boosting.
#' @param weighSimByPrediction (Optional) power to downweigh features that are
#' not well predicted between training and test sets. See details.
#' @param CVfold Non-negative integer specifying cross-validation. Zero means
#' no cross-validation will be performed. values above zero specify the number
#' of samples to be considered test data for each step of cross-validation.
#' @param returnFactor Logical: should a factor be returned?
#' @param randomSeed Integere specifying the seed for the random number
#' generator. If \code{NULL}, the seed will not be set. See
#' \code{\link{set.seed}}.
#' @param verbose Integer controling how verbose the diagnostic messages should
#' be. Zero means silent.
#' @param indent Indentation for the diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#' 
#' \item{predicted}{The back-substitution prediction in the training set.}
#' 
#' \item{predictedTest}{Prediction in the test set.}
#' 
#' \item{featureSignificance}{A vector of feature significance calculated by
#' \code{assocFnc} or a copy of the input \code{featureSignificance} if the
#' latter is non-NULL.}
#' 
#' \item{selectedFeatures}{A vector giving the indices of the features that
#' were selected for the predictor.}
#' 
#' \item{centroidProfile}{The representative profiles of each class (or
#' cluster). Only returned in \code{useQuntile} is \code{NULL}. }
#' 
#' \item{testSample2centroidSimilarities}{A matrix of calculated similarities
#' between the test samples and class/cluster centroids.}
#' 
#' \item{featureValidationWeights}{A vector of validation weights (see Details)
#' for the selected features. If \code{weighFeaturesByValidation} is 0, a unit
#' vector is used and returned.}
#' 
#' \item{CVpredicted}{Cross-validation prediction on the training data. Present
#' only if \code{CVfold} is non-zero.}
#' 
#' \item{sampleClusterLabels}{A list with two components (one per class). Each
#' component is a vector of sample cluster labels for samples in the class.}
#' @author Peter Langfelder
#' @seealso \code{\link{votingLinearPredictor}}
#' @keywords misc
nearestCentroidPredictor = function(
    # Input training and test data
    x, y, 
    xtest = NULL,

    # Feature weights and selection criteria
    featureSignificance = NULL,
    assocFnc = "cor", assocOptions = "use = 'p'",
    assocCut.hi = NULL, assocCut.lo = NULL,
    nFeatures.hi = 10, nFeatures.lo = 10,
    weighFeaturesByAssociation = 0,
    scaleFeatureMean = TRUE, scaleFeatureVar = TRUE,

    # Predictor options 
    centroidMethod = c("mean", "eigensample"),
    simFnc = "cor", simOptions = "use = 'p'",
    useQuantile = NULL,
    sampleWeights = NULL,
    weighSimByPrediction = 0,

    # What should be returned
    CVfold = 0, returnFactor = FALSE,

    # General options
    randomSeed = 12345,
    verbose = 2, indent = 0)
{

  # For now we do not support regression

  centroidMethod = match.arg(centroidMethod)

  if (simFnc=="dist")
  {
    if (verbose > 0)
      printFlush(paste("NearestCentroidPredictor: 'dist' is changed to a suitable", 
                       "Euclidean distance function.\n",
                       "   Note: simOptions will be disregarded."))
    simFnc = ".euclideanDist.forNCP"
    simOptions = "use = 'p'"
  }

  # Convert factors to numeric variables.
  
  ySaved = y
  #if (classify) 
  #{ 
    originalYLevels = sort(unique(y)) 
    y = as.numeric(as.factor(y))
  #}

  x = as.matrix(x)
  doTest = !is.null(xtest)

  if (doTest) 
  {
    xtest = as.matrix(xtest)
    nTestSamples = nrow(xtest)
    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.")
  } else {
    if (weighSimByPrediction > 0)
       stop("weighting similarity by prediction is not possible when xtest = NULL.")
    nTestSamples = 0
  }

  numYLevels = sort(unique(y))

  minY = min(y)
  maxY = max(y)

  nSamples = length(y)
  nVars = ncol(x)

  if (!is.null(assocCut.hi))
  {
    if (is.null(assocCut.lo)) assocCut.lo = -assocCut.hi
  }

  spaces = indentSpaces(indent)

  if (!is.null(useQuantile)) 
  {
     if ( (useQuantile < 0) | (useQuantile > 1) )
      stop("If 'useQuantile' is given, it must be between 0 and 1.")
  }

  if (is.null(sampleWeights)) sampleWeights = rep(1, nSamples)

  # If cross-validation is requested, change the whole flow and use a recursive call.
  if (CVfold > 0)
  {
    if (CVfold > nSamples )
    {
      printFlush("'CVfold' is larger than number of samples. Will perform leave-one-out cross-validation.")
      CVfold = nSamples
    }
    ratio = nSamples/CVfold
    if (floor(ratio)!=ratio)
    {
      smaller = floor(ratio)
      nLarger = nSamples - CVfold * smaller
      binSizes = c(rep(smaller, CVfold-nLarger), rep(smaller +1, nLarger))
    } else
      binSizes = rep(ratio, CVfold)
    if (!is.null(randomSeed))
    {
      if (exists(".Random.seed"))
      {
        saved.seed = .Random.seed
        seedSaved = TRUE
      } else
        seedSaved = FALSE
      set.seed(randomSeed)
    }

    sampleOrder = sample(1:nSamples)
      
    CVpredicted = rep(NA, nSamples)
    CVbin = rep(0, nSamples)

    if (verbose > 0) 
    {
      cat(paste(spaces, "Running cross-validation: "))
      if (verbose==1) pind = initProgInd() else printFlush("")
    }

    if (!is.null(featureSignificance)) 
      printFlush(paste("Warning in nearestCentroidPredictor: \n", 
                       "   cross-validation will be biased if featureSignificance was derived", 
                       "from training data."))

    ind = 1
    for (cv in 1:CVfold)
    {
      if (verbose > 1) printFlush(paste("..cross validation bin", cv, "of", CVfold))
      end = ind + binSizes[cv] - 1
      samples = sampleOrder[ind:end]
      CVbin[samples] = cv
      xCVtrain = x[-samples, , drop = FALSE]
      xCVtest = x[samples, , drop = FALSE]
      yCVtrain = y[-samples]
      yCVtest = y[samples]
      CVsampleWeights = sampleWeights[-samples]
      pr = nearestCentroidPredictor(xCVtrain, yCVtrain, xCVtest,
                        #classify = classify,
                        featureSignificance = featureSignificance, 
                        assocCut.hi = assocCut.hi,
                        assocCut.lo = assocCut.lo,
                        nFeatures.hi = nFeatures.hi,
                        nFeatures.lo = nFeatures.lo,
                        useQuantile = useQuantile, 
                        sampleWeights = CVsampleWeights,

                        CVfold = 0, returnFactor = FALSE,
                        randomSeed = randomSeed,
                        centroidMethod = centroidMethod,
                        assocFnc = assocFnc, assocOptions = assocOptions,
                        scaleFeatureMean = scaleFeatureMean,
                        scaleFeatureVar = scaleFeatureVar,
                        simFnc = simFnc, simOptions = simOptions,
                        weighFeaturesByAssociation = weighFeaturesByAssociation,
                        weighSimByPrediction = weighSimByPrediction,
                        verbose = verbose - 2, indent = indent + 1)

      CVpredicted[samples] = pr$predictedTest
      ind = end + 1
      if (verbose==1) pind = updateProgInd(cv/CVfold, pind)
    }
    if (verbose==1) printFlush("")
  }

  if (nrow(x)!=length(y))
    stop("Number of observations in x and y must equal.")

  # Feature selection:

  xWeighted = x * sampleWeights
  yWeighted = y * sampleWeights

  if (is.null(featureSignificance))
  {
    corEval = parse(text = paste(assocFnc, "(xWeighted, yWeighted, ", assocOptions, ")"))
    featureSignificance = as.vector(eval(corEval))
  } else {
    if (length(featureSignificance)!=nVars)
      stop("Given 'featureSignificance' has incorrect length (must be nFeatures).")
  }
  nGood = nVars
  nNA = sum(is.na(featureSignificance))
  testCentroidSimilarities = list()
  xSD = apply(x, 2, sd, na.rm = TRUE)
  keep = is.finite(featureSignificance) & (xSD>0)
  nKeep = sum(keep)
  keepInd = c(1:nVars)[keep]
  order = order(featureSignificance[keep])
  levels = sort(unique(y))
  nLevels = length(levels)
  if (is.null(assocCut.hi))
  {
    nf = c(nFeatures.hi, nFeatures.lo)
    if (nf[2] > 0) ind1 = c(1:nf[2]) else ind1 = c()
    if (nf[1] > 0) ind2 = c((nKeep-nf[1] + 1):nKeep) else ind2 = c()
    indexSelect = unique(c(ind1, ind2))
    if (length(indexSelect) < 1) 
      stop("No features were selected. At least one of 'nFeatures.hi', 'nFeatures.lo' must be nonzero.")
    indexSelect = indexSelect[indexSelect > 0]
    select = keepInd[order[indexSelect]]
  } else {
    indexSelect = (1:nKeep)[ featureSignificance[keep] >= assocCut.hi |
                               featureSignificance[keep] <= assocCut.lo ]
    if (length(indexSelect)<2) 
       stop(paste("'assocCut.hi'", assocCut.hi, "and assocCut.lo", assocCut.lo, 
                  "are too stringent, less than 3 features were selected.\n",
                  "Please relax the cutoffs.")); 
    select = keepInd[indexSelect]
  }
  if ((length(select) < 3) && (simFnc!='dist')) 
  {
       stop(paste("Less than 3 features were selected. Please either relax", 
                  "the selection criteria of use simFnc = 'dist'."))
  }

  selectedFeatures = select
  nSelect = length(select)

  xSel = x[, select]
  selectSignif = featureSignificance[select]

  if (scaleFeatureMean)
  {
    if (scaleFeatureVar)
    {
      xSD = apply(xSel, 2, sd, na.rm = TRUE)
    } else 
      xSD = rep(1, nSelect)
    xMean = apply(xSel, 2, mean, na.rm = TRUE)
  } else {
    if (scaleFeatureVar)
    {
      xSD = sqrt(apply(xSel^2, 2, sum, na.rm = TRUE)) / pmax(apply(!is.na(xSel), 2, sum) - 1, rep(1, nSelect))
    } else
      xSD = rep(1, nSelect)
    xMean = rep(0, nSelect)
  }
  xSel = scale(xSel, center = scaleFeatureMean, scale = scaleFeatureVar)
  if (doTest)
  {
    xtestSel = xtest[, select]
    xtestSel = (xtestSel - matrix(xMean, nTestSamples, nSelect, byrow = TRUE) ) / 
             matrix(xSD, nTestSamples, nSelect, byrow = TRUE)
  } else
    xtestSel = NULL

  xWeighted = xSel * sampleWeights

  if (weighSimByPrediction > 0)
  {
    pr = .quickGeneVotingPredictor.CV(xSel, xtestSel, c(1:nSelect))
    dCV = sqrt(colMeans( (pr$CVpredicted - xSel)^2, na.rm = TRUE))
    dTS = sqrt(colMeans( (pr$predictedTest - xtestSel)^2, na.rm = TRUE))
    dTS[dTS==0] = min(dTS[dTS>0])
    validationWeight = (dCV/dTS)^weighSimByPrediction
    validationWeight[validationWeight > 1] = 1
  } else
    validationWeight = rep(1, nSelect)

  nTestSamples = if (doTest) nrow(xtest) else 0

  predicted = rep(0, nSamples)
  predictedTest = rep(0, nTestSamples)

  clusterLabels = list()
  clusterNumbers = list()

  if ( (centroidMethod=="eigensample") )
  {
    if (sum(is.na(xSel)) > 0)
    { 
       xImp = t(impute.knn(t(xSel), k = min(10, nSelect - 1))$data)
    } else 
       xImp = xSel
    if (doTest && sum(is.na(xtestSel))>0)
    {
       xtestImp = t(impute.knn(t(xtestSel), k = min(10, nSelect - 1))$data)
    } else
       xtestImp = xtestSel
  }

  clusterNumbers = rep(1, nLevels)
  sampleModules = list()
  
  # Trivial cluster labels: clusters equal case classes
  for (l in 1:nLevels)
    clusterLabels[[l]] = rep(l, sum(y==levels[l]))

  nClusters = sum(clusterNumbers)
  centroidSimilarities = array(NA, dim = c(nSamples, nClusters))
  testCentroidSimilarities = array(NA, dim = c(nTestSamples, nClusters))
  #if (classify)
  #{
    cluster2level = rep(c(1:nLevels), clusterNumbers)
    featureWeight = validationWeight; 
    if (is.null(useQuantile))
    {
      # Form centroid profiles for each cluster and class
      centroidProfiles = array(0, dim = c(nSelect, nClusters))
      for (cl in 1:nClusters)
      {
        l = cluster2level[cl]
        clusterSamples = c(1:nSamples)[ y==l ] [ clusterLabels[[l]]==cl ]
        if (centroidMethod=="mean")
        {
          centroidProfiles[, cl] = apply(xSel[clusterSamples, , drop = FALSE], 
                                                    2, mean, na.rm = TRUE)
        } else if (centroidMethod=="eigensample")
        {
          cp = svd(xSel[clusterSamples,], nu = 0, nv = 1)$v[, 1]
          cor = cor(t(xSel[clusterSamples,]), cp)
          if (sum(cor, na.rm = TRUE) < 0) cp = -cp
          centroidProfiles[, cl] = cp
        }
      }
      if (weighFeaturesByAssociation > 0)
        featureWeight = featureWeight * sqrt(abs(selectSignif)^weighFeaturesByAssociation)
      # Back-substitution prediction: 
      wcps = centroidProfiles * featureWeight
      wxSel = t(xSel) * featureWeight
      distExpr = paste0( simFnc, "( wcps, wxSel, ", simOptions, ")")
      sample.centroidSim = eval(parse(text = distExpr))

      # Actual prediction: for each sample, calculate distances to centroid profiles
      if (doTest)
      {
        wxtestSel = t(xtestSel) * featureWeight
        distExpr = paste0( simFnc, "( wcps, wxtestSel, ", simOptions, ")")
        testSample.centroidSim = eval(parse(text = distExpr))
      }
    } else {
      labelVector = y
      for (l in 1:nLevels)
        labelVector[y==l] = clusterLabels[[l]]
      keepSamples = labelVector!=0
      nKeepSamples = sum(keepSamples)
      keepLabels = labelVector[keepSamples]
      if (weighFeaturesByAssociation > 0)
        featureWeight = featureWeight * sqrt(abs(selectSignif)^weighFeaturesByAssociation)
      wxSel = t(xSel) * featureWeight
      wxSel.keepSamples = t(xSel[keepSamples, ]) * featureWeight
  
      # Back-substitution prediction: 
      
      distExpr = paste0( simFnc, "( wxSel.keepSamples, wxSel, ", simOptions, ")")
      dst = eval(parse(text = distExpr))
      # Test prediction:
      if (doTest)
      {
        wxtestSel = t(xtestSel) * featureWeight
        distExpr = paste0( simFnc, "( wxSel.keepSamples, wxtestSel, ", simOptions, ")")
        dst.test = eval(parse(text = distExpr))
        sample.centroidSim = matrix(0, nClusters, nSamples)
        testSample.centroidSim = matrix(0, nClusters, nTestSamples)
      }
      for (l in 1:nClusters)
      {
         #x = try ( {
         lSamples = c(1:nKeepSamples)[keepLabels==l]
         sample.centroidSim[l, ] = colQuantileC(dst[lSamples, ], 1-useQuantile)
         testSample.centroidSim[l, ] = colQuantileC(dst.test[lSamples, ], 1-useQuantile)
         #} )
         #if (class(x) == 'try-error') browser(text = "zastavka.")
      }
    }
   centroidSimilarities = t(sample.centroidSim)
   prediction = cluster2level[apply(sample.centroidSim, 2, which.max)]
   # Save predictions
   predicted = prediction
   if (doTest)
   {
     testCentroidSimilarities = t(testSample.centroidSim)
     testprediction = cluster2level[apply(testSample.centroidSim, 2, which.max)]
     predictedTest = testprediction
   }
  #} else 
  #  stop("Prediction for continouos variables is not implemented yet. Sorry!")

  # Reformat output if factors are to be returned
    
  if (returnFactor)
  {
      predicted.out = factor(originalYLevels[[t]][predicted])
      if (doTest) predictedTest.out = factor(originalYLevels[[t]][predictedTest])
      if (CVfold > 0)
        CVpredicted.out = factor(originalYLevels[[t]][CVpredicted])
  } else {
    # Turn ordinal predictions into levels of input traits
    predicted.out = originalYLevels[predicted]
    if (doTest) predictedTest.out = originalYLevels[predictedTest]
    if (CVfold > 0)
      CVpredicted.out = originalYLevels[CVpredicted]
  }
  out = list(predicted = predicted.out,
             predictedTest = if (doTest) predictedTest.out else NULL,
             featureSignificance = featureSignificance, 
             selectedFeatures = selectedFeatures,
             centroidProfiles = if (is.null(useQuantile)) centroidProfiles else NULL,
             testSample2centroidSimilarities = if (doTest) testCentroidSimilarities else NULL,
             featureValidationWeights = validationWeight
             )
  if (CVfold > 0)
    out$CVpredicted = CVpredicted.out

  out
}


