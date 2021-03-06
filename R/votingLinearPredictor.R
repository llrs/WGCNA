.devianceResidual = function(y)
{
  event = y[, ncol(y) ]
  fit = summary(coxph(y~1, na.action = na.exclude))
  CumHazard = predict(fit, type = "expected")
  martingale1 = event - CumHazard
  deviance0 = ifelse(event == 0, 2 * CumHazard, -2 * log(CumHazard) +
                     2 * CumHazard - 2)
  sign(martingale1) * sqrt(deviance0)
}


.dropThirdDim = function(x)
{
  d = dim(x)
  dim(x) = c(d[1], d[2]*d[3])
  x
}
 

#-----------------------------------------------------------------------------------------------------
#
# Voting linear predictor: given a set of vectors y and a set of vectors x, 
# predictor is given by the correlations of y with x[,i] 
# 
#-----------------------------------------------------------------------------------------------------





#' Voting linear predictor
#' 
#' Predictor based on univariate regression on all or selected given features
#' that pools all predictions using weights derived from the univariate linear
#' models.
#' 
#' The predictor calculates the association of each (selected) feature with the
#' response and uses the association to calculate the weight of the feature as
#' \code{sign(association) * (association)^featureWeightPower}. Optionally,
#' this weight is multiplied by \code{priorWeights}. Further, a feature
#' prediction weight can be used to downweigh features that are not well
#' predicted by other features (see below).
#' 
#' For classification, the (continuous) result of the above calculation is
#' turned into ordinal values essentially by rounding.
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
#' low (for example, due to excessive noise or outliers).  Such features can be
#' downweighed using the argument \code{weighByPrediction}. The extra factor is
#' min(1, (root mean square prediction error in test set)/(root mean square
#' cross-validation prediction error in the trainig data)^weighByPrediction),
#' that is it is never bigger than 1.
#' 
#' @param x Training features (predictive variables). Each column corresponds
#' to a feature and each row to an observation.
#' @param y The response variable. Can be a single vector or a matrix with
#' arbitrary many columns. Number of rows (observations) must equal to the
#' number of rows (observations) in x.
#' @param xtest Optional test set data. A matrix of the same number of columns
#' (i.e., features) as \code{x}.  If test set data are not given, only the
#' prediction on training data will be returned.
#' @param classify Should the response be treated as a categorical variable?
#' Classification really only works with two classes. (The function will run
#' for multiclass problems as well, but the results will be sub-optimal.)
#' @param CVfold Optional specification of cross-validation fold. If 0 (the
#' default), no cross-validation is performed.
#' @param randomSeed Random seed, used for observation selection for
#' cross-validation. If \code{NULL}, the random generator is not reset.
#' @param assocFnc Function to measure association. Usually a measure of
#' correlation, for example Pearson correlation or \code{\link{bicor}}.
#' @param assocOptions Character string specifying the options to be passed to
#' the association function.
#' @param featureWeightPowers Powers to which to raise the result of
#' \code{assocFnc} to obtain weights. Can be a single number or a vector of
#' arbitrary length; the returned value will contain one prediction per power.
#' @param priorWeights Prior weights for the features. If given, must be either
#' (1) a vector of the same length as the number of features (columns in
#' \code{x}); (2) a matrix of dimensions length(featureWeightPowers)x(number of
#' features); or (3) array of dimensions (number of response
#' variables)xlength(featureWeightPowers)x(number of features).
#' @param weighByPrediction (Optional) power to downweigh features that are not
#' well predicted between training and test sets. See details.
#' @param nFeatures.hi Optional restriction of the number of features to use.
#' If given, this many features with the highest association and lowest
#' association (if \code{nFeatures.lo} is not given) will be used for
#' prediction.
#' @param nFeatures.lo Optional restriction of the number of lowest (i.e., most
#' negatively) associated features to use.  Only used if \code{nFeatures.hi} is
#' also non-NULL.
#' @param dropUnusedDimensions Logical: should unused dimensions be dropped
#' from the result?
#' @param verbose Integer controling how verbose the diagnostic messages should
#' be. Zero means silent.
#' @param indent Indentation for the diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{predicted}{The
#' back-substitution prediction on the training data. Normally an array of
#' dimensions (number of observations) x (number of response variables) x
#' length(featureWeightPowers), but unused are dropped unless
#' \code{dropUnusedDimensions = FALSE}.}
#' 
#' \item{weightBase}{Absolute value of the associations of each feature with
#' each response.}
#' 
#' \item{variableImportance}{The weight of each feature in the prediction
#' (including the sign).}
#' 
#' \item{predictedTest}{If input \code{xtest} is non-NULL, the predicted test
#' response, in format analogous to \code{predicted} above.}
#' 
#' \item{CVpredicted}{If input \code{CVfold} is non-zero, cross-validation
#' prediction on the training data.}
#' @note It makes little practical sense to supply neither \code{xtest} nor
#' \code{CVfold} since the prediction accuracy on training data will be highly
#' biased.
#' @author Peter Langfelder
#' @seealso \code{\link{bicor}} for robust correlation that can be used as an
#' association measure
#' @keywords misc
votingLinearPredictor = function(x, y, xtest = NULL, 
                    classify = FALSE, 
                    CVfold = 0,
                    randomSeed = 12345,
                    assocFnc = "cor", assocOptions = "use = 'p'",
                    featureWeightPowers = NULL, priorWeights = NULL,
                    weighByPrediction = 0,
                    nFeatures.hi = NULL,
                    nFeatures.lo = NULL,
                    dropUnusedDimensions = TRUE,
                    verbose = 2, indent = 0)
{

  # Special handling of a survival response

  if (is.Surv(y))
  {
    pr = votingLinearPredictor(x, .devianceResidual(y), xtest, classify = FALSE,
                             CVfold = CVfold, randomSeed = randomSeed,
                             assocFnc = assocFnc, assocOptions = assocOptions,
                             featureWeightPowers = featureWeightPowers, priorWeights = priorWeights,
                             nFeatures.hi = nFeatures.hi, nFeatures.lo = nFeatures.lo,
                             dropUnusedDimensions = dropUnusedDimensions,
                             verbose = verbose, indent = indent)
    # may possibly need completion here. 
    return(pr)
  }       

  # Standard code for a numeric response

  ySaved = y

  if (classify)
  { 
    # Convert factors to numeric variables.
    if (is.null(dim(y)))
    {
      ySaved = list(y)
      if (classify)
      {
         originalYLevels = list( sort(unique(y)) )
         y = as.factor(y)
      }
      y = as.matrix(as.numeric(y))
    } else {
      ySaved = y
      if (classify)
      {
         y = as.data.frame(lapply(as.data.frame(y), as.factor))
         originalYLevels = lapply(as.data.frame(apply(ySaved, 2, unique)), sort)
      } else
         y = as.data.frame(y)
      y = as.matrix(sapply(lapply(y, as.numeric), I))
    }
    numYLevels = lapply(as.data.frame(apply(y, 2, unique)), sort)
    minY = apply(y, 2, min)
    maxY = apply(y, 2, max)
  } else 
    y = as.matrix(y)

  doTest = !is.null(xtest)

  x = as.matrix(x)
  nSamples = nrow(y)

  nTraits = ncol(y)
  nVars = ncol(x)

  if (is.null(featureWeightPowers)) featureWeightPowers = 0
  nPredWPowers = length(featureWeightPowers)
    
  if (is.null(rownames(x)))
  {
    sampleNames = paste0("Sample.", c(1:nSamples))
  } else
    sampleNames = rownames(x)
   
  if (doTest) 
  { 
    xtest = as.matrix(xtest)
    nTestSamples = nrow(xtest)
    if (is.null(rownames(xtest)))
    {
      testSampleNames = paste0("testSample.", c(1:nTestSamples))
    } else
      testSampleNames = rownames(xtest)
    if (ncol(x)!=ncol(xtest))
      stop("Number of learning and testing predictors (columns of x, xtest) must equal.")
  }

  if (is.null(colnames(y)))
  {
    traitNames = paste0("Response.", c(1:nTraits))
  } else
    traitNames = colnames(y)

  if (is.null(colnames(x)))
  {
    featureNames = paste0("Feature.", c(1:nVars))
  } else
    featureNames = colnames(x)

  powerNames = paste0("Power.", featureWeightPowers)
   

  spaces = indentSpaces(indent)
 
  # If cross-validation is requested, change the whole flow and use a recursive call.
  if (CVfold > 0)
  {
    if (CVfold > nSamples )
    {
      printFlush("CVfold is larger than nSamples. Will perform leave-one-out cross-validation.")
      CVfold = nSamples
    }
    ratio = nSamples/CVfold
    if (floor(ratio)!=ratio)
    {
      smaller =floor(ratio)
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
      
    CVpredicted = array(NA, dim = c(nSamples, nTraits, nPredWPowers))
    CVbin = rep(0, nSamples)

    if (verbose > 0) 
    {
      cat(paste(spaces, "Running cross-validation: "))
      if (verbose==1) pind = initProgInd() else printFlush("")
    }

    ind = 1
    for (cv in 1:CVfold)
    {
      if (verbose > 1) printFlush(paste("..cross validation bin", cv, "of", CVfold))
      end = ind + binSizes[cv] - 1
      samples = sampleOrder[ind:end]
      CVbin[samples] = cv
      xCVtrain = x[-samples, , drop = FALSE]
      xCVtest = x[samples, , drop = FALSE]
      yCVtrain = y[-samples, , drop = FALSE]
      yCVtest = y[samples,, drop = FALSE ]
      pr = votingLinearPredictor(xCVtrain, yCVtrain, xCVtest,
                        classify = FALSE,
                        CVfold = 0,
                        assocFnc = assocFnc, assocOptions = assocOptions, 
                        featureWeightPowers = featureWeightPowers,
                        priorWeights = priorWeights, 
                        nFeatures.hi = nFeatures.hi, nFeatures.lo = nFeatures.lo,
                        dropUnusedDimensions = dropUnusedDimensions,
                        verbose = verbose - 1, indent = indent + 1)
      CVpredicted[samples, , ] = pr$predictedTest
      ind = end + 1
      if (verbose==1) pind = updateProgInd(cv/CVfold, pind)
    }
    if (verbose==1) printFlush("")

    collectGarbage()
  }

  if (nrow(x)!=nrow(y))
    stop("Number of observations in x and y must equal.")

  xSD = apply(x, 2, sd, na.rm = TRUE)

  validFeatures = xSD > 0

  xMean = apply(x, 2, mean, na.rm = TRUE)
  x = scale(x)
  if (doTest) 
  {
     xtest = (xtest - matrix(xMean, nTestSamples, nVars, byrow = TRUE) ) / 
               matrix(xSD, nTestSamples, nVars, byrow = TRUE)
     xtest[, !validFeatures] = 0
  }

  # This prevents NA's generated from zero xSD to contaminate the results
  x[, !validFeatures] = 0

  xSD[!validFeatures] = 1

  obsMean = apply(y, 2, mean, na.rm = TRUE)
  obsSD = apply(y, 2, sd, na.rm = TRUE)

  if (sum(!is.finite(obsSD)) > 0)
    stop("Something is wrong with given trait: not all standard deviations are finite.")

  if (sum(obsSD==0) > 0)
    stop("Some of given traits have variance zero. Prediction on such traits will not work.")

  y = scale(y)

  if (is.null(priorWeights))
  {
    priorWeights = array(1, dim = c(nTraits, nPredWPowers, nVars))
  } else {
    dimPW = dim(priorWeights)
    if (length(dimPW)<=1) 
    {
      if (length(priorWeights!=nVars))
         stop ("priorWeights are a vector - must have length = number of variables (ncol(x)).")
      priorWeights = matrix(priorWeights, nrow = nSamples * nPredWPowers, ncol = nVars, 
                            byrow = TRUE)
      dim(priorWeights) = c(nTraits, nPredWPowers, nVars)
    } else if (length(dimPW)==2)
    {
      if ((dimPW[1]!=nPredWPowers) | (dimPW[2]!=nVars) )
        stop(paste("If priorWeights is two-dimensional, 1st dimension must equal",
                   "number of predictor weight powers, 2nd must equal number of variables"))
      # if (verbose>0) printFlush("..converting dimensions of priorWeights..")
      newWeights = array(0, dim = c(nTraits, nPredWPowers, nVars))
      for (trait in 1:nTraits)
        newWeights[trait, , ] = priorWeights[,]
      priorWeights = newWeights
      collectGarbage()
    } else if ( (dimPW[1] != nTraits) | (dimPW[3] != nVars) | (dimPW[2] != nPredWPowers) )
        stop(paste("priorWeights have incorrect dimensions. Dimensions must be",
                   "(ncol(y), length(featureWeightPowers), ncol(x))."))
  }
      
  varImportance = array(0, dim = c(nTraits, nPredWPowers, nVars))
  predictMat = array(0, dim = c(nSamples, nTraits, nPredWPowers))
  if (doTest) predictTestMat = array(0, dim = c(nTestSamples, nTraits, nPredWPowers))

  corEval = parse(text = paste(assocFnc, "(x, y ", prepComma(assocOptions), ")"))
  r = eval(corEval)
  if (!is.null(nFeatures.hi))
  {
    if (is.null(nFeatures.lo)) nFeatures.lo = nFeatures.hi
    # Zero out associations (and therefore weights) for features that do not make the cut
    rank = apply(r, 2, rank, na.last = TRUE)
    nFinite = colSums(!is.na(r))
    for (t in 1:nTraits)
      r[ rank[, t]>nFeatures.lo & rank[, t] <= nFinite[t]-nFeatures.hi, t] = 0
  }
  r[is.na(r)] = 0; 

  validationWeights = rep(1, nVars)
  if (weighByPrediction > 0)
  {
    select = c(1:nVars)[rowSums(r!=0) > 0]
    nSelect = length(select)
    pr = .quickGeneVotingPredictor.CV(x[, select], xtest[, select], c(1:nSelect))
    dCV = sqrt(colMeans( (pr$CVpredicted - x[, select])^2, na.rm = TRUE))
    dTS = sqrt(colMeans( (pr$predictedTest - xtest[, select])^2, na.rm = TRUE))
    dTS[dTS==0] = min(dTS[dTS>0])
    w = (dCV/dTS)^weighByPrediction
    w[w>1] = 1
    validationWeights[select] = w; 
  } 

  finiteX = (is.finite(x) + 1)-1
  x.fin = x
  x.fin[!is.finite(x)] = 0

  if (doTest)
  {
    finiteXTest =  (is.finite(xtest) + 1)-1
    xtest.fin = xtest
    xtest.fin[!is.finite(xtest)] = 0
  }

  for (power in 1:nPredWPowers)
  {
    prWM = priorWeights[, power, ]
    dim(prWM) = c(nTraits, nVars)

    weights = abs(r)^featureWeightPowers[power] * t(prWM) * validationWeights
    #weightSum = apply(weights, 2, sum)
    weightSum = finiteX %*% weights
    RWeights = sign(r) * weights

    predictMat[ , , power] = 
         x.fin %*% RWeights / weightSum
    predMean = apply(.dropThirdDim(predictMat[ , , power, drop = FALSE]), 2, mean, na.rm = TRUE)
    predSD = apply(.dropThirdDim(predictMat[, , power, drop = FALSE]), 2, sd, na.rm = TRUE)
    predictMat[, , power] = scale(predictMat[, , power]) *
                      matrix(obsSD, nrow = nSamples, ncol = nTraits, byrow = TRUE) +
                      matrix(obsMean, nrow = nSamples, ncol = nTraits, byrow = TRUE)
    if (doTest) 
    {
      weightSum.test = finiteXTest %*% weights
      predictTestMat[, , power] = 
         (xtest.fin %*% RWeights /  weightSum.test -
                  matrix(predMean, nrow = nTestSamples, ncol = nTraits, byrow = TRUE) ) /
                  matrix(predSD, nrow = nTestSamples, ncol = nTraits, byrow = TRUE) * 
                  matrix(obsSD, nrow = nTestSamples, ncol = nTraits, byrow = TRUE) +
                  matrix(obsMean, nrow = nTestSamples, ncol = nTraits, byrow = TRUE); 
    }
    varImportance[ , power, ] = RWeights

  }

  dimnames(predictMat) = list(sampleNames, traitNames, powerNames)
  if (doTest) dimnames(predictTestMat) = list(testSampleNames, traitNames, powerNames)

  dimnames(r) = list(featureNames, traitNames)
  dimnames(varImportance) = list(traitNames, powerNames, featureNames)

  discretize = function(x, min, max)
  {
      fac = round(x)
      fac[fac < min] = min
      fac[fac > max] = max
      # dim(fac) = dim(x)
      fac
  }

  trafo = function(x, drop)
  {
    if (classify)
    {
      for (power in 1:nPredWPowers) for (t in 1:nTraits)
      {
        disc = discretize(x[, t, power], minY[t], maxY[t])
        x[, t, power] = originalYLevels[[t]] [disc]
      }
    } 
    x = x[ , , , drop = drop]
    x
  }

  out = list(predicted = trafo(predictMat, drop = dropUnusedDimensions), 
             weightBase = abs(r), 
             variableImportance = varImportance)
  if (doTest) out$predictedTest = trafo(predictTestMat, drop = dropUnusedDimensions)
  if (CVfold > 0)
  {
    dimnames(CVpredicted) = list(sampleNames, traitNames, powerNames)
    CVpredFac = trafo(CVpredicted, drop = dropUnusedDimensions)
    out$CVpredicted = CVpredFac
  } 
  out
}




#=======================================================================================================
#
# Quick linear predictor for genes. 
#
#=======================================================================================================

# Assume we have expression data (training and test). Run the prediction
# in a CV-like way. Keep the predictions for CV samples and for the test samples; at the end average the
# test predictions to get one final test prediction. 
 
# In each CV run: 
# scale the training data and keep scale
# scale the test data using the training scale
# get correlations of predicted vectors and predictors (note: number of predicted vectors should be
# relatively small, but all genes may be predictors - however, here we can also make some restrictions)
# Form the ensemble of predictors for each gene
# Predict the CV samples and test samples

.quickGeneVotingPredictor = function(x, xtest, predictedIndex, nPredictorGenes = 20, power = 3,
                         corFnc = "bicor", corOptions = "use = 'p'", verbose = 0)
{
  nSamples = nrow(x)
  nTestSamples = nrow(xtest)
  nGenes = ncol(x)
  nPredicted = length(predictedIndex)
  if (nPredictorGenes >=nGenes) nPredictorGenes = nGenes -1

  geneMeans = as.numeric(colMeans(x, na.rm = TRUE))
  geneSD = as.numeric(sqrt(colMeans(x^2, na.rm = TRUE) - geneMeans^2))
  
  xScaled = (x-matrix(geneMeans, nSamples, nGenes, byrow = TRUE))/
                 matrix(geneSD, nSamples, nGenes, byrow = TRUE)

  xTestScaled = (xtest-matrix(geneMeans, nTestSamples, nGenes, byrow = TRUE))/
                 matrix(geneSD, nTestSamples, nGenes, byrow = TRUE)
  corExpr = parse(text = paste(corFnc, " (xScaled, xScaled[, predictedIndex] ", prepComma(corOptions), ")"))
  significance = eval(corExpr)

  predictedTest = matrix(NA, nTestSamples, nPredicted)
  
  if (verbose > 0) pind = initProgInd()
  for (g in 1:nPredicted)
  {
    gg = predictedIndex[g]
    sigOrder = order(-significance[, g])
    useGenes = sigOrder[2:(nPredictorGenes+1)]
    w.base = significance[useGenes, g]
    w = w.base * abs(w.base)^(power -1)
    bsub = rowSums(xScaled[, useGenes, drop = FALSE] * matrix(w, nSamples, nPredictorGenes, byrow = TRUE))
    bsub.scale = sqrt(mean(bsub^2))
    predictedTest[, g] = rowSums(xTestScaled[, useGenes, drop = FALSE] * 
                               matrix(w, nTestSamples, nPredictorGenes, byrow = TRUE)) / bsub.scale *
                             geneSD[gg] + geneMeans[gg]
    if (verbose > 0) pind = updateProgInd(g/nPredicted, pind)
  }

  predictedTest
}
      
.quickGeneVotingPredictor.CV = function(x, xtest = NULL, predictedIndex, nPredictorGenes = 20, 
                            power = 3, CVfold = 10,
                            corFnc = "bicor", corOptions = "use = 'p'")
{
  nSamples = nrow(x)
  nTestSamples = if (is.null(xtest)) 0 else nrow(xtest)
  nGenes = ncol(x)
  nPredicted = length(predictedIndex)
  
  ratio = nSamples/CVfold
  if (floor(ratio)!=ratio)
  {
    smaller =floor(ratio)
    nLarger = nSamples - CVfold * smaller
    binSizes = c(rep(smaller, CVfold-nLarger), rep(smaller +1, nLarger))
  } else
    binSizes = rep(ratio, CVfold)

  sampleOrder = sample(1:nSamples)

  CVpredicted = matrix(NA, nSamples, nPredicted)
  if (!is.null(xtest)) predictedTest = matrix(0, nTestSamples, nPredicted) else predictedTest = NULL

  cvStart = 1
  for (cv in 1:CVfold)
  {
    end = cvStart + binSizes[cv] - 1
    oob = sampleOrder[cvStart:end]
    CVx = x[-oob, , drop = FALSE]
    if (is.null(xtest)) CVxTest = x[oob, , drop = FALSE] else CVxTest = rbind( x[oob, , drop = FALSE], xtest)
    pred = .quickGeneVotingPredictor(CVx, CVxTest, predictedIndex, nPredictorGenes, power, corFnc,
                                    corOptions)
    CVpredicted[oob, ] = pred[c(1:binSizes[cv]), ]
    if (!is.null(xtest)) predictedTest = predictedTest + pred[-c(1:binSizes[cv]), ]
    cvStart = end + 1
  }

  if (!is.null(xtest)) predictedTest = predictedTest/CVfold

  list(CVpredicted = CVpredicted, predictedTest= predictedTest)
}





#' Remove leading principal components from data
#' 
#' This function calculates a fixed number of the first principal components of
#' the given data and returns the residuals of a linear regression of each
#' column on the principal components.
#' 
#' 
#' @param x Input data, a numeric matrix. All entries must be non-missing and
#' finite.
#' @param n Number of principal components to remove. This must be smaller than
#' the smaller of the number of rows and columns in \code{x}.
#' @return A matrix of residuals of the same dimensions as \code{x}.
#' @author Peter Langfelder
#' @seealso \code{\link{svd}} for singular value decomposition,
#' \code{\link{lm}} for linear regression
#' @keywords misc
removePrincipalComponents = function(x, n)
{
  if (sum(is.na(x)) > 0) x = t(impute.knn(t(x))$data)
  svd = svd(x, nu = n, nv = 0)
  PCs = as.data.frame(svd$u)
  names(PCs) = paste0("PC", c(1:n))
  fit = lm(x~., data = PCs)
  res = residuals(fit)
  res
}

#========================================================================================================
#
# Lin's network screening correlation functions
#
#========================================================================================================

.corWeighted = function(expr, y, ...) {
  modules = blockwiseModules(expr, ...)
  MEs = modules$MEs
  MEs = MEs[, colnames(MEs)!="MEgrey"]
  ns = networkScreening(y, MEs, expr, getQValues=F)
  ns$cor.Weighted
}

.corWeighted.new = function(expr, y, ...) {
  modules = blockwiseModules(expr, ...)
  MEs = modules$MEs
  MEs = MEs[, colnames(MEs)!="MEgrey"]
  scaledMEs = scale(MEs)
  ES = t(as.matrix(cor( y, MEs, use="p")))
  weightedAverageME = as.numeric(as.matrix(scaledMEs)%*%ES)/ncol(MEs)

  w=0.25
  y.new = w*scale(y) + (1-w)*weightedAverageME
  GS.new = as.numeric(cor(y.new, expr, use="p"))
  GS.new
}



