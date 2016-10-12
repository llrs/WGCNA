#===================================================================================================
#
# Multi-variate empirical Bayes with proper accounting for sigma
#
#===================================================================================================

.colWeightedMeans.x = function(data, weights, na.rm)
{
  nc = ncol(data)
  means = rep(NA, nc)

  for (c in 1:nc)
    means[c] = weighted.mean(data[, c], weights[, c], na.rm = na.rm)

  names(means) = colnames(data)
  means
}

.weightedScale = function(data, weights)
{
  weightSums = colSums(weights)
  means = .colWeightedMeans.x(data, weights, na.rm = TRUE)
  V1 = colSums(weights)
  V2 = colSums(weights^2)
  centered = data - matrix(means, nrow(data), ncol(data), byrow = TRUE)
  scale = sqrt( colSums(centered^2* weights, na.rm = TRUE) / ( V1 - V2/V1))

  scaled = centered/ matrix(scale,  nrow(data), ncol(data), byrow = TRUE)

  attr(scaled, "scaled:center") = means
  attr(scaled, "scaled:scale") = scale
  scaled
}

.weightedVar = function(x, weights)
{
  V1 = sum(weights)
  V2 = sum(weights^2)
  mean = sum(x * weights)/V1
  centered = x-mean
  sum(centered^2 * weights)/(V1-V2/V1)
}



#' Empirical Bayes-moderated adjustment for unwanted covariates
#' 
#' This functions removes variation in high-dimensional data due to unwanted
#' covariates while preserving variation due to retained covariates. To prevent
#' numerical instability, it uses Empirical bayes-moderated linear regression,
#' optionally in a robust (outlier-resistant) form.
#' 
#' This function uses Empirical Bayes-moderated (EB) linear regression to
#' remove variation in \code{data} due to the variables in
#' \code{removedCovariates} while retaining variation due to variables in
#' \code{retainedCovariates}, if any are given. The EB step uses simple normal
#' priors on the regression coefficients and inverse gamma priors on the
#' variances. The procedure starts with multivariate ordinary linear regression
#' of individual columns in \code{data} on \code{retainedCovariates} and
#' \code{removedCovariates}. To make the coefficients comparable, columns of
#' \code{data} are scaled to (weighted if weights are given) mean 0 and
#' variance 1. The resulting regression coefficients are used to determine the
#' parameters of the normal prior (mean, covariance, and inverse gamma or
#' median and biweight mid-covariance if robust priors are used), and the
#' variances are used to determine the parameters of the inverse gamma prior.
#' The EB step then essentially shrinks the coefficients toward their means,
#' with the amount of shrinkage determined by the prior covariance.
#' 
#' Using appropriate weights can make the data adjustment robust to outliers.
#' This can be achieved automatically by using the argument
#' \code{automaticWeights = "bicov"}. When bicov weights are used, we also
#' recommend setting the argument \code{maxPOutliers} to a maximum proportion
#' of samples that could be outliers. This is especially important if some of
#' the design variables are binary and can be expected to have a strong effect
#' on some of the columns in \code{data}, since standard biweight
#' midcorrelation (and its weights) do not work well on bimodal data.
#' 
#' The automatic bicov weights are determined from \code{data} only. It is
#' implicitly assumed that there are no outliers in the retained and removed
#' covariates. Outliers in the covariates are more difficult to work with
#' since, even if the regression is made robust to them, they can influence the
#' adjusted values for the sample in which they appear. Unless the the
#' covariate outliers can be attributed to a relevant vriation in experimental
#' conditions, samples with covariate outliers are best removed entirely before
#' calling this function.
#' 
#' @param data A 2-dimensional matrix or data frame of numeric data to be
#' adjusted. Variables (for example, genes or methylation profiles) should be
#' in columns and observations (samples) should be in rows.
#' @param retainedCovariates A vector or two-dimensional object (matrix or data
#' frame) giving the covariates whose effect on the data is to be retained. May
#' be \code{NULL} if there are no such "retained" covariates.
#' @param removedCovariates A vector or two-dimensional object (matrix or data
#' frame) giving the covariates whose effect on the data is to be removed. At
#' least one such covariate must be given.
#' @param weights Optional 2-dimensional matrix or data frame of the same
#' dimensions as \code{data} giving weights for each entry in \code{data}.
#' @param weightType One of (unique abbreviations of) \code{"apriori"} or
#' \code{"empirical"}. Determines whether a standard (\code{"apriori"}) or a
#' modified (\code{"empirical"}) weighted regression is used. The
#' \code{"apriori"} choice is suitable for weights that have been determined
#' without knowledge of the actual \code{data}, while \code{"empirical"} is
#' appropriate for situations where one wants to down-weigh cartain entries of
#' \code{data} because they may be outliers. In either case, the weights should
#' be determined in a way that is independent of the covariates (both retained
#' and removed).
#' @param stopOnSmallWeights Logical: should presence of small \code{"apriori"}
#' weights trigger an error? Because standard weighted regression assumes that
#' all weights are non-zero (otherwise estimates of standard errors will be
#' biased), this function will by default complain about the presence of too
#' small \code{"apriori"} weights.
#' @param tol Convergence criterion used in the numerical equation solver. When
#' the relative change in coefficients falls below this threshold, the system
#' will be considered to have converged.
#' @param maxIterations Maximum number of iterations to use.
#' @param scaleMeanToSamples Optional specification of samples (given as a
#' vector of indices) to whose means the resulting adjusted data should be
#' scaled (more precisely, shifted). If not given, the mean of all samples will
#' be used.
#' @param robustPriors Logical: should robust priors be used? This essentially
#' means replacing mean by median and covariance by biweight mid-covariance.
#' @param automaticWeights One of (unique abrreviations of) \code{"none"} or
#' \code{"bicov"}, instructing the function to calculate weights from the given
#' \code{data}. Value \code{"none"} will result in trivial weights; value
#' \code{"bicov"} will result in biweight midcovariance weights being used.
#' @param aw.maxPOutliers If \code{automaticWeights} above is \code{"bicov"},
#' this argument gets passed to the function \code{\link{bicovWeights}} and
#' determines the maximum proportion of outliers in calculating the weights.
#' See \code{\link{bicovWeights}} for more details.
#' @return A list with the following components: \item{adjustedData}{A matrix
#' of the same dimensions as the input \code{data}, giving the adjusted data.
#' If input \code{data} has non-NULL \code{dimnames}, these are copied.}
#' 
#' \item{residuals}{A matrix of the same dimensions as the input \code{data},
#' giving the residuals, that is, adjusted data with zero means.}
#' 
#' \item{coefficients}{A matrix of regression coefficients. Rows correspond to
#' the design matrix variables (mean, retained and removed covariates) and
#' columns correspond to the variables (columns) in \code{data}.}
#' 
#' \item{coefficiens.scaled}{A matrix of regression coefficients corresponding
#' to columns in \code{data} scaled to mean 0 and variance 1.}
#' 
#' \item{sigmaSq}{Estimated error variances (one for each column of input
#' \code{data}.}
#' 
#' \item{sigmaSq.scaled}{Estimated error variances corresponding to columns in
#' \code{data} scaled to mean 0 and variance 1.}
#' 
#' \item{fittedValues}{Fitted values calculated from the means and coefficients
#' corresponding to the removed covariates, i.e., roughly the values that are
#' subtracted out of the data.}
#' 
#' \item{adjustedData.OLS}{A matrix of the same dimensions as the input
#' \code{data}, giving the data adjusted by ordinary least squares. This
#' component should only be used for diagnostic purposes, not as input for
#' further downstream analyses, as the OLS adjustment is inferior to EB
#' adjustment. }
#' 
#' \item{residuals.OLS}{A matrix of the same dimensions as the input
#' \code{data}, giving the residuals obtained from ordinary least squares
#' regression, that is, OLS-adjusted data with zero means.}
#' 
#' \item{coefficients.OLS}{A matrix of ordinary least squares regression
#' coefficients.  Rows correspond to the design matrix variables (mean,
#' retained and removed covariates) and columns correspond to the variables
#' (columns) in \code{data}.}
#' 
#' \item{coefficiens.OLS.scaled}{A matrix of ordinary least squares regression
#' coefficients corresponding to columns in \code{data} scaled to mean 0 and
#' variance 1.  These coefficients are used to calculate priors for the EB
#' step.}
#' 
#' \item{sigmaSq.OLS}{Estimated OLS error variances (one for each column of
#' input \code{data}.}
#' 
#' \item{sigmaSq.OLS.scaled}{Estimated OLS error variances corresponding to
#' columns in \code{data} scaled to mean 0 and variance 1. These are used to
#' calculate variance priors for the EB step.}
#' 
#' \item{fittedValues.OLS}{OLS fitted values calculated from the means and
#' coefficients corresponding to the removed covariates.}
#' 
#' \item{weights}{A matrix of weights used in the regression models. The matrix
#' has the same dimension as the input \code{data}.}
#' 
#' \item{dataColumnValid}{Logical vector with one element per column of input
#' \code{data}, indicating whether the column was adjusted. Columns with zero
#' variance or too many missing data cannot be adjusted.}
#' 
#' \item{dataColumnWithZeroVariance}{Logical vector with one element per column
#' of input \code{data}, indicating whether the column had zero variance.}
#' 
#' \item{coefficientValid}{Logical matrix of the dimension (number of
#' covariates +1) times (number of variables in \code{data}), indicating
#' whether the corresponding regression coefficient is valid. Invalid
#' regression coefficients may be returned as missing values or as zeroes.}
#' @author Peter Langfelder
#' @seealso \code{\link{bicovWeights}} for suitable weights that make the
#' adjustment robust to outliers.
#' @keywords models regression
empiricalBayesLM = function(
  data, 
  removedCovariates,
  retainedCovariates = NULL, 
  weights = NULL,
  weightType = c("apriori", "empirical"),
  stopOnSmallWeights = TRUE,
  tol = 1e-4, maxIterations = 1000,
  scaleMeanToSamples = NULL,
  robustPriors = FALSE,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1)
{

  nSamples = nrow(data)
  designMat = NULL
  #mean.x = NULL
  #scale.x = NULL

  automaticWeights = match.arg(automaticWeights)
  if (automaticWeights=="bicov")
  {
    weightType = "empirical"
    weights = bicovWeights(data, maxPOutliers = aw.maxPOutliers)
  }

  weightType = match.arg(weightType)
  wtype = match(weightType, c("apriori", "empirical"))

  if (!is.null(retainedCovariates))
  {
     if (is.null(dim(retainedCovariates))) 
        retainedCovariates = data.frame(retainedCovariate = retainedCovariates)
     if (any(is.na(retainedCovariates))) stop("All elements of 'retainedCovariates' must be finite.")
     if (nrow(retainedCovariates)!=nSamples)
       stop("Numbers of rows in 'data' and 'retainedCovariates' differ.")
     retainedCovariates = as.data.frame(retainedCovariates)
     mm = model.matrix(~., data = retainedCovariates)[, -1, drop = FALSE]
     colSDs = colSds(mm)
     if (any(colSDs==0))
       stop("Some columns in 'retainedCovariates' have zero variance.")
     designMat = mm
  }

  if (is.null(removedCovariates))
    stop("'removedCovariates' must be supplied.")
  if (is.null(dim(removedCovariates))) 
      removedCovariates = data.frame(removedCovariate = removedCovariates)
  if (any(is.na(removedCovariates))) stop("All elements of 'removedCovariates' must be finite.")
  if (nrow(removedCovariates)!=nSamples)
    stop("Numbers of rows in 'data' and 'removedCovariates' differ.")
  removedCovariates = as.data.frame(removedCovariates)
  mm = model.matrix(~., data = removedCovariates)[, -1, drop = FALSE]
  colSDs = colSds(mm)
  if (any(colSDs==0))
    stop("Some columns in 'removedCovariates' have zero variance.")
  designMat = cbind(designMat, mm)
  removedColumns = (ncol(designMat)-ncol(mm) + 1):ncol(designMat)

  y.original = as.matrix(data)
  N.original = ncol(y.original)
  if (any(!is.finite(y.original))) 
  {
    warning(immediate. = TRUE,
            "Found missing and non-finite data. These will be removed.")
  }
    
  if (is.null(weights))
    weights = matrix(1, nSamples, ncol(y.original))

  if (any(!is.finite(weights)))
    stop("Given 'weights' contain some infinite or missing entries. All weights must be present and finite.")

  if (any(weights<0))
    stop("Given 'weights' contain negative entries. All weights must be non-negative.")

  originalWeights = weights

  dimnamesY = dimnames(y.original)

  varY = colVars(y.original, na.rm = TRUE)
  varYMissing = is.na(varY)
  varYZero = varY==0
  varYZero[is.na(varYZero)] = FALSE
  keepY = !(varYZero | varYMissing)

  y = y.original[, keepY]
  weights = weights[, keepY]
  yFinite = is.finite(y)
  weights[!yFinite] = 0
  nSamples.y = colSums(yFinite)

  if (weightType == "apriori")
  { 
    # Check weights
    if (any(weights[yFinite]<1))
    {
      if (stopOnSmallWeights)
      {
         stop("When weights are determined 'apriori', small weights are not allowed.\n",
              "Weights must be at least 1. Use weightType='empirical' if weights were determined from data.")
      } else
         warning(immediate. = TRUE,
                 "Small weights found. This can lead to unreliable fit with weight type 'apriori'.\n",
                 "Proceed with caution.")
    }
  }

  N = ncol(y)
  nc = ncol(designMat)

  # Scale y to mean zero and variance 1. This is needed for the prior to make sense.

  if (is.null(scaleMeanToSamples)) scaleMeanToSamples = c(1:nSamples)

  mean.y.target = .colWeightedMeans.x(y[scaleMeanToSamples, ], weights[scaleMeanToSamples, ], na.rm = TRUE)
  y = .weightedScale(y, weights)
  mean.y = attr(y, "scaled:center")
  mean.y.mat = matrix(mean.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
  scale.y = attr(y, "scaled:scale")
  scale.y.mat = matrix(scale.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
  y[!yFinite] = 0

  # Get the means of the design matrix with respect to all weight vectors.
  V1 = colSums(weights)
  V2 = colSums(weights^2)
  means.dm = t(designMat) %*% weights / matrix(V1, nrow = nc, ncol = N, byrow = TRUE)

  # Ordinary regression to get starting point for beta

  beta.OLS = matrix(NA, nc, N)
  betaValid = matrix(TRUE, nc, N)
  sigma.OLS = rep(NA, N)
  regressionValid = rep(TRUE, N)
  oldWeights = rep(-1, nSamples)
  for (i in 1:N)
  {
    w1 = weights[, i]
    y1 = y[, i]
    if (any(w1!=oldWeights))
    {
      centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
      #dmVar = colSds(centeredDM)
      dmVar.w = apply(centeredDM, 2, .weightedVar, weights = w1)
      #keepDM = dmVar > 0 & dmVar.w > 0
      keepDM = dmVar.w > 0
      centeredDM.keep = centeredDM[, keepDM]
      xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1)
      xtxInv = try(solve(xtx), silent = TRUE)
      oldWeights = w1
    }

    if (!inherits(xtx, "try-error"))
    {
      # The following is really (xtxInv %*% xy1), where xy1 = t(centeredDM) %*% (y[,i]*weights[,i])
      beta.OLS[keepDM, i] = colSums(xtxInv * colSums(centeredDM.keep * y1 * w1))
      betaValid[!keepDM, i] = FALSE
      y.pred = centeredDM.keep %*% beta.OLS[keepDM, i]
      if (weightType=="apriori")
      {
        # Standard calculation of sigma^2 in weighted refression
        sigma.OLS[i] = sum(w1 * (y1 - y.pred)^2)/(sum(yFinite[, i])-nc-1)
      } else {
        xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2)
        sigma.OLS[i] = sum( w1* (y1-y.pred)^2) / (V1[i] - V2[i]/V1[i] - sum(xtxw2 * xtxInv))
      }
    } else {
      regressionValid[i] = FALSE
      betaValid[, i] = FALSE
    }
  }
  if (any(!regressionValid))
     warning(immediate. = TRUE,
             "empiricalBayesLM: OLS regression failed in ", sum(!regressionValid), " variables.")

  # beta.OLS has columns corresponding to variables in data, and rows corresponding to columns in x.

  # Debugging...
  if (FALSE)
  {
    fit = lm(y~., data = data.frame(designMat))
    fit2 = lm(y~., data = as.data.frame(centeredDM))

    max(abs(fit2$coefficients[1,]))

    all.equal(c(fit$coefficients[-1, ]), c(fit2$coefficients[-1, ]))
    all.equal(c(fit$coefficients[-1, ]), c(beta.OLS))
    sigma.fit = apply(fit$residuals, 2, var)*(nSamples-1)/(nSamples-1-nc)
    all.equal(as.numeric(sigma.fit), as.numeric(sigma.OLS))

  }

  # Priors on beta : mean and variance
  if (robustPriors)
  {
    if (is.na(beta.OLS[, regressionValid]))
      stop("Some of OLS coefficients are missing. Please use non-robust priors.")
    prior.means = rowMedians(beta.OLS[, regressionValid], na.rm = TRUE)
    prior.covar = .bicov(t(beta.OLS[, regressionValid]))
  } else {
    prior.means = rowMeans(beta.OLS[, regressionValid], na.rm = TRUE)
    prior.covar = cov(t(beta.OLS[, regressionValid]), use = "complete.obs")
  }
  prior.inverse = solve(prior.covar)

  # Prior on sigma: mean and variance (median and MAD are bad estimators since the distribution is skewed)
  sigma.m = mean(sigma.OLS[regressionValid], na.rm = TRUE)
  sigma.v = var(sigma.OLS[regressionValid], na.rm = TRUE)

  # Turn the sigma mean and variance into the parameters of the inverse gamma distribution
  prior.a = sigma.m^2/sigma.v + 2
  prior.b = sigma.m * (prior.a-1)

  # Calculate the EB estimates
  beta.EB = beta.OLS
  sigma.EB = sigma.OLS
  
  for (i in which(regressionValid))
  {
    # Iterate to solve for EB regression coefficients (betas) and the residual variances (sigma)
    # It appears that this has to be done individually for each variable.

    difference = 1
    iteration = 1

    keepDM = betaValid[, i]
    beta.old = as.matrix(beta.OLS[keepDM, i])
    sigma.old = sigma.OLS[i]
    y1 = y[, i]
    w1 = weights[, i]

    centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
    centeredDM.keep = centeredDM[, keepDM]
    xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1)
    xtxInv = solve(xtx)

    if (all(keepDM))
    {
      prior.inverse.keep = prior.inverse
    } else
      prior.inverse.keep = solve(prior.covar[keepDM, keepDM])

    while (difference > tol && iteration <= maxIterations)
    {
      y.pred = centeredDM.keep %*% beta.old
      if (wtype==1)
      {
        # Apriori weights.
        fin1 = yFinite[, i]
        nSamples1 = sum(fin1)
        sigma.new = (sum(w1*(y1-y.pred)^2) + 2*prior.b)/ (nSamples1-nc + 2 * prior.a + 1)
      } else {
        # Empirical weights
        V1 = sum(w1)
        V2 = sum(w1^2)
        xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2)
        sigma.new = (sum( w1* (y1-y.pred)^2) + 2*prior.b) / (V1 - V2/V1 - sum(xtxw2 * xtxInv) + 2*prior.a + 2)
      }

      A = (prior.inverse.keep + xtx/sigma.new)/2
      A.inv = solve(A)
      #B = as.numeric(t(y1*w1) %*% centeredDM/sigma.new) + as.numeric(prior.inverse %*% prior.means)
      B = colSums(centeredDM.keep * y1 * w1/sigma.new) + colSums(prior.inverse.keep * prior.means[keepDM])

      #beta.new = A.inv %*% as.matrix(B)/2
      beta.new = colSums(A.inv * B)/2  # ...a different and hopefully faster way of writing the above

      difference = max( abs(sigma.new-sigma.old)/(sigma.new + sigma.old),
                        abs(beta.new-beta.old)/(beta.new + beta.old))

      beta.old = beta.new
      sigma.old = sigma.new
      iteration = iteration + 1
    }
    if (iteration > maxIterations) warning(immediate. = TRUE, 
                                           "Exceeded maximum number of iterations for variable ", i, ".")
    beta.EB[keepDM, i] = beta.old
    sigma.EB[i] = sigma.old
  }

  # Put output together. Will return the coefficients for lm and EB-lm, and the residuals with added mean.

  fitAndCoeffs = function(beta, sigma)
  {
      #fitted.removed = fitted = matrix(NA, nSamples, N)
      fitted.removed = fitted = y
      beta.fin = beta
      beta.fin[is.na(beta)] = 0
      for (i in which(regressionValid))
      {
         centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
         fitted.removed[, i] = centeredDM[, removedColumns, drop = FALSE] %*% 
                                   beta.fin[removedColumns, i, drop = FALSE]
         fitted[, i] = centeredDM %*% beta.fin[, i, drop = FALSE]
      }
      #browser()
      
      residuals = (y - fitted.removed) * scale.y.mat
      # Residuals now have weighted column means equal zero.

      meanShift =  matrix(mean.y.target, nSamples, N, byrow = TRUE)

      residuals.all = residualsWithMean.all = fitted.all = matrix(NA, nSamples, N.original)

      residuals.all[, keepY] = residuals
      residuals.all[, varYZero] = 0
      residuals.all[is.na(y.original)] = NA

      residualsWithMean.all[, keepY] = residuals.all[, keepY] + meanShift
      residualsWithMean.all[is.na(y.original)] = NA

      fitted.all[, keepY] = fitted * scale.y.mat + meanShift

      beta.all = beta.all.scaled = matrix(NA, nc+1, N.original)
      sigma.all = sigma.all.scaled = rep(NA, N.original)
      sigma.all.scaled[keepY] = sigma
      sigma.all[keepY] = sigma * scale.y^2

      beta.all[-1, keepY] = beta.fin * matrix(scale.y, nrow = nc, ncol = N, byrow = TRUE)
      beta.all.scaled[-1, keepY] = beta.fin
      beta.all[varYZero] = beta.all.scaled[-1, varYZero] = 0
      #alpha = mean.y - t(as.matrix(mean.x)) %*% beta * scale.y
      alpha = mean.y - colSums(beta * means.dm, na.rm = TRUE) * scale.y
      beta.all[1, keepY] = alpha
      beta.all.scaled[1, keepY] = 0

      dimnames(residualsWithMean.all) = dimnamesY
      colnames(beta.all) = colnames(beta.all.scaled) = colnames(y.original)
      rownames(beta.all) = rownames(beta.all.scaled) = c("(Intercept)", colnames(designMat))

      list(residuals = residuals.all,
           residualsWithMean = residualsWithMean.all,
           beta = beta.all,
           beta.scaled = beta.all.scaled,
           sigmaSq = sigma.all,
           sigmaSq.scaled = sigma.all.scaled,
           fittedValues = fitted.all)
  }

  fc.OLS = fitAndCoeffs(beta.OLS, sigma.OLS)
  fc.EB = fitAndCoeffs(beta.EB, sigma.EB)

  betaValid.all = matrix(FALSE, nc+1, N.original)
  betaValid.all[-1, keepY] = betaValid
  betaValid.all[1, keepY] = TRUE
  dimnames(betaValid.all) = dimnames(fc.OLS$beta)

  list( adjustedData = fc.EB$residualsWithMean,
       residuals = fc.EB$residuals,
       coefficients = fc.EB$beta,
       coefficients.scaled = fc.EB$beta.scaled,
       sigmaSq = fc.EB$sigmaSq,
       sigmaSq.scaled = fc.EB$sigmaSq.scaled,
       fittedValues = fc.EB$fittedValues,

       # OLS results
       adjustedData.OLS = fc.OLS$residualsWithMean,
       residuals.OLS = fc.OLS$residuals,
       coefficients.OLS = fc.OLS$beta,
       coefficients.OLS.scaled = fc.OLS$beta.scaled,
       sigmaSq.OLS = fc.OLS$sigmaSq,
       sigmaSq.OLS.scaled = fc.OLS$sigmaSq.scaled,
       fittedValues.OLS = fc.OLS$fittedValues,


       # Weights used in the model
       weights = originalWeights,

       # indices of valid fits
       dataColumnValid = keepY,
       dataColumnWithZeroVariance = varYZero,
       coefficientValid = betaValid.all)
}


#===================================================================================================
#
# Robust functions
#
#===================================================================================================

.bicov = function(x)
{
  x = as.matrix(x)
  if (any(!is.finite(x))) 
    stop("All entries of 'x' must be finite.")
  nc = ncol(x)
  nr = nrow(x)
  mx = colMedians(x)
  mx.mat = matrix(mx, nr, nc, byrow = TRUE)
  mads = colMads(x, constant = 1)
  mad.mat = matrix(mads, nr, nc, byrow = TRUE)
  u = abs((x - mx.mat)/(9 * mad.mat))
  a = matrix(as.numeric(u<1), nr, nc)

  topMat = a * (x-mx.mat) * (1-u^2)^2
  top = nr * t(topMat) %*% topMat
  botVec = colSums(a * (1-u^2) * (1-5*u^2))
  bot = botVec %o% botVec

  out = top/bot
  dimnames(out) = list(colnames(x), colnames(x))
  out
}


# Argument refWeight:
# w = (1-u^2)^2
# u^2 = 1-sqrt(w)
# referenceU = sqrt(1-sqrt(referenceW))



#' Weights used in biweight midcovariance
#' 
#' The function calculates weights used in the calculation of biweight
#' midcovariance and midcorrelation. The weights are designed such that
#' outliers get smaller weights; the weights become zero for data points more
#' than 9 median absolute deviations from the median.
#' 
#' 
#' @param x A vector or a two-dimensional array (matrix or data frame).  If
#' two-dimensional, the weights will be calculated separately on each column.
#' @param pearsonFallback Logical: if the median absolute deviation is zero,
#' should standard deviation be substituted?
#' @param maxPOutliers Optional specification of the maximum proportion of
#' outliers, i.e., data with weights equal to \code{outlierReferenceWeight}
#' below.
#' @param outlierReferenceWeight A number between 0 and 1 specifying what is to
#' be considered an outlier when calculating the proportion of outliers.
#' @param defaultWeight Value used for weights that correspond to a finite
#' \code{x} but the weights themselves would not be finite, for example, when a
#' column in \code{x} is constant.
#' @return A vector or matrix of the same dimensions as the input \code{x}
#' giving the weights.
#' @author Peter Langfelder
#' @seealso \code{\link{bicor}}
#' @references This function is based on Equation (3) in
#' 
#' Langfelder P, Horvath S (2012) Fast R Functions for Robust Correlations and
#' Hierarchical Clustering Journal of Statistical Software 46(11) 1-17 PMID:
#' 23050260 PMCID: PMC3465711
#' 
#' That article also describes the Pearson fallback and maximum proportion of
#' outliers in detail. For a full discussion of the biweight midcovariance and
#' midcorrelation, see
#' 
#' Wilcox RR (2005). Introduction to Robust Estimation and Hypothesis Testing.
#' 2nd edition. Academic Press, Section 9.3.8, page 399 as well as Section
#' 3.12.1, page 83.
#' @keywords misc
#' @examples
#' 
#' x = rnorm(100);
#' x[1] = 10;
#' plot(x, bicovWeights(x));
#' 
bicovWeights = function(x, pearsonFallback = TRUE, maxPOutliers = 1,
                        outlierReferenceWeight = 0.5625,
                        defaultWeight = 0)
{
  referenceU = sqrt(1-sqrt(outlierReferenceWeight^2))
  dimX = dim(x)
  dimnamesX = dimnames(x)
  x = as.matrix(x)
  nc = ncol(x)
  nr = nrow(x)
  mx = colMedians(x, na.rm = TRUE)
  mx.mat = matrix(mx, nr, nc, byrow = TRUE)
  mads = colMads(x, constant = 1, na.rm = TRUE)
  madZero = mads==0
  if (any(madZero)) 
  {
     warning(immediate. = TRUE,
             "MAD is zero in some columns of 'x'.")
     if (pearsonFallback)
     {
       sds = colSds(x[, madZero, drop = FALSE], na.rm = TRUE)
       mads[madZero] = sds * qnorm(0.75)
     }
  }
  mad.mat = matrix(mads, nr, nc, byrow = TRUE)
  u = (x - mx.mat)/(9 * mad.mat)
  if (maxPOutliers < 0.5)
  {
    lq = colQuantileC(u, p = maxPOutliers)
    uq = colQuantileC(u, p = 1-maxPOutliers)
    lq[is.na(lq)] = 0
    uq[is.na(uq)] = 0

    lq[lq>-referenceU] = -referenceU
    uq[uq < referenceU] = referenceU
    lq = abs(lq)
    changeNeg = which(lq>referenceU)
    changePos = which(uq > referenceU)

    for (c in changeNeg)
    {
      neg1 = u[, c] < 0
      neg1[is.na(neg1)] = FALSE
      u[neg1, c] = u[neg1, c] * referenceU/lq[c]
    }

    for (c in changePos)
    {
      pos1 = u[, c] > 0
      pos1[is.na(pos1)] = FALSE
      u[pos1, c] = u[pos1, c] * referenceU/uq[c]
    }
  }

  a = matrix(as.numeric(abs(u)<1), nr, nc)
  weights = a * (1-u^2)^2
  weights[is.na(x)] = 0
  weights[!is.finite(weights)] = defaultWeight
  dim(weights) = dimX
  if (!is.null(dimX)) dimnames(weights) = dimnamesX
  weights
}

#=====================================================================================================
#
# Diagnostic plot for coefficients
#
#=====================================================================================================




