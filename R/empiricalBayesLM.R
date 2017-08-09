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

# Defaults for certain fitting functions

.initialFit.defaultWeightName = function(fitFnc)
{
  wNames = c(rlm = "w", lmrob = "rweights")
  if (fitFnc %in% names(wNames)) return(wNames[ match(fitFnc, names(wNames))])
  NULL
}

.initialFit.defaultOptions = function(fitFnc)
{
  defOpt = list(
    rlm = list(),
    lmrob = list(model = FALSE, x = FALSE, y = FALSE, control = lmRob.control()))
  if (fitFnc %in% names(defOpt)) return(defOpt[[ match(fitFnc, names(defOpt))]])
  list()
}

.initialFit.requiresFormula = function(fitFnc)
{
  reqForm = c(rlm = FALSE, lmrob = TRUE)
  if (fitFnc %in% names(reqForm)) return(reqForm[ match(fitFnc, names(reqForm))])
  FALSE
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
#' \code{removedCovariates}. Alternatively, the user may specify an intial fit
#' function (e.g., robust linear regression). To make the coefficients
#' comparable, columns of \code{data} are scaled to (weighted if weights are
#' given) mean 0 and variance 1. The resulting regression coefficients are used
#' to determine the parameters of the normal prior (mean, covariance, and
#' inverse gamma or median and biweight mid-covariance if robust priors are
#' used), and the variances are used to determine the parameters of the inverse
#' gamma prior. The EB step then essentially shrinks the coefficients toward
#' their means, with the amount of shrinkage determined by the prior
#' covariance.
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
#' @param removedCovariates A vector or two-dimensional object (matrix or data
#' frame) giving the covariates whose effect on the data is to be removed. At
#' least one such covariate must be given.
#' @param retainedCovariates A vector or two-dimensional object (matrix or data
#' frame) giving the covariates whose effect on the data is to be retained. May
#' be \code{NULL} if there are no such "retained" covariates.
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
#' ======= \item{data}{ A 2-dimensional matrix or data frame of numeric data to
#' be adjusted. Variables (for example, genes or methylation profiles) should
#' be in columns and observations (samples) should be in rows. }
#' \item{removedCovariates}{ A vector or two-dimensional object (matrix or data
#' frame) giving the covariates whose effect on the data is to be removed. At
#' least one such covariate must be given. } \item{retainedCovariates}{ A
#' vector or two-dimensional object (matrix or data frame) giving the
#' covariates whose effect on the data is to be retained. May be \code{NULL} if
#' there are no such "retained" covariates. }
#'
#' \item{initialFitFunction}{ Function name to perform the initial fit. The
#' default is to use the internal implementation of linear model fitting. The
#' function must take arguments \code{formula} and \code{data} or \code{x} and
#' \code{y}, plus possibly additional arguments. The return value must be a
#' list with component \code{coefficients}, either \code{scale} or
#' \code{residuals}, and weights must be returned in component specified by
#' \code{initialFit.returnWeightName}.  See \code{\link{lm}}, \code{\link{rlm}}
#' and other standard fit functions for examples of suitable functions.}
#'
#' \item{initialFitOptions}{ Optional specifications of extra arguments for
#' \code{initialFitFunction}, apart from \code{formula} and \code{data} or
#' \code{x} and \code{y}.  Defaults are provided for functions
#' \code{\link[MASS]{rlm}} and \code{\link[robust]{lmRob}}, i.e., if either of
#' these two functions are used as \code{initialFitFunction}, suitable initial
#' fit options will be chosen automatically.}
#'
#' \item{initialFitRequiresFormula}{ Logical: does the initial fit function
#' need \code{formula} and \code{data} arguments? If \code{TRUE},
#' \code{initialFitFunction} will be called with arguments \code{formula} and
#' \code{data}, otherwise with arguments \code{x} and \code{y}.}
#'
#' \item{initialFit.returnWeightName}{ Name of the component of the return
#' value of \code{initialFitFunction} that contains the weights used in the
#' fit. Suitable default value will be chosen automatically for
#' \code{\link[MASS]{rlm}} and \code{\link[robust]{lmRob}}. } \item{weights}{
#' Optional 2-dimensional matrix or data frame of the same dimensions as
#' \code{data} giving weights for each entry in \code{data}. These weights will
#' be used in the initial fit and are are separate from the ones returned by
#' \code{initialFitFunction} if it is specified. }
#'
#' \item{weightType}{ One of (unique abbreviations of) \code{"apriori"} or
#' \code{"empirical"}. Determines whether a standard (\code{"apriori"}) or a
#' modified (\code{"empirical"}) weighted regression is used. The
#' \code{"apriori"} choice is suitable for weights that have been determined
#' without knowledge of the actual \code{data}, while \code{"empirical"} is
#' appropriate for situations where one wants to down-weigh cartain entries of
#' \code{data} because they may be outliers. In either case, the weights should
#' be determined in a way that is independent of the covariates (both retained
#' and removed). }
#'
#' \item{stopOnSmallWeights}{ Logical: should presence of small
#' \code{"apriori"} weights trigger an error? Because standard weighted
#' regression assumes that all weights are non-zero (otherwise estimates of
#' standard errors will be biased), this function will by default complain
#' about the presence of too small \code{"apriori"} weights. } \item{tol}{
#' Convergence criterion used in the numerical equation solver. When the
#' relative change in coefficients falls below this threshold, the system will
#' be considered to have converged. } \item{maxIterations}{ Maximum number of
#' iterations to use. } \item{scaleMeanToSamples}{ Optional specification of
#' samples (given as a vector of indices) to whose means the resulting adjusted
#' data should be scaled (more precisely, shifted). If not given, the mean of
#' all samples will be used. } \item{robustPriors}{ Logical: should robust
#' priors be used? This essentially means replacing mean by median and
#' covariance by biweight mid-covariance. } \item{automaticWeights}{ One of
#' (unique abrreviations of) \code{"none"} or \code{"bicov"}, instructing the
#' function to calculate weights from the given \code{data}. Value
#' \code{"none"} will result in trivial weights; value \code{"bicov"} will
#' result in biweight midcovariance weights being used. }
#' \item{aw.maxPOutliers}{ If \code{automaticWeights} above is \code{"bicov"},
#' this argument gets passed to the function \code{\link{bicovWeights}} and
#' determines the maximum proportion of outliers in calculating the weights.
#' See \code{\link{bicovWeights}} for more details. }
#'
#' \item{minDesignDeviation}{ Minimum standard deviation for columns of the
#' design matrix to be retained. Columns with standard deviations below this
#' number will be removed (effectively removing the corresponding terms from
#' the design). }
#'
#' \item{garbageCollectInterval}{ Number of variables after which to call
#' garbage collection. }
#'
#' \item{getOLSAdjustedData}{Logical: should data adjusted by ordinary least
#' squares or by \code{initialFitFunction}, if specified, be returned?}
#'
#' \item{getResiduals}{Logical: should the residuals (adjusted values without
#' the means) be returned?}
#'
#' \item{getFittedValues}{Logical: should fitted values be returned?}
#'
#' \item{getWeights}{Logical: should the final weights be returned?}
#'
#' \item{getEBadjustedData}{Logical: should the EB step be performed and the
#' adjusted data returned? If this is \code{FALSE}, the function acts as a
#' rather slow but still potentially useful adjustment using standard fit
#' functions.}
#'
#' \item{verbose}{Level of verbosity. Zero means silent, higher values result
#' in more diagnostic messages being printed.} \item{indent}{Indentation of
#' diagnostic messages. Each unit adds two spaces.}
#' @author Peter Langfelder
#' @seealso \code{\link{bicovWeights}} for suitable weights that make the
#' adjustment robust to outliers.
#' @keywords models regression
#' @importFrom robust glmRob.control
empiricalBayesLM = function(
  data,
  removedCovariates,
  retainedCovariates = NULL,
  initialFitFunction = NULL,
  initialFitOptions = NULL,
  initialFitRequiresFormula = NULL,
  initialFit.returnWeightName = NULL,
  weights = NULL,
  weightType = c("apriori", "empirical"),
  stopOnSmallWeights = TRUE,
  tol = 1e-4, maxIterations = 1000,
  scaleMeanToSamples = NULL,
  robustPriors = FALSE,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1,
  minDesignDeviation = 1e-10,
  garbageCollectInterval = 50000,
  getOLSAdjustedData = TRUE,
  getResiduals = TRUE,
  getFittedValues = TRUE,
  getWeights = TRUE,
  getEBadjustedData = TRUE,
  verbose = 0, indent = 0)
{

  spaces = indentSpaces(indent)

  nSamples = nrow(data)
  designMat = NULL
  #mean.x = NULL
  #scale.x = NULL

  automaticWeights = match.arg(automaticWeights)
  if (automaticWeights=="bicov")
  {
    weightType = "empirical"
    if (verbose > 0) printFlush(paste(spaces, "..calculating weights.."))
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

  if (verbose > 0) {
    printFlush(paste(spaces, "..checking for non-varying responses.."))
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

    if (verbose > 0) printFlush(paste(spaces, "..standardizing responses.."))
    if (is.null(scaleMeanToSamples)) scaleMeanToSamples = c(1:nSamples)

    mean.y.target = .colWeightedMeans.x(y[scaleMeanToSamples, ], weights[scaleMeanToSamples, ], na.rm = TRUE)
    y = .weightedScale(y, weights)
    mean.y = attr(y, "scaled:center")
    mean.y.mat = matrix(mean.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
    scale.y = attr(y, "scaled:scale")
    scale.y.mat = matrix(scale.y, nrow = nrow(y), ncol = ncol(y), byrow = TRUE)
    y[!yFinite] = 0

    if (verbose > 0) printFlush(paste(spaces, "..initial model fitting.."))
    # Prepaer ordinary regression to get starting points for beta and sigma
    beta.OLS = matrix(NA, nc, N)
    betaValid = matrix(TRUE, nc, N)
    sigma.OLS = rep(NA, N)
    regressionValid = rep(TRUE, N)

    # Get the means of the design matrix with respect to all weight vectors.
    V1 = colSums(weights)
    V2 = colSums(weights^2)
    means.dm = t(designMat) %*% weights / matrix(V1, nrow = nc, ncol = N, byrow = TRUE)
    i = 0
    on.exit(printFlush(paste0("Error occurred at i = ", i)))
    #on.exit(browser())
    pindStep = max(1, floor(N/100))
    if (is.null(initialFitFunction)) {
      # Ordinary weighted least squares, in two varieties.
      oldWeights = rep(-1, nSamples)
      if (verbose > 1) pind = initProgInd()
      for (i in 1:N) {
        w1 = weights[, i]
        y1 = y[, i]
        if (any(w1!=oldWeights)) {
          centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
          #dmVar = colSds(centeredDM)
          dmVar.w = apply(centeredDM, 2, .weightedVar, weights = w1)
          #keepDM = dmVar > 0 & dmVar.w > 0
          keepDM = dmVar.w > minDesignDeviation^2
          xtxInv = try(
            {
              centeredDM.keep = centeredDM[, keepDM, drop = FALSE]
              xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1)
              solve(xtx)
            }, silent = TRUE)
        }

        if (!inherits(xtxInv, "try-error")) {
          oldWeights = w1
          # The following is really (xtxInv %*% xy1), where xy1 = t(centeredDM) %*% (y[,i]*weights[,i])
          beta.OLS[keepDM, i] = colSums(xtxInv * colSums(centeredDM.keep * y1 * w1))
          betaValid[!keepDM, i] = FALSE
          y.pred = centeredDM.keep %*% beta.OLS[keepDM, i, drop = FALSE]
          if (weightType=="apriori")
          {
            # Standard calculation of sigma^2 in weighted refression
            sigma.OLS[i] = sum((w1>0) * (y1 - y.pred)^2)/(sum(w1>0)-nc-1)
          } else {
            xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2)
            sigma.OLS[i] = sum( w1* (y1-y.pred)^2) / (V1[i] - V2[i]/V1[i] - sum(xtxw2 * xtxInv))
          }
        } else {
          regressionValid[i] = FALSE
          betaValid[, i] = FALSE
        }
        if (i%%garbageCollectInterval ==0) gc()
        if (verbose > 1 && i%%pindStep==0) pind = updateProgInd(i/N, pind)
      }
      if (verbose > 1) {
        pind = updateProgInd(i/N, pind)
        printFlush()
      }
      rweights = weights
    } else {
      fitFnc = match.fun(initialFitFunction)
      if (is.null(initialFit.returnWeightName))
        initialFit.returnWeightName = .initialFit.defaultWeightName(initialFitFunction)

      if (is.null(initialFitOptions))
        initialFitOptions = .initialFit.defaultOptions(initialFitFunction)

      if (is.null(initialFitRequiresFormula))
        initialFitRequiresFormula = .initialFit.requiresFormula(initialFitFunction)

      rweights = weights
      if (verbose > 1) pind = initProgInd()
      for (i in 1:N) {
        y1 = y[, i]
        w1 = weights[, i]
        dmVar.w = apply(designMat, 2, .weightedVar, weights = w1)
        keepDM = dmVar.w > 0
        if (initialFitRequiresFormula) {
          modelData = data.frame(designMat[, keepDM, drop = FALSE])
          fit = try(do.call(fitFnc, c(list(formula = "y1~.", data = data.frame(y1 = y1, modelData), w = w1),
                                      initialFitOptions)))
        } else {
          fit = try(do.call(fitFnc, c(list(x = cbind(intercept = rep(1, nSamples),
                                                     designMat[, keepDM, drop = FALSE]),
                                           y = y1, w = w1),
                                      initialFitOptions)))
        }
        if (!inherits(fit, "try-error")) {
          beta.OLS[keepDM, i] = fit$coefficients[-1]
          betaValid[!keepDM, i] = FALSE
          rw1 = getElement(fit, initialFit.returnWeightName)
          if (length(rw1)!=nSamples)
            stop("Length of component '", initialFit.returnWeightName,
                 "' does not match the number of samples.\n",
                 "Please check that the name of the component containing robustness weights\n",
                 "is specified correctly.")
          rweights[, i] = rw1 * w1
          if (!is.null(fit$scale)) {
            sigma.OLS[i] = fit$scale
          } else {
            # This is not the greatest estimate since it is non-robust and does not take the final weights into
            # account. The hope is that this will never be used.
            sigma.OLS[i] = sum((w1>0) * fit$residuals^2)/(sum(w1>0)-nc-1)
          }
        } else {
          regressionValid[i] = FALSE
          betaValid[, i] = FALSE
        }
        if (i%%garbageCollectInterval ==0) gc()
        if (verbose > 1 && i%%pindStep==0) pind = updateProgInd(i/N, pind)
      }
      if (verbose > 1) {pind = updateProgInd(i/N, pind); printFlush()}
    }

    if (any(!regressionValid))
      warning(immediate. = TRUE,
              "empiricalBayesLM: initial regression failed in ", sum(!regressionValid), " variables.")

    if (all(!regressionValid))
      stop("Initial regression model failed for all columns in 'data'.\n",
           "Last model returned the following error:\n\n",
           if (is.null(initialFitFunction)) xtx else fit,
           "\n\nPlease check that the model is correctly specified.")

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

    if (getEBadjustedData) {
      # Priors on beta : mean and variance
      if (verbose > 0) printFlush(paste0(spaces, "..calculating priors.."))
      if (robustPriors) {
        if (is.na(beta.OLS[, regressionValid]))
          stop("Some of OLS coefficients are missing. Please use non-robust priors.")
        prior.means = rowMedians(beta.OLS[, regressionValid, drop = FALSE], na.rm = TRUE)
        prior.covar = .bicov(t(beta.OLS[, regressionValid, drop = FALSE]))
      } else {
        prior.means = rowMeans(beta.OLS[, regressionValid, drop = FALSE], na.rm = TRUE)
        prior.covar = cov(t(beta.OLS[, regressionValid, drop = FALSE]), use = "complete.obs")
      }
      prior.inverse = solve(prior.covar)

      # Prior on sigma: mean and variance (median and MAD are bad estimators since the distribution is skewed)
      sigma.m = mean(sigma.OLS[regressionValid], na.rm = TRUE)
      sigma.v = var(sigma.OLS[regressionValid], na.rm = TRUE)

      # Turn the sigma mean and variance into the parameters of the inverse gamma distribution
      prior.a = sigma.m^2/sigma.v + 2
      prior.b = sigma.m * (prior.a-1)

      # Calculate the EB estimates
      if (verbose > 0) printFlush(paste0(spaces, "..calculating final coefficients.."))
      beta.EB = beta.OLS
      sigma.EB = sigma.OLS

      # Get the means of the design matrix with respect to all rweight vectors.
      rV1 = colSums(rweights)
      rV2 = colSums(rweights^2)
      rmeans.dm = t(designMat) %*% weights / matrix(V1, nrow = nc, ncol = N, byrow = TRUE)

      if (verbose > 1) pind = initProgInd()
      for (i in which(regressionValid)) {
        # Iterate to solve for EB regression coefficients (betas) and the residual variances (sigma)
        # It appears that this has to be done individually for each variable.

        difference = 1
        iteration = 1

        y1 = y[, i]
        w1 = rweights[, i]

        centeredDM = designMat - matrix(rmeans.dm[, i], nSamples, nc, byrow = TRUE)
        dmVar.w = apply(centeredDM, 2, .weightedVar, weights = w1)
        keepDM = dmVar.w > 0 & betaValid[, i]

        beta.old = as.matrix(beta.OLS[keepDM, i, drop = FALSE])
        sigma.old = sigma.OLS[i]

        centeredDM.keep = centeredDM[, keepDM, drop = FALSE]
        xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1)
        xtxInv = solve(xtx)

        if (all(keepDM)) {
          prior.inverse.keep = prior.inverse
        } else
          prior.inverse.keep = solve(prior.covar[keepDM, keepDM, drop = FALSE])

        while (difference > tol && iteration <= maxIterations) {
          y.pred = centeredDM.keep %*% beta.old
          if (wtype==1) {
            # Apriori weights.
            fin1 = yFinite[, i]
            nSamples1 = sum(fin1)
            sigma.new = (sum(fin1 * (y1-y.pred)^2) + 2*prior.b)/ (nSamples1-nc + 2 * prior.a + 1)
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
        if (i%%garbageCollectInterval ==0) gc()
        if (verbose > 1 && i%%pindStep==0) pind = updateProgInd(i/N, pind)
      }
      if (verbose > 1) {pind = updateProgInd(i/N, pind); printFlush()}
    } ### if (getEBAdjustedData)

    on.exit(NULL)


    # Put output together. Will return the coefficients for lm and EB-lm, and the residuals with added mean.

    if (verbose > 0) printFlush(paste(spaces, "..calculating adjusted data.."))

    fitAndCoeffs = function(beta, sigma) {
      #fitted.removed = fitted = matrix(NA, nSamples, N)
      fitted.removed = y
      if (getFittedValues) fitted = y
      beta.fin = beta
      beta.fin[is.na(beta)] = 0
      for (i in which(regressionValid)) {
        centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
        fitted.removed[, i] = centeredDM[, removedColumns, drop = FALSE] %*%
          beta.fin[removedColumns, i, drop = FALSE]
        if (getFittedValues) fitted[, i] = centeredDM %*% beta.fin[, i, drop = FALSE]
        if (i%%garbageCollectInterval ==0) gc()
      }
      #browser()

      residuals = (y - fitted.removed) * scale.y.mat
      # Residuals now have weighted column means equal zero.

      meanShift =  matrix(mean.y.target, nSamples, N, byrow = TRUE)

      if (getFittedValues)  {
        fitted.all = matrix(NA, nSamples, N.original)
        fitted.all[, keepY] = fitted * scale.y.mat + meanShift
        dimnames(fitted.all) = dimnamesY
      }

      residualsWithMean.all = matrix(NA, nSamples, N.original)
      if (getResiduals) {
        residuals.all = residualsWithMean.all
        residuals.all[, keepY] = residuals
        residuals.all[, varYZero] = 0
        residuals.all[is.na(y.original)] = NA
      }

      residualsWithMean.all[, keepY] = residuals + meanShift
      residualsWithMean.all[, varYZero] = 0
      residualsWithMean.all[is.na(y.original)] = NA


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

      list(residuals = if (getResiduals) residuals.all else NULL,
           residualsWithMean = residualsWithMean.all,
           beta = beta.all,
           beta.scaled = beta.all.scaled,
           sigmaSq = sigma.all,
           sigmaSq.scaled = sigma.all.scaled,
           fittedValues = if (getFittedValues) fitted.all else NULL)
    }

    fc.OLS = fitAndCoeffs(beta.OLS, sigma.OLS)
    if (getEBadjustedData) fc.EB = fitAndCoeffs(beta.EB, sigma.EB)

    betaValid.all = matrix(FALSE, nc+1, N.original)
    betaValid.all[-1, keepY] = betaValid
    betaValid.all[1, keepY] = TRUE
    dimnames(betaValid.all) = dimnames(fc.OLS$beta)

    finalWeights = originalWeights
    finalWeights[, keepY] = rweights


    list( adjustedData = if (getEBadjustedData) fc.EB$residualsWithMean else NULL,
          residuals = if (getResiduals & getEBadjustedData) fc.EB$residuals else NULL,
          coefficients = if (getEBadjustedData) fc.EB$beta else NULL,
          coefficients.scaled = if (getEBadjustedData) fc.EB$beta.scaled else NULL,
          sigmaSq = if (getEBadjustedData) fc.EB$sigmaSq else NULL,
          sigmaSq.scaled = if (getEBadjustedData) fc.EB$sigmaSq.scaled else NULL,
          fittedValues = if (getFittedValues && getEBadjustedData) fc.EB$fittedValues else NULL,

          # OLS results
          adjustedData.OLS = if (getOLSAdjustedData) fc.OLS$residualsWithMean else NULL,
          residuals.OLS = if (getResiduals) fc.OLS$residuals else NULL,
          coefficients.OLS = fc.OLS$beta,
          coefficients.OLS.scaled = fc.OLS$beta.scaled,
          sigmaSq.OLS = fc.OLS$sigmaSq,
          sigmaSq.OLS.scaled = fc.OLS$sigmaSq.scaled,
          fittedValues.OLS = if (getFittedValues) fc.OLS$fittedValues else NULL,


          # Weights used in the model
          initialWeights = if (getWeights) originalWeights else NULL,
          finalWeights = if (getWeights) finalWeights else NULL,


          # indices of valid fits
          dataColumnValid = keepY,
          dataColumnWithZeroVariance = varYZero,
          coefficientValid = betaValid.all)
  }
}


#===================================================================================================
#
# Linear model coefficients
#
#===================================================================================================

.linearModelCoefficients = function(
  responses,
  predictors,

  weights = NULL,
  automaticWeights = c("none", "bicov"),
  aw.maxPOutliers = 0.1,

  getWeights = FALSE,

  garbageCollectInterval = 50000,
  minDesignDeviation = 1e-10,

  verbose = 0, indent = 0)
{

  spaces = indentSpaces(indent)

  designMat = NULL

  automaticWeights = match.arg(automaticWeights)
  if (automaticWeights=="bicov")
  {
    if (verbose > 0) printFlush(paste(spaces, "..calculating weights.."))
    weights = bicovWeights(responses, maxPOutliers = aw.maxPOutliers)
  }

  if (is.null(dim(predictors)))
    predictors = data.frame(retainedCovariate = predictors)

  keepSamples = rowSums(is.na(predictors))==0

  responses = responses[keepSamples, , drop = FALSE]
  predictors = predictors[keepSamples, , drop = FALSE]
  nSamples = nrow(responses)

  if (nrow(predictors)!=nSamples)
    stop("Numbers of rows in 'responses' and 'predictors' differ.")
  predictors = as.data.frame(predictors)
  mm = model.matrix(~., data = predictors)
  colSDs = colSds(mm[, -1, drop = FALSE], na.rm = TRUE)
  if (any(colSDs==0))
    stop("Some columns in 'predictors' have zero variance.")
  designMat = mm

  y.original = as.matrix(responses)
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

  if (verbose > 0) printFlush(paste(spaces, "..checking for non-varying responses.."))
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
  y[!yFinite] = 0

  N = ncol(y)
  nc = ncol(designMat)

  if (verbose > 0) printFlush(paste(spaces, "..model fitting.."))
  # Prepaer ordinary regression to get starting points for beta and sigma
  beta.OLS = matrix(NA, nc, N)
  betaValid = matrix(TRUE, nc, N)
  sigma.OLS = rep(NA, N)
  regressionValid = rep(TRUE, N)

  # Get the means of the design matrix with respect to all weight vectors.
  V1 = colSums(weights)
  V2 = colSums(weights^2)
  means.dm = t(designMat) %*% weights / matrix(V1, nrow = nc, ncol = N, byrow = TRUE)
  i = 0
  on.exit(printFlush(paste0("Error occurred at i = ", i)))
  #on.exit(browser())
  pindStep = max(1, floor(N/100))

  # Ordinary weighted least squares
  oldWeights = rep(-1, nSamples)
  if (verbose > 1) pind = initProgInd()
  for (i in 1:N)
  {
    w1 = weights[, i]
    y1 = y[, i]
    if (any(w1!=oldWeights))
    {
      #centeredDM = designMat - matrix(means.dm[, i], nSamples, nc, byrow = TRUE)
      centeredDM = designMat; # Do not center the design mat.
      #dmVar = colSds(centeredDM)
      dmVar.w = apply(centeredDM, 2, .weightedVar, weights = w1)
      #keepDM = dmVar > 0 & dmVar.w > 0
      keepDM = dmVar.w > minDesignDeviation^2
      keepDM[1] = TRUE; # For the intercept
      centeredDM.keep = centeredDM[, keepDM, drop = FALSE]
      xtx = t(centeredDM.keep) %*% (centeredDM.keep * w1)
      xtxInv = try(solve(xtx), silent = TRUE)
      oldWeights = w1
    }

    if (!inherits(xtxInv, "try-error"))
    {
      # The following is really (xtxInv %*% xy1), where xy1 = t(centeredDM) %*% (y[,i]*weights[,i])
      beta.OLS[keepDM, i] = colSums(xtxInv * colSums(centeredDM.keep * y1 * w1))
      betaValid[!keepDM, i] = FALSE
      y.pred = centeredDM.keep %*% beta.OLS[keepDM, i, drop = FALSE]
      #if (weightType=="apriori")
      #{
      # Standard calculation of sigma^2 in weighted refression
      sigma.OLS[i] = sum((w1>0) * (y1 - y.pred)^2)/(sum(w1>0)-nc-1)
      #} else {
      #  xtxw2 =  t(centeredDM.keep) %*% (centeredDM.keep *w1^2)
      #  sigma.OLS[i] = sum( w1* (y1-y.pred)^2) / (V1[i] - V2[i]/V1[i] - sum(xtxw2 * xtxInv))
      #}
    } else {
      regressionValid[i] = FALSE
      betaValid[, i] = FALSE
    }
    if (i%%garbageCollectInterval ==0) gc()
    if (verbose > 1 && i%%pindStep==0) pind = updateProgInd(i/N, pind)
  }
  if (verbose > 1) {pind = updateProgInd(i/N, pind); printFlush()}

  on.exit(NULL)

  if (any(!regressionValid))
    warning(immediate. = TRUE,
            "linearModelCoefficients: initial regression failed in ", sum(!regressionValid), " variables.")

  if (all(!regressionValid))
    stop("Initial regression model failed for all columns in 'data'.\n",
         "Last model returned the following error:\n\n",
         xtx,
         "\n\nPlease check that the model is correctly specified.")

  # beta.OLS has columns corresponding to variables in responses, and rows corresponding to columns in x.

  # Extend results to all variables.

  beta.all = matrix(NA, nc, N.original)
  beta.all[, keepY] = beta.OLS
  dimnames(beta.all) = list(colnames(designMat), dimnamesY[[2]])

  sigma.all = rep(NA, N.original)
  sigma.all[keepY] = sigma.OLS


  betaValid.all = matrix(FALSE, nc, N.original)
  betaValid.all[, keepY] = betaValid
  dimnames(betaValid.all) = dimnames(beta.all)


  list(coefficients = beta.all,
       sigmaSq = sigma.all,
       # Weights used in the model
       weights = if (getWeights) weights else NULL,

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
#' Calculation of weights and the intermediate weight factors used in the
#' calculation of biweight midcovariance and midcorrelation. The weights are
#' designed such that outliers get smaller weights; the weights become zero for
#' data points more than 9 median absolute deviations from the median.
#'
#' These functions are based on Equations (1) and (3) in Langfelder and Horvath
#' (2012). The weight factor is denoted \code{u} in that article.
#'
#' Langfelder and Horvath (2012) also describe the Pearson fallback and maximum
#' proportion of outliers in detail. For a full discussion of the biweight
#' midcovariance and midcorrelation, see Wilcox (2005).
#'
#' @aliases bicovWeights bicovWeightFactors bicovWeightsFromFactors
#' @param x A vector or a two-dimensional array (matrix or data frame).  If
#' two-dimensional, the weights will be calculated separately on each column.
#' @param u A vector or matrix of weight factors, usually calculated by
#' \code{bicovWeightFactors}.
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
#' @param defaultFactor Value used for factors that correspond to a finite
#' \code{x} but the weights themselves would not be finite, for example, when a
#' column in \code{x} is constant.
#' @return A vector or matrix of the same dimensions as the input \code{x}
#' giving the bisquare weights (\code{bicovWeights} and
#' \code{bicovWeightsFromFactors}) or the bisquare factors
#' (\code{bicovWeightFactors}).
#' @author Peter Langfelder
#' @seealso \code{\link{bicor}}
#' @references Langfelder P, Horvath S (2012) Fast R Functions for Robust
#' Correlations and Hierarchical Clustering Journal of Statistical Software
#' 46(11) 1-17 PMID: 23050260 PMCID: PMC3465711 Wilcox RR (2005). Introduction
#' to Robust Estimation and Hypothesis Testing. 2nd edition. Academic Press,
#' Section 9.3.8, page 399 as well as Section 3.12.1, page 83.
#' @keywords misc
#' @examples
#' x <-  rnorm(100)
#' x[1] <- 10
#' \dontrund{plot(x, bicovWeights(x))}
#' @export
bicovWeights = function(x, pearsonFallback = TRUE, maxPOutliers = 1,
                        outlierReferenceWeight = 0.5625,
                        defaultFactor = NA)
{
  referenceU = sqrt(1-sqrt(outlierReferenceWeight))
  dimX = dim(x)
  dimnamesX = dimnames(x)
  x = as.matrix(x)
  nc = ncol(x)
  nr = nrow(x)
  mx = colMedians(x, na.rm = TRUE)
  mx.mat = matrix(mx, nr, nc, byrow = TRUE)
  mads = colMads(x, constant = 1, na.rm = TRUE)
  madZero = replaceMissing(mads==0)
  if (any(madZero, na.rm = TRUE)) {
    warning(immediate. = TRUE,
            "MAD is zero in some columns of 'x'.")
    if (pearsonFallback) {
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
  if (!is.null(defaultFactor)) u[!is.finite(u)] = defaultFactor
  u
}


bicovWeights = function(x, pearsonFallback = TRUE, maxPOutliers = 1,
                        outlierReferenceWeight = 0.5625,
                        defaultWeight = 0)
{
  dimX = dim(x)
  dimnamesX = dimnames(x)
  x = as.matrix(x)
  nc = ncol(x)
  nr = nrow(x)

  u = bicovWeightFactors(x, pearsonFallback = pearsonFallback,
                         maxPOutliers = maxPOutliers,
                         outlierReferenceWeight = outlierReferenceWeight,
                         defaultFactor = NA)

  a = matrix(as.numeric(abs(u)<1), nr, nc)
  weights = a * (1-u^2)^2
  weights[is.na(x)] = defaultWeight
  weights[!is.finite(weights)] = defaultWeight
  dim(weights) = dimX
  if (!is.null(dimX)) dimnames(weights) = dimnamesX
  weights
}

#' @rdname bicovWeightFactors
#' @export
bicovWeightsFromFactors = function(u, defaultWeight = 0) {
  dimU = dim(u)
  u = as.matrix(u)
  a = matrix(as.numeric(abs(u)<1), nrow(u), ncol(u))
  weights = a * (1-u^2)^2
  weights[is.na(u)] = defaultWeight
  weights[!is.finite(weights)] = defaultWeight
  dim(weights) = dimU
  if (!is.null(dimU)) dimnames(weights) = dimnames(u)
  weights
}
