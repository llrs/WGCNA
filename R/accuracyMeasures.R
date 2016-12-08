# Accuracy measures, modified from the WGCNA version.

# Helper function: contingency table of 2 variables that will also include
# rows/columns for levels that do not appear in x or y.

.table2.allLevels <- function(x, y, levels.x = sort(unique(x)),
                             levels.y = sort(unique(y)), setNames = FALSE) {
  nx <- length(levels.x)
  ny <- length(levels.y)
  tab <- table(x, y)

  out <- matrix(0, nx, ny)
  if (setNames) {
    rownames(out) <- levels.x
    colnames(out) <- levels.y
  }
  out[match(rownames(tab), levels.x), match(colnames(tab), levels.y)] <- tab
  out
}

# accuracy measures

#' Accuracy measures for a 2x2 confusion matrix or for vectors of predicted and
#' observed values.
#'
#' The function calculates various prediction accuracy statistics for
#' predictions of binary or quantitative (continuous) responses. For binary
#' classification, the function calculates the error rate, accuracy,
#' sensitivity, specificity, positive predictive value, and other accuracy
#' measures. For quantitative prediction, the function calculates correlation,
#' R-squared, error measures, and the C-index.
#'
#' The rows of the 2x2 table tab must correspond to a test (or predicted)
#' outcome and the columns to a true outcome ("gold standard"). A table that
#' relates a predicted outcome to a true test outcome is also known as
#' confusion matrix. Warning: To correctly calculate sensitivity and
#' specificity, the positive and negative outcome must be properly specified so
#' they can be matched to the appropriate rows and columns in the confusion
#' table.
#'
#' Interchanging the negative and positive levels swaps the estimates of the
#' sensitivity and specificity but has no effect on the error rate or accuracy.
#' Specifically, denote by \code{pos} the index of the positive level in the
#' confusion table, and by \code{neg} th eindex of the negative level in the
#' confusion table.  The function then defines number of true
#' positives=TP=tab[pos, pos], no.false positives =FP=tab[pos, neg], no.false
#' negatives=FN=tab[neg, pos], no.true negatives=TN=tab[neg, neg].  Then
#' Specificity= TN/(FP+TN) Sensitivity= TP/(TP+FN) NegativePredictiveValue=
#' TN/(FN + TN) PositivePredictiveValue= TP/(TP + FP) FalsePositiveRate =
#' 1-Specificity FalseNegativeRate = 1-Sensitivity Power = Sensitivity
#' LikelihoodRatioPositive = Sensitivity / (1-Specificity)
#' LikelihoodRatioNegative = (1-Sensitivity)/Specificity. The naive error rate
#' is the error rate of a constant (naive) predictor that assigns the same
#' outcome to all samples. The prediction of the naive predictor equals the
#' most frequenly observed outcome. Example: Assume you want to predict disease
#' status and 70 percent of the observed samples have the disease. Then the
#' naive predictor has an error rate of 30 percent (since it only misclassifies
#' 30 percent of the healthy individuals).
#'
#' @param predicted either a a 2x2 confusion matrix (table) whose entries
#' contain non-negative integers, or a vector of predicted values. Predicted
#' values can be binary or quantitative (see \code{type} below). If a 2x2
#' matrix is given, it must have valid column and row names that specify the
#' levels of the predicted and observed variables whose counts the matrix is
#' giving (e.g., the function \code{\link{table}} sets the names
#' appropriately.) If it is a 2x2 table and the table contains non-negative
#' real (non-integer) numbers the function outputs a warning.
#' @param observed if \code{predicted} is a vector of predicted values, this
#' (\code{observed}) must be a vector of the same length giving the "gold
#' standard" (or observed) values. Ignored if \code{predicted} is a 2x2 table.
#' @param type character string specifying the type of the prediction problem
#' (i.e., values in the \code{predicted} and \code{observed} vectors). The
#' default \code{"auto"} decides type automatically: if \code{predicted} is a
#' 2x2 table or if the number of unique values in the concatenation of
#' \code{predicted} and \code{observed} is 2, the prediction problem (type) is
#' assumed to be binary, otherwise it is assumed to be quantitative.
#' Inconsistent specification (for example, when \code{predicted} is a 2x2
#' matrix and \code{type} is \code{"quantitative"}) trigger errors.
#' @param levels a 2-element vector specifying the two levels of binary
#' variables. Only used if \code{type} is \code{"binary"} (or \code{"auto"}
#' that results in the binary type). Defaults to either the column names of the
#' confusion matrix (if the matrix is specified) or to the sorted unique values
#' of \code{observed} and \code{opredicted}.
#' @param negativeLevel the binary value (level) that corresponds to the
#' negative outcome. Note that the default is the second of the sorted levels
#' (for example, if levels are 1,2, the default negative level is 2). Only used
#' if \code{type} is \code{"binary"} (or \code{"auto"} that results in the
#' binary type).
#' @param positiveLevel the binary value (level) that corresponds to the
#' positive outcome. Note that the default is the second of the sorted levels
#' (for example, if levels are 1,2, the default negative level is 2). Only used
#' if \code{type} is \code{"binary"} (or \code{"auto"} that results in the
#' binary type).
#' @return Data frame with two columns: \item{Measure}{this column contais
#' character strings that specify name of the accuracy measure.}
#' \item{Value}{this column contains the numeric estimates of the corresponding
#' accuracy measures.}
#' @author Steve Horvath and Peter Langfelder
#' @references See wikipedia \url{http://en.wikipedia.org/wiki/Sensitivity_and_specificity}
#' @keywords misc
#' @examples
#'
#' m <- 100
#' trueOutcome <- sample( c(1,2), m, replace=TRUE)
#' predictedOutcome <- trueOutcome
#' # now we noise half of the entries of the predicted outcome
#' predictedOutcome[ 1:(m/2)] <- sample(predictedOutcome[ 1:(m/2)])
#' tab <- table(predictedOutcome, trueOutcome)
#' accuracyMeasures(tab)
#'
#' # Should get the same result:
#' accuracyMeasures(predictedOutcome, trueOutcome)
accuracyMeasures <- function(predicted, observed = NULL,
                             type = c("auto", "binary", "quantitative"),
                             levels = if (isTRUE(all.equal(dim(predicted),
                                                    c(2,2)))) {
                                 colnames(predicted)
                             } else if (is.factor(predicted)) {
                                 sort(unique(c(as.character(predicted),
                                               as.character(observed))))
                             } else {sort(unique(c(observed, predicted)))},
                             negativeLevel = levels[2],
                             positiveLevel = levels[1]) {
    type <- match.arg(type)
    if (type == "auto") {
        if (!is.null(dim(predicted))) {
            if (all.equal(dim(predicted), c(2,2))) {
                type <- "binary"
            } else {
                stop("If supplying a matrix in 'predicted', it must be a 2x2 ",
                     "contingency table.")
            }
        } else {
            if (is.null(observed)) {
                stop("When 'predicted' is a vector, 'observed' must be given ",
                     "and have the same length as 'predicted'.")
            }
            if (length(levels) == 2) {
                type <- "binary"
            } else {
                type <- "quantitative"
            }
        }
    }

    if (type == "binary") {
        if (is.null(dim(predicted))) {
            if (is.null(observed)) {
                stop("When 'predicted' is a vector, 'observed' must be given ",
                     "and have the same length as 'predicted'.")
            }
            if ( length(predicted) != length(observed)) {
                stop("When both 'predicted' and 'observed' are given, they ",
                     "must be vectors of the same length.")
            }
            if (length(levels) != 2) {
                stop("'levels' must contain 2 entries (the possible values of ",
                     "the binary variables\n   'predicted' and 'observed').")
            }
            tab <- table(predicted, observed)
        } else {
            tab <- predicted
            if (is.null(colnames(tab)) | is.null(rownames(tab))) {
                stop("When 'predicted' is a contingency table, it must have ",
                     "valid colnames and rownames.")
            }
        }

        if (  ncol(tab) != 2 |  nrow(tab) != 2 ) {
            stop("The input table must be a 2x2 table. ")
        }
        if (negativeLevel == positiveLevel) {
            stop("'negativeLevel' and 'positiveLevel' cannot be the same.")
        }
        neg <- match(negativeLevel, colnames(tab))
        if (is.na(neg)){
            stop("Cannot find the negative level ", negativeLevel,
                 " among the colnames of the contingency table.\n   Please ",
                 "check the input and try again.")
        }
        pos <- match(positiveLevel, colnames(tab))
        if (is.na(pos)) {
            stop("Cannot find the positive level ", positiveLevel,
                 " among the colnames of the contingency table.\n   Please ",
                 "check the input and try again.")
        }

        if (sum(is.na(tab))) {
            warning("Missing data should not be present in input.\n",
                    "Suggestion: check whether NA should be coded as 0.")
        }
        is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {
            abs(x - round(x)) < tol
        }

        if (sum(!is.wholenumber(tab), na.rm = T) > 0){
            stop("The input table contains non-integers, which does not make ",
                 "sense.")
        }
        if (sum(tab < 0, na.rm = T) > 0) {
            stop("The input table cannot contain negative numbers.")
        }
        num1 <- sum(diag(tab), na.rm = T)
        denom1 <- sum(tab, na.rm = T)
        if (denom1 == 0) {
            warning("The input table has zero observations (sum of all cells ",
                    "is zero).")
        }
        TP <- tab[pos, pos]
        FP <- tab[pos, neg]
        FN <- tab[neg, pos]
        TN <- tab[neg, neg]

        error.rate <- ifelse(denom1 == 0, NA, 1 - num1/denom1)
        Accuracy <- ifelse(denom1 == 0, NA, num1/denom1 )
        Specificity <- ifelse(FP + TN == 0, NA, TN / (FP + TN))
        Sensitivity <- ifelse(TP + FN == 0, NA, TP / (TP + FN))
        NegativePredictiveValue <- ifelse(FN + TN == 0, NA, TN / (FN + TN))
        PositivePredictiveValue <- ifelse(TP + FP == 0, NA, TP / (TP + FP))
        FalsePositiveRate <- 1 - Specificity
        FalseNegativeRate <- 1 - Sensitivity
        Power <- Sensitivity
        LikelihoodRatioPositive <- ifelse(1 - Specificity == 0, NA,
                                          Sensitivity / (1 - Specificity))
        LikelihoodRatioNegative <- ifelse(Specificity == 0, NA,
                                          (1 - Sensitivity) / Specificity)
        NaiveErrorRate <- ifelse(denom1 == 0, NA,
                                 min(c(tab[pos, pos] + tab[neg, pos],
                                       tab[pos, neg]+ tab[neg, neg]))/denom1)
        out <- data.frame(
            Measure = c("Error.Rate","Accuracy", "Specificity","Sensitivity",
                        "NegativePredictiveValue", "PositivePredictiveValue",
                        "FalsePositiveRate","FalseNegativeRate","Power",
                        "LikelihoodRatioPositive","LikelihoodRatioNegative",
                        "NaiveErrorRate", "NegativeLevel", "PositiveLevel"),
            Value = c(error.rate,Accuracy, Specificity, Sensitivity,
                      NegativePredictiveValue, PositivePredictiveValue,
                      FalsePositiveRate, FalseNegativeRate, Power,
                      LikelihoodRatioPositive, LikelihoodRatioNegative,
                      NaiveErrorRate, negativeLevel, positiveLevel))
    } else if (type == "quantitative") {
        if (!is.null(dim(predicted))) {
            stop("When 'type' is \"quantitative\", 'predicted' cannot be a ",
                 "2-dimensional matrix.")
        }
        if (length(predicted)!=length(observed)) {
            stop("'predicted' and 'observed' must be vectors of the ",
                 "same length.")
        }
        cr <- cor(predicted, observed, use = 'p')
        # TODO convert the output from a data.frame to a character vector
        out <- data.frame(
            Measure = c("Cor", "R.squared", "MeanSquareError",
                        "MedianAbsoluteError", "Cindex"),
            Value = c(cr,
                      cr^2,
                      mean( (predicted-observed)^2,na.rm=TRUE),
                      median((predicted-observed)^2,na.rm=TRUE),
                      rcorr.cens(predicted,observed,outx=TRUE)[[1]]))
    }
    out
}
