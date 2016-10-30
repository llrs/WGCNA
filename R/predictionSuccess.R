
# corPredictionSuccess ####
# The function corPredictionSuccess can be used to determine which method is
# best for predicting correlations in a new test set. corTestSet should be a
# vector of correlations in the test set.
# The parameter topNumber specifies that the top number most positive and the
# top most negative predicted correlations
#  TopNumber is a vector of integers.
# corPrediction should be a data frame of predictions for the correlations.
# Output a list with the following components:
# meancorTestSetPositive = mean test set correlation among the topNumber of
# genes which are predicted to have positive correlations.
# meancorTestSetNegative = mean test set correlation among the topNumber of
# genes which are predicted to have negative correlations.
# meancorTestSetOverall = (meancorTestSetPositive - meancorTestSetNegative)/2
#' Quantification of success of gene screening
#'
#' This function calculates the success of gene screening.
#'
#' For each column in \code{corPrediction}, the function evaluates the mean
#' \code{corTestSet} for the number of top genes (ranked by the column in
#' \code{corPrediction}) given in \code{topNumber}. The higher the mean
#' \code{corTestSet} (for positive \code{corPrediction}) or negative (for
#' negative \code{corPrediction}), the more successful the prediction.
#'
#' @param corPrediction a vector or a matrix of prediction statistics
#' @param corTestSet correlation or other statistics on test set
#' @param topNumber a vector of the number of top genes to consider
#' @return \item{meancorTestSetOverall }{ difference of
#' \code{meancorTestSetPositive} and \code{meancorTestSetNegative} below }
#' \item{meancorTestSetPositive}{ mean \code{corTestSet} on top genes with
#' positive \code{corPrediction} } \item{meancorTestSetNegative}{ mean
#' \code{corTestSet} on top genes with negative \code{corPrediction} } ...
#' @author Steve Horvath
#' @seealso \code{\link{relativeCorPredictionSuccess}}
#' @keywords misc
corPredictionSuccess <- function(corPrediction, corTestSet, topNumber = 100) {
    nPredictors = dim(as.matrix(corPrediction))[[2]]
    nGenes = dim(as.matrix(corPrediction))[[1]]
    if (length(as.numeric(corTestSet)) != nGenes) {
        stop("non - compatible dimensions of 'corPrediction' and ",
             "'corTestSet'")
    }
    out1 = rep(NA, nPredictors)
    meancorTestSetPositive = matrix(ncol = nPredictors,
                                    nrow = length(topNumber))
    meancorTestSetNegative = matrix(ncol = nPredictors,
                                    nrow = length(topNumber))
    for (i in c(1:nPredictors)) {
        rankpositive = rank(-as.matrix(corPrediction)[, i],
                            ties.method = "first")
        ranknegative = rank(as.matrix(corPrediction)[, i],
                            ties.method = "first")
        for (j in c(1:length(topNumber))) {
            meancorTestSetPositive[j, i] = mean(corTestSet[rankpositive <= topNumber[j]], na.rm = TRUE)
            meancorTestSetNegative[j, i] = mean(corTestSet[ranknegative <= topNumber[j]], na.rm = TRUE)
        } # end of j loop over topNumber
    } # end of i loop over predictors
    meancorTestSetOverall = data.frame((meancorTestSetPositive - meancorTestSetNegative) /
                                           2)
    dimnames(meancorTestSetOverall)[[2]] = names(data.frame(corPrediction))
    meancorTestSetOverall = data.frame(topNumber = topNumber,
                                       meancorTestSetOverall)
    meancorTestSetPositive = data.frame(meancorTestSetPositive)
    dimnames(meancorTestSetPositive)[[2]] = names(data.frame(corPrediction))
    meancorTestSetPositive = data.frame(topNumber = topNumber,
                                        meancorTestSetPositive)
    meancorTestSetNegative = data.frame(meancorTestSetNegative)
    dimnames(meancorTestSetNegative)[[2]] = names(data.frame(corPrediction))
    meancorTestSetNegative = data.frame(topNumber = topNumber,
                                        meancorTestSetNegative)
    datout = list(
        meancorTestSetOverall = meancorTestSetOverall,
        meancorTestSetPositive = meancorTestSetPositive,
        meancorTestSetNegative  = meancorTestSetNegative
    )
    datout
} # end of function corPredictionSuccess

# relativeCorPredictionSuccess ####
# The function relativeCorPredictionSuccess can be used to test whether a gene
# screening method is significantly better than a standard method.
# For each gene screening method (column of corPredictionNew) it provides a
# Kruskal Wallis test p-value for comparison with the vector
# corPredictionStandard, TopNumber is a vector of integers.
# corTestSet should be a vector of correlations in the test set.
# corPredictionNew should be a data frame of predictions for the correlations.
# corPredictionStandard should be the standard prediction (correlation in the
# training data).
# The function outputs a p-value for the Kruskal test that the new correlation
# prediction methods outperform the standard correlation prediction method.
#' Compare prediction success
#'
#' Compare prediction success of several gene screening methods.
#'
#'
#' @param corPredictionNew Matrix of predictor statistics
#' @param corPredictionStandard Reference presdictor statistics
#' @param corTestSet Correlations of predictor variables with trait in test set
#' @param topNumber A vector giving the numbers of top genes to consider
#' @return Data frame with components \item{topNumber}{copy of the input
#' \code{topNumber}} \item{kruskalp}{Kruskal-Wallis p-values}
#' @author Steve Horvath
#' @seealso \code{\link{corPredictionSuccess}}
#' @keywords misc
relativeCorPredictionSuccess <- function(corPredictionNew,
                                         corPredictionStandard,
                                         corTestSet,
                                         topNumber = 100) {
    nPredictors = dim(as.matrix(corPredictionNew))[[2]]
    nGenes = dim(as.matrix(corPredictionNew))[[1]]
    if (length(as.numeric(corTestSet)) != nGenes)
        stop("non - compatible dimensions of 'corPrediction' and 'corTestSet'.")
    if (length(as.numeric(corTestSet)) != length(corPredictionStandard))
        stop("non - compatible dimensions of 'corTestSet' and
             'corPredictionStandard'.")
    kruskalp = matrix(nrow = length(topNumber), ncol = nPredictors)
    for (i in c(1:nPredictors)) {
        rankhighNew = rank(-as.matrix(corPredictionNew)[, i],
                           ties.method = "first")
        ranklowNew = rank(as.matrix(corPredictionNew)[, i],
                          ties.method = "first")
        for (j in c(1:length(topNumber))) {
            highCorNew = as.numeric(corTestSet[rankhighNew <= topNumber[j]])
            lowCorNew = as.numeric(corTestSet[ranklowNew <= topNumber[j]])
            highCorStandard = as.numeric(corTestSet[rank(-as.numeric(corPredictionStandard),
                                                         ties.method = "first") <= topNumber[j]])
            lowCorStandard = as.numeric(corTestSet[rank(as.numeric(corPredictionStandard), ties.method = "first")
                                                   <= topNumber[j]])
            signedCorNew = c(highCorNew,-lowCorNew)
            signedCorStandard = c(highCorStandard,-lowCorStandard)
            x1 = c(signedCorNew, signedCorStandard)
            Grouping = rep(c(2, 1), c(length(signedCorNew),
                                      length(signedCorStandard)))
            sign1 = sign(cor(Grouping, x1, use = "p"))
            if (sign1 == 0)
                sign1 = 1
            kruskalp[j, i] = kruskal.test(x = x1, g = Grouping)$p.value * sign1
            #print(names(data.frame(corPredictionNew))[[i]])
            #print(paste("This correlation is positive if the new method is
            #better than the old method" ,
            # signif (cor(Grouping, x1, use = "p"), 3)))
        } # end of j loop
    } # end of i loop
    kruskalp[kruskalp < 0] = 1
    kruskalp = data.frame(kruskalp)
    dimnames(kruskalp)[[2]] = paste0(names(data.frame(corPredictionNew)),
                                     ".kruskalP")
    kruskalp = data.frame(topNumber = topNumber, kruskalp)
    kruskalp
} # end of function relativeCorPredictionSuccess
