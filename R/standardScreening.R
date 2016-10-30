
# standardScreeningCensoredTime ####
#' Standard Screening with regard to a Censored Time Variable
#'
#' The function standardScreeningCensoredTime computes association measures
#' between the columns of the input data datE and a censored time variable
#' (e.g. survival time). The censored time is specified using two input
#' variables "time" and "event". The event variable is binary where 1 indicates
#' that the event took place (e.g. the person died) and 0 indicates censored
#' (i.e. lost to follow up).  The function fits univariate Cox regression
#' models (one for each column of datE) and outputs a Wald test p-value, a
#' logrank p-value, corresponding local false discovery rates (known as
#' q-values, Storey et al 2004), hazard ratios. Further it reports the
#' concordance index (also know as area under the ROC curve) and optionally
#' results from dichotomizing the columns of datE.
#'
#'
#' If input option fastCalculation=TRUE, then the function outputs correlation
#' test p-values (and q-values) for correlating the columns of datE with the
#' expected hazard (if no covariate is fit). Specifically, the expected hazard
#' is defined as the deviance residual of an intercept only Cox regression
#' model. The results are very similar to those resulting from a univariate Cox
#' model where the censored time is regressed on the columns of dat.
#' Specifically, this computational speed up is facilitated by the insight that
#' the p-values resulting from a univariate Cox regression
#' coxph(Surv(time,event)~datE[,i]) are very similar to those from
#' corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples)
#'
#' @param time numeric variable showing time to event or time to last follow
#' up.
#' @param event Input variable \code{time} specifies the time to event or time
#' to last follow up. Input variable \code{event} indicates whether the event
#' happend (=1) or whether there was censoring (=0).
#' @param datExpr a data frame or matrix whose columns will be related to the
#' censored time.
#' @param percentiles numeric vector which is only used when
#' dichotomizationResults=T. Each value should lie between 0 and 1. For each
#' value specified in the vector percentiles, a binary vector will be defined
#' by dichotomizing the column value according to the corresponding quantile.
#' Next a corresponding p-value will be calculated.
#' @param dichotomizationResults logical. If this option is set to TRUE then
#' the values of the columns of datE will be dichotomized and corresponding Cox
#' regression p-values will be calculated.
#' @param qValues logical. If this option is set to TRUE (default) then
#' q-values will be calculated for the Cox regression p-values.
#' @param fastCalculation logical. If set to TRUE, the function outputs
#' correlation test p-values (and q-values) for correlating the columns of datE
#' with the expected hazard (if no covariate is fit). Specifically, the
#' expected hazard is defined as the deviance residual of an intercept only Cox
#' regression model. The results are very similar to those resulting from a
#' univariate Cox model where the censored time is regressed on the columns of
#' dat. Specifically, this computational speed up is facilitated by the insight
#' that the p-values resulting from a univariate Cox regression
#' coxph(Surv(time,event)~datE[,i]) are very similar to those from
#' corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples).
#' @return If \code{fastCalculation} is \code{FALSE}, the function outputs a
#' data frame whose rows correspond to the columns of datE and whose columns
#' report \item{ID}{column names of the input data datExpr.}
#' \item{pvalueWald}{Wald test p-value from fitting a univariate Cox regression
#' model where the censored time is regressed on each column of datExpr.}
#' \item{qValueWald}{local false discovery rate (q-value) corresponding to the
#' Wald test p-value. } \item{pvalueLogrank}{Logrank p-value resulting from the
#' Cox regression model. Also known as score test p-value. For large sample
#' sizes this sould be similar to the Wald test p-value. }
#' \item{qValueLogrank}{local false discovery rate (q-value) corresponding to
#' the Logrank test p-value. } \item{HazardRatio}{hazard ratio resulting from
#' the Cox model. If the value is larger than 1, then high values of the column
#' are associated with shorter time, e.g. increased hazard of death. A hazard
#' ratio equal to 1 means no relationship between the column and time. HR<1
#' means that high values are associated with longer time, i.e. lower hazard.}
#' \item{CI.LowerLimitHR}{Lower bound of the 95 percent confidence interval of
#' the hazard ratio. } \item{CI.UpperLimitHR}{Upper bound of the 95 percent
#' confidence interval of the hazard ratio. } \item{C.index}{concordance index,
#' also known as C-index or area under the ROC curve. Calculated with the
#' rcorr.cens option outx=TRUE (ties are ignored).}
#' \item{MinimumDichotPvalue}{This is the smallest p-value from the
#' dichotomization results. To see which dichotomized variable (and percentile)
#' corresponds to the minimum, study the following columns. }
#' \item{pValueDichot0.1}{This columns report the p-value when the column is
#' dichotomized according to the specified percentile (here 0.1). The
#' percentiles are specified in the input option percentiles. }
#' \item{pvalueDeviance}{The p-value resulting from using a correlation test to
#' relate the expected hazard (deviance residual) with each (undichotomized)
#' column of datE. Specifically, the Fisher transformation is used to calculate
#' the p-value for the Pearson correlation. The resulting p-value should be
#' very similar to that of a univariate Cox regression model.}
#' \item{qvalueDeviance}{Local false discovery rate (q-value) corresponding to
#' pvalueDeviance.} \item{corDeviance}{Pearson correlation between the expected
#' hazard (deviance residual) with each (undichotomized) column of datExpr.}
#' @author Steve Horvath
#' @keywords misc
standardScreeningCensoredTime <- function(time,
                                          event,
                                          datExpr,
                                          percentiles = seq(from = 0.1, to = 0.9, by = 0.2),
                                          dichotomizationResults = FALSE,
                                          qValues = TRUE,
                                          fastCalculation = TRUE) {
    datExpr = data.frame(datExpr)
    no.Columns = dim(as.matrix(datExpr))[[2]]
    m = dim(as.matrix(datExpr))[[1]]
    if (length(time)  != m)
        stop(
            "The length of the time variable does not equal the number of ",
            "rows of datExpr.\nConsider transposing datExpr."
        )
    if (length(event)  != m)
        stop(
            "The length of the event variable does not equal the number of ",
            "rows of datExpr.\nConsider transposing datExpr."
        )
    if (fastCalculation) {
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0,
                           2 * CumHazard,
                           -2 * log(CumHazard)  +
                               2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr,
                                     use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance,
                                         nSamples = no.nonMissing)
        qvalueDeviance = rep(NA, length(pvalueDeviance))
        rest1 = !is.na(pvalueDeviance)
        qvalueDeviance [rest1] = qvalue(pvalueDeviance [rest1])$qvalues

        datout = data.frame(ID = dimnames(datExpr)[[2]],
                            pvalueDeviance,
                            qvalueDeviance,
                            corDeviance)
    }
    if (!fastCalculation) {
        pvalueWald = rep(NA, no.Columns)
        HazardRatio = rep(NA, no.Columns)
        CI.UpperLimitHR = rep(NA, no.Columns)
        CI.LowerLimitHR = rep(NA, no.Columns)
        C.index = rep(NA, no.Columns)
        pvalueLogrank = rep(NA, no.Columns)
        pValuesDichotomized = data.frame(matrix(nrow = no.Columns,
                                                ncol = length(percentiles)))
        names(pValuesDichotomized) = paste0("pValueDichotPercentile",
                                            as.character(percentiles))
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0,
                           2 * CumHazard,
                           -2 * log(CumHazard)  +
                               2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr,
                                     use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance,
                                         nSamples = no.nonMissing)


        for (i in 1:no.Columns) {
            Column = as.numeric(as.matrix(datExpr[, i]))
            var1 = var(Column, na.rm = TRUE)
            if (var1 == 0 | is.na(var1)) {
                pvalueWald[i] = NA
                pvalueLogrank[i] = NA
                HazardRatio[i] = NA
                CI.UpperLimitHR[i] = NA
                CI.LowerLimitHR[i] = NA
                C.index[i] = NA
            }  # end of              if (var1 == 0 | is.na(var1))
            if (var1  != 0 & !is.na(var1)) {
                cox1 = summary(coxph(Surv(time, event) ~ Column,
                                     na.action = na.exclude))
                pvalueWald[i] = cox1$coef[5]
                pvalueLogrank[i] = cox1$sctest[[3]]
                HazardRatio[i] = exp(cox1$coef[1])
                CI.UpperLimitHR[i] = exp(cox1$coef[1] + 1.96  *
                                             cox1$coef[3])
                CI.LowerLimitHR[i] = exp(cox1$coef[1] - 1.96  *
                                             cox1$coef[3])
                C.index[i] = rcorr.cens(Column, Surv(time, event),
                                        outx = TRUE)[[1]]
            } # end of   if (var1  != 0 & !is.na(var1))


            if (dichotomizationResults) {
                quantilesE = as.numeric(quantile(Column, prob = percentiles))
                for (j in 1:length(quantilesE)) {
                    ColumnDichot = I(Column > quantilesE[j])
                    var1 = var(ColumnDichot, na.rm = TRUE)
                    if (var1 == 0 | is.na(var1)) {
                        pValuesDichotomized[i, j] = NA
                    } # end of if
                    if (var1  != 0 & !is.na(var1)) {
                        coxh = summary(coxph(
                            Surv(time, event) ~
                                ColumnDichot,
                            na.action = na.exclude
                        ))
                        pValuesDichotomized[i, j] = coxh$coef[5]
                    } # end of if
                } # end of for (j
                MinimumDichotPvalue = apply(pValuesDichotomized,
                                            1, min, na.rm = TRUE)
            } # end of if (dichotomizationResults)



            if (!qValues) {
                datout = data.frame(
                    ID = dimnames(datExpr)[[2]],
                    pvalueWald,
                    pvalueLogrank,
                    pvalueDeviance,
                    corDeviance,
                    HazardRatio,
                    CI.LowerLimitHR,
                    CI.UpperLimitHR,
                    C.index
                )
            }      # end of      if (!qValues) {

        } # end of for (i in 1:no.Columns)


        if (qValues) {
            qvalueWald = rep(NA, length(pvalueWald))
            rest1 = !is.na(pvalueWald)
            qvalueWald [rest1] = qvalue(pvalueWald[rest1])$qvalues

            qvalueLogrank = rep(NA, length(pvalueLogrank))
            rest1 = !is.na(pvalueLogrank)
            qvalueLogrank [rest1] = qvalue(pvalueLogrank[rest1])$qvalues

            qvalueDeviance = rep(NA, length(pvalueDeviance))
            rest1 = !is.na(pvalueDeviance)
            qvalueDeviance [rest1] = qvalue(pvalueDeviance[rest1])$qvalues

            datout = data.frame(
                ID = dimnames(datExpr)[[2]],
                pvalueWald,
                qvalueWald,
                pvalueLogrank,
                qvalueLogrank,
                pvalueDeviance,
                qvalueDeviance,
                corDeviance,
                HazardRatio,
                CI.LowerLimitHR,
                CI.UpperLimitHR,
                C.index
            )
        } # end of  if (qValues)


        if (dichotomizationResults) {
            datout = data.frame(datout, MinimumDichotPvalue,
                                pValuesDichotomized)
        }
    }
    datout
} # end of function standardScreeningCensoredTime

# standardScreeningNumericTrait ####
#' Standard screening for numeric traits
#'
#' Standard screening for numeric traits based on Pearson correlation.
#'
#' The function calculates the correlations, associated p-values, area under
#' the ROC, and q-values
#'
#' @param datExpr data frame containing expression data (or more generally
#' variables to be screened), with rows corresponding to samples and columns to
#' genes (variables)
#' @param yNumeric a numeric vector giving the trait measurements for each
#' sample
#' @param corFnc correlation function.  Defaults to Pearson correlation but can
#' also be \code{\link{bicor}}.
#' @param corOptions list specifying additional arguments to be passed to the
#' correlation function given by \code{corFnc}.
#' @param alternative alternative hypothesis for the correlation test
#' @param qValues logical: should q-values be calculated?
#' @param areaUnderROC logical: should are under the receiver-operating curve
#' be calculated?
#' @return
#'
#' Data frame with the following components:
#'
#' \item{ID }{Gene (or variable) identifiers copied from
#' \code{colnames(datExpr)}} \item{cor}{correlations of all genes with the
#' trait} \item{Z}{Fisher Z statistics corresponding to the correlations}
#' \item{pvalueStudent }{Student p-values of the correlations}
#' \item{qvalueStudent }{(if input \code{qValues==TRUE}) q-values of the
#' correlations calculated from the p-values} \item{AreaUnderROC }{(if input
#' \code{areaUnderROC==TRUE}) area under the ROC} \item{nPresentSamples}{number
#' of samples present for the calculation of each association. }
#' @author Steve Horvath
#' @seealso \code{\link{standardScreeningBinaryTrait}},
#' \code{\link{standardScreeningCensoredTime}}
#' @keywords misc
standardScreeningNumericTrait <- function(datExpr,
                                          yNumeric,
                                          corFnc = cor,
                                          corOptions = list(use = 'p'),
                                          alternative = c("two.sided", "less", "greater"),
                                          qValues = TRUE,
                                          areaUnderROC = TRUE) {
    datExpr = as.matrix(datExpr)
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    if (length(yNumeric)  != nSamples)
        stop("the length of the sample trait y does not equal the number of ",
             "rows of datExpr")
    corPearson = rep(NA, nGenes)
    pvalueStudent = rep(NA, nGenes)
    AreaUnderROC = rep(NA, nGenes)
    nPresent = Z = rep(NA, nGenes)

    corFnc = match.fun(corFnc)
    corOptions$y = yNumeric
    corOptions$x = as.matrix(datExpr)
    cp = do.call(corFnc, corOptions)
    corPearson = as.numeric(cp)

    finMat = !is.na(datExpr)
    np = t(finMat) %*% (!is.na(as.matrix(yNumeric)))

    nPresent = as.numeric(np)

    ia = match.arg(alternative)
    T = sqrt(np - 2) * corPearson / sqrt(1 - corPearson ^ 2)
    if (ia == "two.sided") {
        p = 2 * pt(abs(T), np - 2, lower.tail = FALSE)
    }
    else if (ia == "less") {
        p = pt(T, np - 2, lower.tail = TRUE)
    }
    else if (ia == "greater") {
        p = pt(T, np - 2, lower.tail = FALSE)
    }
    pvalueStudent = as.numeric(p)

    Z = 0.5 * log((1 + corPearson) / (1 - corPearson)) * sqrt(nPresent  - 2)

    if (areaUnderROC)
        for (i in 1:dim(datExpr)[[2]])
        {
            AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = TRUE)[[1]]
        }

    q.Student = rep(NA, length(pvalueStudent))
    rest1 = !is.na(pvalueStudent)
    if (qValues)
    {
        x = try({
            q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues
        },
        silent = TRUE)
        if (inherits(x, "try - error"))
            printFlush(
                paste(
                    "Warning in standardScreeningNumericTrait: function qvalue ",
                    "returned an error.\n",
                    "The returned qvalues will be invalid. The qvalue error: ",
                    x,
                    "\n"
                )
            )
    }
    if (is.null(colnames(datExpr))) {
        ID = paste0("Variable.", 1:ncol(datExpr))
    } else
        ID = colnames(datExpr)

    output = data.frame(
        ID = ID,
        cor = corPearson,
        Z = Z,
        pvalueStudent = pvalueStudent
    )
    if (qValues)
        output$qvalueStudent = q.Student
    if (areaUnderROC)
        output$AreaUnderROC = AreaUnderROC

    output$nPresentSamples = nPresent

    output
}
