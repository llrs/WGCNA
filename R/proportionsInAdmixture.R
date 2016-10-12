#' Estimate the proportion of pure populations in an admixed population based
#' on marker expression values.
#' 
#' Assume that \code{datE.Admixture} provides the expression values from a
#' mixture of cell types (admixed population) and you want to estimate the
#' proportion of each pure cell type in the mixed samples (rows of
#' \code{datE.Admixture}).  The function allows you to do this as long as you
#' provide a data frame \code{MarkerMeansPure} that reports the mean expression
#' values of markers in each of the pure cell types.
#' 
#' The methods implemented in this function were motivated by the gene
#' expression deconvolution approach described by Abbas et al (2009), Lu et al
#' (2003), Wang et al (2006). This approach can be used to predict the
#' proportions of (pure) cells in a complex tissue, e.g. the proportion of
#' blood cell types in whole blood. To define the markers, you may need to have
#' expression data from pure populations. Then you can define markers based on
#' a significant t-test or ANOVA across the pure populations. Next use the pure
#' population data to estimate corresponding mean expression values. Hopefully,
#' the array platforms and normalization methods for
#' \code{datE.MarkersAdmixtureTranspose} and \code{MarkerMeansPure} are
#' comparable. When dealing with Affymetrix data: we have successfully used it
#' on untransformed MAS5 data. For statisticians: To estimate the proportions,
#' we use the coefficients of a linear model. Specifically: \code{datCoef=
#' t(lm(datE.MarkersAdmixtureTranspose
#' ~MarkerMeansPure[,-1])$coefficients[-1,])} where \code{datCoef} is a matrix
#' whose rows correspond to the mixed samples (rows of \code{datE.Admixture})
#' and the columns correspond to pure populations (e.g. cell types), i.e. the
#' columns of \code{MarkerMeansPure[,-1]}. More details can be found in Abbas
#' et al (2009).
#' 
#' @param MarkerMeansPure is a data frame whose first column reports the name
#' of the marker and the remaining columns report the mean values of the
#' markers in each of the pure populations. The function will estimate the
#' proportion of pure cells which correspond to columns 2 through of
#' \code{dim(MarkerMeansPure)[[2]]} of \code{MarkerMeansPure}. Rows that
#' contain missing values (NA) will be removed.
#' @param datE.Admixture is a data frame of expression data, e.g. the columns
#' of \code{datE.Admixture} could correspond to thousands of genes.  The rows
#' of \code{datE.Admixture} correspond to the admixed samples for which the
#' function estimates the proportions of pure populations.  Some of the markers
#' specified in the first column of \code{MarkerMeansPure} should correspond to
#' column names of \code{datE.Admixture}.
#' @param calculateConditionNumber logical. Default is FALSE. If set to TRUE
#' then it uses the \code{kappa} function to calculates the condition number of
#' the matrix \code{MarkerMeansPure[,-1]}.  This allows one to determine
#' whether the linear model for estimating the proportions is well specified.
#' Type \code{help(kappa)} to learn more.  \code{kappa()} computes by default
#' (an estimate of) the 2-norm condition number of a matrix or of the R matrix
#' of a QR decomposition, perhaps of a linear fit.
#' @param coefToProportion logical. By default, it is set to TRUE. When
#' estimating the proportions the function fits a multivariate linear model.
#' Ideally, the coefficients of the linear model correspond to the proportions
#' in the admixed samples. But sometimes the coefficients take on negative
#' values or do not sum to 1. If \code{coefToProportion=TRUE} then negative
#' coefficients will be set to 0 and the remaining coefficients will be scaled
#' so that they sum to 1.
#' @return A list with the following components
#' \item{PredictedProportions}{data frame that contains the predicted
#' proportions. The rows of \code{PredictedProportions} correspond to the
#' admixed samples, i.e. the rows of \code{datE.Admixture}. The columns of
#' \code{PredictedProportions} correspond to the pure populations, i.e. the
#' columns of \code{MarkerMeansPure[,-1].} } \item{datCoef=datCoef}{data frame
#' of numbers that is analogous to \code{PredictedProportions}. In general,
#' \code{datCoef} will only be different from \code{PredictedProportions} if
#' \code{coefToProportion=TRUE}. See the description of \code{coefToProportion}
#' } \item{conditionNumber}{This is the condition number resulting from the
#' \code{kappa} function. See the description of calculateConditionNumber. }
#' \item{markersUsed}{vector of character strings that contains the subset of
#' marker names (specified in the first column of \code{MarkerMeansPure}) that
#' match column names of \code{datE.Admixture} and that contain non-missing
#' pure mean values. }
#' @note This function can be considered a wrapper of the \code{lm} function.
#' @author Steve Horvath, Chaochao Cai
#' @seealso \code{\link{lm}}, \code{\link{kappa}}
#' @references Abbas AR, Wolslegel K, Seshasayee D, Modrusan Z, Clark HF (2009)
#' Deconvolution of Blood Microarray Data Identifies Cellular Activation
#' Patterns in Systemic Lupus Erythematosus. PLoS ONE 4(7): e6098.
#' doi:10.1371/journal.pone.0006098
#' 
#' Lu P, Nakorchevskiy A, Marcotte EM (2003) Expression deconvolution: a
#' reinterpretation of DNA microarray data reveals dynamic changes in cell
#' populations. Proc Natl Acad Sci U S A 100: 10370-10375.
#' 
#' Wang M, Master SR, Chodosh LA (2006) Computational expression deconvolution
#' in a complex mammalian organ. BMC Bioinformatics 7: 328.
#' @keywords misc
proportionsInAdmixture<-function 
(MarkerMeansPure, 
datE.Admixture, 
calculateConditionNumber = FALSE,
coefToProportion = TRUE)
{
    datE.Admixture = data.frame(datE.Admixture)
    if (sum(is.na(names(datE.Admixture))) > 0) {
        warning("Some of the column names of datE.Admixture are missing. Recommendation: check or assign column names. But for your convenience, we remove the corresponding columns of datE.Admixture from the analysis.")
        datE.Admixture = datE.Admixture[, !is.na(names(datE.Admixture))]
    }
    if (sum(names(datE.Admixture) == "") > 0) {
        warning("Some of the column names of datE.Admixture are missing. Recommendation: check or assign column names.  But for your convenience, we remove the corresponding columns of datE.Admixture from the analysis.")
        datE.Admixture = datE.Admixture[, names(datE.Admixture) != 
            ""]
    }
    MarkerID = MarkerMeansPure[, 1]
    if (sum(is.na(MarkerID)) > 0) {
        warning("Some of the marker are missing (NA). Recommendation: check the first column of the input MarkerMeansPure. It should contain marker names. But for your convenience, we remove the corresponding markers from the analysis.")
        MarkerMeansPure = MarkerMeansPure[!is.na(MarkerID), ]
        MarkerID = MarkerMeansPure[, 1]
    }
    if (sum(MarkerID == "", na.rm = T) > 0) {
        warning("Some of the marker names are empty strings. Recommendation: check the first column of the input MarkerMeansPure. It should contain marker names.  But for your convenience, we remove the corresponding markers from the analysis.")
        MarkerMeansPure = MarkerMeansPure[MarkerID != "", ]
        MarkerID = MarkerMeansPure[, 1]
    }
    noMissingValuesMarker = as.numeric(apply(is.na(MarkerMeansPure[, 
        -1]), 1, sum))
    if (max(noMissingValuesMarker, na.rm = T) > 0) {
        warning("Some of the markers  (rows of MarkerMeansPure) contain missing values. This is problematic.\nFor your convenience, we remove the corresponding markers (rows) from the analysis.")
        MarkerMeansPure = MarkerMeansPure[noMissingValuesMarker == 
            0, ]
        MarkerID = MarkerMeansPure[, 1]
    }
    match1 = match(MarkerID, names(datE.Admixture))
    match1 = match1[!is.na(match1)]
    if (length(match1) == 0) 
        stop("None of the marker names correspond to column names of the input datE.Admixture. Possible solutions: Transpose datE.Admixture or MarkerMeansPure. Or make sure to assign suitable names to the columns of datE.Admixture, e.g. as follows dimnames(datE.Admixture)[[2]]=GeneSymbols.")
    if (length(match1) < dim(MarkerMeansPure)[[1]]) {
        warning(paste("Only", length(match1), "out of ", dim(MarkerMeansPure)[[1]], 
            "rows of MarkerMeansPure correspond to columns of datE.Admixture. \nIf this suprises you, check the the first column of MarkerMeansPure or the column names of datE.Admixture. \nThe output contains a list of markers that could be identified."))
    }
    datE.MarkersAdmixtureTranspose = t(datE.Admixture[, match1])
    match2 = match(names(datE.Admixture)[match1], MarkerID)
    match2 = match2[!is.na(match2)]
    MarkerMeansPure = MarkerMeansPure[match2, ]
    if (sum(as.character(MarkerMeansPure[, 1]) != dimnames(datE.MarkersAdmixtureTranspose)[[1]], 
        na.rm = T) > 0) 
        stop("I am sorry but things do not line up. Maybe you need to look inside the R code. Specifically,\nas.character(MarkerMeansPure) != dimnames(datE.MarkersAdmixtureTranspose)[[1]]")
    conditionNumber = NA
    if (dim(MarkerMeansPure)[[2]] == 2) {
        A = as.matrix(MarkerMeansPure[, -1], ncol = 1)
    }
    else {
        A = as.matrix(MarkerMeansPure[, -1])
    }
    if (dim(as.matrix(A))[[2]] > 1 & dim(as.matrix(A))[[1]] > 
        1 & calculateConditionNumber) {
        conditionNumber = kappa(A)
    }
    datCoef = t(lm(datE.MarkersAdmixtureTranspose ~ A)$coefficients[-1, 
        ])
    coef2prop = function(coef) {
        prop = rep(NA, length(coef))
        coef[coef < 0] = 0
        if (sum(coef, na.rm = T) > 0 & !is.na(sum(coef, na.rm = T))) {
            prop = coef/sum(coef, na.rm = T)
        }
        prop
    }
    if (coefToProportion) {
        PredictedProportions = data.frame(t(apply(datCoef, 1, 
            coef2prop)))
    }
    else {
        PredictedProportions = datCoef
    }
    dimnames(PredictedProportions)[[1]] = dimnames(datE.Admixture)[[1]]
    out = list(PredictedProportions = PredictedProportions, datCoef = datCoef, 
        conditionNumber = conditionNumber, markersUsed = as.character(MarkerMeansPure[, 
            1]))
    out
}
