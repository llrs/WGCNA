
# standardScreeningBinaryTrait ####
#' Standard screening for binatry traits
#'
#' The function standardScreeningBinaryTrait computes widely used statistics
#' for relating the columns of the input data frame (argument datE) to a binary
#' sample trait (argument y). The statistics include Student t-test p-value and
#' the corresponding local false discovery rate (known as q-value, Storey et al
#' 2004), the fold change, the area under the ROC curve (also known as
#' C-index), mean values etc. If the input option KruskalTest is set to TRUE,
#' it also computes the Kruskal Wallist test p-value and corresponding q-value.
#' The Kruskal Wallis test is a non-parametric, rank-based group comparison
#' test.
#'
#'
#' @param datExpr a data frame or matrix whose columns will be related to the
#' binary trait
#' @param y a binary vector whose length (number of components) equals the
#' number of rows of datE
#' @param corFnc correlation function. Defaults to Pearson correlation.
#' @param corOptions a list specifying options to corFnc. An empty list must be
#' specified as \code{list()} (supplying \code{NULL} instead will trigger an
#' error).
#' @param kruskalTest logical: should the Kruskal test be performed?
#' @param qValues logical: should the q-values be calculated?
#' @param var.equal logical input parameter for the Student t-test. It
#' indicates whether to treat the two variances (corresponding to the binary
#' grouping) are being equal. If TRUE then the pooled variance is used to
#' estimate the variance otherwise the Welch (or Satterthwaite) approximation
#' to the degrees of freedom is used. Warning: here the default value is TRUE
#' which is different from the default value of t.test. Type help(t.test) for
#' more details.
#' @param na.action character string for the Student t-test: indicates what
#' should happen when the data contain missing values NAs.
#' @param getAreaUnderROC logical: should area under the ROC curve be
#' calculated? The calculation slows the function down somewhat.
#' @return A data frame whose rows correspond to the columns of datE and whose
#' columns report \item{ID}{column names of the input \code{datExpr}.}
#' \item{corPearson}{pearson correlation with a binary numeric version of the
#' input variable. The numeric variable equals 1 for level 1 and 2 for level 2.
#' The levels are given by levels(factor(y)).} \item{t.Student}{Student's
#' t-test statistic} \item{pvalueStudent}{two-sided Student t-test p-value.}
#' \item{qvalueStudent}{(if input \code{qValues==TRUE}) q-value (local false
#' discovery rate) based on the Student T-test p-value (Storey et al 2004).}
#' \item{foldChange}{a (signed) ratio of mean values. If the mean in the first
#' group (corresponding to level 1) is larger than that of the second group, it
#' equals meanFirstGroup/meanSecondGroup. But if the mean of the second group
#' is larger than that of the first group it equals
#' -meanSecondGroup/meanFirstGroup (notice the minus sign).}
#' \item{meanFirstGroup}{means of columns in input \code{datExpr} across
#' samples in the first group.} \item{meanSecondGroup}{means of columns in
#' input \code{datExpr} across samples in the second group.}
#'
#' \item{SE.FirstGroup}{standard errors of columns in input \code{datExpr}
#' across samples in the first group. Recall that SE(x)=sqrt(var(x)/n) where n
#' is the number of non-missing values of x. } \item{SE.SecondGroup}{standard
#' errors of columns in input \code{datExpr} across samples in the second
#' group.} \item{areaUnderROC}{the area under the ROC, also known as the
#' concordance index or C.index. This is a measure of discriminatory power. The
#' measure lies between 0 and 1 where 0.5 indicates no discriminatory power. 0
#' indicates that the "opposite" predictor has perfect discriminatory power. To
#' compute it we use the function \link[Hmisc]{rcorr.cens} with
#' \code{outx=TRUE} (from Frank Harrel's package Hmisc). Only present if input
#' \code{getAreUnderROC} is \code{TRUE}.} \item{nPresentSamples}{number of
#' samples with finite measurements for each gene.}
#'
#' If input \code{kruskalTest} is \code{TRUE}, the following columns further
#' summarize results of Kruskal-Wallis test: \item{stat.Kruskal}{Kruskal-Wallis
#' test statistic.} \item{stat.Kruskal.signed}{(Warning: experimental)
#' Kruskal-Wallis test statistic including a sign that indicates whether the
#' average rank is higher in second group (positive) or first group (negative).
#' } \item{pvaluekruskal}{Kruskal-Wallis test p-values.}
#' \item{qkruskal}{q-values corresponding to the Kruskal-Wallis test p-value
#' (if input \code{qValues==TRUE}).}
#' @author Steve Horvath
#' @references Storey JD, Taylor JE, and Siegmund D. (2004) Strong control,
#' conservative point estimation, and simultaneous conservative consistency of
#' false discovery rates: A unified approach. Journal of the Royal Statistical
#' Society, Series B, 66: 187-205.
#' @keywords misc
#' @examples
#'
#'
#' require(survival) # For is.Surv in rcorr.cens
#' m <- 50
#' y <- sample(c(1,2), m, replace = TRUE)
#' datExprSignal <- simulateModule(scale(y), 30)
#' datExprNoise <- simulateModule(rnorm(m), 150)
#' datExpr <- data.frame(datExprSignal, datExprNoise)
#'
#' Result1 <- standardScreeningBinaryTrait(datExpr, y)
#' Result1[1:5, ]
#'
#' # use unequal variances and calculate q-values
#' Result2 <- standardScreeningBinaryTrait(datExpr, y,
#'                                        var.equal = FALSE, qValue = TRUE)
#' Result2[1:5, ]
#'
#' # calculate Kruskal Wallis test and q-values
#' Result3 <- standardScreeningBinaryTrait(datExpr, y, kruskalTest = TRUE,
#'                                        qValue=TRUE)
#' Result3[1:5, ]
standardScreeningBinaryTrait <- function(datExpr, y,
                                         corFnc = cor,
                                         corOptions = list(use = 'p'),
                                         kruskalTest=FALSE,
                                         qValues = FALSE,
                                         var.equal=FALSE,
                                         na.action="na.exclude",
                                         getAreaUnderROC = TRUE) {

   datExpr=data.frame(datExpr)
   levelsy=levels(factor(y))
   if (length(levelsy)>2 )
     stop("The sample trait y contains more than 2 levels. Please input a binary variable y")
   if (length(levelsy)==1 )
     stop("The sample trait y is constant. Please input a binary sample trait with some variation.")
   yNumeric=as.numeric(factor(y))
   if (length(yNumeric) !=dim(datExpr)[[1]] )
     stop("the length of the sample trait y does not equal the number of rows of datExpr")
   pvalueStudent = t.Student = Z.Student = rep(NA, dim(datExpr)[[2]] )
   pvaluekruskal = stat.Kruskal = Z.Kruskal = sign.Kruskal = rep(NA, dim(datExpr)[[2]] )
   nPresent = rep(0, dim(datExpr)[[2]] )
   AreaUnderROC=rep(NA, dim(datExpr)[[2]] )

  if (var.equal)
     printFlush(paste("Warning: T-test that assumes equal variances in each group is requested.\n",
                      "This is not the default option for t.test. We recommend to use var.equal=FALSE."))

  corFnc = match.fun(corFnc)
  corOptions$y = yNumeric
  corOptions$x = datExpr
  corPearson=as.numeric(do.call(corFnc, corOptions))
  nGenes = dim(datExpr)[[2]]
  nPresent1 = as.numeric( t(as.matrix(!is.na(yNumeric) & yNumeric==1)) %*% ! is.na(datExpr) )
  nPresent2 = as.numeric( t(as.matrix(!is.na(yNumeric) & yNumeric==2)) %*% ! is.na(datExpr) )
  nPresent = nPresent1 + nPresent2

  for (i in 1:nGenes) {
        no.present1 = nPresent1[i]
        no.present2 = nPresent2[i]
        no.present = nPresent[i]
        if (no.present1<2 | no.present2<2 )
        {
           pvalueStudent[i]= t.Student[i] = NA
        } else {
          tst = try(t.test( as.numeric(datExpr[,i])~yNumeric,var.equal=var.equal,na.action=na.action),
                    silent = TRUE)
          if (!inherits(tst, "try-error"))
          {
            pvalueStudent[i] = tst$p.value
            t.Student[i] = -tst$statistic
            # The - sign above is intentional to make the sign of t consistent with correlation
          } else {
            printFlush(paste("standardScreeningBinaryTrait: An error ocurred in t.test for variable",
                             i, ":\n", tst))
            printFlush(paste("Will return missing value(s) for this variable.\n\n"))
          }

        }
        if (getAreaUnderROC) AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = TRUE)[[1]]
        if (kruskalTest) {
            if (no.present<5 )
            {
               pvaluekruskal[i] = stat.Kruskal[i] = NA
            } else {
               kt = try(kruskal.test(datExpr[, i] ~ factor(yNumeric),  na.action="na.exclude"), silent = TRUE)
               if (!inherits(kt, "try-error"))
               {
                 pvaluekruskal[i] = kt$p.value
                 stat.Kruskal[i] = kt$statistic
                 # Find which side is higher
                 r = rank(datExpr[, i])
                 means = tapply(r, factor(yNumeric), mean, na.rm = TRUE)
                 sign.Kruskal[i] = 2 * ( (means[1] < means[2]) - 0.5)
                 # sign.Kruskal is 1 if the ranks in group 1 are smaller than in group 2
               } else {
                 printFlush(paste("standardScreeningBinaryTrait: An error ocurred in kruskal.test for variable",
                                  i, ":\n", kt))
                 printFlush(paste("Will return missing value(s) for this variable.\n\n"))
               }

           }
        }
    }
    q.Student=rep(NA, length(pvalueStudent) )
    rest1= ! is.na(pvalueStudent)
    if (qValues)
    {
      x = try({ q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues }, silent = TRUE)
      if (inherits(x, "try-error"))
        printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                         "calculated q-values will be invalid. qvalue error:\n\n", x, "\n"))
      if (kruskalTest) {
         q.kruskal=rep(NA, length(pvaluekruskal) )
         rest1= ! is.na(pvaluekruskal)
         xx = try( { q.kruskal[rest1] = qvalue(pvaluekruskal[rest1])$qvalues} , silent = TRUE)
         if (inherits(xx, "try-error"))
           printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                            "calculated q-values will be invalid. qvalue error:\n\n", xx, "\n"))
      }
    }
    meanLevel1 = as.numeric(apply(datExpr[!is.na(y) & y == levelsy[[1]], ], 2, mean, na.rm = TRUE))
    meanLevel2 = as.numeric(apply(datExpr[!is.na(y) & y == levelsy[[2]], ], 2, mean, na.rm = TRUE))

    Z.Student = qnorm(pvalueStudent/2, lower.tail = FALSE) * sign(t.Student)
    if (kruskalTest)
      Z.Kruskal = qnorm(pvaluekruskal/2, lower.tail = FALSE) * sign(stat.Kruskal)


    stderr1=function(x) {no.present=sum(!is.na(x))
    if (no.present<2) out=NA else {out=sqrt(var(x,na.rm=TRUE)/no.present) }
    out } # end of function stderr1

    SE.Level1 = as.numeric(apply(datExpr[y == levelsy[[1]] & !is.na(y), ], 2, stderr1))
    SE.Level2 =as.numeric(apply(datExpr[y == levelsy[[2]] & !is.na(y), ], 2, stderr1))

    FoldChangeLevel1vsLevel2 = ifelse(meanLevel1/meanLevel2 > 1, meanLevel1/meanLevel2, -meanLevel2/meanLevel1)

    output = data.frame(ID = dimnames(datExpr)[[2]], corPearson = corPearson,
        t.Student = t.Student,
        pvalueStudent = pvalueStudent,
        FoldChange = FoldChangeLevel1vsLevel2,
         meanFirstGroup = meanLevel1,
        meanSecondGroup = meanLevel2,
         SE.FirstGroup = SE.Level1,
      SE.SecondGroup = SE.Level2)

    if (getAreaUnderROC)
       output$AreaUnderROC = AreaUnderROC

    if (kruskalTest) {
        output = data.frame(output, stat.Kruskal = stat.Kruskal,
                            stat.Kruskal.signed = sign.Kruskal * stat.Kruskal,
                            pvaluekruskal = pvaluekruskal)
    }

   if (qValues && !inherits(x, "try-error")) output=data.frame(output, q.Student)
   if (qValues &  kruskalTest ) {
      if ( !inherits(xx, "try-error")) output=data.frame(output, q.kruskal)
   }
   names(output)[3:5] = paste(names(output)[3:5], levelsy[[1]], "vs", levelsy[[2]], sep = ".")
   output = data.frame(output, nPresentSamples = nPresent)
   output
}

# metaAnalysis ####
#' Meta-analysis of binary and continuous variables
#'
#' This is a meta-analysis complement to functions
#' \code{\link{standardScreeningBinaryTrait}} and
#' \code{\link{standardScreeningNumericTrait}}. Given expression (or other)
#' data from multiple independent data sets, and the corresponding clinical
#' traits or outcomes, the function calculates multiple screening statistics in
#' each data set, then calculates meta-analysis Z scores, p-values, and
#' optionally q-values (False Discovery Rates). Three different ways of
#' calculating the meta-analysis Z scores are provided: the Stouffer method,
#' weighted Stouffer method, and using user-specified weights.
#'
#' The Stouffer method of combines Z statistics by simply taking a mean of
#' input Z statistics and multiplying it by \code{sqrt(n)}, where \code{n} is
#' the number of input data sets. We refer to this method as
#' \code{Stouffer.equalWeights}. In general, a better (i.e., more powerful)
#' method of combining Z statistics is to weigh them by the number of degrees
#' of freedom (which approximately equals \code{n}). We refer to this method as
#' \code{weightedStouffer}. Finally, the user can also specify custom weights,
#' for example if a data set needs to be downweighted due to technical
#' concerns; however, specifying own weights by hand should be done carefully
#' to avoid possible selection biases.
#'
#' @param multiExpr Expression data (or other data) in multi-set format (see
#' \code{\link{checkSets}}). A vector of lists; in each list there must be a
#' component named \code{data} whose content is a matrix or dataframe or array
#' of dimension 2.
#' @param multiTrait Trait or ourcome data in multi-set format. Only one trait
#' is allowed; consequesntly, the \code{data} component of each component list
#' can be either a vector or a data frame (matrix, array of dimension 2).
#' @param binary Logical: is the trait binary (\code{TRUE}) or continuous
#' (\code{FALSE})? If not given, the decision will be made based on the content
#' of \code{multiTrait}.
#' @param metaAnalysisWeights Optional specification of set weights for
#' meta-analysis. If given, must be a vector of non-negative weights, one entry
#' for each set contained in \code{multiExpr}.
#' @param corFnc Correlation function to be used for screening. Should be
#' either the default \code{\link{cor}} or its robust alternative,
#' \code{\link{bicor}}.
#' @param corOptions A named list giving extra arguments to be passed to the
#' correlation function.
#' @param getQvalues Logical: should q-values (FDRs) be calculated?
#' @param getAreaUnderROC Logical: should area under the ROC be calculated?
#' Caution, enabling the calculation will slow the function down considerably
#' for large data sets.
#' @param useRankPvalue Logical: should the \code{\link{rankPvalue}} function
#' be used to obtain alternative meta-analysis statistics?
#' @param rankPvalueOptions Additional options for function
#' \code{\link{rankPvalue}}. These include \code{na.last} (default
#' \code{"keep"}), \code{ties.method} (default \code{"average"}),
#' \code{calculateQvalue} (default copied from input \code{getQvalues}), and
#' \code{pValueMethod} (default \code{"all"}).  See the help file for
#' \code{\link{rankPvalue}} for full details.
#' @param setNames Optional specification of set names (labels). These are used
#' to label the corresponding components of the output. If not given, will be
#' taken from the \code{names} attribute of \code{multiExpr}. If
#' \code{names(multiExpr)} is \code{NULL}, generic names of the form
#' \code{Set_1, Set2, ...} will be used.
#' @param kruskalTest Logical: should the Kruskal test be performed in addition
#' to t-test? Only applies to binary traits.
#' @param var.equal Logical: should the t-test assume equal variance in both
#' groups? If \code{TRUE}, the function will warn the user that the returned
#' test statistics will be different from the results of the standard
#' \code{\link[stats]{t.test}} function.
#' @param metaKruskal Logical: should the meta-analysis be based on the results
#' of Kruskal test (\code{TRUE}) or Student t-test (\code{FALSE})?
#' @param na.action Specification of what should happen to missing values in
#' \code{\link[stats]{t.test}}.
#' @return Data frame with the following components: \item{ID}{ Identifier of
#' the input genes (or other variables) }
#'
#' \item{Z.equalWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with equal weights} \item{p.equalWeights}{ p-values corresponding to
#' \code{Z.Stouffer.equalWeights} } \item{q.equalWeights}{ q-values
#' corresponding to \code{p.Stouffer.equalWeights}, only present if
#' \code{getQvalues} is \code{TRUE}.}
#'
#' \item{Z.RootDoFWeights}{ Meta-analysis Z statistics obtained using
#' Stouffer's method with weights given by the square root of the number of
#' (non-missing) samples in each data set} \item{p.RootDoFWeights}{ p-values
#' corresponding to \code{Z.DoFWeights} } \item{q.RootDoFWeights}{ q-values
#' corresponding to \code{p.DoFWeights}, only present if \code{getQvalues} is
#' \code{TRUE}. }
#'
#' \item{Z.DoFWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with weights given by the number of (non-missing) samples in each
#' data set} \item{p.DoFWeights}{ p-values corresponding to \code{Z.DoFWeights}
#' } \item{q.DoFWeights}{ q-values corresponding to \code{p.DoFWeights}, only
#' present if \code{getQvalues} is \code{TRUE}. }
#'
#' \item{Z.userWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with user-defined weights. Only present if input
#' \code{metaAnalysisWeights} are present.} \item{p.userWeights}{ p-values
#' corresponding to \code{Z.userWeights} } \item{q.userWeights}{ q-values
#' corresponding to \code{p.userWeights}, only present if \code{getQvalues} is
#' \code{TRUE}. }
#'
#' The next set of columns is present only if input \code{useRankPvalue} is
#' \code{TRUE} and contain the output of the function \code{\link{rankPvalue}}
#' with the same column weights as the above meta-analysis. Depending on the
#' input options \code{calculateQvalue} and \code{pValueMethod} in
#' \code{rankPvalueOptions}, some columns may be missing. The following columns
#' are calculated using equal weights for each data set.
#'
#' \item{pValueExtremeRank.equalWeights}{This is the minimum between
#' pValueLowRank and pValueHighRank, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowRank.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the rank method.}
#'
#' \item{pValueHighRank.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the rank method.}
#'
#' \item{pValueExtremeScale.equalWeights}{This is the minimum between
#' pValueLowScale and pValueHighScale, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowScale.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the Scale
#' method.}
#'
#' \item{pValueHighScale.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the Scale
#' method.}
#'
#' \item{qValueExtremeRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueExtremeRank}
#'
#' \item{qValueLowRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueLowRank}
#'
#' \item{qValueHighRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueHighRank}
#'
#' \item{qValueExtremeScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueExtremeScale}
#'
#' \item{qValueLowScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueLowScale}
#'
#' \item{qValueHighScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueHighScale}
#'
#' \item{...}{Analogous columns calculated by weighting each input set using
#' the square root of the number of samples, number of samples, and user
#' weights (if given). The corresponding column names carry the suffixes
#' \code{RootDofWeights}, \code{DoFWeights}, \code{userWeights}.}
#'
#' The following columns contain results returned by
#' \code{\link{standardScreeningBinaryTrait}} or
#' \code{\link{standardScreeningNumericTrait}} (depending on whether the input
#' trait is binary or continuous).
#'
#' For binary traits, the following information is returned for each set:
#'
#' \item{corPearson.Set_1, corPearson.Set_2,...}{Pearson correlation with a
#' binary numeric version of the input variable. The numeric variable equals 1
#' for level 1 and 2 for level 2. The levels are given by levels(factor(y)).}
#'
#' \item{t.Student.Set_1, t.Student.Set_2, ...}{Student t-test statistic}
#'
#' \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{two-sided Student
#' t-test p-value.}
#'
#' \item{qvalueStudent.Set_1, qvalueStudent.Set_2, ...}{(if input
#' \code{qValues==TRUE}) q-value (local false discovery rate) based on the
#' Student T-test p-value (Storey et al 2004).}
#'
#' \item{foldChange.Set_1, foldChange.Set_2, ...}{a (signed) ratio of mean
#' values. If the mean in the first group (corresponding to level 1) is larger
#' than that of the second group, it equals meanFirstGroup/meanSecondGroup.
#' But if the mean of the second group is larger than that of the first group
#' it equals -meanSecondGroup/meanFirstGroup (notice the minus sign).}
#'
#' \item{meanFirstGroup.Set_1, meanSecondGroup.Set_2, ...}{means of columns in
#' input \code{datExpr} across samples in the second group.}
#'
#' \item{SE.FirstGroup.Set_1, SE.FirstGroup.Set_2, ...}{standard errors of
#' columns in input \code{datExpr} across samples in the first group.  Recall
#' that SE(x)=sqrt(var(x)/n) where n is the number of non-missing values of x.
#' }
#'
#' \item{SE.SecondGroup.Set_1, SE.SecondGroup.Set_2, ...}{standard errors of
#' columns in input \code{datExpr} across samples in the second group.}
#'
#' \item{areaUnderROC.Set_1, areaUnderROC.Set_2, ...}{the area under the ROC,
#' also known as the concordance index or C.index.  This is a measure of
#' discriminatory power. The measure lies between 0 and 1 where 0.5 indicates
#' no discriminatory power. 0 indicates that the "opposite" predictor has
#' perfect discriminatory power. To compute it we use the function
#' \link[Hmisc]{rcorr.cens} with \code{outx=TRUE} (from Frank Harrel's package
#' Hmisc).}
#'
#' \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples
#' with finite measurements for each gene.}
#'
#' If input \code{kruskalTest} is \code{TRUE}, the following columns further
#' summarize results of Kruskal-Wallis test:
#'
#' \item{stat.Kruskal.Set_1, stat.Kruskal.Set_2, ...}{Kruskal-Wallis test
#' statistic.}
#'
#' \item{stat.Kruskal.signed.Set_1, stat.Kruskal.signed.Set_2,...}{(Warning:
#' experimental) Kruskal-Wallis test statistic including a sign that indicates
#' whether the average rank is higher in second group (positive) or first group
#' (negative).  }
#'
#' \item{pvaluekruskal.Set_1, pvaluekruskal.Set_2, ...}{Kruskal-Wallis test
#' p-value.}
#'
#' \item{qkruskal.Set_1, qkruskal.Set_2, ...}{q-values corresponding to the
#' Kruskal-Wallis test p-value (if input \code{qValues==TRUE}).}
#'
#' \item{Z.Set1, Z.Set2, ...}{Z statistics obtained from
#' \code{pvalueStudent.Set1, pvalueStudent.Set2, ...} or from
#' \code{pvaluekruskal.Set1, pvaluekruskal.Set2, ...}, depending on input
#' \code{metaKruskal}.}
#'
#' For numeric traits, the following columns are returned:
#'
#' \item{cor.Set_1, cor.Set_2, ...}{correlations of all genes with the trait}
#'
#' \item{Z.Set1, Z.Set2, ...}{Fisher Z statistics corresponding to the
#' correlations}
#'
#' \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{Student p-values of the
#' correlations}
#'
#' \item{qvalueStudent.Set_1, qvalueStudent.Set_1, ...}{(if input
#' \code{qValues==TRUE}) q-values of the correlations calculated from the
#' p-values}
#'
#' \item{AreaUnderROC.Set_1, AreaUnderROC.Set_2, ...}{area under the ROC}
#'
#' \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples
#' present for the calculation of each association. }
#' @author Peter Langfelder
#' @seealso
#' \code{\link{standardScreeningBinaryTrait}},
#' \code{\link{standardScreeningNumericTrait}} for screening functions for
#' individual data sets
#' @references
#' For Stouffer's method, see
#'
#' Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. & Williams, R.M.
#' Jr. 1949. The American Soldier, Vol. 1: Adjustment during Army Life.
#' Princeton University Press, Princeton.
#'
#' A discussion of weighted Stouffer's method can be found in
#'
#' Whitlock, M. C., Combining probability from independent tests: the weighted
#' Z-method is superior to Fisher's approach, Journal of Evolutionary Biology
#' 18:5 1368 (2005)
#' @keywords misc
metaAnalysis <- function(multiExpr,
                         multiTrait,
                         binary = NULL,
                         #consensusQuantile = 0,
                         metaAnalysisWeights = NULL,
                         corFnc = cor,
                         corOptions = list(use = 'p'),
                         getQvalues = FALSE,
                         getAreaUnderROC = FALSE,
                         useRankPvalue = TRUE,
                         rankPvalueOptions = list(),
                         setNames = NULL,
                         kruskalTest = FALSE,
                         var.equal = FALSE,
                         metaKruskal = kruskalTest,
                         na.action = "na.exclude") {
    size = checkSets(multiExpr)
    nSets = size$nSets

    for (set in 1:nSets)
        multiTrait[[set]]$data = as.matrix(multiTrait[[set]]$data)

    tSize = checkSets(multiTrait)
    if (tSize$nGenes != 1)
        stop("This function only works for a single trait.")

    if (size$nSets != tSize$nSets)
        stop("The number of sets in 'multiExpr' and 'multiTrait' ",
             "must be the same.")

    if (!all.equal(size$nSamples, tSize$nSamples))
        stop("Numbers of samples in each set of 'multiExpr' and ",
             "'multiTrait' must be the same.")

    #if (!is.finite(consensusQuantile) || consensusQuantile < 0 ||
    # consensusQuantile > 1)
    #   stop("'consensusQuantile' must be between 0 and 1.")

    if (is.null(setNames))
        setNames = names(multiExpr)

    if (is.null(setNames))
        setNames = paste0("Set_", c(1:nSets))

    if (metaKruskal && !kruskalTest)
        stop("Kruskal statistic meta - analysis requires kruskal test. ",
             "Use kruskalTest = TRUE.")

    if (is.null(binary))
        binary = .isBinary(multiTrait)

    if (!is.null(metaAnalysisWeights))
    {
        if (length(metaAnalysisWeights) != nSets)
            stop(
                "Length of 'metaAnalysisWeights' must equal the number",
                " of sets in 'multiExpr'."
            )
        if (any(!is.finite(metaAnalysisWeights)) || any(metaAnalysisWeights < 0))
            stop("All weights in 'metaAnalysisWeights' must be positive.")
    }

    setResults = list()

    for (set in 1:size$nSets) {
        if (binary) {
            setResults[[set]] = standardScreeningBinaryTrait(
                multiExpr[[set]]$data,
                as.vector(multiTrait[[set]]$data),
                kruskalTest = kruskalTest,
                qValues = getQvalues,
                var.equal = var.equal,
                na.action = na.action,
                corFnc = corFnc,
                corOptions = corOptions
            )
            trafo = TRUE
            if (metaKruskal) {
                metaStat = "stat.Kruskal.signed"
                metaP = "pvaluekruskal"
            } else {
                metaStat = "t.Student"
                metaP = "pvalueStudent"
            }
        } else {
            setResults[[set]] = standardScreeningNumericTrait(
                multiExpr[[set]]$data,
                as.vector(multiTrait[[set]]$data),
                qValues = getQvalues,
                corFnc = corFnc,
                corOptions = corOptions,
                areaUnderROC = getAreaUnderROC
            )
            metaStat = "Z"
            trafo = FALSE
        }
    }

    comb = NULL
    for (set in 1:nSets) {
        if (set == 1) {
            comb = setResults[[set]] [,-1]
            ID = setResults[[set]] [, 1]
            colNames = colnames(comb)
            nColumns = ncol(comb)
            colnames(comb) = paste0("X", c(1:nColumns))
        } else {
            xx = setResults[[set]][,-1]
            colnames(xx) = paste0("X", c(1:nColumns))
            comb = rbind(comb, xx)
        }
    }

    # Re - arrange comb:

    comb = matrix(as.matrix(as.data.frame(comb)), size$nGenes, nColumns * nSets)

    colnames(comb) = paste0(rep(colNames, rep(nSets, nColumns)), ".",
                            rep(setNames, nColumns))

    # Find the columns from which to do meta - analysis
    statCols = grep(paste0("^", metaStat), colnames(comb))
    if (length(statCols) == 0)
        stop("Internal error: no columns for meta - analysis found. Sorry!")
    setStats = comb[, statCols]

    if (trafo) {
        # transform p-values to Z statistics
        # Find the pvalue columns
        pCols = grep(paste0("^", metaP), colnames(comb))
        if (length(pCols) == 0)
            stop("Internal error: no columns for meta - analysis found. Sorry!")
        setP = comb[, pCols]
        # Caution: I assume here that the returned p-values are two - sided.
        setZ = sign(setStats) * qnorm(setP / 2, lower.tail = FALSE)
    } else {
        setZ = setStats
    }

    colnames(setZ) = paste0("Z.", setNames)
    nObsCols = grep("nPresentSamples", colnames(comb))
    nObs = comb[, nObsCols]

    powers = c(0, 0.5, 1)
    nPowers = 3

    metaNames = c("equalWeights", "RootDoFWeights", "DoFWeights")
    if (is.null(metaAnalysisWeights)) {
        nMeta = nPowers
    } else {
        nMeta = nPowers + 1
        metaNames = c(metaNames, "userWeights")
    }
    metaResults = NULL
    for (m in 1:nMeta) {
        if (m <= nPowers) {
            weights = nObs ^ powers[m]
        } else
            weights = matrix(metaAnalysisWeights, size$nGenes, nSets,
                             byrow = TRUE)

        metaZ = rowSums(setZ * weights, na.rm = TRUE) / sqrt(rowSums(weights ^
                                                                         2, na.rm = TRUE))
        p.meta = 2 * pnorm(abs(metaZ), lower.tail = FALSE)
        if (getQvalues) {
            q.meta = qvalue.restricted(p.meta)
            meta1 = cbind(metaZ, p.meta, q.meta)
        } else {
            q.meta = NULL
            meta1 = cbind(metaZ, p.meta)
        }
        colnames(meta1) = paste0(c("Z.", "p.", "q.")[1:ncol(meta1)],
                                 metaNames[m])
        metaResults = cbind(metaResults, meta1)
    }

    # Use rankPvalue to produce yet another meta - analysis

    rankMetaResults = NULL
    if (useRankPvalue) {
        rankPvalueOptions$datS = as.data.frame(setZ)
        if (is.na(match("calculateQvalue", names(rankPvalueOptions))))
            rankPvalueOptions$calculateQvalue = getQvalues
        for (m in 1:nMeta) {
            if (m <= nPowers) {
                weights = nObs ^ powers[m]
            } else
                weights = matrix(metaAnalysisWeights, size$nGenes, nSets,
                                 byrow = TRUE)

            # rankPvalue requires a vector of weights... so compress the weights to a vector.
            # Output a warning if the compression loses information.
            nDifferent = apply(weights, 2, function(x) {
                length(unique(x))
            })
            if (any(nDifferent) > 1)
                printFlush(
                    paste(
                        "Warning in metaAnalysis: rankPvalue requires compressed ",
                        "weights.\nSome weights may not be entirely accurate."
                    )
                )
            rankPvalueOptions$columnweights = colMeans(weights, na.rm = TRUE)
            rankPvalueOptions$columnweights = rankPvalueOptions$columnweights /
                sum(rankPvalueOptions$columnweights)
            rp = do.call(rankPvalue, rankPvalueOptions)
            colnames(rp) = paste0(colnames(rp), ".", metaNames[m])
            rankMetaResults = cbind(rankMetaResults, as.matrix(rp))
        }
    }

    # Put together the output

    out = list(ID = ID,
               metaResults,
               rankMetaResults,
               comb,
               if (trafo)
                   setZ
               else
                   NULL,
               NULL)
    # The last NULL is necessary so the line below works even if nothing else is
    #  NULL

    out = as.data.frame(out[-(which(sapply(out, is.null), arr.ind = TRUE))])

    out
}
