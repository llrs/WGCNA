# The function standardScreeningBinaryTrait computes widely used statistics for relating the columns of
# the input data frame (argument datExpr) to a binary sample trait (argument y). The statistics include
# Student t-test p-value and the corresponding local false discovery rate (known as q-value, Storey et al
# 2004), the fold change, the area under the ROC curve (also known as C-index), mean values etc. If the
# input option kruskalTest is set to TRUE, it also computes the kruskal Wallist test p-value and
# corresponding q-value. The kruskal Wallis test is a non-parametric, rank-based group comparison test.



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
#' m=50
#' y=sample(c(1,2),m,replace=TRUE)
#' datExprSignal=simulateModule(scale(y),30)
#' datExprNoise=simulateModule(rnorm(m),150)
#' datExpr=data.frame(datExprSignal,datExprNoise)
#' 
#' Result1=standardScreeningBinaryTrait(datExpr,y)
#' Result1[1:5,]
#' 
#' 
#' 
#' # use unequal variances and calculate q-values
#' Result2=standardScreeningBinaryTrait(datExpr,y, var.equal=FALSE,qValue=TRUE)
#' Result2[1:5,]
#' 
#' # calculate Kruskal Wallis test and q-values
#' Result3=standardScreeningBinaryTrait(datExpr,y,kruskalTest=TRUE,qValue=TRUE)
#' Result3[1:5,]
#' 
#' 
standardScreeningBinaryTrait=function(datExpr, y, 
           corFnc = cor, corOptions = list(use = 'p'),
           kruskalTest=FALSE, qValues = FALSE, var.equal=FALSE, na.action="na.exclude",
           getAreaUnderROC = TRUE) 
{

   datExpr=data.frame(datExpr)
   levelsy=levels(factor(y))
   if (length(levelsy)>2 ) 
     stop("The sample trait y contains more than 2 levels. Please input a binary variable y")
   if (length(levelsy)==1 ) 
     stop("The sample trait y is constant. Please input a binary sample trait with some variation.")
   yNumeric=as.numeric(factor(y));
   if (length(yNumeric) !=dim(datExpr)[[1]] ) 
     stop("the length of the sample trait y does not equal the number of rows of datExpr")
   pvalueStudent = t.Student = Z.Student = rep(NA, dim(datExpr)[[2]] ) 
   pvaluekruskal = stat.Kruskal = Z.Kruskal = sign.Kruskal = rep(NA, dim(datExpr)[[2]] ) 
   nPresent = rep(0, dim(datExpr)[[2]] ) 
   AreaUnderROC=rep(NA, dim(datExpr)[[2]] ) 

  if (var.equal) 
     printFlush(paste("Warning: T-test that assumes equal variances in each group is requested.\n",
                      "This is not the default option for t.test. We recommend to use var.equal=FALSE."));

  corFnc = match.fun(corFnc);
  corOptions$y = yNumeric;
  corOptions$x = datExpr;
  corPearson=as.numeric(do.call(corFnc, corOptions));
  nGenes = dim(datExpr)[[2]];
  nPresent1 = as.numeric( t(as.matrix(!is.na(yNumeric) & yNumeric==1)) %*% ! is.na(datExpr) );
  nPresent2 = as.numeric( t(as.matrix(!is.na(yNumeric) & yNumeric==2)) %*% ! is.na(datExpr) );
  nPresent = nPresent1 + nPresent2;

  for (i in 1:nGenes) {
        no.present1 = nPresent1[i];
        no.present2 = nPresent2[i];
        no.present = nPresent[i];
        if (no.present1<2 | no.present2<2 ) 
        {
           pvalueStudent[i]= t.Student[i] = NA 
        } else {
          tst = try(t.test( as.numeric(datExpr[,i])~yNumeric,var.equal=var.equal,na.action=na.action),
                    silent = TRUE)
          if (!inherits(tst, "try-error"))
          {
            pvalueStudent[i] = tst$p.value;
            t.Student[i] = -tst$statistic 
            # The - sign above is intentional to make the sign of t consistent with correlation
          } else {
            printFlush(paste("standardScreeningBinaryTrait: An error ocurred in t.test for variable", 
                             i, ":\n", tst));
            printFlush(paste("Will return missing value(s) for this variable.\n\n"));
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
                 pvaluekruskal[i] = kt$p.value;
                 stat.Kruskal[i] = kt$statistic;
                 # Find which side is higher
                 r = rank(datExpr[, i]);
                 means = tapply(r, factor(yNumeric), mean, na.rm = TRUE);
                 sign.Kruskal[i] = 2 * ( (means[1] < means[2]) - 0.5);
                 # sign.Kruskal is 1 if the ranks in group 1 are smaller than in group 2
               } else {
                 printFlush(paste("standardScreeningBinaryTrait: An error ocurred in kruskal.test for variable",
                                  i, ":\n", kt));
                 printFlush(paste("Will return missing value(s) for this variable.\n\n"));
               }

           }
        } 
    }
    q.Student=rep(NA, length(pvalueStudent) )
    rest1= ! is.na(pvalueStudent) 
    if (qValues) 
    {
      x = try({ q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues }, silent = TRUE);
      if (inherits(x, "try-error"))
        printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                         "calculated q-values will be invalid. qvalue error:\n\n", x, "\n")) 
      if (kruskalTest) {
         q.kruskal=rep(NA, length(pvaluekruskal) )
         rest1= ! is.na(pvaluekruskal) 
         xx = try( { q.kruskal[rest1] = qvalue(pvaluekruskal[rest1])$qvalues} , silent = TRUE);
         if (inherits(xx, "try-error"))
           printFlush(paste("Warning in standardScreeningBinaryTrait: function qvalue returned an error.\n",
                            "calculated q-values will be invalid. qvalue error:\n\n", xx, "\n")) 
      }
    }
    meanLevel1 = as.numeric(apply(datExpr[!is.na(y) & y == levelsy[[1]], ], 2, mean, na.rm = TRUE));
    meanLevel2 = as.numeric(apply(datExpr[!is.na(y) & y == levelsy[[2]], ], 2, mean, na.rm = TRUE));
  
    Z.Student = qnorm(pvalueStudent/2, lower.tail = FALSE) * sign(t.Student);
    if (kruskalTest)
      Z.Kruskal = qnorm(pvaluekruskal/2, lower.tail = FALSE) * sign(stat.Kruskal);

    
    stderr1=function(x) {no.present=sum(!is.na(x));
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
      SE.SecondGroup = SE.Level2);

    if (getAreaUnderROC)
       output$AreaUnderROC = AreaUnderROC;

    if (kruskalTest) {
        output = data.frame(output, stat.Kruskal = stat.Kruskal, 
                            stat.Kruskal.signed = sign.Kruskal * stat.Kruskal,
                            pvaluekruskal = pvaluekruskal);
    }

   if (qValues && !inherits(x, "try-error")) output=data.frame(output, q.Student)
   if (qValues &  kruskalTest ) {
      if ( !inherits(xx, "try-error")) output=data.frame(output, q.kruskal)
   }
   names(output)[3:5] = paste(names(output)[3:5], levelsy[[1]], "vs", levelsy[[2]], sep = ".")
   output = data.frame(output, nPresentSamples = nPresent);
   output
}




