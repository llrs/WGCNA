#' Estimate the population-specific mean values in an admixed population.
#' 
#' Uses the expression values from an admixed population and estimates of the
#' proportions of sub-populations to estimate the population specific mean
#' values. For example, this function can be used to estimate the cell type
#' specific mean gene expression values based on expression values from a
#' mixture of cells. The method is described in Shen-Orr et al (2010) where it
#' was used to estimate cell type specific gene expression levels based on a
#' mixture sample.
#' 
#' The function outputs a matrix of coefficients resulting from fitting a
#' regression model. If the proportions sum to 1, then i-th row of the output
#' matrix reports the coefficients of the following model
#' \code{lm(datE.Admixture[,i]~.-1,data=datProportions)}. Aside, the minus 1 in
#' the formula indicates that no intercept term will be fit. Under certain
#' assumptions, the coefficients can be interpreted as the mean expression
#' values in the sub-populations (Shen-Orr 2010).
#' 
#' @param datProportions a matrix of non-negative numbers (ideally proportions)
#' where the rows correspond to the samples (rows of \code{datE.Admixture}) and
#' the columns correspond to the sub-populations of the mixture. The function
#' calculates a mean expression value for each column of \code{datProportions}.
#' Negative entries in datProportions lead to an error message. But the rows of
#' \code{datProportions} do not have to sum to 1, see the argument
#' \code{scaleProportionsTo1}.
#' @param datE.Admixture a matrix of numbers. The rows correspond to samples
#' (mixtures of populations). The columns contain the variables (e.g. genes)
#' for which the means should be estimated.
#' @param scaleProportionsTo1 logical. If set to TRUE (default) then the
#' proportions in each row of \code{datProportions} are scaled so that they sum
#' to 1, i.e. datProportions[i,]=datProportions[i,]/max(datProportions[i,]). In
#' general, we recommend to set it to TRUE.
#' @param scaleProportionsInCelltype logical. If set to TRUE (default) then the
#' proportions in each cell types are recaled and make the mean to 0.
#' @param setMissingProportionsToZero logical. Default is FALSE. If set to TRUE
#' then it sets missing values in \code{datProportions} to zero.
#' @return a numeric matrix whose rows correspond to the columns of
#' \code{datE.Admixture} (e.g. to genes) and whose columns correspond to the
#' columns of \code{datProportions} (e.g. sub populations or cell types).
#' @note This can be considered a wrapper of the \code{lm} function.
#' @author Steve Horvath, Chaochao Cai
#' @references Shen-Orr SS, Tibshirani R, Khatri P, Bodian DL, Staedtler F,
#' Perry NM, Hastie T, Sarwal MM, Davis MM, Butte AJ (2010) Cell type-specific
#' gene expression differences in complex tissues.  Nature Methods, vol 7 no.4
#' @keywords misc
#' @examples
#' 
#' 
#' set.seed(1)
#' # this is the number of complex (mixed) tissue samples, e.g. arrays
#' m=10
#' # true count data (e.g. pure cells in the mixed sample)
#' datTrueCounts=as.matrix(data.frame(TrueCount1=rpois(m,lambda=16),
#' TrueCount2=rpois(m,lambda=8),TrueCount3=rpois(m,lambda=4),
#' TrueCount4=rpois(m,lambda=2)))
#' no.pure=dim(datTrueCounts)[[2]]
#' 
#' # now we transform the counts into proportions
#' divideBySum=function(x) t(x)/sum(x)
#' datProportions= t(apply(datTrueCounts,1,divideBySum))
#' dimnames(datProportions)[[2]]=paste("TrueProp",1:dim(datTrueCounts)[[2]],sep=".")
#' 
#' # number of genes that are highly expressed in each pure population
#' no.genesPerPure=rep(5, no.pure)
#' no.genes= sum(no.genesPerPure)
#' GeneIndicator=rep(1:no.pure, no.genesPerPure)
#' # true mean values of the genes in the pure populations
#' # in the end we hope to estimate them from the mixed samples
#' datTrueMeans0=matrix( rnorm(no.genes*no.pure,sd=.3), nrow= no.genes,ncol=no.pure)
#' for (i in 1:no.pure ){
#' datTrueMeans0[GeneIndicator==i,i]= datTrueMeans0[GeneIndicator==i,i]+1
#' }
#' dimnames(datTrueMeans0)[[1]]=paste("Gene",1:dim(datTrueMeans0)[[1]],sep="." )
#' dimnames(datTrueMeans0)[[2]]=paste("MeanPureCellType",1:dim(datTrueMeans0)[[2]],
#'                                    sep=".")
#' # plot.mat(datTrueMeans0)
#' # simulate the (expression) values of the admixed population samples
#' 
#' noise=matrix(rnorm(m*no.genes,sd=.1),nrow=m,ncol= no.genes)
#' datE.Admixture= as.matrix(datProportions) %*% t(datTrueMeans0) + noise
#' dimnames(datE.Admixture)[[1]]=paste("MixedTissue",1:m,sep=".")
#' 
#' datPredictedMeans=populationMeansInAdmixture(datProportions,datE.Admixture)
#' 
#' par(mfrow=c(2,2))
#' for (i in 1:4 ){
#' verboseScatterplot(datPredictedMeans[,i],datTrueMeans0[,i],
#' xlab="predicted mean",ylab="true mean",main="all populations")
#' abline(0,1)
#' }
#' 
#' #assume we only study 2 populations (ie we ignore the others)
#' selectPopulations=c(1,2)
#' datPredictedMeansTooFew=populationMeansInAdmixture(datProportions[,selectPopulations],
#'                                                    datE.Admixture)
#' 
#' par(mfrow=c(2,2))
#' for (i in 1:length(selectPopulations) ){
#' verboseScatterplot(datPredictedMeansTooFew[,i],datTrueMeans0[,i],
#' xlab="predicted mean",ylab="true mean",main="too few populations")
#' abline(0,1)
#' }
#' 
#' #assume we erroneously add a population
#' datProportionsTooMany=data.frame(datProportions,WrongProp=sample(datProportions[,1]))
#' datPredictedMeansTooMany=populationMeansInAdmixture(datProportionsTooMany,
#'                                  datE.Admixture)
#' 
#' par(mfrow=c(2,2))
#' for (i in 1:4 ){
#'   verboseScatterplot(datPredictedMeansTooMany[,i],datTrueMeans0[,i],
#'   xlab="predicted mean",ylab="true mean",main="too many populations")
#'   abline(0,1)
#' }
#' 
#' 
#' 
populationMeansInAdmixture<-function (
datProportions, 
datE.Admixture, 
scaleProportionsTo1 = TRUE, 
scaleProportionsInCelltype=TRUE,
setMissingProportionsToZero = FALSE
) 
{
    datProportions = as.matrix(datProportions)
    if (dim(datE.Admixture)[[1]] != dim(as.matrix(datProportions))[[1]]) 
        stop("Input error. The numbers of samples are not congruent: dim(datE.Admixture)[[1]] is unequal to dim(datProportions)[[1]]. Hint: Consider transposing one of the matrices.")
    noMissing = apply(is.na(datProportions), 1, sum)
    if (max(noMissing) > 0) {
        warning(paste("Urgent Warning: datProportions contains missing proportions in the following row(s):", 
            paste(which(noMissing > 0), collapse = ","), "\nCheck these rows in datProportions. But for your convenience, I will proceed"))
    }
    if (setMissingProportionsToZero) {
        datProportions[is.na(datProportions)] = 0
    }
    noNegative = apply((datProportions) < 0, 1, sum, na.rm = T)
    if (max(noNegative) > 0) {
        stop(paste("datProportions contains negative numbers. Negative proportions can be found in the following row(s):", 
            paste(which(noNegative > 0), collapse = ","), "\nCheck these rows in datProportions."))
    }
    sumsTo1 = TRUE
    sumProp = apply(datProportions, 1, sum, na.rm = TRUE)
    if (max(sumProp, na.rm = T) > 1.0000001) {
        sumsTo1 = FALSE
        if (scaleProportionsTo1) {
            warning(paste("The sum of proportions is larger than 1 for some rows including:", 
                paste(which(sumProp > 1.0000001)[1:5], collapse = ","), 
                "\nCheck these rows in datProportions. By default, I will scale them so that they sum to 1.\nBut if you do not want this scaling, please set scaleProportionsTo1=FALSE .\n"))
        }
    }
    if (min(sumProp, na.rm = T) < 0.99999) {
        sumsTo1 = FALSE
        if (scaleProportionsTo1) {
            warning(paste("The sum of proportions is smaller than 1 for some rows including:", 
                paste(which(sumProp < 0.99999)[1:5], collapse = ","), 
                "\n Check these rows in datProportions. By default, I will scale them so that they sum to 1.\nBut if you do not want this scaling, please set scaleProportionsTo1=FALSE .\n"))
        }
    }
    if (scaleProportionsTo1) {
        sumsTo1 = TRUE
        for (i in 1:dim(datProportions)[1]) {
            datProportions[i, ] = datProportions[i, ]/sum(datProportions[i, 
                ], na.rm = T)
        }
    }
    if (sumsTo1) {
	if(scaleProportionsInCelltype) {
		for(ci in dim(datProportions)[2]) 
		datProportions[,ci]=datProportions[,ci]-mean(datProportions[,ci])	
	}

        fit1 = lm(as.matrix(datE.Admixture) ~ . - 1, data = data.frame(datProportions))
        datPredictedMeans = t(as.matrix(fit1$coefficients))
    }
    if (!sumsTo1) {
	if(scaleProportionsInCelltype) {
		for(ci in dim(datProportions)[2]) 
		datProportions[,ci]=datProportions[,ci]-mean(datProportions[,ci])	
	}
        fit1 = lm(as.matrix(datE.Admixture) ~ ., data = data.frame(datProportions))
        if (dim(as.matrix(datProportions))[[2]] == 1) {
            datPredictedMeans = (matrix(fit1$coefficients[-1, 
                ], ncol = 1))
        }
        else {
            datPredictedMeans = t(as.matrix(fit1$coefficients[-1, 
                ]))
        }
    }
    dimnames(datPredictedMeans)[[1]] = dimnames(datE.Admixture)[[2]]
    if (is.null(dimnames(datPredictedMeans)[[2]])) {
        dimnames(datPredictedMeans)[[2]] = paste("Mean", 1:dim(datPredictedMeans)[[2]], 
            sep = ".")
    }
    else {
        dimnames(datPredictedMeans)[[2]] = paste("Mean", dimnames(datPredictedMeans)[[2]], 
            sep = ".")
    }
    datPredictedMeans
}
