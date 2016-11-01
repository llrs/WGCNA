# TrueTrait ####
#' Estimate the true trait underlying a list of surrogate markers.
#'
#' Assume an imprecisely measured trait \code{y} that is related to the true,
#' unobserved trait yTRUE as follows yTRUE=y+noise where noise is assumed to
#' have mean zero and a constant variance. Assume you have 1 or more surrogate
#' markers for yTRUE corresponding to the columns of \code{datX}. The function
#' implements several approaches for estimating yTRUE based on the inputs
#' \code{y} and/or \code{datX}.
#'
#' This R function implements formulas described in Klemera and Doubal (2006).
#' The assumptions underlying these formulas are described in Klemera et al.
#' But briefly, the function provides several estimates of the true underlying
#' trait under the following assumptions: 1) There is a true underlying trait
#' that affects \code{y} and a list of surrogate markers corresponding to the
#' columns of \code{datX}. 2) There is a linear relationship between the true
#' underlying trait and \code{y} and the surrogate markers. 3) yTRUE =y +Noise
#' where the Noise term has a mean of zero and a fixed variance. 4) Weighted
#' least squares estimation is used to relate the surrogate markers to the
#' underlying trait where the weights are proportional to 1/ssq.j where ssq.j
#' is the noise variance of the j-th marker.
#'
#' Specifically, output \code{y.true1} corresponds to formula 31,
#' \code{y.true2} corresponds to formula 25, and \code{y.true3} corresponds to
#' formula 34.
#'
#' Although the true underlying trait yTRUE is not known, one can estimate the
#' standard deviation between the estimate \code{y.true2} and yTRUE using
#' formula 33. Similarly, one can estimate the SD for the estimate
#' \code{y.true3} using formula 42. These estimated SDs correspond to output
#' components 2 and 3, respectively. These SDs are valuable since they provide
#' a sense of how accurate the measure is.
#'
#' To estimate the correlations between \code{y} and the surrogate markers, one
#' can specify different correlation measures. The default method is based on
#' the Person correlation but one can also specify the biweight midcorrelation
#' by choosing "bicor", see help(bicor) to learn more.
#'
#' When the \code{datX} is comprised of observations measured in different
#' strata (e.g. different batches or independent data sets) then one can obtain
#' stratum specific estimates by specifying the strata using the argument
#' \code{Strata}. In this case, the estimation focuses on one stratum at a
#' time.
#'
#' @param datX is a vector or data frame whose columns correspond to the
#' surrogate markers (variables) for the true underlying trait. The number of
#' rows of \code{datX} equals the number of observations, i.e. it should equal
#' the length of \code{y}
#' @param y is a numeric vector which specifies the observed trait.
#' @param datXtest can be set as a matrix or data frame of a second,
#' independent test data set. Its columns should correspond to those of
#' \code{datX}, i.e. the two data sets should have the same number of columns
#' but the number or rows (test set observations) can be different.
#' @param corFnc Character string specifying the correlation function to be
#' used in the calculations.  Recomended values are the default Pearson
#' correlation \code{"cor"} or biweight mid-correlation \code{"bicor"}.
#' Additional arguments to the correlation function can be specified using
#' \code{corOptions}.
#' @param corOptions Character string giving additional arguments to the
#' function specified in \code{corFnc}.
#' @param LeaveOneOut.CV logical. If TRUE then leave one out cross validation
#' estimates will be calculated for \code{y.true1} and \code{y.true2} based on
#' \code{datX}.
#' @param skipMissingVariables logical. If TRUE then variables whose values are
#' missing for a given observation will be skipped when estimating the true
#' trait of that particular observation. Thus, the estimate of a particular
#' observation are determined by all the variables whose values are
#' non-missing.
#' @param addLinearModel logical. If TRUE then the function also estimates the
#' true trait based on the predictions of the linear model \code{lm(y~.,
#' data=datX)}
#' @return A list with the following components. \item{datEstimates}{is a data
#' frame whose columns corresponds to estimates of the true underlying trait.
#' The number of rows equals the number of observations, i.e. the length of
#' \code{y}. The first column \code{y.true1} is the average value of
#' standardized columns of \code{datX} where standardization subtracts out the
#' intercept term and divides by the slope of the linear regression model
#' lm(marker~y). Since this estimate ignores the fact that the surrogate
#' markers have different correlations with \code{y}, it is typically inferior
#' to \code{y.true2}.  The second column \code{y.true2} equals the weighted
#' average value of standardized columns of \code{datX}. The standardization is
#' described in section 2.4 of Klemera et al. The weights are proportional to
#' r^2/(1+r^2) where r denotes the correlation between the surrogate marker and
#' \code{y}. Since this estimate does not include \code{y} as additional
#' surrogate marker, it may be slightly inferior to \code{y.true3}. Having said
#' this, the difference between \code{y.true2} and \code{y.true3} is often
#' negligible.  An additional column called \code{y.lm} is added if
#' codeaddLinearModel=TRUE. In this case, \code{y.lm} reports the linear model
#' predictions. Finally, the column \code{y.true3} is very similar to
#' \code{y.true2} but it includes \code{y} as additional surrogate marker. It
#' is expected to be the best estimate of the underlying true trait (see
#' Klemera et al 2006). }
#'
#' \item{datEstimatestest}{is output only if a test data set has been specified
#' in the argument \code{datXtest}. In this case, it contains a data frame with
#' columns \code{ytrue1} and \code{ytrue2}. The number of rows equals the
#' number of test set observations, i.e the number of rows of \code{datXtest}.
#' Since the value of \code{y} is not known in case of a test data set, one
#' cannot calculate \code{y.true3}. An additional column with linear model
#' predictions \code{y.lm} is added if codeaddLinearModel=TRUE.  }
#'
#' \item{datEstimates.LeaveOneOut.CV}{is output only if the argument
#' \code{LeaveOneOut.CV} has been set to \code{TRUE}. In this case, it contains
#' a data frame with leave-one-out cross validation estimates of \code{ytrue1}
#' and \code{ytrue2}. The number of rows equals the length of \code{y}. Since
#' the value of \code{y} is not known in case of a test data set, one cannot
#' calculate \code{y.true3} }
#'
#' \item{SD.ytrue2}{is a scalar. This is an estimate of the standard deviation
#' between the estimate \code{y.true2} and the true (unobserved) yTRUE. It
#' corresponds to formula 33.}
#'
#' \item{SD.ytrue3}{is a scalar. This is an estimate of the standard deviation
#' between \code{y.true3} and the true (unobserved) yTRUE. It corresponds to
#' formula 42.}
#'
#' \item{datVariableInfo}{is a data frame that reports information for each
#' variable (column of \code{datX}) when it comes to the definition of
#' \code{y.true2}. The rows correspond to the number of variables. Columns
#' report the variable name, the center (intercept that is subtracted to scale
#' each variable), the scale (i.e. the slope that is used in the denominator),
#' and finally the weights used in the weighted sum of the scaled variables.}
#'
#' \item{datEstimatesByStratum}{ a data frame that will only be output if
#' \code{Strata} is different from NULL. In this case, it is has the same
#' dimensions as \code{datEstimates} but the estimates were calculated
#' separately for each level of \code{Strata}.}
#'
#' \item{SD.ytrue2ByStratum}{ a vector of length equal to the different levels
#' of \code{Strata}. Each component reports the estimate of \code{SD.ytrue2}
#' for observations in the stratum specified by unique(Strata).}
#'
#' \item{datVariableInfoByStratum}{ a list whose components are matrices with
#' variable information. Each list component reports the variable information
#' in the stratum specified by unique(Strata). }
#' @author Steve Horvath
#' @references Klemera P, Doubal S (2006) A new approach to the concept and
#' computation of biological age. Mechanisms of Ageing and Development 127
#' (2006) 240-248
#'
#' Choa IH, Parka KS, Limb CJ (2010) An Empirical Comparative Study on
#' Validation of Biological Age Estimation Algorithms with an Application of
#' Work Ability Index. Mechanisms of Ageing and Development Volume 131, Issue
#' 2, February 2010, Pages 69-78
#' @keywords misc
#' @examples
#'
#' # observed trait
#' y=rnorm(1000,mean=50,sd=20)
#' # unobserved, true trait
#' yTRUE =y +rnorm(100,sd=10)
#' # now we simulate surrogate markers around the true trait
#' datX=simulateModule(yTRUE,nGenes=20, minCor=.4,maxCor=.9,geneMeans=rnorm(20,50,30)  )
#' True1=TrueTrait(datX=datX,y=y)
#' datTrue=True1$datEstimates
#' par(mfrow=c(2,2))
#' for (i in 1:dim(datTrue)[[2]] ){
#'   meanAbsDev= mean(abs(yTRUE-datTrue[,i]))
#'   verboseScatterplot(datTrue[,i],yTRUE,xlab=names(datTrue)[i],
#'                      main=paste(i, "MeanAbsDev=", signif(meanAbsDev,3)))
#'   abline(0,1)
#' }
#' #compare the estimated standard deviation of y.true2
#' True1[[2]]
#' # with the true SD
#' sqrt(var(yTRUE-datTrue$y.true2))
#' #compare the estimated standard deviation of y.true3
#' True1[[3]]
#' # with the true SD
#' sqrt(var(yTRUE-datTrue$y.true3))
#'
TrueTrait <- function(datX, y,datXtest=NULL, corFnc = "bicor",
                      corOptions = "use = 'pairwise.complete.obs'",
                      LeaveOneOut.CV=FALSE,skipMissingVariables=TRUE,
                      addLinearModel=FALSE){
    datX=as.matrix(datX)
    no.variables=dim(as.matrix(datX))[[2]]
    datVariableInfo=data.frame(matrix(NA, nrow=no.variables, ncol=4))
    names(datVariableInfo)=c("Variable","center", "scale", "weights.y.true2")
    if ( is.null(colnames(datX))){datVariableInfo$Variable=1:no.variables} else { datVariableInfo$Variable=colnames(datX)}
    no.observations=dim(as.matrix(datX))[[1]]
    if (no.observations !=  length(y) ) {
        stop("The number of rows of datX does not correspond to the length of ",
             "y. Consider transposing your input matrix or use a different ",
             "input.")
        }
    if (no.observations==1 ) {
        warning("Only 1 observations, i.e. the length of y is 1. The function",
                " cannot be used. For your convenience, the estimated true",
                " values will be set to the input value of y.")
        y.true1=y
        y.true2=y
        y.true3=y
        }
    if (no.observations>1 ) {
        y.true1=rep(NA, length(y) )
        y.true2=rep(NA, length(y) )
        y.true3=rep(NA, length(y) )
        y.lm=rep(NA, length(y) )
        restNonMissingY= !is.na(y)
        r.characteristic=NA
        SD.ytrue2=NA
        SD.ytrue3=NA
        SsquaredBE=NA

        if (sum(restNonMissingY,na.rm=TRUE) >3 ){
            corX = parse(text = paste(corFnc, "(datX,y ",prepComma(corOptions), ")"))
            rVector= as.numeric(eval(corX))
            datCoef=t(coef(lm(datX~ y,na.action="na.exclude")))
            datVariableInfo$center=datCoef[,1]	#intercept
            datVariableInfo$scale=datCoef[,2]	#slope
            datXscaled=scale(datX,center=datCoef[,1],scale=datCoef[,2] )
            weights0=rVector^2/((1-rVector^2)*var(y,na.rm=TRUE)) # Steve, this is where I made the one change
            weights=weights0/sum(weights0)
            datVariableInfo$weights.y.true2=weights
            y.true1=as.numeric(apply(as.matrix(datXscaled),1,mean))
            y.true2=as.numeric(as.matrix(datXscaled)%*%weights)

            if (skipMissingVariables ) {
                y.true1= as.numeric(apply(as.matrix(datXscaled),1,mean,na.rm=TRUE))
                weightsMatrix=matrix(weights,byrow=TRUE,nrow=dim(as.matrix(datXscaled))[[1]],ncol=length(weights) )
                weightsMatrix[is.na(datXscaled)]=0
                rowsum.weightsMatrix=apply(as.matrix(weightsMatrix),1,sum)
                weightsMatrix=t(scale(t(as.matrix(weightsMatrix)),center=F,scale= rowsum.weightsMatrix))
                datXscaledweighted= as.matrix(datXscaled* weightsMatrix)
                # this corresponds to formula 25 in Klemera et al 2006
                y.true2=as.numeric(apply(datXscaledweighted,1,sum,na.rm=TRUE))
            } #end of if (skipMissingVariables )


            # the following is different from Klemera in that it has an absolute value
            r.characteristic=sum(rVector^2/sqrt(1-rVector^2) )/sum(abs(rVector)/sqrt(1-rVector^2) )
            no.missing=sum(apply(  as.matrix( is.na(datX)),1,sum))
            if (sum(no.missing)>0) {
                warning("The input datX contains missing values.\nI recommend ",
                        "you impute missing values in datX before running this",
                        " function.")
            } # end of if (sum(no.missing)>0)
            # formula 37 from Klemera
            SsquaredBE=var( y.true2-y,na.rm=TRUE) -(1- r.characteristic^2)/r.characteristic^2*var(y,na.rm=TRUE)/no.variables
            # this corresponds to formula 34 in Klemera
            y.true3=(as.numeric( as.matrix(datXscaled)%*% weights0)+y/SsquaredBE )/( sum(weights0)+ 1/SsquaredBE)
            y.true3[is.na(y.true3) ]=y[is.na(y.true3)]
        } # end of if (no.observations>1 )
        SD.ytrue2=sqrt(1-r.characteristic^2)/r.characteristic*sqrt(var(y,na.rm=TRUE)/no.variables)
        # now formula 42
        SD.ytrue3=SD.ytrue2/sqrt(1+SD.ytrue2^2/SsquaredBE )
    } # end of if (sum(restNonMissingY,na.rm=TRUE) >3 )
    datEstimates=data.frame(y, y.true1,y.true2,y.true3)

    if (!is.null(datXtest)){
        datXtest=as.matrix(datXtest)
        no.variablestest=dim(as.matrix(datXtest))[[2]]
        if (no.variablestest != no.variables) {
            stop("the number of variables in the test data is not the same as in the training data")}

        y.true1test=rep(NA, length(y) )
        y.true2test=rep(NA, length(y) )
        y.true3test=rep(NA, length(y) )
        restNonMissingY= !is.na(y)
        if (sum(restNonMissingY,na.rm=TRUE) >3 ){
            datXtestscaled=scale(datXtest,center=datCoef[,1],scale=datCoef[,2] )
            y.true1test=as.numeric(apply(as.matrix(datXtestscaled),1,mean) )
            y.true2test=as.numeric( as.matrix(datXtestscaled)%*%weights)
            if (skipMissingVariables ) {
                y.true1test= as.numeric(apply(as.matrix(datXtestscaled),1,mean,na.rm=TRUE))
                weightsMatrixtest=matrix(weights,byrow=TRUE,nrow=dim(as.matrix(datXtestscaled))[[1]],ncol=length(weights) )
                weightsMatrixtest[is.na(datXtestscaled)]=0
                rowsum.weightsMatrixtest=apply(as.matrix(weightsMatrixtest),1,sum)
                weightsMatrixtest=t(scale(t(as.matrix(weightsMatrixtest)),center=F,scale= rowsum.weightsMatrixtest))
                datXscaledweightedtest= as.matrix(datXtestscaled* weightsMatrixtest)
                # this corresponds to formula 25 in Klemera et al 2006
                y.true2test=as.numeric(apply(datXscaledweightedtest,1,sum,na.rm=TRUE))
            } #end of if (skipMissingVariables )
        } # end of if (sum(restNonMissingY,na.rm=TRUE) >3 )


        datEstimatestest=data.frame(y.true1= y.true1test,y.true2= y.true2test)

    } # end of if (!is.null(datXtest))

    if ( LeaveOneOut.CV  ) {
        y.true1test.LOO=rep(NA,no.observations)
        y.true2test.LOO=rep(NA,no.observations)
        y.lmLOO= rep(NA,no.observations)
        for ( i in 1:no.observations ){
            rm(datCoef)
            rm(corX)

            datX.LOO=datX[-i,]
            datXtest.LOO= matrix(datX[i,],nrow=1)
            y.LOO=y[-i]
            no.variables=dim(as.matrix(datX.LOO))[[2]]
            no.observations=dim(as.matrix(datX.LOO))[[1]]

            if (no.observations==1 ) {
                warning("When dealing with leave one out cross validation, there is only 1 observations in the training data")}

            if (no.observations>1 ) {

                if (addLinearModel) {
                    lmLOO=lm(y.LOO~., data=data.frame(datX.LOO),na.action=na.exclude)
                    y.lmLOO[i]= sum(datXtest.LOO*lmLOO$coeff[-1])+lmLOO$coeff[[1]]
                }

                corX = parse(text = paste(corFnc, "(datX.LOO,y.LOO ",
                                          prepComma(corOptions), ")"))
                rVector= as.numeric(eval(corX))
                datCoef=t(coef(lm(datX.LOO~ y.LOO,na.action="na.exclude")))
                datX.LOOscaled=scale(datX.LOO,center=datCoef[,1],
                                     scale=datCoef[,2] )
                weights0=rVector^2/(1-rVector^2)
                weights=weights0/sum(weights0)
                datXtest.LOOscaled=(datXtest.LOO-datCoef[,1])/datCoef[,2]
                y.true1test.LOO[i]= mean(datXtest.LOOscaled)
                y.true2test.LOO[i]=sum(datXtest.LOOscaled*weights)
                if (skipMissingVariables ) {
                    y.true1test.LOO[i]= mean(datXtest.LOOscaled,na.rm=TRUE)
                    weightsMatrixLOO= weights
                    weightsMatrixLOO[is.na(datXtest.LOOscaled)]=0
                    rowsum.weightsMatrixLOO=sum(weightsMatrixLOO)
                    weightsMatrixLOO=weightsMatrixLOO/rowsum.weightsMatrixLOO
                    datXscaledweightedLOO= datXtest.LOOscaled* weightsMatrixLOO
                    y.true2test.LOO[i]=sum(datXscaledweightedLOO,na.rm=TRUE)
                } #end of if (skipMissingVariables )
            } # end of for loop
        } # end of if (no.observations>1 )
        datEstimates.LeaveOneOut.CV=data.frame(y.true1= y.true1test.LOO,
                                               y.true2= y.true2test.LOO)
    } # end of if ( LeaveOneOut.CV  )

    if (addLinearModel) {
        y.lmTest=rep(NA, dim(as.matrix(datXtest))[[1]] )
        restNonMissingY= !is.na(y)
        if (sum(restNonMissingY,na.rm=TRUE) >3 ){
            lm1=lm(y~., data=data.frame(datX),na.action=na.exclude)
            y.lmTraining=predict(lm1)
            y.lmTraining=predict(lm1)
            if( !is.null(datXtest)) {
                y.lmTest=predict(lm1,newdata=data.frame(datXtest))}
        }
    }
    if ( !is.null(datXtest) & LeaveOneOut.CV & !addLinearModel ) {
        out=  list( datEstimates=datEstimates,
                    datEstimatestest= datEstimatestest,
                    datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV,
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }
    if ( !is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {
        out=  list( datEstimates=datEstimates,
                    datEstimatestest= datEstimatestest, SD.ytrue2=SD.ytrue2,
                    SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }
    if ( is.null(datXtest) & LeaveOneOut.CV & !addLinearModel) {
        out=  list( datEstimates=datEstimates,
                    datEstimates.LeaveOneOut.CV= datEstimates.LeaveOneOut.CV,
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }
    if ( is.null(datXtest) & !LeaveOneOut.CV & !addLinearModel) {
        out=  list( datEstimates=datEstimates, SD.ytrue2=SD.ytrue2,
                    SD.ytrue3=SD.ytrue3, datVariableInfo= datVariableInfo ) }

    if ( !is.null(datXtest) & LeaveOneOut.CV & addLinearModel ) {
        out=  list( datEstimates=data.frame(datEstimates, y.lm=y.lmTraining),
                    datEstimatestest= data.frame(datEstimatestest,
                                                 y.lm=y.lmTest),
                    datEstimates.LeaveOneOut.CV= data.frame(
                        datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO),
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }

    if ( !is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {
        out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining),
                    datEstimatestest= data.frame(datEstimatestest,
                                                 y.lm=y.lmTest),
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }

    if ( is.null(datXtest) & LeaveOneOut.CV & addLinearModel) {
        out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining),
                    datEstimates.LeaveOneOut.CV= data.frame(datEstimates.LeaveOneOut.CV,y.lm=y.lmLOO),
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }
    if ( is.null(datXtest) & !LeaveOneOut.CV & addLinearModel) {
        out=  list( datEstimates=data.frame(datEstimates,y.lm=y.lmTraining),
                    SD.ytrue2=SD.ytrue2, SD.ytrue3=SD.ytrue3,
                    datVariableInfo= datVariableInfo ) }
    out
} # end of function

# addTraitToPCs ####
#' Add trait information to multi-set module eigengene structure
#'
#' Adds trait information to multi-set module eigengene structure.
#'
#' The function simply \code{cbind}'s the module eigengenes and traits for each
#' set. The number of sets and numbers of samples in each set must be
#' consistent between \code{multiMEs} and \code{multiTraits}.
#'
#' @param multiME Module eigengenes in multi-set format. A vector of lists, one
#' list per set. Each list must contain an element named \code{data} that is a
#' data frame with module eigengenes.
#' @param multiTraits Microarray sample trait(s) in multi-set format. A vector
#' of lists, one list per set. Each list must contain an element named
#' \code{data} that is a data frame in which each column corresponds to a
#' trait, and each row to an individual sample.
#' @return A multi-set structure analogous to the input: a vector of lists, one
#' list per set. Each list will contain a component \code{data} with the merged
#' eigengenes and traits for the corresponding set.
#' @author Peter Langfelder
#' @seealso \code{\link{checkSets}}, \code{\link{moduleEigengenes}}
#' @keywords misc
addTraitToMEs <- function(multiME, multiTraits) {
    nSets = length(multiTraits)
    setsize = checkSets(multiTraits)
    nTraits = setsize$nGenes
    nSamples = setsize$nSamples

    if (length(multiME) != nSets) {
        stop("Numbers of sets in multiME and multiTraits parameters differ -
             must be the same.")
    }

    multiMETs = vector(mode = "list", length = nSets)
    for (set in 1:nSets) {
        trait.subs = multiTraits[[set]]$data
        multiMET = as.data.frame(cbind(multiME[[set]]$data, trait.subs))
        colnames(multiMET) = c(colnames(multiME[[set]]$data),
                               colnames(trait.subs))
        if (!is.null(multiME[[set]]$AET)) {
            AET = as.data.frame(cbind(multiME[[set]]$averageExpr, trait.subs))
            colnames(AET) = c(colnames(multiME[[set]]$averageExpr),
                              colnames(trait.subs))
        }
        multiMETs[[set]] = list(data = multiMET)
    }
    multiMETs
    }
