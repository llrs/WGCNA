# This function barplots data across two splitting parameters


#' Bar plots of data across two splitting parameters
#' 
#' This function takes an expression matrix which can be split using two
#' separate splitting parameters (ie, control vs AD with multiple brain
#' regions), and plots the results as a barplot. Group average, standard
#' deviations, and relevant Kruskal-Wallis p-values are returned.
#' 
#' 
#' @param expAll An expression matrix, with rows as samples and genes/probes as
#' columns.  If genes=NA, then column names must be included.
#' @param groups A character vector corresponding to the samples in expAll,
#' with each element the group name of the relevant sample or NA for samples
#' not in any group.  For, example: NA, NA, NA, Con, Con, Con, Con, AD, AD, AD,
#' AD, NA, NA.  This trait will be plotted as adjacent bars for each split.
#' @param split A character vector corresponding to the samples in expAll, with
#' each element the group splitting name of the relevant sample or NA for
#' samples not in any group.  For, example: NA, NA, NA, Hip, Hip, EC, EC, Hip,
#' Hip, EC, EC, NA, NA.  This trait will be plotted as the same color across
#' each split of the barplot.  For the function to work properly, the same
#' split values should be inputted for each group.
#' @param subset A list of one or more genes to compare the expression with.
#' If the list contains more than one gene, the first element contains the
#' group name. For example, Ribosomes, RPL3, RPL4, RPS3.
#' @param genes If entered, this parameter is a list of gene/probe identifiers
#' corresponding to the columns in expAll.
#' @param scale For subsets of genes that include more than one gene, this
#' parameter determines how the genes are combined into a single value.
#' Currently, there are five options: 1) ("N")o scaling (default); 2) first
#' divide each gene by the ("A")verage across samples; 3) first scale genes to
#' ("Z")-score across samples; 4) only take the top ("H")ub gene (ignore all
#' but the highest-connected gene); and 5) take the ("M")odule eigengene.  Note
#' that these scaling methods have not been sufficiently tested, and should be
#' considered experimental.
#' @param graph If TRUE (default), bar plot is made.  If FALSE, only the
#' results are returned, and no plot is made.
#' @param cex1 Sets the graphing parameters of cex.axis and cex.names
#' (default=1.5)
#' @param las1 Sets the graphing parameter las (default=2).
#' @param \dots Other graphing parameters allowed in the barplot function.
#' Note that the parameters for cex.axis, cex.names, and las are superseded by
#' cex1 and las1 and will therefore be ignored.
#' @return \item{splitGroupMeans}{ The group/split averaged expression across
#' each group and split combination.  This is the height of the bars in the
#' graph. } \item{splitGroupSDs}{ The standard deviation of group/split
#' expression across each group and split combination.  This is the height of
#' the error bars in the graph. } \item{splitPvals}{ Kruskal-Wallis p-values
#' for each splitting parameter across groups. } \item{groupPvals}{
#' Kruskal-Wallis p-values for each group parameter across splits. }
#' @author Jeremy Miller
#' @seealso \code{\link{barplot}}, \code{\link{verboseBarplot}}
#' @keywords misc
#' @examples
#' 
#' # Example: first simulate some data
#' set.seed(100)
#' ME.A = sample(1:100,50);  ME.B = sample(1:100,50)
#' ME.C = sample(1:100,50);  ME.D = sample(1:100,50)  
#' ME1     = data.frame(ME.A, ME.B, ME.C, ME.D)
#' simDatA = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.3), signed=TRUE)
#' datExpr = simDatA$datExpr+5
#' datExpr[1:10,]  = datExpr[1:10,]+2
#' datExpr[41:50,] = datExpr[41:50,]-1
#' 
#' # Now split up the data and plot it!
#' subset  = c("Random Genes", "Gene.1", "Gene.234", "Gene.56", "Gene.789")
#' groups  = rep(c("A","A","A","B","B","B","C","C","C","C"),5)
#' split   = c(rep("ZZ",10), rep("YY",10), rep("XX",10), rep("WW",10), rep("VV",10))
#' par(mfrow = c(1,1))
#' results = stratifiedBarplot(datExpr, groups, split, subset)
#' results
#' 
#' # Now plot it the other way
#' results = stratifiedBarplot(datExpr, split, groups, subset)
#' 
#' 
stratifiedBarplot = function (expAll, groups, split, subset, 
genes=NA, scale="N", graph=TRUE, las1=2, cex1=1.5, ...){

## Code to take care of array formatting
	expAll = t(expAll)
	if (length(subset)>1){
		subsetName = subset[1]; 
		subset = subset[2:length(subset)]
	} else { subsetName=subset }
	groupNames = as.character(names(.tableOrd(groups)))
	splitNames = as.character(names(.tableOrd(split))) 
	if (is.na(genes)[1]) genes = rownames(expAll)
	keep   = !is.na(groups)
	groups = groups[keep]
	split  = split[keep]
	expAll = expAll[,keep]
	scale  = substr(scale,1,1)
	
## Collect and scale the expression data
	expSubset = expAll[is.element(genes,subset),]
	if(length(subset)>1){
		if(scale=="A")  expSubset = t(apply(expSubset,1,function(x) return(x/mean(x))))
		if(scale=="Z")  expSubset = t(apply(expSubset,1,function(x) return((x-mean(x))/sd(x))))
		if(scale=="H")  {
			AdjMat = adjacency(t(expSubset),type="signed",power=2)
			diag(AdjMat) = 0
			Degree = rowSums(AdjMat)
			keep   = which(Degree == max(Degree))
			expSubset = expSubset[keep,]	
		}
		if(scale=="M")  {
			me = moduleEigengenes(as.matrix(t(expSubset)), rep("blue",dim(expSubset)[1]))
			expSubset = me$eigengenes$MEblue
		} 
		expSubset = rbind(expSubset,expSubset)
		expSubset = apply(expSubset,2,mean)
	}
	
## Now average the data and output it in a meaningful way
	exp <- std <- matrix(0,ncol=length(splitNames),nrow=length(groupNames))
	splitPvals = rep(1,length(splitNames))
	names(splitPvals) = splitNames
	groupPvals = rep(1,length(groupNames))
	names(groupPvals) = groupNames
	for (c in 1:length(splitNames)){
		expTmp = expSubset[split==splitNames[c]]
		grpTmp = groups[split==splitNames[c]]
		splitPvals[c] = kruskal.test(expTmp,as.factor(grpTmp))$p.value
		for (r in 1:length(groupNames)){
			exp[r,c] = mean(expSubset[(groups==groupNames[r])&(split==splitNames[c])])
			std[r,c] = sd(expSubset[(groups==groupNames[r])&(split==splitNames[c])])
			if(c==1){
				expTmp = expSubset[groups==groupNames[r]]
				splTmp = split[groups==groupNames[r]]
				groupPvals[r] = kruskal.test(expTmp,as.factor(splTmp))$p.value	
			}
		}
	}
	colnames(exp) <- colnames(std) <- splitNames
	rownames(exp) <- rownames(std) <- groupNames
	
## Now plot the results, if requested
	if(graph){
		ylim = c(min(0,min(min(exp-std))),max(max(exp+std))*(1+0.08*length(groupNames)))
		barplot(exp, beside=TRUE, legend.text=TRUE, main=subsetName, las=las1,
				ylim=ylim, cex.axis=cex1, cex.names=cex1,  ...)
		.err.bp(exp,std,TRUE)
	}
	
# Now collect the output and return it.
	out = list(splitGroupMeans = exp, splitGroupSDs = std, 
			   splitPvals = splitPvals, groupPvals = groupPvals)
	return(out)
}
# --------------------------------


.err.bp<-function(daten,error,two.side=F){
## This function was written by Steve Horvath
# The function err.bp  is used to create error bars in a barplot
# usage: err.bp(as.vector(means), as.vector(stderrs), two.side=F)
	if(!is.numeric(daten)) {
	stop("All arguments must be numeric")}
	if(is.vector(daten)){ 
		xval<-(cumsum(c(0.7,rep(1.2,length(daten)-1)))) 
	}else{
		if (is.matrix(daten)){
			xval<-cumsum(array(c(1,rep(0,dim(daten)[1]-1)),
							   dim=c(1,length(daten))))+0:(length(daten)-1)+.5
		}else{
		stop("First argument must either be a vector or a matrix") }
	}
	MW<-0.25*(max(xval)/length(xval)) 
	ERR1<-daten+error 
	ERR2<-daten-error
	for(i in 1:length(daten)){
		segments(xval[i],daten[i],xval[i],ERR1[i])
		segments(xval[i]-MW,ERR1[i],xval[i]+MW,ERR1[i])
		if(two.side){
			segments(xval[i],daten[i],xval[i],ERR2[i])
			segments(xval[i]-MW,ERR2[i],xval[i]+MW,ERR2[i])
		} 
	} 
} 


.tableOrd = function (input){
## This is the same as the "table" function but retains the order
	
## This internal function collects the order
	tableOrd2 = function(input, output=NULL){
		input  = input[!is.na(input)]
		if (length(input)==0) return (output)
		outTmp = input[1] 
		output = c(output, outTmp)
		input  = input[input!=outTmp]
		output = tableOrd2(input, output)
		return(output)
	}
	
## Get the results
	tableOut   = table(input)
	tableOrder = tableOrd2(input)
	return(tableOut[tableOrder])
}  
