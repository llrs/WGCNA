## Merge modules and reassign module genes based on kME values


#' Merge modules and reassign genes using kME.
#' 
#' This function takes an expression data matrix (and other user-defined
#' parameters), calculates the module membership (kME) values, and adjusts the
#' module assignments, merging modules that are not sufficiently distinct and
#' reassigning modules that were originally assigned suboptimally.
#' 
#' 
#' @param datExpr An expression data matrix, with samples as rows, genes (or
#' probes) as column.
#' @param colorh The color vector (module assignments) corresponding to the
#' columns of datExpr.
#' @param ME Either NULL (default), at which point the module eigengenes will
#' be calculated, or pre-calculated module eigengenes for each of the modules,
#' with samples as rows (corresponding to datExpr), and modules corresponding
#' to columns (column names MUST be module colors or module colors prefixed by
#' "ME" or "PC").
#' @param threshPercent Threshold percent of the number of genes in the module
#' that should be included for the various analyses. For example, in a module
#' with 200 genes, if threshPercent=50 (default), then 50 genes will be checked
#' for reassignment and used to test whether two modules should be merged.  See
#' also threshNumber.
#' @param mergePercent If greater than this percent of the assigned genes are
#' above the threshold are in a module other than the assigned module, then
#' these two modules will be merged.  For example, if mergePercent=25
#' (default), and the 70 out of 200 genes in the blue module were more highly
#' correlated with the black module eigengene, then all genes in the blue
#' module would be reassigned to the black module.
#' @param reassignModules If TRUE (default), genes are resassigned to the
#' module with which they have the highest module membership (kME), but only if
#' their kME is above the threshPercent (or threshNumber) threshold of that
#' module.
#' @param convertGrey If TRUE (default), unassigned (grey) genes are assigned
#' as in "reassignModules"
#' @param omitColors These are all of the module assignments which indicate
#' genes that are not assigned to modules (default="grey").  These genes will
#' all be assigned as "grey" by this function.
#' @param reassignScale A value between 0 and 1 (default) which determines how
#' the threshPercent gets scaled for reassigning genes.  Smaller values
#' reassign more genes, but does not affect the merging process.
#' @param threshNumber Either NULL (default) or, if entered, every module is
#' counted as having exactly threshNumber genes, and threshPercent it ignored.
#' This parameter should have the effect of
#' @return \item{moduleColors}{ The NEW color vector (module assignments)
#' corresponding to the columns of datExpr, after module merging and
#' reassignments. } \item{mergeLog}{ A log of the order in which modules were
#' merged, for reference. }
#' @note Note that this function should be considered "experimental" as it has
#' only been beta tested.  Please e-mail jeremyinla@gmail.com if you have any
#' issues with the function.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#' 
#' 
#' ## First simulate some data and the resulting network dendrogram
#' set.seed(100)
#' MEturquoise = sample(1:100,50)
#' MEblue      = sample(1:100,50)
#' MEbrown     = sample(1:100,50)
#' MEyellow    = sample(1:100,50) 
#' MEgreen     = c(MEyellow[1:30], sample(1:100,20))
#' MEred	    = c(MEbrown [1:20], sample(1:100,30))
#' #MEblack	    = c(MEblue  [1:25], sample(1:100,25))
#' ME   = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)#, MEblack)
#' dat1 = simulateDatExpr(ME, 400, c(0.15,0.13,0.12,0.10,0.09,0.09,0.1), signed=TRUE)
#' TOM1  = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
#' tree1 = fastcluster::hclust(as.dist(1-TOM1),method="average")
#' 
#' ## Here is an example using different mergePercentages, 
#' # setting an inclusive threshPercent (91)
#' colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
#' merges = c(65,40,20,5)
#' for (m in merges)  
#'    colorPlot = cbind(colorPlot, 
#'                      moduleMergeUsingKME(dat1$datExpr,colorh1,
#'                         threshPercent=91, mergePercent=m)$moduleColors)
#' plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)
#' 
#' ## Here is an example using a lower reassignScale (so that more genes get reassigned)
#' colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
#' merges = c(65,40,20,5)
#' for (m in merges)  colorPlot = cbind(colorPlot, 
#'   moduleMergeUsingKME(dat1$datExpr,colorh1,threshPercent=91, 
#'                       reassignScale=0.7, mergePercent=m)$moduleColors)
#' plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)
#' 
#' ## Here is an example using a less-inclusive threshPercent (75), 
#' # little if anything is merged.
#' 
#' colorh1  <- colorPlot <- labels2colors(dat1$allLabels)
#' merges = c(65,40,20,5)
#' for (m in merges)  colorPlot = cbind(colorPlot, 
#'   moduleMergeUsingKME(dat1$datExpr,colorh1,
#'                       threshPercent=75, mergePercent=m)$moduleColors)
#' plotDendroAndColors(tree1, colorPlot, c("ORIG",merges), dendroLabels=FALSE)
#' # (Note that with real data, the default threshPercent=50 usually results 
#' # in some modules being merged)
#' 
#' 
moduleMergeUsingKME <- function (datExpr, colorh, ME=NULL, threshPercent=50, mergePercent = 25, reassignModules=TRUE, 
convertGrey=TRUE, omitColors="grey", reassignScale=1, threshNumber=NULL) {

## First assign all of the variables and put everything in the correct format.
	if (length(colorh)!=dim(datExpr)[2]){
		write("Error: color vector much match inputted datExpr columns.",""); return(0)
	}
	if (is.null(ME))  
		ME = (moduleEigengenes(datExpr, colors=as.character(colorh), excludeGrey=TRUE))$eigengenes
	if (dim(ME)[1]!=dim(datExpr)[1]){
		write("Error: ME rows much match inputted datExpr rows (samples).",""); return(0)
	}
	modules = colnames(ME)
	if (length(grep("ME",modules))==length(modules)) modules = substr(modules,3,nchar(modules))
	if (length(grep("PC",modules))==length(modules)) modules = substr(modules,3,nchar(modules))
	if (length(setdiff(modules,as.character(colorh)))>0){
		write("ME cannot include colors with no genes assigned.",""); return(0)
	}
	names(ME) = modules
	
	datCorrs = as.data.frame(cor(datExpr,ME,use="p"))
	colnames(datCorrs) = modules
	modules  = sort(modules[!is.element(modules,omitColors)])
	modulesI = modules # To test whether merging occurs
	datCorrs = datCorrs[,modules]
	iteration = 1
	if(is.null(colnames(datExpr))) colnames(datExpr) = as.character(1:length(colorh))
	rownames(datCorrs) <- colnames(datExpr)
	colorOut = colorh
	colorOut[is.element(colorOut,omitColors)] = "grey"
	mergeLog = NULL
	datExpr = t(datExpr) # For consistency with how the function was originally written.

## Iteratively run this function until no further changes need to be made
	while (!is.na(iteration)){
		write("",""); write("__________________________________________________","")
		write(paste("This is iteration #",iteration,". There are ",length(modules)," modules.",sep=""),"")
		iteration = iteration+1
		
	## Reassign modules if requested by reassignModules and convertGrey
		colorMax  = NULL
		whichMod  = apply(datCorrs,1,which.max)
		cutNumber = round(table(colorOut)*threshPercent/100)
		cutNumber = cutNumber[names(cutNumber)!="grey"] 
		if(!is.null(threshNumber)) cutNumber = rep(threshNumber,length(cutNumber))
		cutNumber = apply(cbind(cutNumber,10),1,max)
		for (i in 1:length(whichMod)) 
			colorMax = c(colorMax,modules[whichMod[i]])
		for (i in 1:length(modules)){
			corrs    = as.numeric(datCorrs[,i])
			cutValue = sort(corrs[colorOut==modules[i]],decreasing=TRUE)[cutNumber[i]]
			inModule = corrs>(cutValue*reassignScale)		
			if(convertGrey)
				colorOut[inModule&(colorOut=="grey")&(colorMax==modules[i])]=modules[i]
			if(reassignModules)
				colorOut[inModule&(colorOut!="grey")&(colorMax==modules[i])]=modules[i]
		}
		
	## Merge all modules meeting the mergePercent and threshPercent criteria
		for (i in 1:length(modules)){
			cutNumber  = round(table(colorOut)*threshPercent/100)
			cutNumber  = cutNumber[names(cutNumber)!="grey"]
			if(!is.null(threshNumber)) cutNumber = rep(threshNumber,length(cutNumber))
			cutNumber = apply(cbind(cutNumber,10),1,max)
			corrs      = as.numeric(datCorrs[,i])
			# Make sure you do not include more genes than are in the module
			numInMod   = sum(colorOut==modules[i])
			cutValue   = sort(corrs[colorOut==modules[i]],decreasing=TRUE)[min(numInMod,cutNumber[modules[i]])]
			colorMod   = colorOut[corrs>=cutValue]
			colorMod   = colorMod[colorMod!="grey"]
			modPercent = 100*table(colorMod)/length(colorMod)
			modPercent = modPercent[names(modPercent)!=modules[i]]
			if(length(modPercent)>1) if(max(modPercent)>mergePercent){
				whichModuleMerge = names(modPercent)[which.max(modPercent)]
				colorOut[colorOut==modules[i]] = whichModuleMerge
				write(paste(modules[i],"has been merged into",whichModuleMerge,"."),"")
				mergeLog = rbind(mergeLog,c(modules[i],whichModuleMerge))
			}
		}
		
	## If no modules were merged, then set iteration to NA 
		modules  = sort(unique(colorOut))		
		modules  = modules[modules!="grey"]
		if (length(modules)==length(modulesI)) iteration=NA
		modulesI = modules
		
	## Recalculate the new module membership values
		MEs      = (moduleEigengenes(t(datExpr), colors=as.character(colorOut)))$eigengenes
		MEs      = MEs[,colnames(MEs)!="MEgrey"]
		datCorrs = as.data.frame(cor(t(datExpr),MEs,use="p"))
		colnames(datCorrs) = modules
	}
	if(!is.null(dim(mergeLog))){
		colnames(mergeLog) = c("Old Module","Merged into New Module")
		rownames(mergeLog) = paste("Merge #",1:dim(mergeLog)[1],sep="")
	}
	return(list(moduleColors=colorOut,mergeLog=mergeLog))
}

