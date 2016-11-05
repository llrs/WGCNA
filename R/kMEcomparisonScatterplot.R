# Plots the kME values of genes in two groups of expression data for each module in an inputted color vector


#' Function to plot kME values between two comparable data sets.
#'
#' Plots the kME values of genes in two groups of expression data for each
#' module in an inputted color vector.
#'
#'
#' @param datExpr1 The first expression matrix (samples=rows, genes=columns).
#' This can either include only the data for group A (in which case dataExpr2
#' must be entered), or can contain all of the data for groups A and B (in
#' which case inA and inB must be entered).
#' @param datExpr2 The second expression matrix, or set to NULL if all data is
#' from same expression matrix.  If entered, datExpr2 must contain the same
#' genes as datExpr1 in the same order.
#' @param colorh The common color vector (module labels) corresponding to both
#' sets of expression data.
#' @param inA,inB Vectors of TRUE/FALSE indicating whether a sample is in group
#' A/B, or a vector of numeric indices indicating which samples are in group
#' A/B. If datExpr2 is entered, these inputs are ignored (thus default = NULL).
#' For these and all other A/B inputs, "A" corresponds to datExpr1 and "B"
#' corresponds to datExpr2 if datExpr2 is entered; otherwise "A" corresponds to
#' datExpr1[inA,] while "B" corresponds to datExpr1[inB,].
#' @param MEsA,MEsB Either the module eigengenes or NULL (default) in which
#' case the module eigengenes will be calculated.  In inputted, MEs MUST be
#' calculated using "moduleEigengenes(<parameters>)$eigengenes" for function to
#' work properly.
#' @param nameA,nameB The names of these groups (defaults = "A" and "B").  The
#' resulting file name (see below) and x and y axis labels for each scatter
#' plot depend on these names.
#' @param plotAll If TRUE, plot gene-ME correlations for all genes.  If FALSE,
#' plot correlations for only genes in the plotted module (default).  Note that
#' the output file name will be different depending on this parameter, so both
#' can be run without overwriting results.
#' @param noGrey If TRUE (default), the grey module genes are ignored.  This
#' parameter is only used if MEsA and MEsB are calculated.
#' @param maxPlot The maximum number of random genes to include (default=1000).
#' Smaller values lead to smaller and less cluttered plots, usually without
#' significantly affecting the resulting correlations. This parameter is only
#' used if plotAll=TRUE.
#' @param pch See help file for "points". Setting pch=19 (default) produces
#' solid circles.
#' @param fileName Name of the file to hold the plots. Since the output format
#' is pdf, the extension should be .pdf .
#' @param ...  Other plotting parameters that are allowable inputs to
#' verboseScatterplot.
#' @return The default output is a file called
#' "kME_correlations_between_[nameA]_and_[nameB]_[all/inMod].pdf", where
#' [nameA] and [nameB] correspond to the nameA and nameB input parameters, and
#' [all/inMod] depends on whether plotAll=TRUE or FALSE. This output file
#' contains all of the plots as separate pdf images, and will be located in the
#' current working directory.
#' @note The function "pdf", which can be found in the grDevices library, is
#' required to run this function.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#'
#' # Example output file ("kME_correlations_between_A_and_B_inMod.pdf") using simulated data.
#'
#' set.seed = 100
#' ME=matrix(0,50,5)
#' for (i in 1:5) ME[,i]=sample(1:100,50)
#' simData1 = simulateDatExpr5Modules(MEturquoise=ME[,1],MEblue=ME[,2],
#'                           MEbrown=ME[,3],MEyellow=ME[,4], MEgreen=ME[,5])
#' simData2 = simulateDatExpr5Modules(MEturquoise=ME[,1],MEblue=ME[,2],
#'                           MEbrown=ME[,3],MEyellow=ME[,4], MEgreen=ME[,5])
#' kMEcomparisonScatterplot(simData1$datExpr,simData2$datExpr,simData1$truemodule)
#'
#'
kMEcomparisonScatterplot <- function (datExpr1, datExpr2, colorh, inA=NULL,
                                      inB=NULL,
                                      MEsA=NULL, MEsB=NULL, nameA="A",
                                      nameB="B", plotAll=FALSE, noGrey=TRUE,
                                      maxPlot=1000,
                                      pch=19,
                                      fileName = ifelse(plotAll, paste0(
                                          "kME_correlations_between_",
                                          nameA,"_and_", nameB, "_all.pdf"),
                                          paste0("kME_correlations_between_",
                                                 nameA, "_and_", nameB,
                                                 "_inMod.pdf")),
                                      ...) {

# First, get the data
	if (is.null(dim(datExpr1))) {
		write ("Error: datExpr1 must be a matrix",""); return(0)
	}
	if (is.null(datExpr2)){
		datA = datExpr1[inA,]
		datB = datExpr1[inB,]
		if ((is.null(dim(datA)))|(is.null(dim(datB)))) {
			 write ("Error: Check input for inA and inB.",""); return(0)
		}
	} else {
		if (is.null(dim(datExpr2))) {
			write ("Error: datExpr2 must be a matrix",""); return(0)
		}
		datA = datExpr1
		datB = datExpr2
	}
	if ((dim(datA)[2]!=length(colorh))|(dim(datB)[2]!=length(colorh))){
		write ("Error: Both sets of input data and color vector must all have same length.",""); return(0)
	}

	if(is.null(MEsA))
		MEsA = (moduleEigengenes(datA, colors=as.character(colorh), excludeGrey=noGrey))$eigengenes
	if(is.null(MEsB))
		MEsB = (moduleEigengenes(datB, colors=as.character(colorh), excludeGrey=noGrey))$eigengenes
	mods  = substring(names(MEsA),3)
	kMEsA = as.data.frame(cor(datA,MEsA,use="p"))
	kMEsB = as.data.frame(cor(datB,MEsB,use="p"))

# Second, make the plots
	xlab  = paste("kME values in",nameA)
	ylab  = paste("kME values in",nameB)
        printFlush(paste("Plotting kME scatterplots into file", fileName))
	if (plotAll){
		pdf(file=fileName)
		numPlot = min(maxPlot,length(colorh))
		these   = sample(1:length(colorh),numPlot)
		for (i in 1:length(mods)){
			plotCol = mods[i]; if(mods[i]=="white") plotCol="black"
			verboseScatterplot(kMEsA[these,i],kMEsB[these,i],main=mods[i],
							   xlab=xlab,ylab=ylab,pch=pch,col=plotCol,...)
		}
		dev.off()
		return("DONE - Plotted All")
	}
	pdf(file=fileName)
	for (i in 1:length(mods)){
		these   = colorh==mods[i]
		plotCol = mods[i]; if(mods[i]=="white") plotCol="black"
		verboseScatterplot(kMEsA[these,i],kMEsB[these,i],main=mods[i],
						   xlab=xlab,ylab=ylab,pch=pch,col=plotCol,...)
	}
	dev.off()
	return("DONE - Plotted only in module")
}
