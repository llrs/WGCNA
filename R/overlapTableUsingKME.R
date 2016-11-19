# Determines significant overlap between modules in two networks based on kME tables.


#' Determines significant overlap between modules in two networks based on kME
#' tables.
#'
#' Takes two sets of expression data (or kME tables) as input and returns a
#' table listing the significant overlap between each module in each data set,
#' as well as the actual genes in common for every module pair.  Modules can be
#' defined in several ways (generally involving kME) based on user input.
#'
#'
#' @param dat1,dat2 Either expression data sets (with samples as rows and genes
#' as columns) or module membership (kME) tables (with genes as rows and
#' modules as columns).  Function reads these inputs based on whether
#' datIsExpression=TRUE or FALSE.  ***Be sure that these inputs include
#' relevant row and column names, or else the function will not work
#' properly.***
#' @param colorh1,colorh2 Color vector (module assignments) corresponding to
#' the genes from dat1/2.  This vector must be the same length as the Gene
#' dimension from dat1/2.
#' @param MEs1,MEs2 If entered (default=NULL), these are the module eigengenes
#' that will be used to form the kME tables. Rows are samples and columns are
#' module assignments.  Note that if datIsExpression=FALSE, these inputs are
#' ignored.
#' @param name1,name2 The names of the two data sets being compared.  These
#' names affect the output parameters.
#' @param cutoffMethod This variable is used to determine how modules are
#' defined in each data set.  Must be one of four options: (1) "assigned" ->
#' use the module assignments in colorh (default); (2) "kME" -> any gene with
#' kME > cutoff is in the module; (3) "numGenes" -> the top cutoff number of
#' genes based on kME is in the module; and (4) "pvalue" -> any gene with
#' correlation pvalue < cutoff is in the module (this includes both positively
#' and negatively-correlated genes).
#' @param cutoff For all cutoffMethods other than "assigned", this parameter is
#' used as the described cutoff value.
#' @param omitGrey If TRUE the grey modules (non-module genes) for both
#' networks are not returned.
#' @param datIsExpression If TRUE (default), dat1/2 is assumed to be expression
#' data.  If FALSE, dat1/2 is assumed to be a table of kME values.
#' @return \item{PvaluesHypergeo}{ A table of p-values showing significance of
#' module overlap based on the hypergeometric test. Note that these p-values
#' are not corrected for multiple comparisons. } \item{AllCommonGenes}{ A
#' character vector of all genes in common between the two data sets. }
#' \item{Genes<name1/2>}{ A list of character vectors of all genes in each
#' module in both data sets.  All genes in the MOD module in data set MM1 could
#' be found using "<outputVariableName>$GenesMM1$MM1_MOD" }
#' \item{OverlappingGenes}{ A list of character vectors of all genes for each
#' between-set comparison from PvaluesHypergeo.  All genes in MOD.A from MM1
#' that are also in MOD.B from MM2 could be found using
#' "<outputVariableName>$OverlappingGenes$MM1_MOD.A_MM2_MOD.B" }
#' @author Jeremy Miller
#' @seealso \code{\link{overlapTable}}
#' @keywords misc
#' @examples
#'
#' # Example: first generate simulated data.
#'
#' set.seed(100)
#' ME.A = sample(1:100,50);  ME.B = sample(1:100,50)
#' ME.C = sample(1:100,50);  ME.D = sample(1:100,50)
#' ME.E = sample(1:100,50);  ME.F = sample(1:100,50)
#' ME.G = sample(1:100,50);  ME.H = sample(1:100,50)
#' ME1     = data.frame(ME.A, ME.B, ME.C, ME.D, ME.E)
#' ME2     = data.frame(ME.A, ME.C, ME.D, ME.E, ME.F, ME.G, ME.H)
#' simDat1 = simulateDatExpr(ME1,1000,c(0.2,0.1,0.08,0.05,0.04,0.3), signed=TRUE)
#' simDat2 = simulateDatExpr(ME2,1000,c(0.2,0.1,0.08,0.05,0.04,0.03,0.02,0.3),
#'                           signed=TRUE)
#'
#' # Now run the function using assigned genes
#' results = overlapTableUsingKME(simDat1$datExpr, simDat2$datExpr,
#'                    labels2colors(simDat1$allLabels), labels2colors(simDat2$allLabels),
#'                    cutoffMethod="assigned")
#' results$PvaluesHypergeo
#'
#' # Now run the function using a p-value cutoff, and inputting the original MEs
#' colnames(ME1) = standardColors(5);  colnames(ME2) = standardColors(7)
#' results = overlapTableUsingKME(simDat1$datExpr, simDat2$datExpr,
#'                       labels2colors(simDat1$allLabels),
#'                       labels2colors(simDat2$allLabels),
#'                       ME1, ME2, cutoffMethod="pvalue", cutoff=0.05)
#' results$PvaluesHypergeo
#'
#' # Check which genes are in common between the black modules from set 1 and
#' # the green module from set 2
#' results$OverlappingGenes$MM1_green_MM2_black
#'
overlapTableUsingKME <- function(dat1, dat2, colorh1, colorh2, MEs1=NULL, MEs2=NULL,
name1="MM1", name2="MM2", cutoffMethod="assigned", cutoff=0.5, omitGrey=TRUE,
datIsExpression=TRUE){

# Run a few tests on the imput data formatting
	if (is.null(dim(dat1))|is.null(dim(dat2))) {
		write("Error: dat1 and dat2 must be matrices.",""); return(0)
	}
	if ((dim(dat1)[datIsExpression+1]!=length(colorh1))|
		(dim(dat2)[datIsExpression+1]!=length(colorh2))){
		write("Error: Both sets of input data and color vectors must have same length.",""); return(0)
	}
	if ((cutoffMethod=="pvalue")&(datIsExpression==FALSE)){
		write("Error: Pvalues are not calculated if datIsExpression=TRUE.  Choose other cutoffMethod.",
			  ""); return(0)
	}

# Find and format the kME values and other variables for both inputs
	G1 = dimnames(dat1)[[datIsExpression+1]];  G2 = dimnames(dat2)[[datIsExpression+1]]
	if(datIsExpression){
		if(is.null(MEs1))
			MEs1 = (moduleEigengenes(dat1, colors=as.character(colorh1), excludeGrey=omitGrey))$eigengenes
		if(is.null(MEs2))
			MEs2 = (moduleEigengenes(dat2, colors=as.character(colorh2), excludeGrey=omitGrey))$eigengenes
		mods1 = colnames(MEs1);  mods2 = colnames(MEs2)
		if (length(grep("ME",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("PC",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("ME",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		if (length(grep("PC",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		out = corAndPvalue(dat1,MEs1);  MM1 = out$cor;  PV1 = out$p;  rm(out)
		out = corAndPvalue(dat2,MEs2);  MM2 = out$cor;  PV2 = out$p;  rm(out)
		colnames(MM1) <- colnames(PV1) <- mods1
		colnames(MM2) <- colnames(PV2) <- mods2
		rownames(MM1) <- rownames(PV1) <- G1
		rownames(MM2) <- rownames(PV2) <- G2
	} else {
		MM1 = dat1[,sort(colnames(dat1))];  mods1 = colnames(MM1)
		MM2 = dat2[,sort(colnames(dat2))];  mods2 = colnames(MM2)
		if (length(grep("ME",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("PC",mods1))==length(mods1)) mods1 = substr(mods1,3,nchar(mods1))
		if (length(grep("ME",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		if (length(grep("PC",mods2))==length(mods2)) mods2 = substr(mods2,3,nchar(mods2))
		colnames(MM1) = mods1;  colnames(MM2) = mods2
		rownames(MM1) = G1;     rownames(MM2) = G2
		if(omitGrey){
			MM1 = MM1[,!is.element(mods1,"grey")];  mods1 = colnames(MM1)
			MM2 = MM2[,!is.element(mods2,"grey")];  mods2 = colnames(MM2)
		}
	}
	if ((length(setdiff(mods1,as.character(colorh1)))>omitGrey)|
		(length(setdiff(mods2,as.character(colorh2)))>omitGrey)){
		write("MEs cannot include colors with no genes assigned.",""); return(0)
	}
	l1 = length(mods1);	 l2 = length(mods2)
	cutoffMethod = substr(cutoffMethod,1,1)
	names=c(name1,name2)
	comGenes = sort(unique(intersect(G1,G2)));  total   = length(comGenes)
	MM1     = MM1[comGenes,];  MM2 = MM2[comGenes,]
    if (datIsExpression){
		PV1 = PV1[comGenes,];  PV2 = PV2[comGenes,]
	}
	names(colorh1) = G1;  colorh1 = colorh1[comGenes]
	names(colorh2) = G2;  colorh2 = colorh2[comGenes]

# Assign each gene in each module to a vector corresponding to the modules
	genes1 <- genes2 <- list()
	if (cutoffMethod=="a"){
		for (i in 1:l1)  genes1[[i]] = comGenes[colorh1==mods1[i]]
		for (i in 1:l2)  genes2[[i]] = comGenes[colorh2==mods2[i]]
	} else if (cutoffMethod=="p") {
		for (i in 1:l1)  genes1[[i]] = comGenes[PV1[,mods1[i]]<=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[PV2[,mods2[i]]<=cutoff]
	} else if (cutoffMethod=="k") {
		for (i in 1:l1)  genes1[[i]] = comGenes[MM1[,mods1[i]]>=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[MM2[,mods2[i]]>=cutoff]
	} else if (cutoffMethod=="n") {
		for (i in 1:l1)  genes1[[i]] = comGenes[rank(-MM1[,mods1[i]])<=cutoff]
		for (i in 1:l2)  genes2[[i]] = comGenes[rank(-MM2[,mods2[i]])<=cutoff]
	} else {
		write("ERROR: cutoffMethod entered is not supported.",""); return(0)
	}
	names(genes1) = paste(names[1],mods1,sep="_")
	names(genes2) = paste(names[2],mods2,sep="_")

# Determine signficance of each comparison and write out all of the gene lists
	ovGenes = list()
	ovNames = rep("",l1*l2)
	pVals   = matrix(1, nrow=l1, ncol=l2)
	rownames(pVals) = paste(names[1],mods1,sep="_")
	colnames(pVals) = paste(names[2],mods2,sep="_")
	i = 0
	for (m1 in 1:l1) for (m2 in 1:l2) {
		i = i+1
		ovGenes[[i]] = sort(unique(intersect(genes1[[m1]],genes2[[m2]])))
		pVals[m1,m2] = .phyper2(total,length(genes1[[m1]]), length(genes2[[m2]]),length(ovGenes[[i]]))
		if (pVals[m1,m2]>10^(-10))   pVals[m1,m2] =
		.phyper2(total,length(genes1[[m1]]), length(genes2[[m2]]),length(ovGenes[[i]]),FALSE)
		ovNames[i] = paste(names[1],mods1[m1],names[2],mods2[m2],sep="_")
	}
	names(ovGenes) = ovNames
	out = list(pVals,comGenes,genes1,genes2,ovGenes)
	names(out) = c("PvaluesHypergeo", "AllCommonGenes", paste0("Genes", names),
	               "OverlappingGenes")
	return(out)
}

.phyper2 <- function (total, group1, group2, overlap, verySig=TRUE ,lt=TRUE){
# This function is the same is phyper, just allows for more sensible input values
	q = overlap
	m = group1
	n = total-group1
	k = group2
	prob = phyper(q, m, n, k, log.p = verySig, lower.tail=lt)
	if (verySig) return(-prob)
	return(1-prob)
}
