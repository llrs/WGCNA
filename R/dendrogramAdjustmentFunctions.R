## This file contains several functions which can be used to adjust the dendrogram
#   in ways which keep the dendrogram mathematically identical (ie, branch swapping,
#   branch reflection, etc).  The goal is to biologically optimize the dendrogram.

# ----------------------------------------------------------------------------- #



#' Optimize dendrogram using branch swaps and reflections.
#' 
#' This function takes as input the hierarchical clustering tree as well as a
#' subset of genes in the network (generally corresponding to branches in the
#' tree), then returns a semi-optimally ordered tree.  The idea is to maximize
#' the correlations between adjacent branches in the dendrogram, in as much as
#' that is possible by adjusting the arbitrary positionings of the branches by
#' swapping and reflecting branches.
#' 
#' 
#' @param hierTOM A hierarchical clustering object (or gene tree) that is used
#' to plot the dendrogram.  For example, the output object from the function
#' hclust or fastcluster::hclust.  Note that elements of hierTOM$order MUST be
#' named (for example, with the corresponding gene name).
#' @param datExpr Gene expression data with rows as samples and columns as
#' genes, or NULL if a pre-made adjacency is entered.  Column names of datExpr
#' must be a subset of gene names of hierTOM$order.
#' @param colorh The module assignments (color vectors) corresponding to the
#' rows in datExpr, or NULL if a pre-made adjacency is entered.
#' @param type What type of network is being entered.  Common choices are
#' "signed" (default) and "unsigned".  With "signed" negative correlations
#' count against, whereas with "unsigned" negative correlations are treated
#' identically as positive correlations.
#' @param adj Either NULL (default) or an adjacency (or any other square)
#' matrix with rows and columns corresponding to a subset of the genes in
#' hierTOM$order.  If entered, datExpr, colorh, and type are all ignored.
#' Typically, this would be left blank but could include correlations between
#' module eigengenes, with rows and columns renamed as genes in the
#' corresponding modules, for example.
#' @param iter The number of iterations to run the function in search of
#' optimal branch ordering.  The default is the square of the number of modules
#' (or the quare of the number of genes in the adjacency matrix).
#' @param useReflections If TRUE, both reflections and branch swapping will be
#' used to optimize dendrogram.  If FALSE (default) only branch swapping will
#' be used.
#' @param allowNonoptimalSwaps If TRUE, there is chance (that decreases with
#' each iteration) of swapping / reflecting branches whether or not the new
#' correlation between expression of genes in adjacent branches is better or
#' worse.  The idea (which has not been sufficiently tested), is that this
#' would prevent the function from getting stuck at a local maxima of
#' correlation.  If FALSE (default), the swapping / reflection of branches only
#' occurs if it results in a higher correlation between adjacent branches.
#' @return \item{hierTOM}{ A hierarchical clustering object with the
#' hierTOM$order variable properly adjusted, but all other variables identical
#' as the heirTOM input. } \item{changeLog}{ A log of all of the changes that
#' were made to the dendrogram, including what change was made, on what
#' iteration, and the Old and New scores based on correlation.  These scores
#' have arbitrary units, but higher is better. }
#' @note This function is very slow and is still in an *experimental* function.
#' We have not had problems with ~10 modules across ~5000 genes, although
#' theoretically it should work for many more genes and modules, depending upon
#' the speed of the computer running R.  Please address any problems or
#' suggestions to jeremyinla@gmail.com.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#' 
#' ## Example: first simulate some data.
#' 
#' MEturquoise = sample(1:100,50)
#' MEblue      = c(MEturquoise[1:25], sample(1:100,25))
#' MEbrown     = sample(1:100,50)
#' MEyellow    = sample(1:100,50) 
#' MEgreen     = c(MEyellow[1:30], sample(1:100,20))
#' MEred	    = c(MEbrown [1:20], sample(1:100,30))
#' ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)
#' dat1   = simulateDatExpr(ME,400,c(0.16,0.12,0.11,0.10,0.10,0.10,0.1), signed=TRUE)
#' TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
#' colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
#' tree1  = fastcluster::hclust(as.dist(1-TOM1),method="average")
#' colorh = labels2colors(dat1$allLabels)
#' 
#' plotDendroAndColors(tree1,colorh,dendroLabels=FALSE)
#' 
#' ## Reassign modules using the selectBranch and chooseOneHubInEachModule functions
#' 
#' datExpr = dat1$datExpr
#' hubs    = chooseOneHubInEachModule(datExpr, colorh)
#' colorh2 = rep("grey", length(colorh))
#' colorh2 [selectBranch(tree1,hubs["blue"],hubs["turquoise"])] = "blue"
#' colorh2 [selectBranch(tree1,hubs["turquoise"],hubs["blue"])] = "turquoise"
#' colorh2 [selectBranch(tree1,hubs["green"],hubs["yellow"])]   = "green"
#' colorh2 [selectBranch(tree1,hubs["yellow"],hubs["green"])]   = "yellow"
#' colorh2 [selectBranch(tree1,hubs["red"],hubs["brown"])]      = "red"
#' colorh2 [selectBranch(tree1,hubs["brown"],hubs["red"])]      = "brown"
#' plotDendroAndColors(tree1,cbind(colorh,colorh2),c("Old","New"),dendroLabels=FALSE)
#' 
#' ## Now swap and reflect some branches, then optimize the order of the branches 
#' # and output pdf with resulting images
#' 
#' \dontrun{
#' pdf("DENDROGRAM_PLOTS.pdf",width=10,height=5)
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Starting Dendrogram")
#' 
#' tree1 = swapTwoBranches(tree1,hubs["red"],hubs["turquoise"])
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Swap blue/turquoise and red/brown")
#' 
#' tree1 = reflectBranch(tree1,hubs["blue"],hubs["green"])
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Reflect turquoise/blue")
#' 
#' # (This function will take a few minutes)
#' out = orderBranchesUsingHubGenes(tree1,datExpr,colorh2,useReflections=TRUE,iter=100)
#' tree1 = out$geneTree
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Semi-optimal branch order")
#' 
#' out$changeLog
#' 
#' dev.off()
#' }
#' 
orderBranchesUsingHubGenes <- function(hierTOM, datExpr=NULL, colorh=NULL, type="signed", 
adj=NULL, iter=NULL, useReflections=FALSE, allowNonoptimalSwaps=FALSE){

# First read in and format all of the variables
	hGenes = hierTOM$labels
	if(is.null(adj)){
		genes = chooseOneHubInEachModule(datExpr, colorh, type=type)
		adj   = adjacency(datExpr[,genes],type=type, power=(2-(type=="unsigned")))
		colnames(adj) <- rownames(adj) <- genes
	}
	genes = rownames(adj)
	if(length(genes)!=length(intersect(genes,hGenes))){
		write("All genes in the adjacency must also be in the gene tree.  Check to make sure","")
		write("that names(hierTOM$labels) is set to the proper gene or probe names and that","")
		write("these correspond to the expression / adjacency gene names.","")
		return(0)
	}
	genes = hGenes[is.element(hGenes,genes)]
	adj   = adj[genes,genes]
	if (is.null(iter)) iter = length(genes)^2
	iters=(1:iter)/iter
	swapAnyway = rep(0,length(iters))  # Quickly decreasing chance of random swap
	if (allowNonoptimalSwaps)  swapAnyway = ((1-iters)^3)/3+0.001
		
# Iterate random swaps in the branch, only accepting the new result
#  if it produces a higher correlation than the old result OR if the 
#  random variable says to swap (which gets less likely each iteration)
	changes=NULL
	for (i in 1:iter){
		swap = 1; 
		if (useReflections) swap = sample(0:1,1)
		gInd = sample(1:length(genes),2)
		g    = genes[gInd]
		if (swap==1) {
			hierTOMnew = swapTwoBranches(hierTOM, g[1], g[2])
		} else hierTOMnew = reflectBranch(hierTOM, g[1], g[2], TRUE)
		oldSum    = .offDiagonalMatrixSum(adj)
		oGenesNew = hGenes[hierTOMnew$order]
		oGenesNew = oGenesNew[oGenesNew%in%genes]
		adjNew    = adj[oGenesNew,oGenesNew]
		newSum    = .offDiagonalMatrixSum(adjNew)
		if ((newSum>oldSum)|((sample(1:1000,1)/1000)<swapAnyway[i])) {
			hierTOM = hierTOMnew
			changes = rbind(changes,c(i,ifelse(swap==1,"Swap","Reflect"),g,oldSum,newSum))
			adj     = adjNew
		}
		write(paste("Interation",i,"of",iter),"")
		collectGarbage()
	}
	
# Perform all of the suggested swappings on the input network.
	
# Output the results
	colnames(changes)=c("Iter.#","Swap?","Gene1","Gene2","OldScore","NewScore")
    out = list(geneTree = hierTOM, changeLog = changes)
	return(out)
}

# ----------------------------------------------------------------------------- #

selectBranch <- function (hierTOM, g1, g2){
## This function selects of all genes in a branch given a gene in the
##  branch (g1) and a gene in a neighboring branch (g2), returning the 
##  indices for genes in the branch in the hierTOM$labels vector 
	
# Convert genes to UNORDERED indices (if given, indices should be ordered)
	if(is.numeric(g1)) g1 = hierTOM$order[g1]
	if(is.numeric(g2)) g2 = hierTOM$order[g2]
	if(!is.numeric(g1)) g1 = which(hierTOM$labels==g1)
	if(!is.numeric(g2)) g2 = which(hierTOM$labels==g2)
	if((length(g1)==0)|(length(g2)==0)|(max(c(g1,g2))>length(hierTOM$labels))){
		write("Input genes are not both legal indices","")
		return(hierTOM)
	}
	
# Now determine which branch is the correct one, and find the genes
	len = length(hierTOM$height)
	tree1 = which(hierTOM$merge==(-g1))%%len
	continue=length(which(hierTOM$merge==tree1))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree1[length(tree1)])%%len
		tree1 = c(tree1,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
	branchIndex = which(hierTOM$height==.minTreeHeight(hierTOM,g1,g2))
	branch=hierTOM$merge[branchIndex,]
	b1 <- NULL
	if(is.element(branch[1],tree1)){
		b1 = .getBranchMembers(hierTOM,branch[1],b1)
	} else b1 = .getBranchMembers(hierTOM,branch[2],b1)
	collectGarbage()
	return(b1)
}

# ----------------------------------------------------------------------------- #

reflectBranch <- function (hierTOM, g1, g2, both=FALSE){
## This function reverses the ordering of all genes in a branch of the
##  clustering tree defined by the minimal branch possible that contains
##  both g1 and g2 (as either ORDERED index or gene names), or just by
##  the genes in g1
	
	b1 = selectBranch(hierTOM, g1, g2)
	if (both) b1 = c(b1,selectBranch(hierTOM, g2, g1))
	
# Now reorder the hierTOM correctly
	ord = hierTOM$order
	i1 = which(ord%in%b1)
	b=1:(min(i1)-1); 
	if(b[length(b)]<b[1]) b = NULL
	e=(max(i1)+1):length(ord); 
	if((max(i1)+1)>length(ord)) e = NULL
	ord = ord[c(b,i1[order(i1,decreasing=T)],e)]
	hierTOM$order = ord
	return(hierTOM)
}

# ----------------------------------------------------------------------------- #



#' Select, swap, or reflect branches in a dendrogram.
#' 
#' swapTwoBranches takes the a gene tree object and two genes as input, and
#' swaps the branches containing these two genes at the nearest branch point of
#' the dendrogram.
#' 
#' reflectBranch takes the a gene tree object and two genes as input, and
#' reflects the branch containing the first gene at the nearest branch point of
#' the dendrogram.
#' 
#' selectBranch takes the a gene tree object and two genes as input, and
#' outputs indices for all genes in the branch containing the first gene, up to
#' the nearest branch point of the dendrogram.
#' 
#' 
#' @aliases swapTwoBranches reflectBranch selectBranch
#' @param hierTOM A hierarchical clustering object (or gene tree) that is used
#' to plot the dendrogram.  For example, the output object from the function
#' hclust or fastcluster::hclust.  Note that elements of hierTOM$order MUST be
#' named (for example, with the corresponding gene name).
#' @param g1 Any gene in the branch of interest.
#' @param g2 Any gene in a branch directly adjacent to the branch of interest.
#' @param both Logical: should the selection include the branch gene \code{g2}?
#' @return swapTwoBranches and reflectBranch return a hierarchical clustering
#' object with the hierTOM$order variable properly adjusted, but all other
#' variables identical as the heirTOM input.
#' 
#' selectBranch returns a numeric vector corresponding to all genes in the
#' requested branch.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#' 
#' ## Example: first simulate some data.
#' 
#' MEturquoise = sample(1:100,50)
#' MEblue      = c(MEturquoise[1:25], sample(1:100,25))
#' MEbrown     = sample(1:100,50)
#' MEyellow    = sample(1:100,50) 
#' MEgreen     = c(MEyellow[1:30], sample(1:100,20))
#' MEred	    = c(MEbrown [1:20], sample(1:100,30))
#' ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred)
#' dat1   = simulateDatExpr(ME,400 ,c(0.16,0.12,0.11,0.10,0.10,0.09,0.15), 
#'                          signed=TRUE)
#' TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
#' colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
#' tree1  = fastcluster::hclust(as.dist(1-TOM1),method="average")
#' colorh = labels2colors(dat1$allLabels)
#' 
#' plotDendroAndColors(tree1,colorh,dendroLabels=FALSE)
#' 
#' ## Reassign modules using the selectBranch and chooseOneHubInEachModule functions
#' 
#' datExpr = dat1$datExpr
#' hubs    = chooseOneHubInEachModule(datExpr, colorh)
#' colorh2 = rep("grey", length(colorh))
#' colorh2 [selectBranch(tree1,hubs["blue"],hubs["turquoise"])] = "blue"
#' colorh2 [selectBranch(tree1,hubs["turquoise"],hubs["blue"])] = "turquoise"
#' colorh2 [selectBranch(tree1,hubs["green"],hubs["yellow"])]   = "green"
#' colorh2 [selectBranch(tree1,hubs["yellow"],hubs["green"])]   = "yellow"
#' colorh2 [selectBranch(tree1,hubs["red"],hubs["brown"])]      = "red"
#' colorh2 [selectBranch(tree1,hubs["brown"],hubs["red"])]      = "brown"
#' plotDendroAndColors(tree1,cbind(colorh,colorh2),c("Old","New"),dendroLabels=FALSE)
#' 
#' ## Now swap and reflect some branches, then optimize the order of the branches
#' 
#' # Open a suitably sized graphics window
#' 
#' sizeGrWindow(12,9);
#' 
#' # partition the screen for 3 dendrogram + module color plots
#' 
#' layout(matrix(c(1:6), 6, 1), heights = c(0.8, 0.2, 0.8, 0.2, 0.8, 0.2));
#' 
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Starting Dendrogram", 
#'                     setLayout = FALSE)
#' 
#' tree1 = swapTwoBranches(tree1,hubs["red"],hubs["turquoise"])
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Swap blue/turquoise and red/brown", 
#'                     setLayout = FALSE)
#' 
#' tree1 = reflectBranch(tree1,hubs["blue"],hubs["green"])
#' plotDendroAndColors(tree1,colorh2,dendroLabels=FALSE,main="Reflect turquoise/blue", 
#'                     setLayout = FALSE)
#' 
#' 
swapTwoBranches <- function (hierTOM, g1, g2){
## This function re-arranges two branches in a heirarchical clustering tree
##  at the nearest branch point of two given genes (or indices) 
	
# Convert genes to indices (ORDERED AS ON THE PLOT)
	if(is.numeric(g1)) g1 = hierTOM$order[g1]
	if(is.numeric(g2)) g2 = hierTOM$order[g2]
	if(!is.numeric(g1)) g1 = which(hierTOM$labels==g1)
	if(!is.numeric(g2)) g2 = which(hierTOM$labels==g2)
	if((length(g1)==0)|(length(g2)==0)|(max(c(g1,g2))>length(hierTOM$labels))){
		write("Input genes are not both legal indices","")
		return(hierTOM)
	}
	
# Now determine the genes in each branch
	branchIndex = which(hierTOM$height==.minTreeHeight(hierTOM,g1,g2))
	b1 <- b2 <- NULL
	b1 = .getBranchMembers(hierTOM,hierTOM$merge[branchIndex,1],b1)
	b2 = .getBranchMembers(hierTOM,hierTOM$merge[branchIndex,2],b2)
	
# Now reorder the hierTOM correctly
	ord = hierTOM$order
	i1 = which(ord%in%b1)
	i2 = which(ord%in%b2)
	if(min(i1)>min(i2)) {tmp = i1; i1=i2; i2=tmp; rm(tmp)}
	b=1:(min(i1)-1); 
	if(b[length(b)]<b[1]) b = NULL
	e=(max(i2)+1):length(ord); 
	if((max(i2)+1)>length(ord)) e = NULL
	ord = ord[c(b,i2,i1,e)]
	hierTOM$order = ord
	return(hierTOM)
}

# ----------------------------------------------------------------------------- #



#' Chooses a single hub gene in each module
#' 
#' chooseOneHubInEachModule returns one gene in each module with high
#' connectivity, given a number of randomly selected genes to test.
#' 
#' 
#' @param datExpr Gene expression data with rows as samples and columns as
#' genes.
#' @param colorh The module assignments (color vectors) corresponding to the
#' rows in datExpr.
#' @param numGenes Th number of random genes to select per module.  Higher
#' number of genes increases the accuracy of hub selection but slows down the
#' function.
#' @param omitColors All colors in this character vector (default is "grey")
#' are ignored by this function.
#' @param power Power to use for the adjacency network (default = 2).
#' @param type What type of network is being entered.  Common choices are
#' "signed" (default) and "unsigned".  With "signed" negative correlations
#' count against, whereas with "unsigned" negative correlations are treated
#' identically as positive correlations.
#' @param \dots Any other parameters accepted by the *adjacency* function
#' @return Both functions output a character vector of genes, where the genes
#' are the hub gene picked for each module, and the names correspond to the
#' module in which each gene is a hub.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#' 
#' ## Example: first simulate some data.
#' 
#' MEturquoise = sample(1:100,50)
#' MEblue      = sample(1:100,50)
#' MEbrown     = sample(1:100,50)
#' MEyellow    = sample(1:100,50) 
#' MEgreen     = c(MEyellow[1:30], sample(1:100,20))
#' MEred	    = c(MEbrown [1:20], sample(1:100,30))
#' MEblack	    = c(MEblue  [1:25], sample(1:100,25))
#' ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred, MEblack)
#' dat1   = simulateDatExpr(ME,300,c(0.2,0.1,0.08,0.051,0.05,0.042,0.041,0.3), 
#'                          signed=TRUE)
#' TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
#' colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
#' tree1 <- tree2 <- fastcluster::hclust(as.dist(1-TOM1),method="average")
#' colorh = labels2colors(dat1$allLabels)
#' hubs    = chooseOneHubInEachModule(dat1$datExpr, colorh)
#' hubs
#' 
#' 
chooseOneHubInEachModule <- function(datExpr, colorh, numGenes=100, 
omitColors="grey", power=2, type="signed",...){
## This function returns the gene in each module with the highest connectivity, given
#   a number of randomly selected genes to test.
	
	numGenes = max(round(numGenes),2)
	keep     = NULL
	isIndex  = FALSE
	modules  = names(table(colorh))
	numCols  = table(colorh)
	if(!(is.na(omitColors)[1]))  modules = modules[!is.element(modules,omitColors)]
	if(is.null(colnames(datExpr))){
		colnames(datExpr) = 1:dim(datExpr)[2]
		isIndex = TRUE
	}
	
	for (m in modules){
		num   = min(numGenes,numCols[m])
		inMod = which(is.element(colorh,m)) 
		keep  = c(keep, sample(inMod,num))
	}
	colorh  = colorh[keep]
	datExpr = datExpr[,keep]
	return(chooseTopHubInEachModule(datExpr, colorh, omitColors, power, type,...))
}

# ----------------------------------------------------------------------------- #



#' Chooses the top hub gene in each module
#' 
#' chooseTopHubInEachModule returns the gene in each module with the highest
#' connectivity, looking at all genes in the expression file.
#' 
#' 
#' @param datExpr Gene expression data with rows as samples and columns as
#' genes.
#' @param colorh The module assignments (color vectors) corresponding to the
#' rows in datExpr.
#' @param omitColors All colors in this character vector (default is "grey")
#' are ignored by this function.
#' @param power Power to use for the adjacency network (default = 2).
#' @param type What type of network is being entered.  Common choices are
#' "signed" (default) and "unsigned".  With "signed" negative correlations
#' count against, whereas with "unsigned" negative correlations are treated
#' identically as positive correlations.
#' @param \dots Any other parameters accepted by the *adjacency* function
#' @return Both functions output a character vector of genes, where the genes
#' are the hub gene picked for each module, and the names correspond to the
#' module in which each gene is a hub.
#' @author Jeremy Miller
#' @keywords misc
#' @examples
#' 
#' ## Example: first simulate some data.
#' 
#' MEturquoise = sample(1:100,50)
#' MEblue      = sample(1:100,50)
#' MEbrown     = sample(1:100,50)
#' MEyellow    = sample(1:100,50) 
#' MEgreen     = c(MEyellow[1:30], sample(1:100,20))
#' MEred	    = c(MEbrown [1:20], sample(1:100,30))
#' MEblack	    = c(MEblue  [1:25], sample(1:100,25))
#' ME     = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen, MEred, MEblack)
#' dat1   = simulateDatExpr(ME,300,c(0.2,0.1,0.08,0.051,0.05,0.042,0.041,0.3), signed=TRUE)
#' TOM1   = TOMsimilarityFromExpr(dat1$datExpr, networkType="signed")
#' colnames(TOM1) <- rownames(TOM1) <- colnames(dat1$datExpr)
#' tree1 <- tree2 <- fastcluster::hclust(as.dist(1-TOM1),method="average")
#' colorh = labels2colors(dat1$allLabels)
#' hubs    = chooseTopHubInEachModule(dat1$datExpr, colorh)
#' hubs
#' 
#' 
chooseTopHubInEachModule <- function(datExpr, colorh, omitColors="grey", 
power=2, type="signed",...){
## This function returns the gene in each module with the highest connectivity.
	
	isIndex = FALSE
	modules = names(table(colorh))
	if(!(is.na(omitColors)[1]))  modules = modules[!is.element(modules,omitColors)]
	if(is.null(colnames(datExpr))){
		colnames(datExpr) = 1:dim(datExpr)[2]
		isIndex = TRUE
	}
	
	hubs = rep(NA,length(modules))
	names(hubs) = modules
	for (m in modules){
		adj = adjacency(datExpr[,colorh==m],power=power,type=type,...)
		hub = which.max(rowSums(adj))
		hubs[m] = colnames(adj)[hub]
	}
	if (isIndex){
		hubs = as.numeric(hubs)
		names(hubs) = modules
	}
	return(hubs)
}

#################################################################################
# Internal functions.............................................................

options(expressions=50000) # Required for .getBranchMembers

.getBranchMembers <- function(hierTOM, ind, members){
# This is a recursive function that gets all the indices of members of
#  a branch in an hClust tree.
	if(ind<0) return(c(members,-ind))
	m1 = hierTOM$merge[ind,1]
	m2 = hierTOM$merge[ind,2]
	if (m1>0) {
		members = .getBranchMembers(hierTOM,m1,members)
	} else members = c(members,-m1)
	if (m2>0) {
		members = .getBranchMembers(hierTOM,m2,members)
	} else members = c(members,-m2)
	return(members)
}

# ----------------------------------------------------------------------------- #

.minTreeHeight <- function(hierTOM,l1,l2) {
## This function finds the minimum height at which two leafs
##  in a hierarchical clustering tree are connected.  l1 and
##  l2 are the UNORDERED indices for the two leafs.
	
## Return 2 (larger than 1, if l1 or l2 is negative). This represents 
##  positions that are off the edge of the tree.
	if((l1<0)|(l2<0)) return(2)
	
## Get the tree for l1
	len = length(hierTOM$height)
	tree1 = which(hierTOM$merge==(-l1))%%len
	continue=length(which(hierTOM$merge==tree1))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree1[length(tree1)])%%len
		tree1 = c(tree1,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
## Get the tree for l2
	tree2 = which(hierTOM$merge==(-l2))%%len
	continue=length(which(hierTOM$merge==tree2))>0
	while(continue){
		nextInd = which(hierTOM$merge==tree2[length(tree2)])%%len
		tree2 = c(tree2,nextInd)
		continue=length(which(hierTOM$merge==nextInd))>0
	}
	
## Now find the index where the two trees first agree
	minTreeLen = min(c(length(tree1),length(tree2)))
	tree1 = tree1[(length(tree1)-minTreeLen+1):length(tree1)]
	tree2 = tree2[(length(tree2)-minTreeLen+1):length(tree2)]
	treeInd = tree1[min(which(tree1==tree2))]
	
## Now find and return the minimum tree height
	return(hierTOM$height[ifelse(treeInd==0,len,treeInd)])
}

# ----------------------------------------------------------------------------- #

.offDiagonalMatrixSum <- function(adj){
    len = dim(adj)[1]	
	output=sum(diag(adj[1:(len-1),2:len]))
	return(output)	
}
