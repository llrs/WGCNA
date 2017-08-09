# Makes a consensus network using all of the default values in the WGCNA library.



#' Consensus clustering based on topological overlap and hierarchical
#' clustering
#' 
#' This function makes a consensus network using all of the default values in
#' the WGCNA library.  Details regarding how consensus modules are formed can
#' be found here:
#' http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/Tutorials/Consensus-NetworkConstruction-man.pdf
#' 
#' 
#' @param multiExpr Expression data in the multi-set format (see checkSets). A
#' vector of lists, one per set. Each set must contain a component data that
#' contains the expression data.  Rows correspond to samples and columns to
#' genes or probes. Two or more sets of data must be included and adjacencies
#' cannot be used.
#' @param softPower Soft thresholding power used to make each of the networks
#' in multiExpr.
#' @param TOM A LIST of matrices holding the topological overlap corresponding
#' to the sets in multiExpr, if they have already been calculated. Otherwise,
#' keep TOM set as NULL (default), and TOM similarities will be calculated
#' using the WGCNA defaults.  If inputted, this variable must be a list with
#' each entree a TOM corresponding to the same entries in multiExpr.
#' @return \item{consensusTOM}{ The TOM difference matrix (1-TOM similarity)
#' corresponding to the consensus network. } \item{consTree}{ Returned value is
#' the same as that of hclust: An object of class hclust which describes the
#' tree produced by the clustering process.  This tree corresponds to the
#' dissimilarity matrix consensusTOM. }
#' @author Peter Langfelder, Steve Horvath, Jeremy Miller
#' @seealso \code{\link{blockwiseConsensusModules}}
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
#' @examples
#' 
#' 
#' # Example consensus network using two simulated data sets
#' 
#' set.seed <- 100
#' MEturquoise <- sample(1:100, 50)
#' MEblue <- sample(1:100, 50)
#' MEbrown <- sample(1:100, 50)
#' MEyellow <- sample(1:100, 50)
#' MEgreen <- sample(1:100, 50)
#' 
#' ME <- data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen)
#' dat1 <- simulateDatExpr(ME,300,c(0.2,  0.10,  0.10,  0.10,  0.10,  0.2),
#'   signed=TRUE)
#' dat2 <- simulateDatExpr(ME, 300, c(0.18, 0.11, 0.11, 0.09, 0.11, 0.23),
#'   signed = TRUE)
#' multiExpr <- list(S1=list(data=dat1$datExpr), S2=list(data=dat2$datExpr))
#' softPower <- 8
#' 
#' consensusNetwork <- consensusDissTOMandTree(multiExpr, softPower)
#' \dontrun{
#' plotDendroAndColors(consensusNetwork$consTree,
#'      cbind(labels2colors(dat1$allLabels),
#'      labels2colors(dat2$allLabels)), c("S1","S2"), dendroLabels=FALSE)
#' }
#' 
consensusDissTOMandTree <- function(multiExpr, softPower, TOM=NULL){
	nGenes <- ncol(multiExpr[[1]]$data)
	nSets <- length(multiExpr)
	if(is.null(TOM)) {
		adjacencies <- TOM <- list()
		for (set in 1:nSets) {
			adjacencies[[set]] <- adjacency(multiExpr[[set]]$data,
			                                power = softPower, type = "signed")
			diag(adjacencies[[set]]) <- 0
			message("Adjacency, set ", set)
			TOM[[set]] <- TOMsimilarity(adjacencies[[set]], TOMType = "signed")
			message("Similarity, set ", set)
			collectGarbage()
		}
	}
	nSets <- length(TOM)
	set.seed(12345)
	scaleP <- 0.95
	nSamples <- as.integer(1/(1-scaleP) * 1000)
	scaleSample <- sample(nGenes*(nGenes-1)/2, size = nSamples)
	TOMScalingSamples <- list()
	scaleQuant <- scalePowers <- rep(1, nSets)
	for (set in 1:nSets) {
		TOMScalingSamples[[set]] <- as.dist(TOM[[set]])[scaleSample]
		scaleQuant[set] <- quantile(TOMScalingSamples[[set]], probs = scaleP,
		                            type = 8)
		if (set > 1) {
			scalePowers[set] <- log(scaleQuant[1])/log(scaleQuant[set])
			TOM[[set]] <- TOM[[set]]^scalePowers[set]
		}
		message("Scaling, set ",set)
	}

	half <- round(nGenes/2)
	haP1 <- half + 1
	kp <- list(list(c(1:half), c(1:half)), list(c(1:half), c(haP1:nGenes)),
		  	    list(c(haP1:nGenes), c(1:half)), list(
		  	        c(haP1:nGenes), c(haP1:nGenes)))
	consensusTOMi <- list()
	for (i in 1:4) {
		a <- kp[[i]][[1]]
		b <- kp[[i]][[2]]
		consensusTOMi[[i]] <- TOM[[1]][a,b]
		for (j in 2:nSets) {
		    consensusTOMi[[i]] <- pmin(consensusTOMi[[i]], TOM[[j]][a,b])
		}
		message(i," of 4 iterations in pMin")
	}
	consensusTOM <- rbind(cbind(consensusTOMi[[1]],consensusTOMi[[2]]),
						 cbind(consensusTOMi[[3]],consensusTOMi[[4]]))
	rownames(consensusTOM) <- colnames(multiExpr[[1]]$data)
	colnames(consensusTOM) <- rownames(consensusTOM)
	consensusTOM <- 1-consensusTOM
	message("Starting dendrogram tree.")
	consTree <- fastcluster::hclust(as.dist(consensusTOM), method = "average")
	message("DONE!!!!")
	out <- list(consensusTOM,consTree)
	names(out) <- c("consensusTOM","consTree")
	return(out)
}

.collect_garbage <- function(){while (gc()[2,4] != gc()[2,4] | gc()[1,4] != gc()[1,4]){}}
