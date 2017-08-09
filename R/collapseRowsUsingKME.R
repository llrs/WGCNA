# This function chooses a single plobe per gene based on a kME table



#' Selects one representative row per group based on kME
#' 
#' This function selects only the most informative probe for each gene in a kME
#' table, only keeping the probe which has the highest kME with respect to any
#' module in the module membership matrix.  This function is a special case of
#' the function collapseRows.
#' 
#' 
#' @param MM A module membership (kME) table with at least a subset of the
#' columns corresponding to kME values.
#' @param Gin Genes labels in a 1 to 1 correspondence with the rows of MM.
#' @param Pin If NULL (default), rownames of MM are assumed to be probe IDs. If
#' entered, Pin must be the same length as Gin and correspond to probe IDs for
#' MM.
#' @param kMEcols A numeric vector showing which columns in MM correspond to
#' kME values.  The default is all of them.
#' @return A list containing the following items: \item{MMcollapsed}{ A numeric
#' matrix with the same columns as the input matrix MM, but with rows
#' corresponding to the genes rather than the probes. } \item{group2row}{ A
#' matrix whose rows correspond to the unique gene labels and whose 2 columns
#' report which gene label (first column called group) is represented by what
#' probe (second column called selectedRowID) } \item{selectedRow}{ A logical
#' vector whose components are TRUE for probes selected as representatives and
#' FALSE otherwise. It has the same length as the vector Pin. }
#' @author Jeremy Miller
#' @seealso \code{\link{collapseRows}}
#' @keywords misc
#' @examples
#' 
#' # Example: first simulate some data
#' set.seed(100)
#' ME.A <- sample(1:100, 50)
#' ME.B <- sample(1:100, 50)
#' ME.C <- sample(1:100, 50)
#' ME.D <- sample(1:100, 50)
#' ME1 <- data.frame(ME.A, ME.B, ME.C, ME.D)
#' simDatA <- simulateDatExpr(ME1, 1000, c(0.2, 0.1, 0.08, 0.05, 0.3),
#' signed = TRUE, verbose = 0)
#' simDatB <- simulateDatExpr(ME1, 1000, c(0.2, 0.1, 0.08, 0.05, 0.3),
#' signed = TRUE, verbose = 0)
#' Gin <- c(colnames(simDatA$datExpr), colnames(simDatB$datExpr))
#' Pin <- paste("Probe", 1:length(Gin), sep=".")
#' datExpr <- cbind(simDatA$datExpr, simDatB$datExpr)
#' MM <- cor(datExpr, ME1)
#' 
#' # Now run the function and see some example output
#' results <- collapseRowsUsingKME(MM, Gin, Pin)
#' head(results$MMcollapsed)
#' head(results$group2Row)
#' head(results$selectedRow)
#' 
collapseRowsUsingKME <- function (MM, Gin, Pin=NULL, kMEcols = 1:dim(MM)[2]) {

	if (is.null(Pin))  {
	    Pin <- rownames(MM)
	}
	rownames(MM) <- Pin
	Gout <- as.character(sort(unique(Gin)))
	cors <- MM[, kMEcols]
	maxC <- apply(cors, 1, max)

	MMout <- matrix(0, nrow = length(Gout), ncol = dim(MM)[2])
	colnames(MMout) <- colnames(MM)
	rownames(MMout) <- Gout
	MM <- as.matrix(MM)

	keepThese <- NULL
	for (g in 1:length(Gout)){
		maxCg <- maxC
		maxCg[Gin != Gout[g]] <- -1000
		keep <- which(maxCg == max(maxCg))[1]
		MMout[g, ] <- MM[keep, ]
		keepThese <- c(keepThese, keep)
	}
	group2Row <- cbind(Gout, Pin[keepThese])
	colnames(group2Row) <- c("group", "selectedRowID")
	selectedRow <- is.element(1:length(Pin), keepThese)
	out <- list(MMout, group2Row, selectedRow)
	names(out) <- c("MMcollapsed", "group2Row", "selectedRow")
	return(out)
}
