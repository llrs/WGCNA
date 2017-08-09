#' Prediction of Weighted Mutual Information Adjacency Matrix by Correlation
#' 
#' AFcorMI computes a predicted weighted mutual information adjacency matrix
#' from a given correlation matrix.
#' 
#' This function is a one-to-one prediction when we consider correlation as
#' unsigned. The prediction corresponds to the
#' \code{AdjacencyUniversalVersion2} discussed in the help file for the
#' function \code{\link{mutualInfoAdjacency}}. For more information about the
#' generation and features of the predicted mutual information adjacency,
#' please refer to the function \code{\link{mutualInfoAdjacency}}.
#' 
#' @param r a symmetric correlation matrix with values from -1 to 1.
#' @param m number of observations from which the correlation was calcuated.
#' @return A matrix with the same size as the input correlation matrix,
#' containing the predicted mutual information of type
#' \code{AdjacencyUniversalVersion2}.
#' @author Steve Horvath, Lin Song, Peter Langfelder
#' @seealso \code{\link{mutualInfoAdjacency}}
#' @keywords misc
#' @examples
#' 
#' #Simulate a data frame datE which contains 5 columns and 50 observations
#' m <- 50
#' x1 <- rnorm(m)
#' r <- 0.5
#' x2 <- r*x1+sqrt(1 - r ^ 2) * rnorm(m)
#' r <- 0.3
#' x3 <- r * (x1 - 0.5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
#' x4 <- rnorm(m)
#' r <- 0.3
#' x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
#' datE <- data.frame(x1, x2, x3, x4, x5)
#' #calculate predicted AUV2
#' cor.data <- cor(datE, use = "p")
#' \dontrun{
#' AUV2 <- AFcorMI(r = cor.data, m = nrow(datE))
#' }
#' 
AFcorMI <- function(r, m) {
    checkAdjMat(r, min = -1 , max = 1)
    if (!is.numeric(m)) {
        stop("Pease introduce a number with the number of observations used")
    }
	D <- 0.43 * m ^ ( -0.3)
	epsilon <- D ^ 2.2
	out <- log(1 + epsilon - r ^ 2) / log(epsilon) * (1 - D) + D
	checkAdjMat(out)
	out
}
