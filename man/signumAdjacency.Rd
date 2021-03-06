% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjacency.R
\name{signumAdjacency}
\alias{signumAdjacency}
\title{Hard-thresholding adjacency function}
\usage{
signumAdjacency(corMat, threshold)
}
\arguments{
\item{corMat}{a matrix of correlations or other measures of similarity.}

\item{threshold}{threshold for connecting nodes: all nodes whose
\code{corMat} is above the threshold will be connected in the resulting
network.}
}
\value{
An unweighted adjacency matrix of the same dimensions as the input
\code{corMat}.
}
\description{
This function transforms correlations or other measures of similarity into
an unweighted network adjacency.
}
\examples{

datExpr <- matrix(rnorm(150), ncol = 5)
corMat <- cor(datExpr)
signumAdjacency(corMat, 0.1)

}
\references{
Bin Zhang and Steve Horvath (2005) "A General Framework for
Weighted Gene Co-Expression Network Analysis", Statistical Applications in
Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
}
\seealso{
\code{\link{adjacency}} for soft-thresholding and creating weighted
networks.
}
\author{
Steve Horvath
}
\keyword{misc}
