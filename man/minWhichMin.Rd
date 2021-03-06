% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/quantileC.R
\name{minWhichMin}
\alias{minWhichMin}
\title{Fast joint calculation of row- or column-wise minima and indices of minimum
elements}
\usage{
minWhichMin(x, byRow = FALSE, dims = 1)
}
\arguments{
\item{x}{A numeric matrix or array.}

\item{byRow}{Logical: should the minima and indices be found for columns
(\code{FALSE}) or rows (\code{TRUE})?}

\item{dims}{Specifies dimensions for which to find the minima and indices.
For \code{byRow = FALSE}, they are calculated for dimensions \code{dims+1}
to \code{n=length(dim(x))}; for For \code{byRow = TRUE}, they are calculated
for dimensions 1,...,\code{dims}.}
}
\value{
A list with two components, \code{min} and \code{which}; each is a
vector or array with dimensions

\code{dim(x)[(dims+1):n]} (with \code{n=length(dim(x))}) if \code{byRow =
FALSE}, and

\code{dim(x)[1:dims]} if \code{byRow = TRUE}.
}
\description{
Fast joint calculation of row- or column-wise minima and indices of minimum
elements. Missing data are removed.
}
\author{
Peter Langfelder
}
\keyword{stats}
