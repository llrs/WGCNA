% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Functions-multiData.R
\name{nSets}
\alias{nSets}
\title{Number of sets in a multi-set variable}
\usage{
nSets(multiSet, ...)
}
\arguments{
\item{multiSet}{vector of lists; in each list there must be a component
named \code{data} whose content is a matrix or dataframe or array of
dimension 2.}

\item{\dots}{Other arguments to function \code{\link{checkSets}}.}
}
\value{
A single integer that equals the number of sets given in the input
\code{multiSet}.
}
\description{
A convenience function that returns the number of sets in a multi-set
variable.
}
\seealso{
\code{\link{checkSets}}
}
\author{
Peter Langfelder
}
\keyword{misc}
