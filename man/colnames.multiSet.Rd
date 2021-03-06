% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiData.R
\name{colnames<-}
\alias{colnames<-}
\alias{colnames<-.default}
\alias{colnames}
\alias{colnames.default}
\alias{colnames.multiSet}
\title{Extract the colnames of the multiSet}
\usage{
colnames(multiSet) <- value

\method{colnames}{default}(multiSet) <- value

colnames(multiSet)

\method{colnames}{default}(multiSet)

\method{colnames}{multiSet}(multiSet)
}
\arguments{
\item{multiSet}{multiSet object to retrieve the columns of the first data}

\item{\dots}{other arguments passed}
}
\value{
\code{colnames.multiSet} returns the vector of column names of the
\code{data} component. The function assumes the column names in all sets are
the same.
}
\description{
Extract the colnames of the multiSet
}
