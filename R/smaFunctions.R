
###########################################################################
# Statistics for Microarray Analysis for R
# Discriminant analysis
#
# Date : August 21, 2000
# Last update : April 13, 2001
#
# Authors: Sandrine Dudoit, Yee Hwa (Jean) Yang, and Jane Fridlyand.
##########################################################################

#' Red and Green Color Specification
#'
#' This function creates a vector of n ``contiguous'' colors, corresponding to
#' n intensities (between 0 and 1) of the red, green and blue primaries, with
#' the blue intensities set to zero. The values returned by
#' \code{rgcolors.func} can be used with a \code{col} specification in
#' graphics functions or in \code{\link{par}}.
#'
#'
#' @param n the number of colors (>= 1) to be used in the red and green
#' palette.
#' @return a character vector of color names. Colors are specified directly in
#' terms of their RGB components with a string of the form "\#RRGGBB", where
#' each of the pairs RR, GG, BB consist of two hexadecimal digits giving a
#' value in the range 00 to FF.
#' @author
#' Sandrine Dudoit, \email{sandrine@@stat.berkeley.edu} \cr Jane
#' Fridlyand, \email{janef@@stat.berkeley.edu}
#' @seealso
#' \code{\link{plotCor}}, \code{\link{plotMat}}, \code{\link{colors}},
#' \code{\link{rgb}}, \code{\link{image}}.
#' @keywords color
#' @examples
#'
#' rgcolors.func(n=5)
#' ## The following vector is returned:
#' ## "#00FF00" "#40BF00" "#808000" "#BF4000" "#FF0000"
#'
rgcolors.func<-function(n = 50) {
  k <- round(n/2)
  r <- c(rep(0, k), seq(0, 1, length = k))
  g <- c(rev(seq(0, 1, length = k)), rep(0, k))
  res <- rgb(r, g, rep(0, 2 * k))
  res
}

#' Red and Green Color Image of Correlation Matrix
#'
#' This function produces a red and green color image of a correlation matrix
#' using an RGB color specification. Increasingly positive correlations are
#' represented with reds of increasing intensity, and increasingly negative
#' correlations are represented with greens of increasing intensity.
#'
#'
#' @param x a matrix of numerical values.
#' @param new If \code{new=F}, \code{x} must already be a correlation matrix.
#' If \code{new=T}, the correlation matrix for the columns of \code{x} is
#' computed and displayed in the image.
#' @param nrgcols the number of colors (>= 1) to be used in the red and green
#' palette.
#' @param labels vector of character strings to be placed at the tickpoints,
#' labels for the columns of \code{x}.
#' @param labcols colors to be used for the labels of the columns of \code{x}.
#' \code{labcols} can have either length 1, in which case all the labels are
#' displayed using the same color, or the same length as \code{labels}, in
#' which case a color is specified for the label of each column of \code{x}.
#' @param title character string, overall title for the plot.
#' @param \dots graphical parameters may also be supplied as arguments to the
#' function (see \code{\link{par}}). For comparison purposes, it is good to set
#' \code{zlim=c(-1, 1)}.
#' @author Sandrine Dudoit, \email{sandrine@@stat.berkeley.edu}
#' @seealso \code{\link{plotMat}}, \code{\link{rgcolors.func}},
#' \code{\link{cor}}, \code{\link{image}}, \code{\link{rgb}}.
#' @keywords hplot
plotCor <- function(x, new=FALSE, nrgcols=50, labels=FALSE, labcols=1,
                    title="", ...) {
#   X <- x
   n<-ncol(x)
   corr<-x

   if(new)
     corr<-cor(x, use = 'p')

   image(1:n, 1:n, corr[, n:1], col=rgcolors.func(nrgcols), axes=FALSE, xlab="", ylab="", ... )

  if(length(labcols)==1){
    axis(2, at=n:1, labels=labels, las=2, cex.axis=0.6, col.axis=labcols)
    axis(3, at=1:n, labels=labels, las=2, cex.axis=0.6, col.axis=labcols)
  }

  if(length(labcols)==n){
    cols<-unique(labcols)
    for(i in 1:length(cols)){
      which<-(1:n)[labcols==cols[i]]
      axis(2, at=(n:1)[which], labels=labels[which], las=2, cex.axis=0.6, col.axis=cols[i])
      axis(3, at=which, labels=labels[which], las=2, cex.axis=0.6, col.axis=cols[i])
     }
  }

  mtext(title, side=3, line=3)
  box()
}

#' Red and Green Color Image of Data Matrix
#'
#' This function produces a red and green color image of a data matrix using an
#' RGB color specification. Larger entries are represented with reds of
#' increasing intensity, and smaller entries are represented with greens of
#' increasing intensity.
#'
#'
#' @param x a matrix of numbers.
#' @param nrgcols the number of colors (>= 1) to be used in the red and green
#' palette.
#' @param rlabels vector of character strings to be placed at the row
#' tickpoints, labels for the rows of \code{x}.
#' @param clabels vector of character strings to be placed at the column
#' tickpoints, labels for the columns of \code{x}.
#' @param rcols colors to be used for the labels of the rows of \code{x}.
#' \code{rcols} can have either length 1, in which case all the labels are
#' displayed using the same color, or the same length as \code{rlabels}, in
#' which case a color is specified for the label of each row of \code{x}.
#' @param ccols colors to be used for the labels of the columns of \code{x}.
#' \code{ccols} can have either length 1, in which case all the labels are
#' displayed using the same color, or the same length as \code{clabels}, in
#' which case a color is specified for the label of each column of \code{x}.
#' @param title character string, overall title for the plot.
#' @param \dots graphical parameters may also be supplied as arguments to the
#' function (see \code{\link{par}}).  E.g. \code{zlim=c(-3, 3)}
#' @author Sandrine Dudoit, \email{sandrine@@stat.berkeley.edu}
#' @seealso \code{\link{plotCor}}, \code{\link{rgcolors.func}},
#' \code{\link{cor}}, \code{\link{image}}, \code{\link{rgb}}.
#' @keywords hplot
plotMat <- function(x, nrgcols=50, rlabels=FALSE, clabels=FALSE, rcols=1,
                    ccols=1, title="", ...) {
#  X <-x
  n<-nrow(x)
  p<-ncol(x)

  image(1:p, 1:n, t(x[n:1, ]), col=rgcolors.func(nrgcols), axes=FALSE, xlab="",
        ylab="", ... )

  if(length(ccols)==1) {
    axis(3, at=1:p, labels=clabels, las=2, cex.axis=0.6, col.axis=ccols)
      }

  if(length(ccols)==p){
    cols<-unique(ccols)
    for(i in 1:length(cols)){
      which<-(1:p)[ccols==cols[i]]
      axis(3, at=which, labels=clabels[which], las=2, cex.axis=0.6,
           col.axis=cols[i])
     }
  }

  if(length(rcols)==1){
    axis(2, at=n:1, labels=rlabels, las=2, cex.axis=0.6, col.axis=rcols)
      }

  if(length(rcols)==n){
    cols<-unique(rcols)
    for(i in 1:length(cols)){
      which<-(1:n)[rcols==cols[i]]
      axis(2, at=(n:1)[which], labels=rlabels[which], las=2, cex.axis=0.6,
           col.axis=cols[i])
     }
  }

  mtext(title, side=3, line=3)
  box()
}

