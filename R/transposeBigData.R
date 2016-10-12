#' Transpose a big matrix or data frame
#'
#' This transpose command partitions a big matrix (or data frame) into blocks
#' and applies the t() function to each block separately.
#'
#' Assume you have a very large matrix with say 500k columns. In this case, the
#' standard transpose function of R \code{t()} can take a long time. Solution:
#' Split the original matrix into sub-matrices by dividing the columns into
#' blocks. Next apply \code{t()} to each sub-matrix. The same holds if the
#' large matrix contains a large number of rows.  The function
#' \code{transposeBigData} automatically checks whether the large matrix
#' contains more rows or more columns. If the number of columns is larger than
#' or equal to the number of rows then the block wise splitting will be applied
#' to columns otherwise to the rows.
#'
#' @param x a matrix or data frame
#' @param blocksize a positive integer larger than 1, which determines the
#' block size. Default is 20k.
#' @return A matrix or data frame (depending on the input \code{x}) which is
#' the transpose of \code{x}.
#' @note This function can be considered a wrapper of \code{\link{t}()}
#' @author Steve Horvath, UCLA
#' @seealso The standard function \code{\link{t}} .
#' @references Any linear algebra book will explain the transpose.
#' @keywords misc
#' @examples
#'
#' x=data.frame(matrix(1:10000,nrow=4,ncol=2500))
#' dimnames(x)[[2]]=paste("Y",1:2500,sep="")
#' xTranspose=transposeBigData(x)
#' x[1:4,1:4]
#' xTranspose[1:4,1:4]
#'
#'
transposeBigData = function (x, blocksize = 20000) {
    isdataframe = is.data.frame(x)
    ismatrix = is.matrix(x)
    if (!(isdataframe | ismatrix))  {
        stop("Input is neither a data frame nor a matrix")
    }
    if (blocksize < 2) {
        stop("This blocksize makes no sense. It should be a positive integer>1.")
    }
    nrow1 = nrow(x)
    ncol1 = ncol(x)
    xTranspose = matrix(NA, nrow = ncol1, ncol = nrow1)
    if (nrow1 <= ncol1) {
        no.blocks = as.integer(ncol1/blocksize)
        if (no.blocks >= 1) {
            for (i in 1:no.blocks) {
                blockIndex = (i - 1) * blocksize + 1:blocksize
                xTranspose[blockIndex, ] = t(x[, blockIndex])
            }
        }
        if (ncol1 - no.blocks * blocksize == 1) {
            xTranspose[ncol1, ] = t(x[, ncol1])
        }
        if (ncol1 - no.blocks * blocksize > 1) {
            finalindex = (no.blocks * blocksize + 1):ncol1
            xTranspose[finalindex, ] = t(x[, finalindex])
        }
    }
    if (nrow1 > ncol1) {
        no.blocks = as.integer(nrow1/blocksize)
        if (no.blocks >= 1) {
            for (i in 1:no.blocks) {
                blockIndex = (i - 1) * blocksize + 1:blocksize
                xTranspose[, blockIndex] = t(x[blockIndex, ])
            }
        }
        if (nrow1 - no.blocks * blocksize == 1) {
            xTranspose[, nrow1] = t(x[nrow1, ])
        }
        if (nrow1 - no.blocks * blocksize > 1) {
            finalindex = (no.blocks * blocksize + 1):nrow1
            xTranspose[, finalindex] = t(x[finalindex, ])
        }
    }
    if (isdataframe) {
        xTranspose = data.frame(xTranspose)
        dimnames(xTranspose)[[1]] = dimnames(x)[[2]]
        dimnames(xTranspose)[[2]] = dimnames(x)[[1]]
    }
    if (ismatrix) {
        dimnames(xTranspose)[[1]] = dimnames(x)[[2]]
        dimnames(xTranspose)[[2]] = dimnames(x)[[1]]
    }
    xTranspose
}
