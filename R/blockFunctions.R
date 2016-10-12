# Function to calculate an appropriate blocksize


#' Attempt to calculate an appropriate block size to maximize efficiency of
#' block-wise calcualtions.
#' 
#' The function uses a rather primitive way to estimate available memory and
#' use it to suggest a block size appropriate for the many block-by-block
#' calculations in this package.
#' 
#' Multiple functions within the WGCNA package use a divide-and-conquer (also
#' known as block-by-block, or block-wise) approach to handling large data
#' sets. This function is meant to assist in choosing a suitable block size,
#' given the size of the data and the available memory.
#' 
#' If the entire expected result fits into the allowed memory (after taking
#' into account the expected overhead), the returned block size will equal the
#' input \code{matrixSize}.
#' 
#' The internal estimation of available memory works by returning the size of
#' largest successfully allocated block of memory. It is hoped that this will
#' lead to reasonable results but some operating systems may actually allocate
#' more than is available. It is therefore preferable that the user specifies
#' the available memory by hand.
#' 
#' @param matrixSize the relevant dimension (usually the number of columns) of
#' the matrix that is to be operated on block-by-block.
#' @param rectangularBlocks logical indicating whether the bocks of data are
#' rectangular (of size \code{blockSize} times \code{matrixSize}) or square (of
#' size \code{blockSize} times \code{blockSize}).
#' @param maxMemoryAllocation maximum desired memory allocation, in bytes.
#' Should not exceed 2GB or total installed RAM (whichever is greater) on
#' 32-bit systems, while on 64-bit systems it should not exceed the total
#' installed RAM. If not supplied, the available memory will be estimated
#' internally.
#' @param overheadFactor overhead factor for the memory use by R. Recommended
#' values are between 2 (for simple calculations) and 4 or more for complicated
#' calculations where intermediate results (for which R must also allocate
#' memory) take up a lot of space.
#' @return A single integer giving the suggested block size, or
#' \code{matrixSize} if the entire calculation is expected to fit into memory
#' in one piece.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#' 
#' # Suitable blocks for handling 30,000 genes within 2GB (=2^31 bytes) of memory
#' blockSize(30000, rectangularBlocks = TRUE, maxMemoryAllocation = 2^31)
#' 
#' @export blockSize
blockSize = function(matrixSize, rectangularBlocks = TRUE,
                     maxMemoryAllocation = NULL, overheadFactor = 3) {
    if (is.null(maxMemoryAllocation)) {
        maxAlloc = .checkAvailableMemory()
    } else {
        maxAlloc = maxMemoryAllocation/8
    }
    maxAlloc = maxAlloc/overheadFactor

    if (rectangularBlocks){
        blockSz = floor(maxAlloc/matrixSize)
    } else
        blockSz = floor(sqrt(maxAlloc))

    return(min(matrixSize, blockSz))
}
