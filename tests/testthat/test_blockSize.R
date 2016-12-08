library("WGCNA")
context("Testing blockSize")

test_that("blockSize", {
    block1 <- blockSize(30000, rectangularBlocks = TRUE,
                        maxMemoryAllocation = 2^31)
    block2 <- blockSize(30000, rectangularBlocks = TRUE)
    block3 <- blockSize(30000, rectangularBlocks = FALSE)
    block4 <- blockSize(30000, rectangularBlocks = TRUE,
                        maxMemoryAllocation = 2^31, overheadFactor = 3)
    expect_gt(block2, block1)
    expect_equal(block1, block4)
    expect_error(blockSize())
})
