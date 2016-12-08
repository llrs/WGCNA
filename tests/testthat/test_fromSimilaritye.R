library("WGCNA")
context("Testing fromSimilarity")

test_that("matrixToNetwork", {
    mat <- matrix(rnorm(25), ncol = 5)
    adj <- matrixToNetwork(mat, "max")
    expect_equal(ncol(adj), 5L)
    expect_equal(nrow(adj), 5L)

    })

test_that("similarity", {
    mat <- matrix(rnorm(25), ncol = 5)
    adj <- matrixToNetwork(mat, "max")
    expect_true(checkSimilarity(adj))
})
