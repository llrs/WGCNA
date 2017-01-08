library("WGCNA")
context("Testing errors")

test_that("stdErr", {
    x <- c(0.2, 0.5, 0.8)
    test <- stdErr(x)
    expect_lt(test, 0.24494901)
    expect_length(test, 1L)
})
