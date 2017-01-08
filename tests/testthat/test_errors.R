library("WGCNA")
context("Testing errors")

test_that("stdErr", {
    x <- c(0.2, 0.5, 0.8)
    test <- stdErr(x)
    expect_equal(test, 0.1732051)
    expect_length(test, 1L)
})
