library("WGCNA")
context("Testing quantile optimizations in C")

test_that("colQuantileC", {
    set.seed(100)
    data <- rnorm(100)
    p <- 0.1
    test <- colQuantileC(data, p)
    expect_gt(test, -1.2)
    expect_equivalent(test, quantile(data, 0.1))
    })

test_that("rowQuantileC", {
    set.seed(100)
    data <- rnorm(100)
    p <- 0.1
    data <- matrix(data, ncol = 20)
    test <- rowQuantileC(data, p)
    expect_gt(test[1], -0.6)
    expect_equivalent(test, apply(data, 1, quantile, 0.1))
})
