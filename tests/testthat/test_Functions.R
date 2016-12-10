library("WGCNA")
context("Testing Functions")

test_that("alignExpr", {
    set.seed(1)
    datExpr <- matrix(rnorm(10), 5)
    test <- alignExpr(datExpr)
    expect_true(is.data.frame(test))
    expect_equal(ncol(test), ncol(datExpr))
    expect_equal(nrow(test), nrow(datExpr))

    expect_error(alignExpr(datExpr, y = seq(by = 0.2)))
    y <- seq(from=0.2, to = 1, by = 0.2)
    test2 <- alignExpr(datExpr, y)
    expect_true(is.data.frame(test2))
    expect_equal(ncol(test2), ncol(datExpr))
    expect_equal(nrow(test2), nrow(datExpr))
})

test_that("vectorizeMatrix", {
    expect_error(vectorizeMatrix(rnorm(10)))
    M <- matrix(rnorm(10), ncol = 2)
    test <- vectorizeMatrix(M)
    expect_length(test, 10)
    test2 <- vectorizeMatrix(M, TRUE)
    expect_length(test2, 10)
})

test_that("metaZfunction", {
    expect_error(metaZfunction(c(0.5, 0.2)))
    M <- matrix(rnorm(10), 2)
    test <- metaZfunction(M)
    expect_length(test, 2)
    test2 <- metaZfunction(M, c(0.1, 0.2))
    expect_equal(length(test), length(test2))
    expect_true(is.numeric(test))
    expect_true(is.numeric(test2))
})

test_that("prepComma", {
    expect_equal(prepComma(""), "")
    expect_warning(prepComma(c("ha", "he")))
})

test_that("prependZeros", {
    expect_true(is.character(prependZeros(1)))
    expect_equal(nchar(prependZeros(1, 4)), 4L)
    expect_error(prependZeros(10, 1))
})
