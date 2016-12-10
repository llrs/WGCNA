library("WGCNA")
context("Testing transposeBigData")

test_that("transposeBigData", {
    x <- data.frame(matrix(1:10000, nrow = 4, ncol = 2500))
    colnames(x) <- paste0("Y", 1:2500)
    xTranspose <- transposeBigData(x)
    expect_true(all(x[1:4, 1:4] == t(xTranspose[1:4, 1:4])))
    expect_true(all(colnames(x) == rownames(xTranspose)))
    expect_true(all(rownames(x) == colnames(xTranspose)))

    x <- matrix(1:10000, nrow = 4, ncol = 2500)
    colnames(x) <- paste0("Y", 1:2500)
    xTranspose <- transposeBigData(x)
    expect_true(all(x[1:4, 1:4] == t(xTranspose[1:4, 1:4])))
    expect_true(all(colnames(x) == rownames(xTranspose)))
    expect_true(all(rownames(x) == colnames(xTranspose)))

    expect_error(transposeBigData(x, 1))

    x <- matrix(1:10000, nrow = 2500, ncol = 4)
    colnames(x) <- paste0("Y", 1:2500)
    xTranspose <- transposeBigData(x)
    expect_true(all(x[1:4, 1:4] == t(xTranspose[1:4, 1:4])))
    expect_true(all(colnames(x) == rownames(xTranspose)))
    expect_true(all(rownames(x) == colnames(xTranspose)))

    expect_error(transposeBigData("a"))
})
