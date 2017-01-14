library("WGCNA")
context("Testing pquantile")


a <- matrix(c(1:12), 3, 4)
b <- a + 1
c <- a + 2
colnames(a) <- paste0("col_", c(1:4))

test_that("pquantile", {
    test1 <- pquantile(prob = 0.5, a, b, c)
    expect_true(test1[1, 1] == 2L)
})

test_that("pmean", {
    test2 <- pmean(a, b, c)
    expect_true(test1[1, 1] == 2L)
})

test_that("pmedian", {
    test2 <- pmean(a, b, c)
    expect_true(test1[1, 1] == 2L)
})
