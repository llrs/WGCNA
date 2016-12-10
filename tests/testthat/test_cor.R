library("WGCNA")
context("Testing corFunctions")

test_that("cor", {
    set.seed(1)
    x <- matrix(rnorm(100), 10)
    cor.p <- cor(x, method = "pearson")
    cor.k <- cor(x, method = "kendall")
    cor.s <- cor(x, method = "spearman")

    expect_equal(ncol(cor.p), 10L)
    expect_equal(ncol(cor.k), 10L)
    expect_equal(ncol(cor.s), 10L)

    expect_equal(nrow(cor.p), 10L)
    expect_equal(nrow(cor.k), 10L)
    expect_equal(nrow(cor.s), 10L)

    expect_true(checkSimilarity(cor.p))
    expect_true(checkSimilarity(cor.k))
    expect_true(checkSimilarity(cor.s))

    cors.p <- cor(x, method = "pearson", cosine = TRUE)
    cors.k <- cor(x, method = "kendall", cosine = TRUE)
    cors.s <- cor(x, method = "spearman", cosine = TRUE)

    expect_equal(ncol(cors.p), 10L)
    expect_equal(ncol(cors.k), 10L)
    expect_equal(ncol(cors.s), 10L)

    expect_equal(nrow(cors.p), 10L)
    expect_equal(nrow(cors.k), 10L)
    expect_equal(nrow(cors.s), 10L)

    expect_true(checkSimilarity(cors.p))
    expect_true(checkSimilarity(cors.k))
    expect_true(checkSimilarity(cors.s))
})


test_that("bicor", {
    set.seed(1)
    x <- matrix(rnorm(100), 10)
    cor.r <- bicor(x)
    expect_true(checkSimilarity(cor.r))
    expect_equal(ncol(cor.r), 10L)
    expect_equal(nrow(cor.r), 10L)
})

test_that("cor1", {
    set.seed(1)
    x <- matrix(rnorm(100), 10)
    cor.r <- cor1(x)
    expect_true(checkSimilarity(cor.r))
    expect_equal(ncol(cor.r), 10L)
    expect_equal(nrow(cor.r), 10L)
})

test_that("corFast", {
    set.seed(1)
    x <- matrix(rnorm(100), 10)
    cor.r <- corFast(x, quick = 1)
    expect_true(checkSimilarity(cor.r))
    expect_equal(ncol(cor.r), 10L)
    expect_equal(nrow(cor.r), 10L)
})
