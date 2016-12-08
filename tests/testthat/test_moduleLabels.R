library("WGCNA")
context("Testing moduleLabels")

test_that("normalizeLabels", {
    labels <- c(0, 1, rep(2, 2))
    test1 <- normalizeLabels(labels)
    expect_equal(test1, c(0, 2, 1, 1))
    test2 <- normalizeLabels(labels, FALSE)
    expect_equal(test2, c(2, 3, 1, 1))
})

test_that("moduleNumber", {
    datExpr <- matrix(rnorm(150), ncol = 5)
    hc <- hclust(dist(datExpr))
    m1 <- moduleNumber(hc)
    expect_equal(length(m1), 30)
})

test_that("numbers2colors", {
    num <- 0:30
    col1 <- numbers2colors(num)
    col2 <- numbers2colors(num, TRUE)
    expect_true(ncol(col1) == 1)
    expect_true(nrow(col1) == 31)
})
