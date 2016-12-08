library("WGCNA")
context("Testing cutree")

test_that("cutree", {
    datExpr <- matrix(rnorm(150), ncol = 5)
    hc <- hclust(dist(datExpr))
    col <- cutreeStaticColor(hc, minSize = 2)
    num <- cutreeStatic(hc, minSize = 2)
    expect_true(length(col) == 30)
    expect_true(length(num) == 30)
})
