library("WGCNA")
context("Testing coClustering")

test_that("coClustering", {
    set.seed(1)
    nModules <- 10
    nGenes <- 1000
    cl1 <- sample(c(1:nModules), nGenes, replace = TRUE)
    cl2 <- sample(c(1:nModules), nGenes, replace = TRUE)
    coc <- coClustering(cl1, cl2)
    coc3 <- coClustering(cl1, cl2, tupletSize = 3)
    expect_gt(coc[1], coc3[1])
    ref <- coClustering(cl1, cl1)
    expect_length(ref, 10L)
    expect_true(all(ref == 1))
})
