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

test_that("coClustering.permutationTest", {
    set.seed(1)
    nModules <- 5
    nGenes <- 100
    cl1 <- sample(c(1:nModules), nGenes, replace = TRUE)
    cl2 <- sample(c(1:nModules), nGenes, replace = TRUE)
    cc <- coClustering(cl1, cl2)
    ccPerm <- coClustering.permutationTest(cl1, cl2, nPermutations = 20,
                                           verbose = 0)
    expect_length(ccPerm$observed, 5L)
    expect_length(ccPerm$Z, 5L)
    expect_length(ccPerm$permuted.mean, 5L)
    expect_length(ccPerm$permuted.sd, 5L)
    expect_equal(ncol(ccPerm$permuted.cc), 5L)
    expect_equal(nrow(ccPerm$permuted.cc), 20L)
})
