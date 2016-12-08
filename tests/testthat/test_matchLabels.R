library("WGCNA")
context("Testing matchLabels")

test_that("overlapTable", {
    set.seed(1)
    nModules <- 10
    nGenes <- 1000
    labels1 <- sample(c(1:nModules), nGenes, replace = TRUE)
    labels2 <- sample(c(1:nModules), nGenes, replace = TRUE)

    test <- overlapTable(labels1, labels2)
    expect_length(test, 2)
    expect_equal(ncol(test[[1]]), 10)

    test2 <- overlapTable(c(labels1, NA), c(labels2, NA), na.rm = TRUE)
    expect_true(all(test2$countTable == test$countTable))

    test3 <- overlapTable(labels1, labels2, ignore = 3)
    expect_equal(ncol(test3[[1]]), 9)
    expect_equal(ncol(test3[[1]]), ncol(test3[[2]]))
})

test_that("matchLabels", {
    set.seed(1)
    nModules <- 10
    nGenes <- 1000
    labels1 <- sample(c(1:nModules), nGenes, replace = TRUE)
    labels2 <- sample(c(1:nModules), nGenes, replace = TRUE)
    test <- matchLabels(labels1, labels2)

    expect_length(test, 1000)

    test2 <- matchLabels(labels1, labels2, ignoreLabels = 2)
    expect_length(test2, 1000)

    test3 <- matchLabels(labels2, labels1)
    expect_false(all(test == test3))
})
