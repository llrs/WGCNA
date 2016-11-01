library("WGCNA")
context("Testing multiSet objects and functions")

test_that("fixDataStructure fixes data structure", {
    singleSetData <- matrix(1:4, 2, 2)
    testing <- fixDataStructure(singleSetData)
    result <- list(structure(list(data = structure(1:4, .Dim = c(2L, 2L))),
                             .Names = "data"))

    expect_equal(testing, result)
    expect_equal(length(testing), 1L)
    expect_match(names(testing[[1]]), "data")
    expect_equal(dim(testing[[1]]$data), c(2L, 2L))
    expect_equal(testing[[1]]$data, singleSetData)
})

test_that("multiSet works properly", {
    data1 <- matrix(rnorm(100), 20, 5)
    data2 <- matrix(rnorm(50), 10, 5)
    md <- multiSet(Set1 = data1, Set2 = data2)

    expect_is(md, "multiSet")
    expect_equal(length(md), 2L)
    expect_equal(names(md), c("Set1", "Set2"))
    expect_is(md$Set1, "list")
})

test_that("list2multiSet works", {
    data1 <- matrix(rnorm(100), 20, 5)
    data2 <- matrix(rnorm(50), 10, 5)
    md <- list2multiSet(list(Set1 = data1, Set2 = data2))

    expect_is(md, "multiSet")
    expect_equal(length(md), 2L)
    expect_equal(names(md), c("Set1", "Set2"))
    expect_is(md$Set1, "list")
})

test_that("multiSet2list", {
    data1 <- matrix(rnorm(100L), 20L, 5L)
    data2 <- matrix(rnorm(50L), 10L, 5L)
    md <- multiSet(Set1 = data1, Set2 = data2)
    testing <- multiSet2list(md)

    expect_is(testing, "list")
    expect_equal(names(testing), c("Set1", "Set2"))
    expect_equal(dim(testing$Set1), c(20L, 5L))
})

test_that("multiSet is subsetable", {
    data1 <- matrix(rnorm(100L), 20L, 5L)
    data2 <- matrix(rnorm(50L), 10L, 5L)
    colnames(data1) <- LETTERS[1:5]
    colnames(data2) <- LETTERS[2:6]
    md <- multiSet(Set1 = data1, Set2 = data2)

    expect_equal(dim(md[, c(1,3)]), NULL)
    expect_equal(dim(md[, c(1,3)]$Set1$data), c(20L, 2L))
    expect_equal(dim(md[, c(1,3)]$Set2$data), c(10L, 2L))

    expect_equal(dim(md[c(1,3), ]), NULL)
    expect_equal(dim(md[c(1,3), ]$Set1$data), c(2L, 5L))
    expect_equal(dim(md[c(1,3), ]$Set2$data), c(2L, 5L))

    expect_equal(dim(md[, c("B", "C")]$Set1$data), c(20L, 2L))
    expect_equal(dim(md[, c("B", "C")]$Set2$data), c(10L, 2L))
    expect_error(md[, c("B", "A")], "subscript out of bounds")
})

test_that("checkSets works properly", {
    data1 <- matrix(rnorm(100L), 20L, 5L)
    data2 <- matrix(rnorm(50L), 10L, 5L)
    colnames(data1) <- LETTERS[1:5]
    colnames(data2) <- LETTERS[2:6]
    md <- multiSet(Set1 = data1, Set2 = data2)
    testing <- checkSets(md)

    expect_equal(testing$nSets, 2L)
    expect_equal(testing$nGenes, 5L)
    expect_equal(testing$nSamples, c(20L, 10L))
    expect_true(testing$structureOK)
})

test_that("union and intersect work", {
    data1 <- matrix(rnorm(100L), 20L, 5L)
    data2 <- matrix(rnorm(50L), 10L, 5L)
    colnames(data1) <- LETTERS[1:5]
    colnames(data2) <- LETTERS[2:6]
    uni <- multiUnion(list(data1, data2))
    int <- multiIntersect(list(data1, data2))

    expect_equal(length(uni), 150L)
    expect_equal(length(int), 0L)
})
