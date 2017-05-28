library("WGCNA")
context("Testing blockwiseModulesC")

set.seed(123)
datExpr <- matrix(runif(500*25), ncol = 500, nrow = 25)
sim <- cor(datExpr)

test_that("TOMsimilarityFromExpr", {
  TOM <- TOMsimilarityFromExpr(datExpr)
  expect_equal(dim(TOM), c(500, 500))
  expect_true(all(diag(TOM) == 1))
  expect_true(isSymmetric(TOM))
})


test_that("TOMsimilarity", {
  expect_error(TOMsimilarity(datExpr), "is not square")
  expect_error(TOMsimilarity(sim), "not between 0 and 1")

  TOM <- TOMsimilarity(abs(sim))
  expect_equal(dim(TOM), c(500, 500))
  expect_true(all(diag(TOM) == 1))
  expect_true(isSymmetric(TOM))
})

test_that("blockwiseModules", {
  expect_warning(test <- blockwiseModules(datExpr), "Not in multiSet format")
  expect_length(test, 10)
  expect_equal(names(test), c("colors", "unmergedColors", "MEs", "goodSamples",
                           "goodGenes", "dendrograms", "TOMFiles",
                           "blockGenes", "blocks", "MEsOK"))
  expect_equal(dim(test$MEs), c(25, 1))
  expect_length(test$goodSamples, 25)
  expect_true(all(test$goodSamples))
  expect_length(test$colors, 500)
  expect_length(test$unmergedColors, 500)
  expect_equal(class(test$dendrograms[[1]]), "hclust")
  expect_true(test$MEsOK)
  expect_true(all(test$blocks == 1))
  expect_equal(test$blockGenes[[1]], seq_len(500))
})


test_that("projectiveKMeans",  {
  test <- projectiveKMeans(datExpr)
  expect_equal(names(test), c("clusters", "centers"))
  expect_true(all(test$clusters == 1))
  expect_length(test$clusters, 500)
  expect_length(test$centers, 25)

})
