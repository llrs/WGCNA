library("WGCNA")
context("Testing useNThreads")

test_that("simulateModule works", {
  expect_null(disableWGCNAThreads())
})
test_that("WGNCAnTrheads", {
  expect_true(is.numeric(WGCNAnThreads()))
})


test_that("allocateJobs", {
  test <- allocateJobs(100, 2)
  expect_true(is.list(test))
  expect_length(test, 2L)
  expect_equal(lengths(test), c(50, 50))
})
