library("WGCNA")
context("Testing simulateNetworks")

test_that("simulateModule works", {
  m <- 50
  y <- sample(c(1,2), m, replace = TRUE)
  test1 <- simulateModule(scale(y), 30)
  expect_equal(ncol(test1), 30L)
  expect_equal(nrow(test1), 50L)

  test2 <- simulateModule(scale(y), 30, signed = TRUE)
  expect_equal(ncol(test2), 30L)
  expect_equal(nrow(test2), 50L)

  test3 <- simulateModule(scale(y), 30, signed = TRUE,
                          geneMeans = sample(c(4, 5, 10), 30, replace = TRUE))
  expect_equal(ncol(test3), 30L)
  expect_equal(nrow(test3), 50L)
  expect_error(simulateModule(scale(y), 30, signed = TRUE,
                              geneMeans = sample(c(4, 5, 10), 25,
                                                 replace = TRUE)),
               "must equal")
  expect_error(simulateModule(scale(y), 30, signed = TRUE,
                              geneMeans = NA), "finite")

})

test_that("simulateEigengeneNetwork works",{})
test_that("simulateSmallLayer works",{})
test_that("simulateDatExpr works",{})
test_that("simulateMultiExpr works",{})

test_that("simulateDatExpr5Modules works",{
  test <- simulateDatExpr5Modules(MEturquoise = 1, MEblue = 2, MEbrown = 1,
                                  MEyellow = 3, MEgreen = 45)
  expect_equal(ncol(test$datExpr), 2000L)
  expect_equal(nrow(test$datExpr), 1L)
  expect_equal(as.numeric(test$datME), c(1, 2, 1, 3, 45))

  test2 <- simulateDatExpr5Modules(MEturquoise = c(1, 2), MEblue = c(2, 1),
                                   MEbrown = c(1, 1),
                                   MEyellow = c(2, 3), MEgreen = c(2, 45))
  expect_equal(nrow(test2$datExpr), 2L)
  expect_error(simulateDatExpr5Modules(MEturquoise = c(1, 2),
                                                MEblue = c(2, 1),
                                                MEbrown = c(1, 1),
                                                MEyellow = c(2, 3),
                                                MEgreen = 2), "not consistent")
})
