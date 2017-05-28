library("WGCNA")
context("Testing progressIndicator")

test_that("initProgInd", {
  prog <- initProgInd()
  max <- 50
  for (c in 1:max) {
    prog <- updateProgInd(c/max, prog)
  }
})
