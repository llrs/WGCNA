library("WGCNA")
context("Testing smaFunctions")

test_that("rgcolors.func", {
    expect_true(all.equal(rgcolors.func(n = 2), c("#000000", "#000000")))
})

