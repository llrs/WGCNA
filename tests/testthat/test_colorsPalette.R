library("WGCNA")
context("Testing colorsPalette")

test_that("greenBlackRed", {
    expect_length(greenBlackRed(50), 51L)
    expect_length(greenBlackRed(50, 2), 51L)
})

test_that("greenWhiteRed", {
    expect_length(greenWhiteRed(50, warn = FALSE), 51L)
    expect_length(greenWhiteRed(50, 2, warn = FALSE), 51L)
    expect_warning(greenWhiteRed(50, warn = TRUE))
})

test_that("redWhiteGreen", {
    expect_length(redWhiteGreen(50), 51L)
    expect_length(redWhiteGreen(50, 2), 51L)
})

test_that("blueWhiteRed", {
    expect_length(blueWhiteRed(50), 50L)
    expect_length(blueWhiteRed(50, 2), 50L)
    expect_error(blueWhiteRed(50, 2, 2))
    expect_length(blueWhiteRed(50, 2, 0.5), 50L)
})
