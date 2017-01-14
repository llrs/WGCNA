library("WGCNA")
context("Testing qvalues")

test_that("qvalue", {
    expect_error(qvalue(0.5), "use another lambda method")
    expect_error(qvalue(1.2), "valid range")
    expect_error(qvalue(-0.1), "valid range")
    expect_error(qvalue(0.1, 2), "Lambda must be within")
    expect_error(qvalue(0.1, c(1, 2)), "you need at least 4 values")
})
