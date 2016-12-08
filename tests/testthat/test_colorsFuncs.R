library("WGCNA")
context("Testing colorsFuncs")

test_that("labels2colors", {
    labels <- c(0:20)
    col <- labels2colors(labels)
    expect_equal(length(col), 21)

    labels <- matrix(letters[1:9], 3, 3)
    labelc <- labels2colors(labels)
    labelunc <- labels2colors(labels, commonColorCode = FALSE)
    expect_equal(labelc[, 1], labelunc[, 1])
    expect_true(all(labelc[, 1] != labelc[, 2]))
})

test_that("StandardColors", {
    st <- standardColors()
    expect_length(st, 435)
    expect_length(standardColors(10), 10)
})
