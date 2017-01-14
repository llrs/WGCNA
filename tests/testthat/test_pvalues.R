library("WGCNA")
context("Testing pvalues")

test_that("corPvalueFisher", {
    test <- corPvalueFisher(0.5, 10)
    expect_equal(test, 0.146132858882154)
    test <- corPvalueFisher(0.5, 10, FALSE)
    expect_equal(test, 0.926933570558923)
    expect_error(corPvalueFisher(1.2, 10))
})

test_that("corPvalueStudent", {
    test <- corPvalueStudent(0.5, 10)
    expect_equal(test, 0.14111328125)
})

test_that("rankPvalue", {
    z1 <- 1:5
    z2 <- 1:5
    z3 <- c(2, 1, 3:5)
    z4 <- c(2, 1, 3:5)
    z5 <- c(2:3, 1, 4:5)
    datS <- data.frame(z1, z2, z3, z4, z5)
    expect_error(rankPvalue(datS), "valid p-values")
    datS <-  datS/rowSums(datS)

})
