library("WGCNA")
context("Testing accuracyMeasures")

test_that("accuracyMeasures in table format", {
    m <- 100
    trueOutcome <- sample(c(1,2), m, replace = TRUE)
    predictedOutcome <- trueOutcome
    # now we noise half of the entries of the predicted outcome
    predictedOutcome[ 1:(m/2)] <- sample(predictedOutcome[ 1:(m/2)] )
    tab <- table(predictedOutcome, trueOutcome)
    test1 <- accuracyMeasures(tab)
    test2 <- accuracyMeasures(predictedOutcome, trueOutcome)
    expect_equal(nrow(test1), 14L)
    expect_equal(ncol(test1), 2L)
    expect_equal(nrow(test2), 14L)
    expect_equal(ncol(test2), 2L)
    expect_true(all.equal(colnames(test1), c("Measure", "Value")))
    expect_true(all.equal(colnames(test2), c("Measure", "Value")))

    # Test if they are both the same matrix
    expect_true(all.equal(dim(test1), dim(test2)))
    expect_equal(sum(test1 == test2), 28)
})
