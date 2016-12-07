library("WGCNA")
context("Testing accuracyMeasures")

test_that("accuracyMeasures in table format", {
    m <- 100
    trueOutcome <- sample(c(1,2), m, replace = TRUE)
    predictedOutcome <- trueOutcome
    # now we noise half of the entries of the predicted outcome
    predictedOutcome[ 1:(m/2)] <- sample(predictedOutcome[ 1:(m/2)] )
    tab <- table(predictedOutcome, trueOutcome)
    test <- accuracyMeasures(tab)
    expect_equal(nrow(test), 14L)
    expect_equal(ncol(test), 2L)
    expect_true(all.equal(colnames(test), c("Measure", "Value")))
})
