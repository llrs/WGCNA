library("WGCNA")
context("Testing kMEcomparisonScatterplot")

test_that("kMEcomparisonScatterplot", {
    set.seed <- 100
    ME <- matrix(0, 50, 5)
    for (i in 1:5) {
        ME[, i] <- sample(1:100, 50)
    }
    simData1 <- simulateDatExpr5Modules(MEturquoise = ME[, 1], MEblue = ME[, 2],
                              MEbrown = ME[,3], MEyellow = ME[, 4],
                              MEgreen = ME[, 5])
    simData2 <- simulateDatExpr5Modules(MEturquoise = ME[, 1], MEblue = ME[, 2],
                              MEbrown = ME[, 3], MEyellow = ME[, 4],
                              MEgreen = ME[, 5])
    expect_error(kMEcomparisonScatterplot())
    expect_error(kMEcomparisonScatterplot(simData1))
    expect_error(kMEcomparisonScatterplot(simData1$datExpr, simData2$datExpr))
    expect_error(kMEcomparisonScatterplot(1, 2, "grey"))
    expect_error(kMEcomparisonScatterplot(simData1$datExpr, 2, "grey"))
    expect_error(kMEcomparisonScatterplot(simData1$datExpr, simData2$datExpr,
                                          "grey"))
})
