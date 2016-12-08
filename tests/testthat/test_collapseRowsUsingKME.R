library("WGCNA")
context("Testing collapseRowsUsingKME")

test_that("collapseRowsUsingKME", {
    set.seed(100)
    ME.A <- sample(1:100, 50)
    ME.B <- sample(1:100, 50)
    ME.C <- sample(1:100, 50)
    ME.D <- sample(1:100, 50)
    ME1 <- data.frame(ME.A, ME.B, ME.C, ME.D)
    simDatA <- simulateDatExpr(ME1, 1000, c(0.2, 0.1, 0.08, 0.05, 0.3),
    signed = TRUE, verbose = 0)
    simDatB <- simulateDatExpr(ME1, 1000, c(0.2, 0.1, 0.08, 0.05, 0.3),
    signed = TRUE, verbose = 0)
    Gin <- c(colnames(simDatA$datExpr), colnames(simDatB$datExpr))
    Pin <- paste("Probe", 1:length(Gin), sep=".")
    datExpr <- cbind(simDatA$datExpr, simDatB$datExpr)
    MM <- cor(datExpr, ME1)
    results <- collapseRowsUsingKME(MM, Gin, Pin)

    expect_equal(ncol(results$MMcollapsed), 4L)
    expect_equal(nrow(results$MMcollapsed), 1000L)

    expect_equal(ncol(results$group2Row), 2L)
    expect_equal(nrow(results$group2Row), 1000L)

    expect_length(results$selectedRow, 2000L)
    expect_true(is.logical(results$selectedRow))
})
