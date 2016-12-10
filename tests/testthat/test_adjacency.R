library("WGCNA")
context("Testing adjacency related")

test_that("checkAdjMat functionallity", {

    expect_error(checkAdjMat(1))

    adj <- matrix(c(1, 2, 1, 1), ncol = 2, nrow = 2)
    expect_error(checkAdjMat(adj)) # Checking symmetry

    adj <- matrix(c(1, 2, 2, 1), ncol = 2, nrow = 2)
    expect_error(checkAdjMat(adj))
    expect_true(checkAdjMat(adj, max = 2))

    adj <- matrix(c(1, 2, 2, 1), ncol = 1, nrow = 2)
    expect_error(checkAdjMat(adj, max = 2))
    expect_error(checkAdjMat(adj))

    adj <- matrix(c(-1, 1, 1, -1), ncol = 2, nrow = 2)
    expect_error(checkAdjMat(adj))
    expect_true(checkAdjMat(adj, min = -2))
})

test_that("unsignedAdjacency", {
    datExpr <- matrix(rnorm(150), ncol = 5)
    adj <- unsignedAdjacency(datExpr)
    expect_true(all(adj <= 1))
    expect_true(all(adj >= 0))

    adj <- unsignedAdjacency(datExpr, datExpr)
    expect_true(all(adj <= 1))
    expect_true(all(adj >= 0))
})

test_that("adjacency", {
    datExpr <- matrix(seq(-1, 1, length.out = 25), 5)
    datExpr[lower.tri(datExpr)] <- t(datExpr)[lower.tri(datExpr)]
    adj <- adjacency(datExpr, selectCols = c(2, 4))
    expect_equal(ncol(adj), 2)
    expect_equal(nrow(adj), 5)
})

test_that("similarity", {
    datExpr <- matrix(rnorm(100), 10, 10)
    similarity <- similarity(datExpr)
    expect_true(all(similarity <= 1))
    expect_true(all(similarity >= -1))
})

test_that("similarity + power = adjacency", {
    datExpr <- matrix(rnorm(100), 10, 10)
    similarity <- similarity(datExpr)
    adj0 <- similarity^6
    adj <- adjacency(datExpr)
    expect_equal(adj0, adj)
})

test_that("adjacency.polyReg", {
    m <- 50
    x1 <- rnorm(m)
    r <- 0.5
    x2 <- r * x1 + sqrt(1 - r^2) * rnorm(m)
    r <- 0.3
    x3 <- r * (x1 - 0.5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
    x4 <- rnorm(m)
    r <- 0.3
    x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
    datE <- data.frame(x1, x2, x3, x4, x5)

    adj.max <- adjacency.polyReg(datE, symmetrizationMethod = "max")
    adj.mean <- adjacency.polyReg(datE, symmetrizationMethod = "mean")
    adj.none <- adjacency.polyReg(datE, symmetrizationMethod = "none")

    expect_true(all(adj.max <= 1))
    expect_true(all(adj.max >= 0))

    expect_error(checkAdjMat(adj.none), "not symmetric")
})

test_that("adjacency.splineReg", {
    m <- 50
    x1 <- rnorm(m)
    r <- 0.5
    x2 <- r * x1 + sqrt(1 - r^2) * rnorm(m)
    r <- 0.3
    x3 <- r * (x1 - 0.5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
    x4 <- rnorm(m)
    r <- 0.3
    x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
    datE <- data.frame(x1, x2, x3, x4, x5)

    adj.max <- adjacency.splineReg(datE, symmetrizationMethod = "max")
    adj.mean <- adjacency.splineReg(datE, symmetrizationMethod = "mean")
    adj.none <- adjacency.splineReg(datE, symmetrizationMethod = "none")

    expect_true(all(adj.max <= 1))
    expect_true(all(adj.max >= 0))

    expect_error(checkAdjMat(adj.none), "not symmetric")
})

test_that("signumAdjacency", {
    datExpr <- matrix(rnorm(150), ncol = 5)
    corMat <- cor(datExpr)
    adj <- signumAdjacency(corMat, 0.1)
    expect_true(sum(adj == 1, adj == 0) == length(adj))
})


test_that("sigmoidAdjacency", {
    datExpr <- matrix(rnorm(150), ncol = 5)
    corMat <- cor(datExpr)
    adj <- sigmoidAdjacency(abs(corMat), 0.1)
    expect_equal(ncol(adj), 5L)
    expect_equal(nrow(adj), 5L)
})
