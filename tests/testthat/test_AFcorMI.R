library("WGCNA")
context("Testing AFcorMI")

test_that("AFcorMI", {
    m <- 50
    x1 <- rnorm(m)
    r <- 0.5
    x2 <- r*x1+sqrt(1 - r ^ 2) * rnorm(m)
    r <- 0.3
    x3 <- r * (x1 - 0.5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
    x4 <- rnorm(m)
    r <- 0.3
    x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
    datE <- data.frame(x1, x2, x3, x4, x5)
    cor.data <- cor(datE, use = "p")
    AUV2 <- AFcorMI(r = cor.data, m = nrow(datE))
    expect_true(all.equal(dim(AUV2), c(5, 5)))
    })
