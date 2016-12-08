library("WGCNA")
context("Testing coxRegressionResiduals")

test_that("coxRegressionResiduals", {
    library(survival)
    time1 <- sample(1:100)
    event1 <- sample(c(1, 0), 100, replace=TRUE)

    event1[1:5] <- NA
    time1[1:5] <- NA
    datResiduals <- coxRegressionResiduals(time = time1, event = event1)
    cor.1 <- cor(datResiduals, use = "p")
    expect_equal(ncol(datResiduals), 2L)
    expect_equal(nrow(datResiduals), 100L)
    z <- rnorm(100)
    datResiduals <- coxRegressionResiduals(time = time1, event = event1,
                                           datCovariates <- data.frame(z))
    cor.2 <- cor(datResiduals, use = "p")
    expect_lt(cor.2[1, 2], cor.1[1, 2])
})
