
# addErrorBars ####
#' Add error bars to a barplot.
#'
#' This function adds error bars to an existing barplot.
#'
#'
#' @param means vector of means plotted in the barplot
#' @param errors vector of standard errors (signle positive values) to be
#' plotted.
#' @param two.side should the error bars be two-sided?
#' @return None.
#' @author Steve Horvath and Peter Langfelder
#' @keywords hplot
addErrorBars <- function(means, errors, two.side = FALSE) {
    if (!is.numeric(means)) {
        stop("All arguments must be numeric")
    }

    if (is.null(dim(means)) || length(dim(means)) == 1) {
        xval <- cumsum(c(0.7, rep(1.2, length(means) - 1)))
    } else {
        if (length(dim(means)) == 2) {
            xval <- cumsum(array(c(1, rep(0, dim(means)[1] - 1)),
                                 dim = c(1, length(means))))
            xval <- xval + 0:(length(means) - 1) + .5
        } else {
            stop("First argument must either be a vector or a matrix")
        }
    }
    MW <- 0.25 * (max(xval) / length(xval))
    ERR1 <- means + errors
    ERR2 <- means - errors
    for (i in 1:length(means)) {
        segments(xval[i], means[i], xval[i], ERR1[i])
        segments(xval[i] - MW, ERR1[i], xval[i] + MW, ERR1[i])
        if (two.side) {
            segments(xval[i], means[i], xval[i], ERR2[i])
            segments(xval[i] - MW, ERR2[i], xval[i] + MW, ERR2[i])
        }
    }
}

# stdErr ####
#' Standard error of the mean of a given vector.
#'
#' Returns the standard error of the mean of a given vector. Missing values are
#' ignored.
#'
#'
#' @param x a numeric vector
#' @return Standard error of the mean of x.
#' @author Steve Horvath
#' @keywords misc
#' @examples
#' x <- c(0.2, 0.5, 0.8)
#' stdErr(x)
stdErr <- function(x) {
    sqrt(var(x, na.rm = TRUE) / sum(x, na.rm = TRUE))
}

#===============================================================================
# The following two functions are for displaying the pair - wise correlation in
# a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel = panel.smooth,
# lower.panel = .panel.cor, diag.panel = panel.hist)" to
# put the correlation coefficients on the lower panel.

.panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y / max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#===============================================================================
# This function is used in "pairs()" function. The problem of the original
# panel.cor is that when the correlation coefficient is very small, the lower
# panel will have a large font instead of a mini-font in a saved .ps file.
# This new function uses a format for corr = 0.2 when corr<0.2, but it still
# reports the original value of corr, with a minimum format.

.panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    txt1 <- txt
    r1 <- r
    if (r < 0.2) {
        r1 <- 0.2
        txt1 <- format(c(r1, 0.123456789), digits = digits)[1]
        txt1 <- paste0(prefix, txt1)
    }
    if (missing(cex.cor)) {
        cex <- 0.8 / strwidth(txt1)
    }
    cex <- cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    text(0.5, 0.5, txt, cex = cex)
}
