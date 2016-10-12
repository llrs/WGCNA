#Plot functions extracted from Functions.R

#' @name verboseScatterplot
#' @rdname verboseScatterplot
#' @title Scatterplot annotated by regression line and p-value
#' @description
#' Produce a scatterplot annotated by the correlation, p-value, and regression
#' line.
#' @inheritParams verboseIplot
#' @param displayAsZero Correlations whose absolute value is smaller than this
#' number will be displayed as zero. This can result in a more intuitive display
#' (for example, cor=0 instead of cor=2.6e-17).
#' @paramm col color of the plotted symbols. Recycled as necessary.
#' @param bg fill color of the plotted symbols (used for certain symbols).
#' Recycled as necessary.
#' @param lmFnc linear model fit function. Used to calculate the linear model
#' fit line if 'abline' is TRUE. For example, robust linear models are
#' implemented in the function \code{\link{rlm}}.
#' @details
#' Irrespective of the specified correlation function, the p-value is always
#' calculated for pearson correlation.
#' @return
#' If sample above is given, the indices of the plotted points are returned
#' invisibly.
#' @author
#' Steve Horvath and Peter Langfelder
#' @seealso
#' \code{\link{plot.default}} for standard scatterplots
#' @export
verboseScatterplot <- function(x, y,
                               sample = NULL,
                               corFnc = "cor", corOptions = "use = 'p'",
                               main  = "", xlab = NA, ylab = NA, cex = 1,
                               cex.axis = 1.5,
                               cex.lab = 1.5, cex.main = 1.5, abline = FALSE,
                               abline.color = 1, abline.lty = 1,
                               corLabel = corFnc,
                               displayAsZero = 1e-5,
                               col = 1, bg = 0,
                               lmFnc = lm,
                               ...) {
    if (is.na(xlab)) xlab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(ylab)) ylab = as.character(match.call(expand.dots = FALSE)$y)
    x = as.numeric(as.character(x))
    y = as.numeric(as.character(y))
    corExpr = parse(text = paste(corFnc, "(x, y ", prepComma(corOptions), ")"))
    #cor = signif (cor(x, y, use = "p", method = correlationmethod), 2)
    cor = signif (eval(corExpr), 2)
    if (abs(cor) < displayAsZero) cor = 0
    corp = signif (corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))), 2)
    #corpExpr = parse(text = paste("cor.test(x, y, ", corOptions, ")"))
    #corp = signif (cor.test(x, y, use = "p",
    #method = correlationmethod)$p.value, 2)
    #corp = signif (eval(corpExpr)$p.value, 2)
    if (corp < 10^(- 200)) corp = "<1e-200" else corp = paste0(" = ", corp)
    if (!is.na(corLabel))
    {
        mainX = paste0(main, " ", corLabel, " = ", cor, ", p", corp)
    } else
        mainX = main

    if (!is.null(sample)) {
        if (length(sample) == 1) {
            sample = sample(length(x), sample)
        }
        if (length(col)<length(x)) col = rep(col, ceiling(length(x)/length(col)))
        if (length(bg)<length(x))  bg = rep(bg, ceiling(length(x)/length(bg)))
        plot(x[sample], y[sample], main = mainX, xlab = xlab, ylab = ylab,
             cex = cex,
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main,
             col = col[sample], bg = bg[sample], ...)
    } else {
        plot(x, y, main = mainX, xlab = xlab, ylab = ylab, cex = cex,
             cex.axis = cex.axis, cex.lab = cex.lab, cex.main = cex.main, col = col, bg = bg, ...)
    }
    if (abline) {
        lmFnc = match.fun(lmFnc)
        fit = lmFnc(y~x)
        abline(reg = fit, col = abline.color, lty = abline.lty)
    }
    invisible(sample)
}

#' @name verboseBoxplot
#' @rdname verboseBoxplot
#' @title Boxplot annotated by a Kruskal-Wallis p-value
#' @description
#' Plot a boxplot annotated by the Kruskal-Wallis p-value. Uses the function
#' \code{\link{boxplot}} for the actual drawing.
#' @inheritParams verboseIplot
#' @param g a factor or a an object coercible to a factor giving the groups that
#' will go into each box.
#' @param notch logical: should the notches be drawn? See \code{\link{boxplot}} and
#' \code{\link{boxplot.stats}} for details.
#' @param varwidth logical: if TRUE, the boxes are drawn with widths
#' proportional to the square-roots of the number of observations in the groups.
#' @return
#' Returns the value returned by the function boxplot.
#' @author
#' Steve Horvath
#' @seealso
#' \code{\link{boxplot}}
#' @export
verboseBoxplot <- function(x, g,
                           main  = "", xlab = NA, ylab = NA, cex = 1,
                           cex.axis = 1.5, cex.lab = 1.5, cex.main = 1.5,
                           notch = TRUE, varwidth = TRUE, ...) {
    if (is.na(xlab)) xlab = as.character(match.call(expand.dots = FALSE)$g)
    #print(xlab1)
    if (is.na(ylab)) ylab = as.character(match.call(expand.dots = FALSE)$x)
    #print(ylab1)
    p1 = signif (kruskal.test(x, factor(g))$p.value, 2)
    #if (p1< 5.0 * 10^(- 22)) p1 = "< 5e-22"
    boxplot(x~factor(g), notch = notch, varwidth = varwidth,
            main = paste(main, "p  = ", p1),
            xlab = xlab, ylab = ylab, cex = cex, cex.axis = cex.axis,
            cex.lab = cex.lab, cex.main = cex.main, ...)
}

#' @name verboseBarplot
#' @rdname verboseBarplot
#' @title Barplot with error bars, annotated by Kruskal-Wallis or ANOVA p-value
#' @description
#' Produce a barplot with error bars, annotated by Kruskal-Wallis or ANOVA
#' p-value.
#' @inheritParams verboseIplot
#' @param color a vector giving the colors of the bars in the barplot.
#' @param numberStandardErrors size of the error bars in terms of standard
#' errors. See details.
#' @param KruskalTest logical: should Kruskal-Wallis test be performed?
#' See details.
#' @param AnovaTest logical: should ANOVA be performed? See details.
#' @param two.sided logical: should the printed p-value be two-sided?
#' See details.
#' @param addCellCounts logical: should counts be printed above each bar?
#' @param horiz logical:should the bars be drawn horizontally?
#' @param ... other parameters to the function \code{\link{barplot}}
#' @details
#' This function creates a barplot of a numeric variable (input x) across the
#' levels of a grouping variable (input g). The height of the bars equals the
#' mean value of x across the observations with a given level of g. By default,
#' the barplot also shows plus/minus one standard error. If you want only plus
#' one standard error (not minus) choose two.sided=TRUE. But the number of
#' standard errors can be determined with the input numberStandardErrors. For
#' example, if you want a 95% confidence interval around the mean, choose
#' numberStandardErrors=2. If you don't want any standard errors set
#' numberStandardErrors=-1. The function also outputs the p-value of a Kruskal
#' Wallis test (Fisher test for binary input data), which is a non-parametric
#' multi group comparison test. Alternatively, one can use Analysis of Variance
#' (Anova) to compute a p-value by setting AnovaTest=TRUE. Anova is a
#' generalization of the Student t-test to multiple groups. In case of two
#' groups, the Anova p-value equals the Student t-test p-value. Anova should
#' only be used if x follows a normal distribution. Anova also assumes
#' homoscedasticity (equal variances). The Kruskal Wallis test is often
#' advantageous since it makes no distributional assumptions. Since the Kruskal
#' Wallis test is based on the ranks of x, it is more robust with regard to
#' outliers. All p-values are two-sided.
#' @return
#' A boxplot
#' @author
#' Steve Horvath
#' @seealso
#' \code{\link{barplot}}
#' @examples
#' group=sample(c(1,2),100,replace=TRUE)
#'
#' height=rnorm(100,mean=group)
#'
#' par(mfrow=c(2,2))
#' verboseBarplot(height,group, main="1 SE, Kruskal Test")
#'
#' verboseBarplot(height,group,numberStandardErrors=2,
#'                main="2 SE, Kruskal Test")
#'
#' verboseBarplot(height,group,numberStandardErrors=2,AnovaTest=TRUE,
#'                main="2 SE, Anova")
#'
#' verboseBarplot(height,group,numberStandardErrors=2,AnovaTest=TRUE,
#'                main="2 SE, Anova, only plus SE", two.sided=FALSE)
#' @export
verboseBarplot <- function(x, g, main = "",
                           xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5, cex.lab = 1.5,
                           cex.main = 1.5, color = "grey", numberStandardErrors = 1,
                           KruskalTest = TRUE, AnovaTest = FALSE, two.sided = TRUE,
                           addCellCounts = FALSE, horiz = FALSE, ...) {
    stderr1 <- function(x) {
        sqrt(var(x, na.rm = TRUE)/sum(!is.na(x)))
    }
    SE = tapply(x, factor(g), stderr1)
    err.bp <- function(dd, error, two.sided = FALSE, numberStandardErrors,
                       horiz = FALSE) {
        if (!is.numeric(dd)) {
            stop("All arguments must be numeric")
        }
        if (is.vector(dd)) {
            xval = (cumsum(c(0.7, rep(1.2, length(dd) - 1))))
        }
        else {
            if (is.matrix(dd)) {
                xval = cumsum(array(c(1, rep(0, dim(dd)[1]  -
                                                 1)), dim = c(1, length(dd)))) + 0:(length(dd)  -
                                                                                        1) + 0.5
            }
            else {
                stop("First argument must either be a vector or a matrix")
            }
        }
        MW = 0.25 * (max(xval)/length(xval))
        NoStandardErrors = 1
        ERR1 = dd + numberStandardErrors * error
        ERR2 = dd - numberStandardErrors * error
        if (horiz) {
            for (i in 1:length(dd)) {
                segments(dd[i], xval[i], ERR1[i], xval[i])
                segments(ERR1[i], xval[i] - MW, ERR1[i], xval[i]  +
                             MW)
                if (two.sided) {
                    segments(dd[i], xval[i], ERR2[i], xval[i])
                    segments(ERR2[i], xval[i] - MW, ERR2[i], xval[i]  +
                                 MW)
                }
            }
        }
        else {
            for (i in 1:length(dd)) {
                segments(xval[i], dd[i], xval[i], ERR1[i])
                segments(xval[i] - MW, ERR1[i], xval[i] + MW,
                         ERR1[i])
                if (two.sided) {
                    segments(xval[i], dd[i], xval[i], ERR2[i])
                    segments(xval[i] - MW, ERR2[i], xval[i] + MW,
                             ERR2[i])
                }
            }
        }
    }
    if (is.na(ylab))
        ylab = as.character(match.call(expand.dots = FALSE)$x)
    if (is.na(xlab))
        xlab = as.character(match.call(expand.dots = FALSE)$g)
    Means1 = tapply(x, factor(g), mean, na.rm = TRUE)

    if (length(unique(x)) > 2) {
        p1 = signif (kruskal.test(x ~ factor(g))$p.value, 2)
        if (AnovaTest)
            p1 = signif (anova(lm(x ~ factor(g)))$Pr[[1]], 2)
    }
    else {
        p1 = tryCatch(
            signif (fisher.test(x, g, alternative = "two.sided")$p.value,2),
            error <- function(e) {
                NA}
        )
    }
    if (AnovaTest | KruskalTest)
        main = paste(main, "p  = ", p1)
    ret = barplot(Means1, main = main, col = color, xlab = xlab,
                  ylab = ylab, cex = cex, cex.axis = cex.axis, cex.lab = cex.lab,
                  cex.main = cex.main, horiz = horiz, ...)
    if (addCellCounts) {
        cellCountsF <- function(x) {  sum(!is.na(x)) }
        cellCounts = tapply(x, factor(g), cellCountsF)
        mtext(text = cellCounts, side = if (horiz) 2 else 1, outer = FALSE,
              at = ret, col = "darkgrey", las = 2, cex = .8, ...)
    } # end of if (addCellCounts)
    abline(h = 0)
    if (numberStandardErrors > 0) {
        err.bp(as.vector(Means1), as.vector(SE), two.sided = two.sided,
               numberStandardErrors = numberStandardErrors, horiz = horiz)
    }
    attr(ret, "height") = as.vector(Means1)
    attr(ret, "stdErr") = as.vector(SE)
    invisible(ret)
}
