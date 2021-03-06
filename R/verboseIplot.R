#' Scatterplot with density
#' 
#' Produce a scatterplot that shows density with color and is annotated by the
#' correlation, MSE, and regression line.
#' 
#' Irrespective of the specified correlation function, the MSE is always
#' calculated based on the residuals of a linear model.
#' 
#' @param x numerical vector to be plotted along the x axis.
#' @param y numerical vector to be plotted along the y axis.
#' @param xlim define the range in x axis
#' @param ylim define the range in y axis
#' @param nBinsX number of bins along the x axis
#' @param nBinsY number of bins along the y axis
#' @param ztransf Function to transform the number of counts per pixel, which
#' will be mapped by the function in colramp to well defined colors. The user
#' has to make sure that the transformed density lies in the range [0,zmax],
#' where zmax is any positive number (>=2).
#' @param gamma color correction power
#' @param sample either a number of points to be sampled or a vector of indices
#' input \code{x} and \code{y} for points to be plotted.  Useful when the input
#' vectors are large and plotting all points is not practical.
#' @param corFnc character string giving the correlation function to annotate
#' the plot.
#' @param corOptions character string giving further options to the correlation
#' function.
#' @param main main title for the plot.
#' @param xlab label for the x-axis.
#' @param ylab label for the y-axis.
#' @param cex character expansion factor for plot annotations.
#' @param cex.axis character expansion factor for axis annotations.
#' @param cex.lab character expansion factor for axis labels.
#' @param cex.main character expansion factor for the main title.
#' @param abline logical: should the linear regression fit line be plotted?
#' @param abline.color color specification for the fit line.
#' @param abline.lty line type for the fit line.
#' @param corLabel character string to be used as the label for the correlation
#' value printed in the main title.
#' @param \dots other arguments to the function plot.
#' @return If sample above is given, the indices of the plotted points are
#' returned invisibly.
#' @note This funtion is based on verboseScatterplot (Steve Horvath and Peter
#' Langfelder), iplot (Andreas Ruckstuhl, Rene Locher) and greenWhiteRed(Peter
#' Langfelder )
#' @author Chaochao Cai, Steve Horvath
#' @seealso \link{image} for more parameters
#' @keywords graphics
verboseIplot <- function(x, y, xlim = NA, ylim = NA, nBinsX = 150, nBinsY = 150,
                        ztransf = function(x){x}, gamma = 1, sample = NULL,
                        corFnc = "cor", corOptions = "use = 'p'", main = "",
                        xlab = NA, ylab = NA, cex = 1, cex.axis = 1.5,
                        cex.lab = 1.5, cex.main = 1.5, abline = FALSE,
                        abline.color = 1, abline.lty = 1, corLabel = corFnc,
                        ...) {
	if (is.na(xlab)) {
		xlab <- deparse(substitute(x))
	}
	if (is.na(ylab)) {
		ylab <- deparse(substitute(y))
	}
	x <- as.numeric(as.character(x))
	y <- as.numeric(as.character(y))
	xy <- data.frame(x,y)
	xy <- xy[!is.na(x)&!is.na(y),]

	if (sum(is.na(xlim)) != 0)
		xlim <- c(min(xy[, 1])-10^-10*diff(range(xy[, 1])), max(xy[, 1]))
	if (sum(is.na(ylim)) != 0)
		ylim <- c(min(xy[, 2])-10^-10*diff(range(xy[, 2])), max(xy[, 2]))

	corExpr <- parse(text = paste(corFnc, "(x, y ", prepComma(corOptions),")"))
	cor <- signif(eval(corExpr), 2)
	corp <- signif(corPvalueStudent(cor, sum(is.finite(x) & is.finite(y))), 2)
	if (corp < 10^(-200)) {
		corp <- "<1e-200"
	} else {
	    corp <- paste0("=", corp)
	}

	resid <- lm(y~x)$residuals
	MSE <- round(mean(resid^2), 2)
	if (!is.na(corLabel)) {
		mainX <- paste0(main, " ", corLabel, "=", cor, " MSE = ", MSE)
	} else {
	    mainX <- main
	}

	if (!is.null(sample)) {
	    if (length(sample) == 1) {
	        sample <- sample(length(x), sample)
	    }
	    xy <- xy[sample, ]
	}

	sx <- seq(xlim[1], xlim[2], by = diff(xlim)/nBinsX)
	sy <- seq(ylim[1], ylim[2], by = diff(ylim)/nBinsY)
	den <- ztransf(table(cut(xy[, 1], breaks =  sx), cut(xy[, 2], breaks = sy)))

	lsx <- length(sx)
	lsy <- length(sy)
	xx <- 0.5 * (sx[-1] + sx[-lsx])
	yy <- 0.5 * (sy[-1] + sy[-lsy])

	whiteBlueGreenRedBlack <- function(n) {
		quarter <- as.integer(n/5)
		power <- (1/gamma)
		red <- c(seq(1, 0, length.out = quarter)^power,
	         seq(0, 0, length.out = quarter)^power,
		         seq(0, 1, length.out = quarter)^power,
		         seq(1, 1, length.out = quarter)^power,
		         seq(1, 0, length.out = quarter)^power)
		green <- c(seq(1, 1, length.out = quarter)^power,
		           seq(1, 1, length.out = quarter)^power,
		           seq(1, 1, length.out = quarter)^power,
		           seq(1, 0, length.out = quarter)^power,
		           seq(0, 0, length.out = quarter)^power)
		blue <- c(seq( 1, 1, length.out = quarter)^power,
		          seq(1, 0, length.out = quarter)^power,
		          seq(0, 0, length.out = quarter)^power,
		          seq(0, 0, length.out = quarter)^power,
		          seq(0, 0, length.out = quarter)^power)
		col <- rgb(red, green, blue, maxColorValue = 1)
		col
	}


	image(x = xx, y = yy, den,  xaxs = "r", yaxs = "r", xlab = xlab,
	      ylab = ylab, cex = cex, main=mainX, cex.axis = cex.axis,
	      cex.lab = cex.lab, cex.main = cex.main,
	      col=whiteBlueGreenRedBlack(50))

	if (abline) {
	    fit <- lm(y ~ x)
	    abline(reg = fit, col = abline.color, lty = abline.lty)
	}
	invisible(sample)
}
