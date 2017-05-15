
# plotDendroAndColors ####
#' Dendrogram plot with color annotation of objects
#'
#' This function plots a hierarchical clustering dendrogram and color
#' annotation(s) of objects in the dendrogram underneath.
#'
#' The function slits the plotting device into two regions, plots the given
#' dendrogram in the upper region, then plots color rows in the region below
#' the dendrogram.
#'
#' @param dendro a hierarchical clustering dendrogram such as one produced by
#' \code{\link[stats]{hclust}}.
#' @param colors Coloring of objects on the dendrogram. Either a vector (one
#' color per object) or a matrix (can also be an array or a data frame) with
#' each column giving one color per object. Each column will be plotted as a
#' horizontal row of colors under the dendrogram.
#' @param groupLabels Labels for the colorings given in \code{colors}. The
#' labels will be printed to the left of the color rows in the plot. If the
#' argument is given, it must be a vector of length equal to the number of
#' columns in \code{colors}. If not given, \code{names(colors)} will be used if
#' available. If not, sequential numbers starting from 1 will be used.
#' @param rowText Optional labels to identify colors in the color rows.  If
#' given, must be either the same dimensions as \code{colors} or must have the
#' same number of rows and \code{textPositions} must be used to specify which
#' columns of \code{colors} each column of \code{rowText} corresponds to. Each
#' label that occurs will be displayed once, under the largest continuous block
#' of the corresponding \code{colors}.
#' @param rowTextAlignment Character string specifying whether the labels
#' should be left-justified to the start of the largest block of each label,
#' centered in the middle, or right-justified to the end of the largest block.
#' @param rowTextIgnore Optional specifications of labels that should be
#' ignored when displaying them using \code{rowText} above.
#' @param textPositions optional numeric vector of the same length as the
#' number of columns in \code{rowText} giving the color rows under which the
#' text rows should appear.
#' @param setLayout logical: should the plotting device be partitioned into a
#' standard layout? If \code{FALSE}, the user is responsible for partitioning.
#' The function expects two regions of the same width, the first one
#' immediately above the second one.
#' @param autoColorHeight logical: should the height of the color area below
#' the dendrogram be automatically adjusted for the number of traits? Only
#' effective if \code{setLayout} is \code{TRUE}.
#' @param colorHeight specifies the height of the color area under dendrogram
#' as a fraction of the height of the dendrogram area. Only effective when
#' \code{autoColorHeight} above is \code{FALSE}.
#' @param rowWidths optional specification of relative row widths for the color
#' and text (if given) rows. Need not sum to 1.
#' @param dendroLabels dendrogram labels. Set to \code{FALSE} to disable
#' dendrogram labels altogether; set to \code{NULL} to use row labels of
#' \code{datExpr}.
#' @param addGuide logical: should vertical "guide lines" be added to the
#' dendrogram plot? The lines make it easier to identify color codes with
#' individual samples.
#' @param guideAll logical: add a guide line for every sample? Only effective
#' for \code{addGuide} set \code{TRUE}.
#' @param guideCount number of guide lines to be plotted. Only effective when
#' \code{addGuide} is \code{TRUE} and \code{guideAll} is \code{FALSE}.
#' @param guideHang fraction of the dendrogram height to leave between the top
#' end of the guide line and the dendrogram merge height. If the guide lines
#' overlap with dendrogram labels, increase \code{guideHang} to leave more
#' space for the labels.
#' @param addTextGuide logical: should guide lines be added for the text rows
#' (if given)?
#' @param cex.colorLabels character expansion factor for trait labels.
#' @param cex.dendroLabels character expansion factor for dendrogram (sample)
#' labels.
#' @param cex.rowText character expansion factor for text rows (if given).
#' @param marAll a vector of length 4 giving the bottom, left, top and right
#' margins of the combined plot. There is no margin between the dendrogram and
#' the color plot underneath.
#' @param saveMar logical: save margins setting before starting the plot and
#' restore on exit?
#' @param abHeight optional specification of the height for a horizontal line
#' in the dendrogram, see \code{\link{abline}}.
#' @param abCol color for plotting the horizontal line.
#' @param \dots other graphical parameters to \code{\link{plot.hclust}}.
#' @return None.
#' @author Peter Langfelder
#' @seealso \code{\link{plotColorUnderTree}}
#' @keywords hplot
plotDendroAndColors <- function(dendro,
                                colors,
                                groupLabels = NULL,
                                rowText = NULL,
                                rowTextAlignment = c("left", "center", "right"),
                                rowTextIgnore = NULL,
                                textPositions = NULL,
                                setLayout = TRUE,
                                autoColorHeight = TRUE,
                                colorHeight = 0.2,
                                rowWidths = NULL,
                                dendroLabels = NULL,
                                addGuide = FALSE,
                                guideAll = FALSE,
                                guideCount = 50,
                                guideHang = 0.20,
                                addTextGuide = FALSE,
                                cex.colorLabels = 0.8,
                                cex.dendroLabels = 0.9,
                                cex.rowText = 0.8,
                                marAll = c(1, 5, 3, 1),
                                saveMar = TRUE,
                                abHeight = NULL,
                                abCol = "red",
                                ...) {
    oldMar = par("mar")
    on.exit(par(oldMar))
    if (!is.null(dim(colors))) {
        nRows = dim(colors)[2]
    } else
        nRows = 1
    if (!is.null(rowText)) {
        nRows = nRows + if (is.null(textPositions)) {
            nRows
        } else {
            length(textPositions)
        }
    }

    if (autoColorHeight)
        colorHeight = 0.2 + 0.3 * (1 - exp(-(nRows - 1) / 6))
    if (setLayout)
        layout(matrix(c(1:2), 2, 1),
               heights = c(1 - colorHeight, colorHeight))
    par(mar = c(0, marAll[2], marAll[3], marAll[4]))
    plot(dendro, labels = dendroLabels, cex = cex.dendroLabels, ...)
    if (addGuide)
        addGuideLines(
            dendro,
            count =
                if (guideAll)
                    length(dendro$height) + 1
            else
                guideCount,
            hang = guideHang
        )
    if (!is.null(abHeight))
        abline(h = abHeight, col = abCol)
    par(mar = c(marAll[1], marAll[2], 0, marAll[4]))
    plotColorUnderTree(
        dendro,
        colors,
        groupLabels,
        cex.rowLabels = cex.colorLabels,
        rowText = rowText,
        rowTextAlignment = rowTextAlignment,
        rowTextIgnore = rowTextIgnore,
        textPositions = textPositions,
        cex.rowText = cex.rowText,
        rowWidths = rowWidths,
        addTextGuide = addTextGuide
    )
    if (saveMar)
        par(mar = oldMar)
}

# plotMEpairs ####
#' Pairwise scatterplots of eigengenes
#'
#' The function produces a matrix of plots containing pairwise scatterplots of
#' given eigengenes, the distribution of their values and their pairwise
#' correlations.
#'
#' The function produces an NxN matrix of plots, where N is the number of
#' eigengenes. In the upper traingle it plots pairwise scatterplots of module
#' eigengenes (plus the trait \code{y}, if given). On the diagonal it plots
#' histograms of sample values for each eigengene. Below the diagonal, it
#' displays the pairwise correlations of the eigengenes.
#'
#' @param datME a data frame containing expression data, with rows
#' corresponding to samples and columns to genes. Missing values are allowed
#' and will be ignored.
#' @param y optional microarray sample trait vector. Will be treated as an
#' additional eigengene.
#' @param main main title for the plot.
#' @param clusterMEs logical: should the module eigengenes be ordered by their
#' dendrogram?
#' @param \dots additional graphical parameters to the function
#' \code{\link{pairs}}
#' @return None.
#' @author Steve Horvath
#' @seealso \code{\link{pairs}}
#' @keywords hplot
plotMEpairs <- function(datME,
                        y = NULL,
                        main = "Relationship between module eigengenes",
                        clusterMEs = TRUE,
                        ...) {
    if (dim(as.matrix(datME))[[2]] == 1 & is.null(y)) {
        hist(datME, ...)
    } else {
        datMEordered = datME
        if (clusterMEs & dim(as.matrix(datME))[[1]]  > 1)
        {
            dissimME = (1 - t(cor(
                datME, method = "p", use = "p"
            ))) / 2
            hclustdatME = fastcluster::hclust(as.dist(dissimME),
                                              method = "average")
            datMEordered = datME[, hclustdatME$order]
        } # end of if
        if (!is.null(y)) {
            if (length(y)   != dim(as.matrix(datMEordered))[[1]])
                stop(
                    paste(
                        "The length of the outcome vector 'y' does not match
                        the number of rows of 'datME'.\n",
                        "     The columns of datME should correspond to the
                        module eigengenes.\n",
                        "     The rows correspond to the array samples.
                        Hint: consider transposing datME."
                    )
                    )
            datMEordered = data.frame(y, datMEordered)
        } # end of if
        pairs(
            datMEordered,
            upper.panel = panel.smooth,
            lower.panel = .panel.cor,
            diag.panel = .panel.hist,
            main = main,
            ...
        )
    } # end if
} # end of function

# ScaleFreePlot ####
# The function ScaleFreePlot creates a plot for checking scale free topology
# when truncated1 = TRUE is specificed, it provides the R^2 measures for the
# following degree distributions:
# a) scale free topology,
# b) log - log R^2 and
# c) truncated exponential R^2
#' Visual check of scale-free topology
#'
#' A simple visula check of scale-free network ropology.
#'
#' The function plots a log-log plot of a histogram of the given
#' \code{connectivities}, and fits a linear model plus optionally a truncated
#' exponential model. The \eqn{R^2} of the fit can be considered an index of
#' the scale freedom of the network topology.
#'
#' @param connectivity vector containing network connectivities.
#' @param nBreaks number of breaks in the connectivity dendrogram.
#' @param truncated logical: should a truncated exponential fit be calculated
#' and plotted in addition to the linear one?
#' @param removeFirst logical: should the first bin be removed from the fit?
#' @param main main title for the plot.
#' @param \dots other graphical parameter to the \code{plot} function.
#' @return None.
#' @author Steve Horvath
#' @seealso
#' \code{\link{softConnectivity}} for connectivity calculation in
#' weigheted networks.
#' @references
#' Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
scaleFreePlot <- function(connectivity, nBreaks = 10, truncated = FALSE,
                          removeFirst = FALSE, main = "", ...) {
    k = connectivity
    discretized.k = cut(k, nBreaks)
    dk = tapply(k, discretized.k, mean)
    p.dk = as.vector(tapply(k, discretized.k, length) / length(k))
    breaks1 = seq(from = min(k),
                  to = max(k),
                  length = nBreaks + 1)
    hist1 = suppressWarnings(hist(
        k,
        breaks = breaks1,
        equidist = FALSE,
        plot = FALSE,
        right = TRUE,
        ...
    ))
    dk2 = hist1$mids
    dk = ifelse(is.na(dk), dk2, dk)
    dk = ifelse(dk == 0, dk2, dk)
    p.dk = ifelse(is.na(p.dk), 0, p.dk)
    log.dk = as.vector(log10(dk))
    if (removeFirst) {
        p.dk = p.dk[-1]
        log.dk = log.dk[-1]
    }
    log.p.dk = as.numeric(log10(p.dk + 1e-09))
    lm1 = lm(log.p.dk ~ log.dk)
    if (truncated == TRUE) {
        lm2 = lm(log.p.dk ~ log.dk + I(10 ^ log.dk))
        OUTPUT = data.frame(
            scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 2),
            slope = round(lm1$coefficients[[2]], 2),
            TruncatedRsquared = round(summary(lm2)$adj.r.squared, 2)
        )
        printFlush("the red line corresponds to the truncated exponential fit")
        title = paste(
            main,
            " scale free R^2 = ",
            as.character(round(summary(lm1)$adj.r.squared, 2)),
            ", slope = ",
            round(lm1$coefficients[[2]], 2),
            ", trunc.R^2 = ",
            as.character(round(summary(lm2)$adj.r.squared, 2))
        )
    } else {
        title = paste(
            main,
            " scale R^2 = ",
            as.character(round(summary(lm1)$adj.r.squared, 2)),
            ", slope = ",
            round(lm1$coefficients[[2]], 2)
        )
        OUTPUT = data.frame(
            scaleFreeRsquared = round(summary(lm1)$adj.r.squared, 2),
            slope = round(lm1$coefficients[[2]], 2)
        )
    }

    suppressWarnings(plot(
        log.dk,
        log.p.dk,
        xlab = "log10(k)",
        ylab = "log10(p(k))",
        main = title,
        ...
    ))
    lines(log.dk, predict(lm1), col = 1)
    if (truncated)
        lines(log.dk, predict(lm2), col = 2)
    OUTPUT
} # end of function

# plotColorUnderTree ####
#' Plot color rows in a given order, for example under a dendrogram
#'
#' Plot color rows encoding information about objects in a given order, for
#' example the order of a clustering dendrogram, usually below the dendrogram
#' or a barplot.
#'
#' It is often useful to plot dendrograms or other plots (e.g., barplots) of
#' objects together with additional information about the objects, for example
#' module assignment (by color) that was obtained by cutting a hierarchical
#' dendrogram or external color-coded measures such as gene significance. This
#' function provides a way to do so. The calling code should section the screen
#' into two (or more) parts, plot the dendrogram (via \code{plot(hclust)}) or
#' other information in the upper section and use this function to plot color
#' annotation in the order corresponding to the dendrogram in the lower
#' section.
#'
#' @aliases plotColorUnderTree plotOrderedColors
#' @param dendro A hierarchical clustering dendrogram such one returned by
#' \code{\link{hclust}}.
#' @param colors Coloring of objects on the dendrogram. Either a vector (one
#' color per object) or a matrix (can also be an array or a data frame) with
#' each column giving one color per object. Each column will be plotted as a
#' horizontal row of colors under the dendrogram.
#' @param rowLabels Labels for the colorings given in \code{colors}. The labels
#' will be printed to the left of the color rows in the plot. If the argument
#' is given, it must be a vector of length equal to the number of columns in
#' \code{colors}. If not given, \code{names(colors)} will be used if available.
#' If not, sequential numbers starting from 1 will be used.
#' @param rowWidths Optional specification of relative row widths for the color
#' and text (if given) rows. Need not sum to 1.
#' @param rowText Optional labels to identify colors in the color rows.  If
#' given, must be of the same dimensions as \code{colors}. Each label that
#' occurs will be displayed once.
#' @param rowTextAlignment Character string specifying whether the labels
#' should be left-justified to the start of the largest block of each label,
#' centered in the middle, or right-justified to the end of the largest block.
#' @param rowTextIgnore Optional specifications of labels that should be
#' ignored when displaying them using \code{rowText} above.
#' @param textPositions optional numeric vector of the same length as the
#' number of columns in \code{rowText} giving the color rows under which the
#' text rows should appear.
#' @param addTextGuide logical: should guide lines be added for the text rows
#' (if given)?
#' @param cex.rowLabels Font size scale factor for the row labels. See
#' \code{\link[graphics]{par}}.
#' @param cex.rowText character expansion factor for text rows (if given).
#' @param \dots Other parameters to be passed on to the plotting method (such
#' as \code{main} for the main title etc).
#' @return None.
#' @note This function replaces \code{plotHclustColors} in package
#' \code{moduleColor}.
#' @author
#' Steve Horvath \email{SHorvath@@mednet.ucla.edu} and Peter Langfelder
#' \email{Peter.Langfelder@@gmail.com}
#' @seealso
#' \code{\link[dynamicTreeCut]{cutreeDynamic}} for module detection in
#' a dendrogram
#'
#' \code{\link{plotDendroAndColors}} for automated plotting of dendrograms and
#' colors in one step.
#' @keywords hplot
plotColorUnderTree <- function(dendro,
                               colors,
                               rowLabels = NULL,
                               rowWidths = NULL,
                               rowText = NULL,
                               rowTextAlignment = c("left", "center", "right"),
                               rowTextIgnore = NULL,
                               textPositions = NULL,
                               addTextGuide = TRUE,
                               cex.rowLabels = 1,
                               cex.rowText = 0.8,
                               ...) {
    plotOrderedColors(
        dendro$order,
        colors = colors,
        rowLabels = rowLabels,
        rowWidths = rowWidths,
        rowText = rowText,
        rowTextAlignment = rowTextAlignment,
        rowTextIgnore = rowTextIgnore,
        textPositions = textPositions,
        addTextGuide = addTextGuide,
        cex.rowLabels = cex.rowLabels,
        cex.rowText = cex.rowText,
        startAt = 0,
        ...
    )
}


plotOrderedColors <- function(order,
                              colors,
                              rowLabels = NULL,
                              rowWidths = NULL,
                              rowText = NULL,
                              rowTextAlignment = c("left", "center", "right"),
                              rowTextIgnore = NULL,
                              textPositions = NULL,
                              addTextGuide = TRUE,
                              cex.rowLabels = 1,
                              cex.rowText = 0.8,
                              startAt = 0,
                              ...) {
    colors = as.matrix(colors)
    dimC = dim(colors)

    if (is.null(rowLabels) &
        (length(dimnames(colors)[[2]]) == dimC[2]))
        rowLabels = colnames(colors)


    sAF = options("stringsAsFactors")
    options(stringsAsFactors = FALSE)
    on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

    nColorRows = dimC[2]
    if (length(order)  != dimC[1])
        stop("ERROR: length of colors vector not compatible with number of
             objects in 'order'.")
    C = colors[order, , drop = FALSE]
    step = 1 / (dimC[1] - 1 + 2 * startAt)
    #barplot(height = 1, col = "white", border = FALSE, space = 0,
    #axes = FALSE, ...)
    barplot(
        height = 1,
        col = "white",
        border = FALSE,
        space = 0,
        axes = FALSE
    )
    charWidth = strwidth("W") / 2
    if (!is.null(rowText)) {
        if (is.null(textPositions)) {
            textPositions = c(1:nColorRows)
        } else if (is.logical(textPositions)) {
            textPositions = c(1:nColorRows)[textPositions]
        }
        nTextRows = length(textPositions)
    } else
        nTextRows = 0
    nRows = nColorRows + nTextRows
    ystep = 1 / nRows
    if (is.null(rowWidths)) {
        rowWidths = rep(ystep, nColorRows + nTextRows)
    } else {
        if (length(rowWidths) != nRows)
            stop("plotOrderedColors: Length of 'rowWidths' must equal
                 the total number of rows.")
        rowWidths = rowWidths / sum(rowWidths)
    }

    hasText = rep(0, nColorRows)
    hasText[textPositions] = 1
    csPosition = cumsum(c(0, hasText[-nColorRows]))

    colorRows = c(1:nColorRows) + csPosition
    rowType = rep(2, nRows)
    rowType[colorRows] = 1

    physicalTextRow = c(1:nRows)[rowType == 2]

    # Has one extra entry but that shouldn't hurt
    yBottom = c(0, cumsum(rowWidths[nRows:1]))
    yTop = cumsum(rowWidths[nRows:1])

    if (!is.null(rowText)) {
        rowTextAlignment = match.arg(rowTextAlignment)
        rowText = as.matrix(rowText)
        textPos = list()
        textPosY = list()
        textLevs = list()
        for (tr in 1:nTextRows) {
            charHeight = max(strheight(rowText[, tr], cex = cex.rowText))
            width1 = rowWidths[physicalTextRow[tr]]
            nCharFit = floor(width1 / charHeight / 1.7 / par("lheight"))
            if (nCharFit < 1) {
                stop("Rows are too narrow to fit text. Consider decreasing
                     cex.rowText.")
            }
            set = textPositions[tr]
            #colLevs = sort(unique(colors[, set]))
            #textLevs[[tr]] = rowText[match(colLevs, colors[, set]), tr]
            textLevs[[tr]] = sort(unique(rowText[, tr]))
            textLevs[[tr]] = textLevs[[tr]][!textLevs[[tr]] %in% rowTextIgnore]
            nLevs = length(textLevs[[tr]])
            textPos[[tr]] = rep(0, nLevs)
            orderedText = rowText[order, tr]
            for (cl in 1:nLevs) {
                ind = orderedText == textLevs[[tr]][cl]
                sind = ind[-1]
                ind1 = ind[-length(ind)]
                starts = c(if (ind[1])
                    1
                    else
                        NULL, which(!ind1 & sind) + 1)
                ends = which(c(ind1 & !sind, ind[length(ind)]))
                if (length(starts) == 0) {
                    starts = 1
                }
                if (length(ends) == 0) {
                    ends = length(ind)
                }
                if (ends[1] < starts[1]) {
                    starts = c(1, starts)
                }
                if (ends[length(ends)] < starts[length(starts)]) {
                    ends = c(ends, length(ind))
                }
                lengths = ends - starts
                long = which.max(lengths)
                textPos[[tr]][cl] = switch(
                    rowTextAlignment,
                    left = starts[long],
                    center = (starts[long] +
                                  ends[long]) / 2 + 0.5,
                    right = ends[long] + 1
                )
            }
            if (rowTextAlignment == "left") {
                yPos = seq(from = 1,
                           to = nCharFit,
                           by = 1) / (nCharFit + 1)
            } else {
                yPos = seq(from = nCharFit,
                           to = 1,
                           by = -1) / (nCharFit + 1)
            }
            textPosY[[tr]] = rep(yPos, ceiling(nLevs / nCharFit) + 5)[1:nLevs][rank(textPos[[tr]])]
            }
    }

    jIndex = nRows

    if (is.null(rowLabels))
        rowLabels = c(1:nColorRows)
    C[is.na(C)] = "grey"
    for (j in 1:nColorRows) {
        jj = jIndex
        ind = (1:dimC[1])
        xl = (ind - 1.5 + startAt) * step
        xr = (ind - 0.5 + startAt) * step
        yb = rep(yBottom[jj], dimC[1])
        yt = rep(yTop[jj], dimC[1])
        if (is.null(dim(C))) {
            rect(xl,
                 yb,
                 xr,
                 yt,
                 col = as.character(C),
                 border = as.character(C))
        } else {
            rect(xl,
                 yb,
                 xr,
                 yt,
                 col = as.character(C[, j]),
                 border = as.character(C[, j]))
        }
        text(
            rowLabels[j],
            pos = 2,
            x = -charWidth / 2 + xl[1],
            y = (yBottom[jj] + yTop[jj]) / 2,
            cex = cex.rowLabels,
            xpd = TRUE
        )
        textRow = match(j, textPositions)
        if (is.finite(textRow)) {
            jIndex = jIndex - 1
            xt = (textPos[[textRow]] - 1.5) * step

            xt[xt < par("usr")[1]] = par("usr")[1]
            xt[xt > par("usr")[2]] = par("usr")[2]

            #printFlush(paste0("jIndex: ", jIndex, ", yBottom: ",
            #yBottom[jIndex],
            #                  ", yTop: ", yTop[jIndex], ", min(textPosY): ",
            #                  min(textPosY[[textRow]]),
            #                  ", max(textPosY): ", max(textPosY[[textRow]])))
            yt = yBottom[jIndex] + (yTop[jIndex] - yBottom[jIndex]) * (textPosY[[textRow]] + 1 /
                                                                           (2 * nCharFit + 2))
            nt = length(textLevs[[textRow]])
            # Add guide lines
            if (addTextGuide)
                for (l in 1:nt)
                    lines(c(xt[l], xt[l]),
                          c(yt[l], yTop[jIndex]),
                          col = "darkgrey",
                          lty = 3)

            textAdj = c(0, 0.5, 1)[match(rowTextAlignment,
                                         c("left", "center", "right"))]
            text(
                textLevs[[textRow]],
                x = xt,
                y = yt,
                adj = c(textAdj, 1),
                xpd = TRUE,
                cex = cex.rowText
            )
            # printFlush("ok")
        }
        jIndex = jIndex - 1
    }
    for (j in 0:(nColorRows + nTextRows)) {
        lines(x = c(0, 1), y = c(yBottom[j + 1], yBottom[j + 1]))
    }
}

# plotClusterTreeSample ####
# This function can be used to create an average linkage hierarchical
# clustering tree or the microarray samples. The rows of datExpr correspond to
# the samples and the columns to the genes.
# You can optionally input a quantitative microarray sample trait.
#' Annotated clustering dendrogram of microarray samples
#'
#' This function plots an annotated clustering dendorgram of microarray
#' samples.
#'
#' The function generates an average linkage hierarchical clustering dendrogram
#' (see \code{\link[stats]{hclust}}) of samples from the given expression data,
#' using Eclidean distance of samples. The dendrogram is plotted together with
#' color annotation for the samples.
#'
#' The trait \code{y} must be numeric. If \code{y} is integer, the colors will
#' correspond to values. If \code{y} is continouos, it will be dichotomized to
#' two classes, below and above median.
#'
#' @param datExpr a data frame containing expression data, with rows
#' corresponding to samples and columns to genes. Missing values are allowed
#' and will be ignored.
#' @param y microarray sample trait. Either a vector with one entry per sample,
#' or a matrix in which each column corresponds to a (different) trait and each
#' row to a sample.
#' @param traitLabels labels to be printed next to the color rows depicting
#' sample traits. Defaults to column names of \code{y}.
#' @param yLabels Optional labels to identify colors in the row identifying the
#' sample classes. If given, must be of the same dimensions as \code{y}. Each
#' label that occurs will be displayed once.
#' @param main title for the plot.
#' @param setLayout logical: should the plotting device be partitioned into a
#' standard layout?  If \code{FALSE}, the user is responsible for partitioning.
#' The function expects two regions of the same width, the first one
#' immediately above the second one.
#' @param autoColorHeight logical: should the height of the color area below
#' the dendrogram be automatically adjusted for the number of traits? Only
#' effective if \code{setLayout} is \code{TRUE}.
#' @param colorHeight Specifies the height of the color area under dendrogram
#' as a fraction of the height of the dendrogram area. Only effective when
#' \code{autoColorHeight} above is \code{FALSE}.
#' @param dendroLabels dendrogram labels. Set to \code{FALSE} to disable
#' dendrogram labels altogether; set to \code{NULL} to use row labels of
#' \code{datExpr}.
#' @param addGuide logical: should vertical "guide lines" be added to the
#' dendrogram plot? The lines make it easier to identify color codes with
#' individual samples.
#' @param guideAll logical: add a guide line for every sample? Only effective
#' for \code{addGuide} set \code{TRUE}.
#' @param guideCount number of guide lines to be plotted. Only effective when
#' \code{addGuide} is \code{TRUE} and \code{guideAll} is \code{FALSE}.
#' @param guideHang fraction of the dendrogram height to leave between the top
#' end of the guide line and the dendrogram merge height. If the guide lines
#' overlap with dendrogram labels, increase \code{guideHang} to leave more
#' space for the labels.
#' @param cex.traitLabels character expansion factor for trait labels.
#' @param cex.dendroLabels character expansion factor for dendrogram (sample)
#' labels.
#' @param marAll a 4-element vector giving the bottom, left, top and right
#' margins around the combined plot. Note that this is not the same as setting
#' the margins via a call to \code{\link{par}}, because the bottom margin of
#' the dendrogram and the top margin of the color underneath are always zero.
#' @param saveMar logical: save margins setting before starting the plot and
#' restore on exit?
#' @param abHeight optional specification of the height for a horizontal line
#' in the dendrogram, see \code{\link{abline}}.
#' @param abCol color for plotting the horizontal line.
#' @param \dots other graphical parameters to \code{\link{plot.hclust}}.
#' @return None.
#' @author Steve Horvath and Peter Langfelder
#' @seealso
#' \code{\link[stats]{dist}}, \code{\link[stats]{hclust}},
#' \code{\link{plotDendroAndColors}}
#' @keywords hplot misc
plotClusterTreeSamples <- function(datExpr,
                                   y = NULL,
                                   traitLabels = NULL,
                                   yLabels = NULL,
                                   main = if (is.null(y)) {
                                       "Sample dendrogram"
                                   } else {
                                       "Sample dendrogram and trait indicator"
                                   },
                                   setLayout = TRUE,
                                   autoColorHeight = TRUE,
                                   colorHeight = 0.3,
                                   dendroLabels = NULL,
                                   addGuide = FALSE,
                                   guideAll = TRUE,
                                   guideCount = NULL,
                                   guideHang = 0.20,
                                   cex.traitLabels = 0.8,
                                   cex.dendroLabels = 0.9,
                                   marAll = c(1, 5, 3, 1),
                                   saveMar = TRUE,
                                   abHeight = NULL,
                                   abCol = "red",
                                   ...) {
    dendro = fastcluster::hclust(dist(datExpr), method = "average")
    if (is.null(y)) {
        oldMar = par("mar")
        par(mar = marAll)
        plot(dendro, main = main, sub = "", xlab = "",
             labels = dendroLabels, cex = cex.dendroLabels)
        if (saveMar)
            par(oldMar)
    } else {
        if (is.null(traitLabels))
            traitLabels = names(as.data.frame(y))
        y = as.matrix(y)
        if (!is.numeric(y)) {
            warning("The microarray sample trait y will be transformed to ",
                    "numeric.")
            dimy = dim(y)
            y = as.numeric(y)
            dim(y) = dimy
        } # end of if (!is.numeric(y))
        if (nrow(as.matrix(datExpr))  != nrow(y)) {
            stop("dim(as.matrix(datExpr))[[1]]  != length(y)\n",
                 "The number of microarray sample arrays does
                 not match the number of samples for the trait.\n")
        }

        if (is.integer(y)) {
            y = y - min(0, min(y, na.rm = TRUE)) + 1
        } else {
            y = (y >= median(y, na.rm = TRUE)) + 1
        }
        plotDendroAndColors(dendro, colors = y, groupLabels = traitLabels,
                            rowText = yLabels, setLayout = setLayout,
                            autoColorHeight = autoColorHeight,
                            colorHeight = colorHeight, addGuide = addGuide,
                            guideAll = guideAll, guideCount = guideCount,
                            guideHang = guideHang,
                            cex.colorLabels = cex.traitLabels,
                            cex.dendroLabels = cex.dendroLabels,
                            marAll = marAll, saveMar = saveMar,
                            abHeight = abHeight, abCol = abCol,
                            main = main, ...)
        }
}

# TOMplot ####
#' Graphical representation of the Topological Overlap Matrix
#'
#' Graphical representation of the Topological Overlap Matrix using a heatmap
#' plot combined with the corresponding hierarchical clustering dendrogram and
#' module colors.
#'
#' The standard \code{heatmap} function uses the \code{\link{layout}} function
#' to set the following layout (when \code{Colors} is given): \preformatted{ 0
#' 0 5 0 0 2 4 1 3 } To get a meaningful heatmap plot, user-set layout must
#' respect this geometry.
#'
#' @param dissim a matrix containing the topological overlap-based
#' dissimilarity
#' @param dendro the corresponding hierarchical clustering dendrogram
#' @param Colors optional specification of module colors to be plotted on top
#' @param ColorsLeft optional specification of module colors on the left side.
#' If \code{NULL}, \code{Colors} will be used.
#' @param terrainColors logical: should terrain colors be used?
#' @param setLayout logical: should layout be set? If \code{TRUE}, standard
#' layout for one plot will be used. Note that this precludes multiple plots on
#' one page. If \code{FALSE}, the user is responsible for setting the correct
#' layout.
#' @param \dots other graphical parameters to \code{\link{heatmap}}.
#' @return None.
#' @author Steve Horvath and Peter Langfelder
#' @seealso \code{\link{heatmap}}, the workhorse function doing the plotting.
#' @keywords misc
TOMplot <- function(dissim, dendro, Colors = NULL, ColorsLeft = Colors,
                    terrainColors = FALSE, setLayout = TRUE, ...) {
    if (is.null(Colors))
        Colors = rep("white", dim(as.matrix(dissim))[[1]])
    if (is.null(ColorsLeft))
        ColorsLeft = Colors
    nNodes = length(Colors)
    if (nNodes < 2) {
        warning("You have only 1 or 2 genes in TOMplot.
                No plot will be produced")
    } else {
        if (nNodes  != length(ColorsLeft))
            stop("ERROR: number of (top) color labels does not equal number of
                 left color labels")
        if (nNodes  != dim(dissim)[[1]])
            stop(
                paste(
                    "ERROR: number of color labels does not equal number of
                    nodes in dissim.\n",
                    "     nNodes  != dim(dissim)[[1]] "
                )
            )
        labeltree = as.character(Colors)
        labelrow  = as.character(ColorsLeft)
        #labelrow[dendro$order[length(labeltree):1]] = labelrow[dendro$order]
        options(expressions = 10000)
        dendro$height = (dendro$height - min(dendro$height)) / (1.15 * (max(dendro$height) - min(dendro$height)))
        if (terrainColors) {
            .heatmap(
                as.matrix(dissim),
                Rowv = as.dendrogram(dendro, hang = 0.1),
                Colv = as.dendrogram(dendro, hang = 0.1),
                scale = "none",
                revC = TRUE,
                ColSideColors = as.character(labeltree),
                RowSideColors = as.character(labelrow),
                labRow = FALSE,
                labCol = FALSE,
                col = terrain.colors(100),
                setLayout = setLayout,
                ...
            )
        } else {
            .heatmap(
                as.matrix(dissim),
                Rowv = as.dendrogram(dendro, hang = 0.1),
                Colv = as.dendrogram(dendro, hang = 0.1),
                scale = "none",
                revC = TRUE,
                ColSideColors = as.character(labeltree),
                RowSideColors = as.character(labelrow),
                labRow = FALSE,
                labCol = FALSE,
                setLayout = setLayout,
                ...
            )
        } #end of if
    }
}

# plotNetworkHeatmap ####
#' Network heatmap plot
#'
#' Network heatmap plot.
#'
#' The function constructs a network from the given expression data (selected
#' by \code{plotGenes}) using the soft-thresholding procedure, optionally
#' calculates Topological Overlap (TOM) and plots a heatmap of the network.
#'
#' Note that all network calculations are done in one block and may fail due to
#' memory allocation issues for large numbers of genes.
#'
#' @param datExpr a data frame containing expression data, with rows
#' corresponding to samples and columns to genes. Missing values are allowed
#' and will be ignored.
#' @param plotGenes a character vector giving the names of genes to be included
#' in the plot. The names will be matched against \code{names(datExpr)}.
#' @param useTOM logical: should TOM be plotted (\code{TRUE}), or
#' correlation-based adjacency (\code{FALSE})?
#' @param power soft-thresholding power for network construction.
#' @param networkType a character string giving the newtork type. Recognized
#' values are (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param main main title for the plot.
#' @return None.
#' @author Steve Horvath
#' @seealso \code{\link{adjacency}}, \code{\link{TOMsimilarity}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords hplot
plotNetworkHeatmap <- function(datExpr, plotGenes, useTOM = TRUE, power = 6 ,
                               networkType = "unsigned",
                               main = "Heatmap of the network") {
    match1 = match(plotGenes, names(data.frame(datExpr)))
    match1 = match1[!is.na(match1)]
    nGenes = length(match1)
    if (sum(!is.na(match1)) != length(plotGenes)) {
        printFlush(
            paste(
                "Warning: Not all gene names were recognized.",
                "Only the following genes were recognized."
            )
        )
        printFlush(paste("   ", names(data.frame(datExpr))[match1],
                         collapse = ", "))
    }
    if (nGenes < 3) {
        warning(
            paste(
                "Since you have fewer than 3 genes, the network will not be
                visualized.\n",
                "   Hint: please input more genes."
            )
        )
        plot(1, 1)
    } else {
        datErest = datExpr[, match1]
        ADJ1 = adjacency(datErest, power = power, type = networkType)
        if (useTOM) {
            diss1 = 1 - TOMsimilarity(ADJ1)
        } else {
            diss1 = 1 - ADJ1
        }
        diag(diss1) = NA
        hier1 = fastcluster::hclust(as.dist(diss1), method = "average")
        colors1 = rep("white", nGenes)
        labeltree = names(data.frame(datErest))
        labelrow  = names(data.frame(datErest))
        labelrow[hier1$order[length(labeltree):1]] = labelrow[hier1$order]
        options(expressions = 10000)
        heatmap(
            as.matrix(diss1),
            Rowv = as.dendrogram(hier1),
            Colv = as.dendrogram(hier1),
            scale = "none",
            revC = TRUE,
            labRow = labeltree,
            labCol = labeltree,
            main = main
        )
    } # end of if (nGenes >  2)
} # end of function

# plotModuleSigninficance ####
#' Barplot of module significance
#'
#' Plot a barplot of gene significance.
#'
#' Given individual gene significances and their module assigment, the function
#' calculates the module significance for each module as the average gene
#' significance of the genes within the module. The result is plotted in a
#' barplot or boxplot form. Each bar or box is labeled by the corresponding
#' module color.
#'
#' @param geneSignificance a numeric vector giving gene significances.
#' @param colors a character vector specifying module assignment for the genes
#' whose significance is given in \code{geneSignificance }. The modules should
#' be labeled by colors.
#' @param boxplot logical: should a boxplot be produced instead of a barplot?
#' @param main main title for the plot.
#' @param ylab y axis label for the plot.
#' @param \dots other graphical parameters to \code{\link{plot}}.
#' @return None.
#' @author Steve Horvath
#' @seealso \code{\link{barplot}}, \code{\link{boxplot}}
#' @references
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24
#' @keywords hplot misc
plotModuleSignificance <- function(geneSignificance, colors, boxplot = FALSE,
                                   main = "Gene significance across modules, ",
                                   ylab = "Gene Significance", ...) {
    if (length(geneSignificance)  != length(colors))
        stop("Error: 'geneSignificance' and 'colors' do not have the same lengths")
    no.colors = length(names(table(colors)))
    if (no.colors == 1)
        pp = NA
    if (no.colors > 1)
    {
        pp = try(kruskal.test(geneSignificance, factor(colors))$p.value)
        if (class(pp) == 'try - error')
            pp = NA
    }
    title = paste0(main, " p-value = ", signif (pp, 2))
    if (boxplot  != TRUE) {
        means1 = as.vector(tapply(geneSignificance, colors, mean, na.rm = TRUE))
        se1 = as.vector(tapply(geneSignificance, colors, stdErr))
        # par(mfrow = c(1, 1))
        barplot(
            means1,
            names.arg = names(table(colors)),
            col = names(table(colors)),
            ylab = ylab,
            main = title,
            ...
        )
        addErrorBars(as.vector(means1), as.vector(1.96 * se1), two.side = TRUE)
    } else {
        boxplot(
            split(geneSignificance, colors),
            notch = TRUE,
            varwidth = TRUE,
            col = names(table(colors)),
            ylab = ylab,
            main = title,
            ...
        )
    }
} # end of function

# plotEigengeneNetworks ####
#' Eigengene network plot
#'
#' This function plots dendrogram and eigengene representations of (consensus)
#' eigengenes networks.  In the case of conensus eigengene networks the
#' function also plots pairwise preservation measures between consensus
#' networks in different sets.
#'
#'
#' Consensus eigengene networks consist of a fixed set of eigengenes
#' "expressed" in several different sets. Network connection strengths are
#' given by eigengene correlations. This function aims to visualize the
#' networks as well as their similarities and differences across sets.
#'
#' The function partitions the screen appropriately and plots eigengene
#' dendrograms in the top row, then a square matrix of plots: heatmap plots of
#' eigengene networks in each set on the diagonal, heatmap plots of pairwise
#' preservation networks below the diagonal, and barplots of aggregate network
#' preservation of individual eigengenes above the diagonal. A preservation
#' plot or barplot in the row i and column j of the square matrix represents
#' the preservation between sets i and j.
#'
#' Individual eigengenes are labeled by their name in the dendrograms; in the
#' heatmaps and barplots they can optionally be labeled by color squares. For
#' compatibility with other functions, the color labels are encoded in the
#' eigengene names by prefixing the color with two letters, such as
#' \code{"MEturquoise"}.
#'
#' Two types of network preservation can be plotted: the \code{"standard"} is
#' simply the difference between adjacencies in the two compared sets. The
#' \code{"hyperbolic"} difference de-emphasizes the preservation of low
#' adjacencies. When \code{"both"} is specified, standard preservation is
#' plotted in the lower triangle and hyperbolic in the upper triangle of each
#' preservation heatmap.
#'
#' If the eigengenes are labeled by color, the bars in the barplot can be split
#' into segments representing the contribution of each eigengene and labeled by
#' the contribution. For example, a yellow segment in a bar labeled by a
#' turquoise square represents the preservation of the adjacency between the
#' yellow and turquoise eigengenes in the two networks compared by the barplot.
#'
#' For large numbers of eigengenes and/or sets, it may be difficult to get a
#' meaningful plot fit a standard computer screen. In such cases we recommend
#' using a device such as \code{\link{postscript}} or \code{\link{pdf}} where
#' the user can specify large dimensions; such plots can be conveniently viewed
#' in standard pdf or postscript viewers.
#'
#' @param multiME either a single data frame containing the module eigengenes,
#' or module eigengenes in the multi-set format (see \code{\link{checkSets}}).
#' The multi-set format is a vector of lists, one per set. Each set must
#' contain a component \code{data} whose rows correspond to samples and columns
#' to eigengenes.
#' @param setLabels A vector of character strings that label sets in
#' \code{multiME}.
#' @param letterSubPlots logical: should subplots be lettered?
#' @param Letters optional specification of a sequence of letters for
#' lettering. Defaults to "ABCD"...
#' @param excludeGrey logical: should the grey module eigengene be excluded
#' from the plots?
#' @param greyLabel label for the grey module. Usually either "grey" or the
#' number 0.
#' @param plotDendrograms logical: should eigengene dendrograms be plotted?
#' @param plotHeatmaps logical: should eigengene network heatmaps be plotted?
#' @param setMargins logical: should margins be set? See
#' \code{\link[graphics]{par}}.
#' @param marDendro a vector of length 4 giving the margin setting for
#' dendrogram plots. See \code{\link[graphics]{par}}. If \code{setMargins} is
#' \code{TRUE} and \code{marDendro} is not given, the function will provide
#' reasonable default values.
#' @param marHeatmap a vector of length 4 giving the margin setting for heatmap
#' plots. See \code{\link[graphics]{par}}. If \code{setMargins} is \code{TRUE}
#' and \code{marDendro} is not given, the function will provide reasonable
#' default values.
#' @param colorLabels logical: should module eigengene names be interpreted as
#' color names and the colors used to label heatmap plots and barplots?
#' @param signed logical: should eigengene networks be constructed as signed?
#' @param heatmapColors color palette for heatmaps. Defaults to
#' \code{\link{heat.colors}} when \code{signed} is \code{FALSE}, and to
#' \code{\link{redWhiteGreen}} when \code{signed} is \code{TRUE}.
#' @param plotAdjacency logical: should module eigengene heatmaps plot
#' adjacency (ranging from 0 to 1), or correlation (ranging from -1 to 1)?
#' @param printAdjacency logical: should the numerical values be printed into
#' the adjacency or correlation heatmap?
#' @param cex.adjacency character expansion factor for printing of numerical
#' values into the adjacency or correlation heatmap
#' @param coloredBarplot logical: should the barplot of eigengene adjacency
#' preservation distinguish individual contributions by color? This is possible
#' only if \code{colorLabels} is \code{TRUE} and module eigengene names encode
#' valid colors.
#' @param barplotMeans logical: plot mean preservation in the barplot? This
#' option effectively rescales the preservation by the number of eigengenes in
#' the network. If means are plotted, the barplot is not colored.
#' @param barplotErrors logical: should standard errors of the mean
#' preservation be plotted?
#' @param plotPreservation a character string specifying which type of
#' preservation measure to plot. Allowed values are (unique abbreviations of)
#' \code{"standard"}, \code{"hyperbolic"}, \code{"both"}.
#' @param zlimPreservation a vector of length 2 giving the value limits for the
#' preservation heatmaps.
#' @param printPreservation logical: should preservation values be printed
#' within the heatmap?
#' @param cex.preservation character expansion factor for preservation display.
#' @param \dots other graphical arguments to function
#' \code{\link{labeledHeatmap}}.
#' @return None.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{labeledHeatmap}}, \code{\link{labeledBarplot}} for annotated
#' heatmaps and barplots
#'
#' \code{\link[stats]{hclust}} for hierarchical clustering and dendrogram plots
#' @references
#'
#' For theory and applications of consensus eigengene networks, see
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co-expression modules. BMC Systems Biology 2007, 1:54
#' @keywords hplot
plotEigengeneNetworks <- function(multiME,
                                  setLabels,
                                  letterSubPlots = FALSE,
                                  Letters = NULL,
                                  excludeGrey = TRUE,
                                  greyLabel = "grey",
                                  plotDendrograms = TRUE,
                                  plotHeatmaps = TRUE,
                                  setMargins = TRUE,
                                  marDendro = NULL,
                                  marHeatmap = NULL,
                                  colorLabels = TRUE,
                                  signed = TRUE,
                                  heatmapColors = NULL,
                                  plotAdjacency = TRUE,
                                  printAdjacency = FALSE,
                                  cex.adjacency = 0.9,
                                  coloredBarplot = TRUE,
                                  barplotMeans = TRUE,
                                  barplotErrors = FALSE,
                                  plotPreservation = "standard",
                                  zlimPreservation = c(0, 1),
                                  printPreservation = FALSE,
                                  cex.preservation = 0.9,
                                  ...) {
    # invertColors = FALSE
    size = checkSets(multiME, checkStructure = TRUE)
    if (!size$structureOK) {
        #printFlush(paste(
        #  "plotEigengeneNetworks: Given multiME does not appear to be a multi -
        #  set structure.\n",
        #  "Will attempt to convert it into a multi - set structure containing
        #  1 set."))
        multiME = fixDataStructure(multiME)
    }

    if (is.null(Letters))
        Letters = LETTERS

    if (is.null(heatmapColors))
        if (signed) {
            heatmapColors = blueWhiteRed(50)
        } else {
            heatmapColors = heat.colors(30)
        }
    nSets = length(multiME)
    cex = par("cex")
    mar = par("mar")
    nPlotCols = nSets
    nPlotRows = as.numeric(plotDendrograms) + nSets * as.numeric(plotHeatmaps)
    if (nPlotRows == 0)
        stop("Nothing to plot: neither dendrograms not heatmaps requested.")
    oldpar <- par(mfrow = c(nPlotRows, nPlotCols))
    on.exit(par(oldpar))
    par(cex = cex)
    if (excludeGrey)
        for (set in 1:nSets)
            multiME[[set]]$data  =
        multiME[[set]]$data[, substring(names(multiME[[set]]$data), 3) != greyLabel]

    plotPresTypes = c("standard", "hyperbolic", "both")
    ipp = pmatch(plotPreservation, plotPresTypes)
    if (is.na(ipp))
        stop(paste(
            "Invalid 'plotPreservation'. Available choices are",
            paste(plotPresTypes, sep = ", ")
        ))

    letter.ind = 1
    if (plotDendrograms)
        for (set in 1:nSets) {
            #par(cex = StandardCex/1.4)
            par(mar = marDendro)
            labels = names(multiME[[set]]$data)
            uselabels = labels[substring(labels, 3) != greyLabel]
            corME = cor(multiME[[set]]$data[substring(labels, 3) != greyLabel,
                                            substring(labels, 3) != greyLabel],
                        use = "p")
            disME = as.dist(1 - corME)
            clust = fastcluster::hclust(disME, method = "average")
            if (letterSubPlots) {
                main = paste0(substring(Letters, letter.ind, letter.ind),
                              ".",
                              setLabels[set])
            } else {
                main = setLabels[set]
            }
            #validColors = is.na(match(uselabels, colors()))
            #plotLabels = ifelse(validColors, substring(uselabels[validColors], 3),
            #uselabels[!validColors])
            plotLabels = uselabels
            plot(
                clust,
                main = main,
                sub = "",
                xlab = "",
                labels = plotLabels,
                ylab = "",
                ylim = c(0, 1)
            )
            letter.ind = letter.ind + 1
        }

    if (plotHeatmaps)
        for (i.row in (1:nSets))
            for (i.col in (1:nSets)) {
                letter.ind = i.row * nSets + i.col
                if (letterSubPlots) {
                    #letter = paste("(", substring(Letters, first = letter.ind,
                    #last = letter.ind), ")", sep = "")
                    letter = paste0(substring(Letters, first = letter.ind,
                                              last = letter.ind),
                                    ".  ")
                } else {
                    letter = NULL
                }
                par(cex = cex)
                if (setMargins) {
                    if (is.null(marHeatmap)) {
                        if (colorLabels) {
                            par(mar = c(1, 2, 3, 4) + 0.2)
                        } else {
                            par(mar = c(6, 7, 3, 5) + 0.2)
                        }
                    } else {
                        par(mar = marHeatmap)
                    }
                }
                nModules = dim(multiME[[i.col]]$data)[2]
                textMat = NULL
                if (i.row == i.col) {
                    corME = cor(multiME[[i.col]]$data, use = "p")
                    pME = corPvalueFisher(corME, nrow(multiME[[i.col]]$data))
                    if (printAdjacency) {
                        textMat = paste(signif (corME, 2), "\n", signif (pME, 1))
                        dim(textMat) = dim(corME)
                    }
                    if (signed) {
                        if (plotAdjacency) {
                            if (printAdjacency) {
                                textMat = paste(signif ((1 + corME) / 2, 2), "\n",
                                                signif (pME, 1))
                                dim(textMat) = dim(corME)
                            }

                            labeledHeatmap((1 + corME) / 2,
                                           names(multiME[[i.col]]$data),
                                           names(multiME[[i.col]]$data),
                                           main = paste(letter, setLabels[[i.col]]),
                                           invertColors = FALSE,
                                           zlim = c(0, 1.0),
                                           colorLabels = colorLabels,
                                           colors = heatmapColors,
                                           setStdMargins = FALSE,
                                           textMatrix = textMat,
                                           cex.text = cex.adjacency,
                                           ...
                            )
                        } else {
                            labeledHeatmap(
                                corME,
                                names(multiME[[i.col]]$data),
                                names(multiME[[i.col]]$data),
                                main = paste(letter, setLabels[[i.col]]),
                                invertColors = FALSE,
                                zlim = c(-1, 1.0),
                                colorLabels = colorLabels,
                                colors = heatmapColors,
                                setStdMargins = FALSE,
                                textMatrix = textMat,
                                cex.text = cex.adjacency,
                                ...
                            )
                        }
                    } else {
                        labeledHeatmap(
                            abs(corME),
                            names(multiME[[i.col]]$data),
                            names(multiME[[i.col]]$data),
                            main = paste(letter, setLabels[[i.col]]),
                            invertColors = FALSE,
                            zlim = c(0, 1.0),
                            colorLabels = colorLabels,
                            colors = heatmapColors,
                            setStdMargins = FALSE,
                            textMatrix = textMat,
                            cex.text = cex.adjacency,
                            ...
                        )
                    }
                } else {
                    corME1 = cor(multiME[[i.col]]$data, use = "p")
                    corME2 = cor(multiME[[i.row]]$data, use = "p")
                    cor.dif = (corME1 - corME2) / 2
                    d = tanh((corME1 - corME2) / (abs(corME1) + abs(corME2)) ^ 2)
                    # d = abs(corME1 - corME2) / (abs(corME1) + abs(corME2))
                    if (ipp == 1 | ipp == 3) {
                        dispd = cor.dif
                        main = paste(letter, "Preservation")
                        if (ipp == 3) {
                            dispd[upper.tri(d)] = d[upper.tri(d)]
                            main = paste(
                                letter,
                                "Hyperbolic preservation (UT)\nStandard preservation
                                (LT)"
                            )
                        }
                        } else {
                            dispd = d
                            main = paste(letter, "Hyperbolic preservation")
                    }
                    if (i.row > i.col) {
                        if (signed) {
                            half = as.integer(length(heatmapColors) / 2)
                            range = c(half:length(heatmapColors))
                            halfColors = heatmapColors[range]
                        } else {
                            halfColors = heatmapColors
                        }
                        if (printPreservation) {
                            printMtx = matrix(
                                paste0(".", as.integer((
                                    1 - abs(dispd)
                                ) * 100)),
                                nrow = nrow(dispd),
                                ncol = ncol(dispd)
                            )
                            printMtx[printMtx == ".100"] = "1"
                        } else {
                            printMtx = NULL
                        }
                        if ((sum((1 - abs(
                            dispd
                        )) < zlimPreservation[1]) || ((1 - abs(
                            dispd
                        )) > zlimPreservation[2]))  > 0)
                            warning(
                                "plotEigengeneNetworks: Correlation preservation data
                                out of zlim range."
                            )
                        labeledHeatmap(
                            1 - abs(dispd),
                            names(multiME[[i.col]]$data),
                            names(multiME[[i.col]]$data),
                            main = main,
                            invertColors = FALSE,
                            colorLabels = colorLabels,
                            zlim = zlimPreservation,
                            colors = halfColors,
                            setStdMargins = FALSE,
                            textMatrix = printMtx,
                            cex.text = cex.preservation,
                            ...
                        )
                    } else {
                        if (ipp == 2) {
                            dp = 1 - abs(d)
                            method = "Hyperbolic:"
                        } else {
                            dp = 1 - abs(cor.dif)
                            method = "Preservation:"
                        }
                        diag(dp) = 0
                        if (barplotMeans) {
                            sum_dp = mean(dp[upper.tri(dp)])
                            means = apply(dp, 2, sum) / (ncol(dp) - 1)
                            if (barplotErrors) {
                                errors = sqrt((
                                    apply(dp ^ 2, 2, sum) / (ncol(dp) - 1) - means ^ 2
                                ) / (ncol(dp) - 2))
                            } else {
                                errors = NULL
                            }
                            labeledBarplot(
                                means,
                                names(multiME[[i.col]]$data),
                                main = paste(letter, "D = ", signif (sum_dp, 2)),
                                ylim = c(0, 1),
                                colorLabels = colorLabels,
                                colored = coloredBarplot,
                                setStdMargins = FALSE,
                                stdErrors = errors,
                                ...
                            )
                        } else {
                            sum_dp = sum(dp[upper.tri(dp)])
                            labeledBarplot(
                                dp,
                                names(multiME[[i.col]]$data),
                                main = paste(letter, method, "sum = ", signif (sum_dp, 3)),
                                ylim = c(0, dim(dp)[[1]]),
                                colorLabels = colorLabels,
                                colored = coloredBarplot,
                                setStdMargins = FALSE,
                                ...
                            )
                        }
                    }
                }
            }
            }

# addGrid ####
#' Add grid lines to an existing plot.
#'
#' This function adds horizontal and/or vertical grid lines to an existing
#' plot. The grid lines are aligned with tick marks.
#'
#' If \code{linesPerTick} is not specified, it is set to 5 if number of tick s
#' is 5 or less, and it is set to 2 if number of ticks is greater than 5.
#'
#' @param linesPerTick Number of lines between successive tick marks (including
#' the line on the tickmarks themselves)
#' @param horiz Draw horizontal grid lines?
#' @param vert Draw vertical tick lines?
#' @param col Specifies color of the grid lines
#' @param lty Specifies line type of grid lines. See \code{\link{par}}.
#' @note The function does not work whenever logarithmic scales are in use.
#' @author Peter Langfelder
#' @keywords hplot
#' @examples
#'
#'   plot(c(1:10), c(1:10))
#'   addGrid()
#'
addGrid <- function(linesPerTick = NULL, horiz = TRUE, vert = FALSE,
                    col = "grey30", lty = 3) {
    box = par("usr")
    if (horiz) {
        ticks = par("yaxp")
        nTicks = ticks[3]
        if (is.null(linesPerTick)) {
            linesPerTick <- ifelse(nTicks < 6, 5, 2)
        }
        spacing = (ticks[2] - ticks[1]) / (linesPerTick * nTicks)
        first = ceiling((box[3] - ticks[1]) / spacing)
        last = floor((box[4] - ticks[1]) / spacing)
        #print(paste("addGrid: first = ", first, ", last  = ", last, "box = ",
        #paste(signif (box, 2), collapse = ", "),
        #"ticks = ", paste(signif (ticks, 2), collapse = ", "),
        #"spacing  = ", spacing))
        for (k in first:last) {
            lines(x = box[c(1, 2)],
                  y = rep(ticks[1] + spacing * k, 2),
                  col = col,
                  lty = lty)
        }
    }
    if (vert) {
        ticks = par("xaxp")
        nTicks = ticks[3]
        if (is.null(linesPerTick)) {
            linesPerTick <- ifelse(nTicks < 6, 5, 2)
        }
        spacing = (ticks[2] - ticks[1]) / (linesPerTick * ticks[3])
        first = ceiling((box[1] - ticks[1]) / spacing)
        last = floor((box[2] - ticks[1]) / spacing)
        #print(paste("addGrid: first = ", first, ", last  = ", last, "box = ",
        #paste(signif (box, 2), collapse = ", "),
        #            "ticks = ", paste(signif (ticks, 2), collapse = ", "),
        #            "spacing  = ", spacing))
        for (l in first:last) {
            lines(x = rep(ticks[1] + spacing * l, 2),
                  y = box[c(3, 4)],
                  col = col,
                  lty = lty)
        }
    }
}

# addGuidesLines ####
#' Add vertical ``guide lines'' to a dendrogram plot
#'
#' Adds vertical ``guide lines'' to a dendrogram plot.
#'
#'
#' @param dendro The dendrogram (see \code{\link{hclust}}) to which the guide
#' lines are to be added.
#' @param all Add a guide line to every object on the dendrogram? Useful if the
#' number of objects is relatively low.
#' @param count Number of guide lines to be plotted. The lines will be
#' equidistantly spaced.
#' @param positions Horizontal positions of the added guide lines. If given,
#' overrides \code{count}.
#' @param col Color of the guide lines
#' @param lty Line type of the guide lines. See \code{\link{par}}.
#' @param hang Fraction of the figure height that will separate top ends of
#' guide lines and the merge heights of the corresponding objects.
#' @author Peter Langfelder
#' @keywords hplot
addGuideLines <- function(dendro, all = FALSE, count = 50, positions = NULL,
                          col = "grey30", lty = 3, hang = 0) {
    if (all) {
        positions = 1:(length(dendro$height) + 1)
    } else {
        if (is.null(positions)) {
            lineSpacing = (length(dendro$height) + 1) / count
            positions = (1:count) *  lineSpacing
        }
    }
    objHeights = rep(0, length(dendro$height + 1))
    objHeights[-dendro$merge[dendro$merge[, 1] < 0, 1]] <- dendro$height[
        dendro$merge[, 1] < 0]
    objHeights[-dendro$merge[dendro$merge[, 2] <  0, 2]] <- dendro$height[
        dendro$merge[, 2] < 0]
    box = par("usr")
    ymin = box[3]
    ymax = box[4]
    objHeights = objHeights - hang * (ymax - ymin)
    objHeights[objHeights < ymin] = ymin
    posHeights = pmin(objHeights[dendro$order][floor(positions)],
                      objHeights[dendro$order][ceiling(positions)])
    for (line in 1:length(positions)) {
        # The last guide line is superfluous
        lines(x = rep(positions[line], 2),
              y = c(ymin, posHeights[line]),
              lty = 3,
              col = col
        )
    }
}

# sizeGrWindow ####
#' Opens a graphics window with specified dimensions
#'
#' If a graphic device window is already open, it is closed and re-opened with
#' specified dimensions (in inches); otherwise a new window is opened.
#'
#'
#' @param width desired width of the window, in inches.
#' @param height desired heigh of the window, in inches.
#' @return None.
#' @author Peter Langfelder
#' @keywords misc
sizeGrWindow <- function(width, height) {
    din = par("din")
    if ((din[1] != width) | (din[2] != height)) {
        dev.off()
        dev.new(width = width, height = height)
    }
}

# scaleFreeFitIndex ####
#' Calculation of fitting statistics for evaluating scale free topology fit.
#'
#' The function scaleFreeFitIndex calculates several indices (fitting
#' statistics) for evaluating scale free topology fit.  The input is a vector
#' (of connectivities) k. Next k is discretized into nBreaks number of
#' equal-width bins.  Let's denote the resulting vector dk.  The relative
#' frequency for each bin is denoted p.dk.
#'
#'
#' @param k numeric vector whose components contain non-negative values
#' @param nBreaks positive integer. This determines the number of equal width
#' bins.
#' @param removeFirst logical. If TRUE then the first bin will be removed.
#' @return Data frame with columns \item{Rsquared.SFT}{the model fitting index
#' (R.squared) from the following model lm(log.p.dk ~ log.dk)}
#' \item{slope.SFT}{the slope estimate from model lm(log(p(k))~log(k))}
#' \item{truncatedExponentialAdjRsquared}{the adjusted R.squared measure from
#' the truncated exponential model given by lm2 = lm(log.p.dk ~ log.dk + dk).}
#' @author Steve Horvath
#' @keywords misc
scaleFreeFitIndex <- function(k, nBreaks = 10, removeFirst = FALSE) {
    # TODO: What does this do?
    discretized.k <- cut(k, nBreaks)
    dk <- tapply(k, discretized.k, mean)
    p.dk <- as.vector(tapply(k, discretized.k, length) / length(k))
    breaks1 <- seq(from = min(k),
                   to = max(k),
                   length = nBreaks + 1)
    hist1 <- hist(k,
                  breaks = breaks1,
                  plot = FALSE,
                  right = TRUE)
    dk2 <- hist1$mids
    dk <- ifelse(is.na(dk), dk2, dk)
    dk <- ifelse(dk == 0, dk2, dk)
    p.dk <- ifelse(is.na(p.dk), 0, p.dk)
    log.dk <- as.vector(log10(dk))
    if (removeFirst) {
        p.dk <- p.dk[-1]
        log.dk <- log.dk[-1]
    }
    log.p.dk <- as.numeric(log10(p.dk + 1e-09))
    lm1 <- lm(log.p.dk ~ log.dk)
    lm2 <- lm(log.p.dk ~ log.dk + I(10 ^ log.dk))
    datout <- data.frame(
        Rsquared.SFT <- summary(lm1)$r.squared,
        slope.SFT <- summary(lm1)$coefficients[2, 1],
        truncatedExponentialAdjRsquared <-
            summary(lm2)$adj.r.squared
    )
    datout
} # end of function scaleFreeFitIndex
