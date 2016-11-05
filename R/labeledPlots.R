
# formatLabels  ####
#' Break long character strings into multiple lines
#'
#' This function attempts to break lomg character strings into multiple lines
#' by replacing a given pattern by a newline character.
#'
#' Each given element of \code{labels} is processed independently. The
#' character string is split using \code{strsplit}, with \code{split} as the
#' splitting pattern. The resulting shorter character strings are then
#' concatenated together with \code{newsplit} as the separator. Whenever the
#' length of the combined result from the start or the previous newline
#' character exceeds \code{maxCharPerLine}, a newline character is inserted (at
#' the previous split).
#'
#' Note that individual segements (i.e., sections of the input between
#' occurrences of \code{split}) whose number of characters exceeds
#' \code{maxCharPerLine} will not be split.
#'
#' @param labels Character strings to be formatted.
#' @param maxCharPerLine Integer giving the maximum number of characters per
#' line.
#' @param split Pattern to be replaced by newline (\code{'\n'}) characters.
#' @param fixed Logical: Should the pattern be interpreted literally
#' (\code{TRUE}) or as a regular expression (\code{FALSE})? See
#' \code{\link{strsplit}} and its argument \code{fixed}.
#' @param newsplit Character string to replace the occurrences of \code{split}
#' above with.
#' @param keepSplitAtEOL When replacing an occurrence of \code{split} with a
#' newline character, should the \code{newsplit} be added before the newline as
#' well?
#' @return A character vector of the same length as input \code{labels}.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' s = "A quick hare jumps over the brown fox"
#' formatLabels(s)
#'
formatLabels <- function(labels, maxCharPerLine = 14, split = " ", fixed = TRUE,
                         newsplit = split, keepSplitAtEOL = TRUE) {
    n = length(labels)
    splitX = strsplit(labels, split = split, fixed = fixed)
    newLabels = rep("", n)
    for (l in 1:n) {
        nl = ""
        line = ""
        if (nchar(labels[l]) > 0) {
            for (s in 1:length(splitX[[l]])) {
                newLen = nchar(line) + nchar(splitX [[l]] [s])
                if (nchar(line) < 5 | newLen <= maxCharPerLine) {
                    nl = paste(nl, splitX[[l]] [s], sep = newsplit)
                    line = paste(line, splitX[[l]] [s], sep = newsplit)
                } else {
                    nl = paste(nl, splitX[[l]] [s],
                               sep = paste0(
                                   ifelse(keepSplitAtEOL, newsplit, ""),
                                   "\n"))
                    line = splitX[[l]] [s]
                }
            }
        }
        newLabels[l] = nl
    }
    substring(newLabels, nchar(newsplit) + 1)
}

# labeledHeatmap  ####
# This function plots a heatmap of the specified matrix and labels the x and y
# axes wit the given labels.
# It is assumed that the number of entries in xLabels and yLabels is consistent
# with the dimensions in.
# If colorLabels == TRUE, the labels are not printed and instead interpreted as
# colors - -  - -  a simple symbol with the appropriate color is printed instead
#  of the label.
# The x, yLabels are expected to have the form "..color" as in "MEgrey" or
# "PCturquoise". xSymbol, ySymbols are additional markers that can be placed
# next to color labels
#' Produce a labeled heatmap plot
#'
#' Plots a heatmap plot with color legend, row and column annotation, and
#' optional text within th heatmap.
#'
#' The function basically plots a standard heatmap plot of the given
#' \code{Matrix} and embellishes it with row and column labels and/or with text
#' within the heatmap entries. Row and column labels can be either character
#' strings or color squares, or both.
#'
#' To get simple text labels, use \code{colorLabels=FALSE} and pass the desired
#' row and column labels in \code{yLabels} and \code{xLabels}, respectively.
#'
#' To label rows and columns by color squares, use \code{colorLabels=TRUE}
#' \code{yLabels} and \code{xLabels} are then expected to represent valid
#' colors. For reasons of compatibility with other functions, each entry in
#' \code{yLabels} and \code{xLabels} is expected to consist of a color
#' designation preceded by 2 characters: an example would be
#' \code{MEturquoise}. The first two characters can be arbitrary, they are
#' stripped.  Any labels that do not represent valid colors will be considered
#' text labels and printed in full, allowing the user to mix text and color
#' labels.
#'
#' It is also possible to label rows and columns by both color squares and
#' additional text annotation. To achieve this, use the above technique to get
#' color labels and, additionally, pass the desired text annotation in the
#' \code{xSymbols} and \code{ySymbols} arguments.
#'
#' @param Matrix numerical matrix to be plotted in the heatmap.
#' @param xLabels labels for the columns. See Details.
#' @param yLabels labels for the rows. See Details.
#' @param xSymbols additional labels used when \code{xLabels} are interpreted
#' as colors. See Details.
#' @param ySymbols additional labels used when \code{yLabels} are interpreted
#' as colors. See Details.
#' @param colorLabels logical: should \code{xLabels} and \code{yLabels} be
#' interpreted as colors? If given, overrides \code{xColorLabels} and
#' \code{yColorLabels} below.
#' @param xColorLabels logical: should \code{xLabels} be interpreted as colors?
#' @param yColorLabels logical: should \code{yLabels} be interpreted as colors?
#' @param checkColorsValid logical: should given colors be checked for validity
#' against the output of \code{colors()} ? If this argument is \code{FALSE},
#' invalid color specification will trigger an error.
#' @param invertColors logical: should the color order be inverted?
#' @param setStdMargins logical: should standard margins be set before calling
#' the plot function? Standard margins depend on \code{colorLabels}: they are
#' wider for text labels and narrower for color labels. The defaults are
#' static, that is the function does not attempt to guess the optimal margins.
#' @param xLabelsPosition a character string specifying the position of labels
#' for the columns. Recognized values are (unique abbreviations of)
#' \code{"top", "bottom"}.
#' @param xLabelsAngle angle by which the column labels should be rotated.
#' @param xLabelsAdj justification parameter for column labels. See
#' \code{\link{par}} and the description of parameter \code{"adj"}.
#' @param xColorWidth width of the color labels for the x axis expressed as a
#' fraction of the smaller of the range of the x and y axis.
#' @param yColorWidth width of the color labels for the y axis expressed as a
#' fraction of the smaller of the range of the x and y axis.
#' @param xColorOffset gap between the y axis and color labels as a fraction of
#' the range of x axis.
#' @param yColorOffset gap between the x axis and color labels as a fraction of
#' the range of y axis.
#' @param colors color pallette to be used in the heatmap. Defaults to
#' \code{\link{heat.colors}}.
#' @param naColor color to be used for encoding missing data.
#' @param textMatrix optional text entries for each cell. Either a matrix of
#' the same dimensions as \code{Matrix} or a vector of the same length as the
#' number of entries in \code{Matrix}.
#' @param cex.text character expansion factor for \code{textMatrix}.
#' @param textAdj Adjustment for the entries in the text matrix. See the
#' \code{adj} argument to \code{\link{text}}.
#' @param cex.lab character expansion factor for text labels labeling the axes.
#' @param cex.lab.x character expansion factor for text labels labeling the x
#' axis. Overrides \code{cex.lab} above.
#' @param cex.lab.y character expansion factor for text labels labeling the y
#' axis. Overrides \code{cex.lab} above.
#' @param colors.lab.x colors for character labels or symbols along x axis.
#' @param colors.lab.y colors for character labels or symbols along y axis.
#' @param bg.lab.x background color for the margin along the x axis.
#' @param bg.lab.y background color for the margin along the y axs.
#' @param plotLegend logical: should a color legend be plotted?
#' @param keepLegendSpace logical: if the color legend is not drawn, should the
#' space be left empty (\code{TRUE}), or should the heatmap fill the space
#' (\code{FALSE})?
#' @param verticalSeparator.x indices of columns after which separator lines
#' (vertical lines between columns) should be drawn. \code{NULL} means no lines
#' will be drawn.
#' @param verticalSeparator.col color(s) of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.lty line type of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.lwd line width of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.ext number giving the extension of the separator
#' line into the margin as a fraction of the margin width. 0 means no
#' extension, 1 means extend all the way through the margin.
#' @param horizontalSeparator.y indices of columns after which separator lines
#' (horizontal lines between columns) should be drawn. \code{NULL} means no
#' lines will be drawn.
#' @param horizontalSeparator.col color(s) of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.lty line type of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.lwd line width of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.ext number giving the extension of the separator
#' line into the margin as a fraction of the margin width. 0 means no
#' extension, 1 means extend all the way through the margin.
#' @param \dots other arguments to function \code{\link{heatmap}}.
#' @return None.
#' @author Peter Langfelder
#' @seealso \code{\link{heatmap}}, \code{\link{colors}}
#' @keywords hplot
#' @examples
#'
#'
#' # This example illustrates 4 main ways of annotating columns and rows of a heatmap.
#' # Copy and paste the whole example into an R session with an interactive plot window
#' # alternatively, you may replace the command sizeGrWindow below by opening
#' # another graphical device such as pdf.
#'
#' # Generate a matrix to be plotted
#'
#' nCol = 8; nRow = 7
#' mat = matrix(runif(nCol*nRow, min = -1, max = 1), nRow, nCol)
#'
#' rowColors = standardColors(nRow)
#' colColors = standardColors(nRow + nCol)[(nRow+1):(nRow + nCol)]
#'
#' rowColors
#' colColors
#'
#' sizeGrWindow(9,7)
#' par(mfrow = c(2,2))
#' par(mar = c(4, 5, 4, 6))
#'
#' # Label rows and columns by text:
#'
#' labeledHeatmap(mat, xLabels = colColors, yLabels = rowColors,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Text-labeled heatmap")
#'
#' # Label rows and columns by colors:
#'
#' rowLabels = paste0("ME", rowColors)
#' colLabels = paste0("ME", colColors)
#'
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Color-labeled heatmap")
#'
#' # Mix text and color labels:
#'
#' rowLabels[3] = "Row 3"
#' colLabels[1] = "Column 1"
#'
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Mix-labeled heatmap")
#'
#' # Color labels and additional text labels
#'
#' rowLabels = paste0("ME", rowColors)
#' colLabels = paste0("ME", colColors)
#'
#' extraRowLabels = paste("Row", c(1:nRow))
#' extraColLabels = paste("Column", c(1:nCol))
#'
#' # Extend margins to fit all labels
#' par(mar = c(6, 6, 4, 6))
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                xSymbols = extraColLabels,
#'                ySymbols = extraRowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Text- + color-labeled heatmap")
#'
#'
labeledHeatmap <- function(Matrix,
                           xLabels,
                           yLabels = NULL,
                           xSymbols = NULL,
                           ySymbols = NULL,
                           colorLabels = NULL,
                           xColorLabels = FALSE,
                           yColorLabels = FALSE,
                           checkColorsValid = TRUE,
                           invertColors = FALSE,
                           setStdMargins = TRUE,
                           xLabelsPosition = "bottom",
                           xLabelsAngle = 45,
                           xLabelsAdj = 1,
                           xColorWidth = 0.05,
                           yColorWidth = 0.05,
                           # FIXME: For offsetting text, these two seem to be switched
                           xColorOffset = par("cxy")[1] / 3,
                           yColorOffset = par("cxy")[2] / 3,
                           # Content of heatmap
                           colors = NULL,
                           naColor = "grey",
                           textMatrix = NULL,
                           cex.text = NULL,
                           textAdj = c(0.5, 0.5),
                           cex.lab = NULL,
                           cex.lab.x = cex.lab,
                           cex.lab.y = cex.lab,
                           colors.lab.x = 1,
                           colors.lab.y = 1,
                           bg.lab.x = NULL,
                           bg.lab.y = NULL,
                           plotLegend = TRUE,
                           keepLegendSpace = plotLegend,
                           # Separator line specification
                           verticalSeparator.x = NULL,
                           verticalSeparator.col = 1,
                           verticalSeparator.lty = 1,
                           verticalSeparator.lwd = 1,
                           verticalSeparator.ext = 0,

                           horizontalSeparator.y = NULL,
                           horizontalSeparator.col = 1,
                           horizontalSeparator.lty = 1,
                           horizontalSeparator.lwd = 1,
                           horizontalSeparator.ext = 0,
                           ...) {
    if (!is.null(colorLabels)) {
        xColorLabels = colorLabels
        yColorLabels = colorLabels
    }

    if (is.null(yLabels) & (!is.null(xLabels)) &
        (dim(Matrix)[1] == dim(Matrix)[2]))
        yLabels = xLabels

    nCols = ncol(Matrix)
    nRows = nrow(Matrix)

    if (checkColorsValid) {
        xValidColors = !is.na(match(substring(xLabels, 3), colors()))
        yValidColors = !is.na(match(substring(yLabels, 3), colors()))
    } else {
        xValidColors = rep(TRUE, length(xLabels))
        yValidColors = rep(TRUE, length(yLabels))
    }

    if (sum(xValidColors) > 0)
        xColorLabInd = c(1:length(xLabels))[xValidColors]
    if (sum(!xValidColors) > 0)
        xTextLabInd = c(1:length(xLabels))[!xValidColors]

    if (sum(yValidColors) > 0)
        yColorLabInd = c(1:length(yLabels))[yValidColors]
    if (sum(!yValidColors) > 0)
        yTextLabInd = c(1:length(yLabels))[!yValidColors]

    if (setStdMargins) {
        if (xColorLabels & yColorLabels)
        {
            par(mar = c(2, 2, 3, 5) + 0.2)
        } else {
            par(mar = c(7, 7, 3, 5) + 0.2)
        }
    }

    xLabPos = charmatch(xLabelsPosition, c("bottom", "top"))
    if (is.na(xLabPos))
        stop("Argument 'xLabelsPosition' must be (a unique abbreviation of)
             'bottom', 'top'")

    if (is.null(colors))
        colors = heat.colors(30)
    if (invertColors)
        colors = .reverseVector(colors)
    labPos = .heatmapWithLegend(
        Matrix,
        signed = FALSE,
        colors = colors,
        naColor = naColor,
        cex.legend = cex.lab,
        plotLegend = plotLegend,
        keepLegendSpace = keepLegendSpace,
        ...
    )
    #if (plotLegend)
    #{
    #  image.plot(t(.reverseRows(Matrix)), xaxt = "n", xlab = "", yaxt = "n",
    #  ylab = "", col = colors, ...)
    #} else {
    #  image(z = t(.reverseRows(Matrix)), xaxt = "n", xlab = "", yaxt = "n",
    #  ylab = "", col = colors, ...)
    #}
    nxlabels = length(xLabels)
    plotbox = labPos$box
    xmin = plotbox[1]
    xmax = plotbox[2]
    ymin = plotbox[3]
    yrange = plotbox[4] - ymin
    ymax = plotbox[4]
    xrange = xmax - xmin
    xLeft = labPos$xLeft
    xRight = labPos$xRight
    yTop = labPos$yTop
    yBot = labPos$yBot

    xspacing = labPos$xMid[2] - labPos$xMid[1]
    yspacing = abs(labPos$yMid[2] - labPos$yMid[1])

    nylabels = length(yLabels)
    offsetx = yColorOffset
    offsety = xColorOffset
    # Transform fractional widths into coordinate widths
    xColW = min(xmax - xmin, ymax - ymin) * xColorWidth
    yColW = min(xmax - xmin, ymax - ymin) * yColorWidth

    if (any(xValidColors))
        offsety = offsety + xColW
    if (any(yValidColors))
        offsetx = offsetx + yColW

    # Create the background for column and row labels.

    extension.left = par("mai")[2] * # left margin width in inches
        # charcter size in user corrdinates/character size in inches
        par("cxy")[1] / par("cin")[1]

    extension.bottom = par("mai")[1]  *
        # charcter size in user corrdinates/character size in inches
        par("cxy")[2] / par("cin")[2] - offsety

    extension.top = par("mai")[3]  *
        # charcter size in user corrdinates/character size in inches
        par("cxy")[2] / par("cin")[2] - offsety

    figureBox = par("usr")
    figXrange = figureBox[2] - figureBox[1]
    figYrange = figureBox[4] - figureBox[3]
    if (!is.null(bg.lab.x)) {
        bg.lab.x = .extend(bg.lab.x, nCols)
        if (xLabPos == 1) {
            y0 = ymin
            ext = extension.bottom
            sign = 1
        } else {
            y0 = ymax
            ext = extension.top
            sign = -1
        }
        figureDims = par("pin")
        angle = xLabelsAngle / 180 * pi
        ratio = figureDims[1] / figureDims[2] * figYrange / figXrange
        ext.x = -sign * ext * 1 / tan(angle) / ratio
        ext.y = sign * ext * sign(sin(angle))
        for (c in 1:nCols)
            polygon(
                x = c(
                    xLeft[c],
                    xLeft[c],
                    xLeft[c] + ext.x,
                    xRight[c] + ext.x,
                    xRight[c],
                    xRight[c]
                ),
                y = c(
                    y0,
                    y0 - sign * offsety,
                    y0 - sign * offsety - ext.y,
                    y0 - sign * offsety - ext.y,
                    y0 - sign * offsety,
                    y0
                ),
                border = bg.lab.x[c],
                col = bg.lab.x[c],
                xpd = TRUE
            )
    }

    if (!is.null(bg.lab.y)) {
        bg.lab.y = .extend(bg.lab.y, nRows)
        reverseRows = TRUE
        if (reverseRows) {
            bg.lab.y = rev(bg.lab.y)
        }
        for (r in 1:nRows)
            rect(
                xmin - extension.left,
                yBot[r],
                xmin,
                yTop[r],
                col = bg.lab.y[r],
                border = bg.lab.y[r],
                xpd = TRUE
            )
    }

    # Write out labels
    if (sum(!xValidColors) > 0) {
        xLabYPos = ifelse(xLabPos == 1, ymin - offsety, ymax + offsety)
        if (is.null(cex.lab))
            cex.lab = 1
        mapply(
            text,
            x = labPos$xMid[xTextLabInd],
            labels = xLabels[xTextLabInd],
            MoreArgs = list(
                y = xLabYPos,
                srt = xLabelsAngle,
                adj = xLabelsAdj,
                xpd = TRUE,
                cex = cex.lab.x,
                col = colors.lab.x
            )
        )
    }
    if (sum(xValidColors) > 0) {
        baseY = ifelse(xLabPos == 1, ymin - offsety, ymax + offsety)
        deltaY = ifelse(xLabPos == 1, xColW,-xColW)
        rect(
            xleft = labPos$xMid[xColorLabInd] - xspacing / 2,
            ybottom = baseY,
            xright = labPos$xMid[xColorLabInd] + xspacing / 2,
            ytop = baseY + deltaY,
            density = -1,
            col = substring(xLabels[xColorLabInd], 3),
            border = substring(xLabels[xColorLabInd], 3),
            xpd = TRUE
        )
        if (!is.null(xSymbols))
            mapply(
                text,
                x = labPos$xMid[xColorLabInd],
                labels = xSymbols[xColorLabInd],
                MoreArgs = list(
                    baseY - sign(deltaY) *  offsety,
                    adj = xLabelsAdj,
                    xpd = TRUE,
                    srt = xLabelsAngle,
                    cex = cex.lab.x,
                    col = colors.lab.x
                )
            )
    }
    if (sum(!yValidColors) > 0) {
        if (is.null(cex.lab))
            cex.lab = 1
        mapply(
            text,
            y = labPos$yMid[yTextLabInd],
            labels = yLabels[yTextLabInd],
            MoreArgs = list(
                x = xmin - offsetx,
                srt = 0,
                adj = c(1, 0.5),
                xpd = TRUE,
                cex = cex.lab.y,
                col = colors.lab.y
            )
        )
    }
    if (sum(yValidColors) > 0) {
        rect(
            xleft = xmin -  offsetx,
            ybottom = rev(labPos$yMid[yColorLabInd]) - yspacing / 2,
            xright = xmin -  offsetx + yColW,
            ytop = rev(labPos$yMid[yColorLabInd]) + yspacing / 2,
            density = -1,
            col = substring(rev(yLabels[yColorLabInd]), 3),
            border = substring(rev(yLabels[yColorLabInd]), 3),
            xpd = TRUE
        )
        #for (i in yColorLabInd)
        #{
        #  lines(c(xmin -  offsetx, xmin -  offsetx + yColW),
        #  y = rep(labPos$yMid[i] - yspacing/2, 2), col = i, xpd = TRUE)
        #  lines(c(xmin -  offsetx, xmin -  offsetx + yColW),
        #  y = rep(labPos$yMid[i] + yspacing/2, 2), col = i, xpd = TRUE)
        #}
        if (!is.null(ySymbols))
            mapply(
                text,
                y = labPos$yMid[yColorLabInd],
                labels = ySymbols[yColorLabInd],
                MoreArgs = list(
                    xmin + yColW - 2 * offsetx,
                    adj = c(1, 0.5),
                    xpd = TRUE,
                    cex = cex.lab.y,
                    col = colors.lab.y
                )
            )
    }

    # Draw separator lines, if requested

    if (length(verticalSeparator.x) > 0) {
        nLines = length(verticalSeparator.x)
        vs.col = .extend(verticalSeparator.col, nLines)
        vs.lty = .extend(verticalSeparator.lty, nLines)
        vs.lwd = .extend(verticalSeparator.lwd, nLines)
        vs.ext = .extend(verticalSeparator.ext, nLines)
        if (any(verticalSeparator.x < 0 |
                verticalSeparator.x > nCols))
            stop("If given. 'verticalSeparator.x' must all be between 0 and the
                 number of columns.")
        x.lines = ifelse(verticalSeparator.x > 0,
                         labPos$xRight[verticalSeparator.x],
                         labPos$xLeft[1])
        for (l in 1:nLines)
            lines(
                rep(x.lines[l], 2),
                c(ymin, ymax),
                col = vs.col[l],
                lty = vs.lty[l],
                lwd = vs.lwd[l]
            )

        angle = xLabelsAngle / 180 * pi
        if (xLabelsPosition == "bottom") {
            sign = 1
            y0 = ymin
        } else {
            sign = -1
            y0 = ymax
        }
        figureDims = par("pin")
        ratio = figureDims[1] / figureDims[2] * figYrange / figXrange
        ext.x = -sign * extension.bottom * 1 / tan(angle) / ratio
        ext.y = sign * extension.bottom * sign(sin(angle))
        for (l in 1:nLines)
            lines(
                c(x.lines[l], x.lines[l], x.lines[l] + vs.ext * ext.x),
                c(y0, y0 - sign * offsety,
                  y0 - sign * offsety - vs.ext * ext.y),
                col = vs.col[l],
                lty = vs.lty[l],
                lwd = vs.lwd[l],
                xpd = TRUE
            )
    }

    if (length(horizontalSeparator.y) > 0) {
        if (any(horizontalSeparator.y < 0 | horizontalSeparator.y > nRows))
            stop("If given. 'horizontalSeparator.y' must all be between 0 and
                 the number of rows.")
        reverseRows = TRUE
        if (reverseRows)
        {
            horizontalSeparator.y = nRows - horizontalSeparator.y + 1
            y.lines = ifelse(horizontalSeparator.y <= nRows,
                             labPos$yBot[horizontalSeparator.y],
                             labPos$yTop[nRows])
        } else {
            y.lines = ifelse(horizontalSeparator.y > 0,
                             labPos$yBot[horizontalSeparator.y],
                             labPos$yTop[1])
        }
        nLines = length(horizontalSeparator.y)
        vs.col = .extend(horizontalSeparator.col, nLines)
        vs.lty = .extend(horizontalSeparator.lty, nLines)
        vs.lwd = .extend(horizontalSeparator.lwd, nLines)
        vs.ext = .extend(horizontalSeparator.ext, nLines)
        for (l in 1:nLines)
            lines(
                c(xmin - vs.ext[l] * extension.left, xmax),
                rep(y.lines[l], 2),
                col = vs.col[l],
                lty = vs.lty[l],
                lwd = vs.lwd[l],
                xpd = TRUE
            )
    }

    if (!is.null(textMatrix)) {
        if (is.null(cex.text))
            cex.text = par("cex")
        if (is.null(dim(textMatrix)))
            if (length(textMatrix) == prod(dim(Matrix))) {
                dim(textMatrix) = dim(Matrix)
            }
        if (!isTRUE(all.equal(dim(textMatrix), dim(Matrix))))
            stop(
                "labeledHeatmap: textMatrix was given, but has dimensions
                incompatible with Matrix."
            )
        for (rw in 1:dim(Matrix)[1])
            for (cl in 1:dim(Matrix)[2]) {
                text(
                    labPos$xMid[cl],
                    labPos$yMid[rw],
                    as.character(textMatrix[rw, cl]),
                    xpd = TRUE,
                    cex = cex.text,
                    adj = textAdj
                )
            }
    }
    axis(1, labels = FALSE, tick = FALSE)
    axis(2, labels = FALSE, tick = FALSE)
    axis(3, labels = FALSE, tick = FALSE)
    axis(4, labels = FALSE, tick = FALSE)
    invisible(labPos)
}

# labeledHeatmap.multipage ####
#' Labeled heatmap divided into several separate plots.
#'
#' This function produces labaled heatmaps divided into several plots. This is
#' useful for large heatmaps where labels on individual columns and rows may
#' become unreadably small (or overlap).
#'
#' The function \code{\link{labeledHeatmap}} is used to produce each plot/page
#' most arguments are described in more detail in the help file for that
#' function.
#'
#' In each plot/page \code{\link{labeledHeatmap}} plots a standard heatmap plot
#' of an appropriate sub-rectangle of \code{Matrix} and embellishes it with row
#' and column labels and/or with text within the heatmap entries. Row and
#' column labels can be either character strings or color squares, or both.
#'
#' To get simple text labels, use \code{colorLabels=FALSE} and pass the desired
#' row and column labels in \code{yLabels} and \code{xLabels}, respectively.
#'
#' To label rows and columns by color squares, use \code{colorLabels=TRUE}
#' \code{yLabels} and \code{xLabels} are then expected to represent valid
#' colors. For reasons of compatibility with other functions, each entry in
#' \code{yLabels} and \code{xLabels} is expected to consist of a color
#' designation preceded by 2 characters: an example would be
#' \code{MEturquoise}. The first two characters can be arbitrary, they are
#' stripped. Any labels that do not represent valid colors will be considered
#' text labels and printed in full, allowing the user to mix text and color
#' labels.
#'
#' It is also possible to label rows and columns by both color squares and
#' additional text annotation. To achieve this, use the above technique to get
#' color labels and, additionally, pass the desired text annotation in the
#' \code{xSymbols} and \code{ySymbols} arguments.
#'
#' If \code{rowsPerPage} (\code{colsPerPage}) is not given, rows (columns) are
#' allocated automatically as uniformly as possible, in contiguous blocks of
#' size at most \code{maxRowsPerPage} (\code{maxColsPerPage}).  The allocation
#' is performed by the function \code{\link{allocateJobs}}.
#'
#' @param Matrix numerical matrix to be plotted in the heatmap.
#' @param xLabels labels for the columns. See Details.
#' @param yLabels labels for the rows. See Details.
#' @param xSymbols additional labels used when \code{xLabels} are interpreted
#' as colors. See Details.
#' @param ySymbols additional labels used when \code{yLabels} are interpreted
#' as colors. See Details.
#' @param textMatrix optional text entries for each cell. Either a matrix of
#' the same dimensions as \code{Matrix} or a vector of the same length as the
#' number of entries in \code{Matrix}.
#' @param rowsPerPage optional list in which each component is a vector
#' specifying which rows should appear together in each plot. If not given,
#' will be generated automatically based on \code{maxRowsPerPage} below and the
#' number of rows in \code{Matrix}.
#' @param maxRowsPerPage integer giving maximum number of rows appearing on
#' each plot (page).
#' @param colsPerPage optional list in which each component is a vector
#' specifying which columns should appear together in each plot. If not given,
#' will be generated automatically based on \code{maxColsPerPage} below and the
#' number of rows in \code{Matrix}.
#' @param maxColsPerPage integer giving maximum number of columns appearing on
#' each plot (page).
#' @param addPageNumberToMain logical: should plot/page number be added to the
#' \code{main} title of each plot?
#' @param zlim Optional specification of the extreme values for the color
#' scale. If not given, will be determined from the input \code{Matrix}.
#' @param main Main title for each plot/page, optionally with the plot/page
#' number added.
#' @param signed logical: should the input \code{Matrix} be converted to colors
#' using a scale centered at zero?
#' @param verticalSeparator.x indices of columns after which separator lines
#' (vertical lines between columns) should be drawn. \code{NULL} means no lines
#' will be drawn.
#' @param verticalSeparator.col color(s) of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.lty line type of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.lwd line width of the vertical separator lines.
#' Recycled if need be.
#' @param verticalSeparator.ext number giving the extension of the separator
#' line into the margin as a fraction of the margin width. 0 means no
#' extension, 1 means extend all the way through the margin.
#' @param horizontalSeparator.y indices of columns after which separator lines
#' (horizontal lines between columns) should be drawn. \code{NULL} means no
#' lines will be drawn.
#' @param horizontalSeparator.col color(s) of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.lty line type of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.lwd line width of the horizontal separator lines.
#' Recycled if need be.
#' @param horizontalSeparator.ext number giving the extension of the separator
#' line into the margin as a fraction of the margin width. 0 means no
#' extension, 1 means extend all the way through the margin.
#' @param \dots other arguments to function \code{\link{labeledHeatmap}}.
#' @return None.
#' @author Peter Langfelder
#' @seealso
#' The workhorse function \code{\link{labeledHeatmap}} for the actual
#' heatmap plot
#'
#' function \code{\link{allocateJobs}} for the allocation of rows/columns to
#' each plot.
#' @keywords misc
labeledHeatmap.multiPage <- function(# Input data and ornaments
    Matrix,
    xLabels,
    yLabels = NULL,
    xSymbols = NULL,
    ySymbols = NULL,
    textMatrix = NULL,

    # Paging options
    rowsPerPage = NULL,
    maxRowsPerPage = 20,
    colsPerPage = NULL,
    maxColsPerPage = 10,
    addPageNumberToMain = TRUE,

    # Further arguments to labeledHeatmap
    zlim = NULL,
    signed = TRUE,
    main = "",

    verticalSeparator.x = NULL,
    verticalSeparator.col = 1,
    verticalSeparator.lty = 1,
    verticalSeparator.lwd = 1,
    verticalSeparator.ext = 0,

    horizontalSeparator.y = NULL,
    horizontalSeparator.col = 1,
    horizontalSeparator.lty = 1,
    horizontalSeparator.lwd = 1,
    horizontalSeparator.ext = 0,

    ...) {
    nr = nrow(Matrix)
    nc = ncol(Matrix)

    if (is.null(rowsPerPage)) {
        nPages.rows = ceiling(nr / maxRowsPerPage)
        rowsPerPage = allocateJobs(nr, nPages.rows)
    } else
        nPages.rows = length(rowsPerPage)

    if (is.null(colsPerPage)) {
        nPages.cols = ceiling(nc / maxColsPerPage)
        colsPerPage = allocateJobs(nc, nPages.cols)
    } else
        nPages.cols = length(colsPerPage)

    if (is.null(zlim)) {
        zlim = range(Matrix, na.rm = TRUE)
        if (signed)
            zlim = c(-max(abs(zlim)), max(abs(zlim)))
    }

    if (!is.null(verticalSeparator.x)) {
        nvs = length(verticalSeparator.x)
        verticalSeparator.col = .extend(verticalSeparator.col, nvs)
        verticalSeparator.lty = .extend(verticalSeparator.lty, nvs)
        verticalSeparator.lwd = .extend(verticalSeparator.lwd, nvs)
        verticalSeparator.ext = .extend(verticalSeparator.ext, nvs)
    }

    if (!is.null(horizontalSeparator.y)) {
        nhs = length(horizontalSeparator.y)
        horizontalSeparator.col = .extend(horizontalSeparator.col, nhs)
        horizontalSeparator.lty = .extend(horizontalSeparator.lty, nhs)
        horizontalSeparator.lwd = .extend(horizontalSeparator.lwd, nhs)
        horizontalSeparator.ext = .extend(horizontalSeparator.ext, nhs)
    }


    page = 1
    multiPage = (nPages.cols > 1 | nPages.rows > 1)

    for (page.col in 1:nPages.cols)
        for (page.row in 1:nPages.rows) {
            rows = rowsPerPage[[page.row]]
            cols = colsPerPage[[page.col]]
            if (!is.null(verticalSeparator.x)) {
                keep.vs = verticalSeparator.x %in% cols
            } else
                keep.vs = numeric(0)
            if (!is.null(horizontalSeparator.y)) {
                keep.hs = horizontalSeparator.y %in% cols
            } else
                keep.hs = numeric(0)

            main.1 = main
            if (addPageNumberToMain & multiPage)
                main.1 = paste0(main, "(page ", page, ")")
            labeledHeatmap(
                Matrix = Matrix[rows, cols, drop = FALSE],
                xLabels = xLabels[cols],
                xSymbols = xSymbols[cols],
                yLabels = yLabels[rows],
                ySymbols = ySymbols[rows],
                textMatrix = textMatrix[rows, cols, drop = FALSE],
                zlim = zlim,
                main = main.1,
                verticalSeparator.x = verticalSeparator.x[keep.vs] - min(cols) + 1,
                verticalSeparator.col = verticalSeparator.col[keep.vs],
                verticalSeparator.lty = verticalSeparator.lty[keep.vs],
                verticalSeparator.lwd = verticalSeparator.lwd[keep.vs],
                verticalSeparator.ext = verticalSeparator.ext[keep.vs],

                horizontalSeparator.y = horizontalSeparator.y[keep.hs] -
                    min(rows) + 1,
                horizontalSeparator.col = horizontalSeparator.col[keep.hs],
                horizontalSeparator.lty = horizontalSeparator.lty[keep.hs],
                horizontalSeparator.lwd = horizontalSeparator.lwd[keep.hs],
                horizontalSeparator.ext = horizontalSeparator.ext[keep.hs],
                ...
            )
            page = page + 1
        }
}

# labeledBarplot  ####
#' Barplot with text or color labels.
#'
#' Produce a barplot with extra annotation.
#'
#'
#' Individual bars in the barplot can be identified either by printing the text
#' of the corresponding entry in \code{labels} underneath the bar at the angle
#' specified by \code{xLabelsAngle}, or by interpreting the \code{labels} entry
#' as a color (see below) and drawing a correspondingly colored square
#' underneath the bar.
#'
#' For reasons of compatibility with other functions, \code{labels} are
#' interpreted as colors after stripping the first two characters from each
#' label. For example, the label \code{"MEturquoise"} is interpreted as the
#' color turquoise.
#'
#' If \code{colored} is set, the code assumes that \code{labels} can be
#' interpreted as colors, and the input \code{Matrix} is square and the rows
#' have the same labels as the columns. Each bar in the barplot is then
#' sectioned into contributions from each row entry in \code{Matrix} and is
#' colored by the color given by the entry in \code{labels} that corresponds to
#' the row.
#'
#' @param Matrix vector or a matrix to be plotted.
#' @param labels labels to annotate the bars underneath the barplot.
#' @param colorLabels logical: should the labels be interpreted as colors? If
#' \code{TRUE}, the bars will be labeled by colored squares instead of text.
#' See details.
#' @param colored logical: should the bars be divided into segments and
#' colored? If \code{TRUE}, assumes the \code{labels} can be interpreted as
#' colors, and the input \code{Matrix} is square and the rows have the same
#' labels as the columns. See details.
#' @param setStdMargins if \code{TRUE}, the function wil set margins \code{c(3,
#' 3, 2, 2)+0.2}.
#' @param stdErrors if given, error bars corresponding to \code{1.96*stdErrors}
#' will be plotted on top of the bars.
#' @param cex.lab character expansion factor for axis labels, including the
#' text labels underneath the barplot.
#' @param xLabelsAngle angle at which text labels under the barplot will be
#' printed.
#' @param \dots other parameters for the function \code{\link{barplot}}.
#' @return None.
#' @author Peter Langfelder
#' @keywords hplot
labeledBarplot <- function(Matrix, labels, colorLabels = FALSE, colored = TRUE,
                           setStdMargins = TRUE, stdErrors = NULL,
                           cex.lab = NULL, xLabelsAngle = 45, ...) {
    if (setStdMargins)
        par(mar = c(3, 3, 2, 2) + 0.2)

    if (colored) {
        colors = substring(labels, 3)
    } else {
        colors = rep("grey", times = ifelse(length(dim(Matrix)) < 2,
                                            length(Matrix), dim(Matrix)[[2]]))
    }

    ValidColors = !is.na(match(substring(labels, 3), colors()))

    if (sum(ValidColors) > 0)
        ColorLabInd = c(1:length(labels))[ValidColors]
    if (sum(!ValidColors) > 0)
        TextLabInd = c(1:length(labels))[!ValidColors]

    colors[!ValidColors] = "grey"

    mp = barplot(
        Matrix,
        col = colors,
        xaxt = "n",
        xlab = "",
        yaxt = "n",
        ...
    )

    if (length(dim(Matrix)) == 2) {
        means = apply(Matrix, 2, sum)
    } else {
        means = Matrix
    }

    if (!is.null(stdErrors))
        addErrorBars(means, 1.96 * stdErrors,
                     two.side = TRUE)

    # axis(1, labels = FALSE)
    nlabels = length(labels)
    plotbox = par("usr")
    xmin = plotbox[1]
    xmax = plotbox[2]
    ymin = plotbox[3]
    yrange = plotbox[4] - ymin
    ymax = plotbox[4]
    # print(paste("yrange:", yrange))
    if (nlabels > 1) {
        spacing = (mp[length(mp)] - mp[1]) / (nlabels - 1)
    } else {
        spacing = (xmax - xmin)
    }
    yoffset = yrange / 30
    xshift = spacing / 2
    xrange = spacing * nlabels
    if (is.null(cex.lab))
        cex.lab = 1
    if (colorLabels) {
        #rect(xshift + ((1:nlabels) - 1) * spacing - spacing/2.1,
        #ymin - spacing/2.1 - spacing/8,
        #     xshift + ((1:nlabels) - 1) * spacing + spacing/2.1,
        #     ymin - spacing/8,
        #     density = - 1, col = substring(labels, 3),
        #     border = substring(labels, 3), xpd = TRUE)
        if (sum(!ValidColors) > 0) {
            text(
                mp[!ValidColors],
                ymin - 0.02,
                srt = 45,
                adj = 1,
                labels = labels[TextLabInd],
                xpd = TRUE,
                cex = cex.lab,
                srt = xLabelsAngle
            )
        }
        if (sum(ValidColors) > 0) {
            rect(
                mp[ValidColors] - spacing / 2.1,
                ymin - 2 * spacing / 2.1 * yrange / xrange - yoffset,
                mp[ValidColors] + spacing / 2.1,
                ymin - yoffset,
                density = -1,
                col = substring(labels[ValidColors], 3),
                border = substring(labels[ValidColors], 3),
                xpd = TRUE
            )
        }
    } else {
        text(((1:nlabels) - 1) * spacing + spacing / 2,
             ymin - 0.02 * yrange,
             srt = 45,
             adj = 1,
             labels = labels,
             xpd = TRUE,
             cex = cex.lab,
             srt = xLabelsAngle
        )
    }
    axis(2, labels = TRUE)
}

# shortenStrings ####

.listRep <- function(data, n) {
    out = list()
    if (n >  0) {
        for (i in 1:n) {
            out[[i]] = data
        }
    }
    out
}

#' Shorten given character strings by truncating at a suitable separator.
#'
#' This function shortens given character strings so they are not longer than a
#' given maximum length.
#'
#' Strings whose length (number of characters) is at most \code{maxLength} are
#' returned unchanged. For those that are longer, the function uses
#' \code{\link{gregexpr}} to search for the occurrences of \code{split} in each
#' given character string. If such occurrences are found at positions between
#' \code{minLength} and \code{maxLength}, the string will be truncated at the
#' last such \code{split}; otherwise, the string will be truncated at
#' \code{maxLength}. The \code{ellipsis} is appended to each truncated string.
#'
#' @param strings Character strings to be shortened.
#' @param maxLength Maximum length (number of characters) in the strings to be
#' retained. See details for when the returned strings can exceed this length.
#' @param minLength Minimum length of the returned strings. See details.
#' @param split Character string giving the split at which the strings can be
#' truncated. This can be a literal string or a regular expression (if the
#' latter, \code{fixed} below must be set to \code{FALSE}).
#' @param fixed Logical: should \code{split} be interpreted as a literal
#' specification (\code{TRUE}) or as a regular expression (\code{FALSE})?
#' @param ellipsis Character string that will be appended to every shorten
#' string, to indicate that the string has been shortened.
#' @param countEllipsisInLength Logical: should the length of the ellipsis
#' count toward the minimum and maximum length?
#' @return A character vector of strings, shortened as necessary. If the input
#' \code{strings} had non-NULL dimensions and dimnames, these are copied to the
#' output.
#' @author Peter Langfelder
#' @seealso
#' \code{\link{gregexpr}}, the workhorse pattern matching function
#' \code{\link{formatLabels}} for splitting strings into multiple lines
#' @keywords misc
shortenStrings <- function(strings, maxLength = 25, minLength = 10,
                           split = " ", fixed = TRUE, ellipsis = "...",
                           countEllipsisInLength = FALSE) {
    dims = dim(strings)
    dnames = dimnames(strings)
    if (is.data.frame(strings)) {
        strings = as.matrix(strings)
        outputDF = TRUE
    } else {
        outputDF = FALSE
    }
    strings = as.character(strings)
    n = length(strings)
    if (n == 0)
        return(character(0))

    newLabels = rep("", n)
    if (length(split) > 0) {
        splitPositions = gregexpr(pattern = split,
                                  text = strings,
                                  fixed = fixed)
    } else {
        splitPositions = .listRep(numeric(0), n)
    }
    if (countEllipsisInLength) {
        maxLength = maxLength - nchar(ellipsis)
        minLength = minLength - nchar(ellipsis)
    }
    for (l in 1:n) {
        if (nchar(strings[l]) <= maxLength) {
            newLabels[l] = strings[l]
        } else {
            splits.1 = splitPositions[[l]]
            suitableSplits = which(splits.1 > minLength &
                                       splits.1 <= maxLength)
            if (length(suitableSplits) > 0) {
                splitPosition = max(splits.1[suitableSplits])
            } else {
                splitPosition = maxLength + 1
            }
            newLabels[l] = paste0(substring(strings[l], 1, splitPosition - 1),
                                  ellipsis)
        }
    }

    dim(newLabels) = dims
    dimnames(newLabels) = dnames
    if (outputDF) {
        as.data.frame(newLabels)
    } else {
        newLabels
    }
}
