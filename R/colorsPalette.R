# Function of palette of colors

# greenBlackRed ####
#' Green-black-red color sequence
#'
#' Generate a green-black-red color sequence of a given length.
#'
#' The function returns a color vector that starts with pure green, gradually
#' turns into black and then to red. The power \code{gamma} can be used to
#' control the behaviour of the quarter- and three quarter-values (between
#' green and black, and black and red, respectively). Higher powers will make
#' the mid-colors more green and red, respectively.
#'
#' @param n number of colors to be returned
#' @param gamma color correction power
#' @return A vector of colors of length \code{n}.
#' @author Peter Langfelder
#' @keywords color
#' @examples
#'
#'   par(mfrow = c(3, 1))
#'   displayColors(greenBlackRed(50));
#'   displayColors(greenBlackRed(50, 2));
#'   displayColors(greenBlackRed(50, 0.5));
#'
greenBlackRed <- function(n, gamma = 1) {
    half = as.integer(n / 2)
    red = c(rep(0, times = half),
            0,
            seq(
                from = 0,
                to = 1,
                length.out = half
            ) ^ (1 / gamma))
    green = c(seq(
        from = 1,
        to = 0,
        length.out = half
    ) ^ (1 / gamma),
    rep(0, times = half + 1))
    blue = rep(0, times = 2 * half + 1)
    col = rgb(red, green, blue, maxColorValue = 1)
    col
}

# greenWhiteRed ####
#' Green-white-red color sequence
#'
#' Generate a green-white-red color sequence of a given length.
#'
#' The function returns a color vector that starts with green, gradually turns
#' into white and then to red. The power \code{gamma} can be used to control
#' the behaviour of the quarter- and three quarter-values (between green and
#' white, and white and red, respectively). Higher powers will make the
#' mid-colors more white, while lower powers will make the colors more
#' saturated, respectively.
#'
#' Typical use of this function is to produce (via function
#' \code{\link{numbers2colors}}) a color representation of numbers within a
#' symmetric interval around 0, for example, the interval [-1, 1]. Note though
#' that since green and red are not distinguishable by people with the most
#' common type of color blindness, we recommend using the analogous palette
#' returned by the function \code{\link{blueWhiteRed}}.
#'
#' @param n number of colors to be returned
#' @param gamma color change power
#' @param warn logical: should the user be warned that this function produces a
#' palette unsuitable for people with most common color blindness?
#' @return A vector of colors of length \code{n}.
#' @author Peter Langfelder
#' @seealso \code{\link{blueWhiteRed}} for a color sequence more friendly to
#' people with the most common type of color blindness;
#'
#' \code{\link{numbers2colors}} for a function that produces a color
#' representation for continuous numbers.
#' @keywords color
#' @examples
#' \dontrun{
#'   par(mfrow = c(3, 1))
#'   displayColors(greenWhiteRed(50));
#'   title("gamma = 1")
#'   displayColors(greenWhiteRed(50, 3));
#'   title("gamma = 3")
#'   displayColors(greenWhiteRed(50, 0.5));
#'   title("gamma = 0.5")
#' }
greenWhiteRed <- function(n, gamma = 1, warn = TRUE) {
    if (warn)
        warning("WGCNA::greenWhiteRed: this palette is not suitable for people\n",
                "with green - red color blindness (the most common kind of color
                blindness).\n",
                "Consider using the function blueWhiteRed instead.")
    half = as.integer(n / 2)
    red = c(seq(from = 0, to = 1, length.out = half) ^ (1 / gamma),
            rep(1, times = half + 1))
    green = c(rep(1, times = half + 1),
              seq(from = 1, to = 0,length.out = half ) ^ (1 / gamma))
    blue = c(seq(from = 0, to = 1, length.out = half) ^ (1 / gamma), 1,
             seq(from = 1, to = 0, length.out = half) ^ (1 / gamma))
    col = rgb(red, green, blue, maxColorValue = 1)
    col
}

# redWhiteGreen  ####
#' Red-white-green color sequence
#'
#' Generate a red-white-green color sequence of a given length.
#'
#' The function returns a color vector that starts with pure green, gradually
#' turns into white and then to red. The power \code{gamma} can be used to
#' control the behaviour of the quarter- and three quarter-values (between red
#' and white, and white and green, respectively). Higher powers will make the
#' mid-colors more white, while lower powers will make the colors more
#' saturated, respectively.
#'
#' @param n number of colors to be returned
#' @param gamma color correction power
#' @return A vector of colors of length \code{n}.
#' @author Peter Langfelder
#' @keywords color
#' @examples
#'
#'   par(mfrow = c(3, 1))
#'   displayColors(redWhiteGreen(50));
#'   displayColors(redWhiteGreen(50, 3));
#'   displayColors(redWhiteGreen(50, 0.5));
#'
redWhiteGreen <- function(n, gamma = 1) {
    half = as.integer(n / 2)
    green = c(seq(
        from = 0,
        to = 1,
        length.out = half
    ) ^ (1 / gamma),
    rep(1, times = half + 1))
    red = c(rep(1, times = half + 1),
            seq(
                from = 1,
                to = 0,
                length.out = half
            ) ^ (1 / gamma))
    blue = c(
        seq(
            from = 0,
            to = 1,
            length.out = half
        ) ^ (1 / gamma),
        1,
        seq(
            from = 1,
            to = 0,
            length.out = half
        ) ^ (1 / gamma)
    )
    col = rgb(red, green, blue, maxColorValue = 1)
    col
}

# blueWhiteRed ####
#' Blue-white-red color sequence
#'
#' Generate a blue-white-red color sequence of a given length.
#'
#' The function returns a color vector that starts with blue, gradually turns
#' into white and then to red. The power \code{gamma} can be used to control
#' the behaviour of the quarter- and three quarter-values (between blue and
#' white, and white and red, respectively). Higher powers will make the
#' mid-colors more white, while lower powers will make the colors more
#' saturated, respectively.
#'
#' @param n number of colors to be returned.
#' @param gamma color change power.
#' @param endSaturation a number between 0 and 1 giving the saturation of the
#' colors that will represent the ends of the scale. Lower numbers mean less
#' saturation (lighter colors).
#' @return A vector of colors of length \code{n}.
#' @author Peter Langfelder
#' @seealso \code{\link{numbers2colors}} for a function that produces a color
#' representation for continuous numbers.
#' @keywords color
#' @examples
#'
#'   par(mfrow = c(3, 1))
#'   displayColors(blueWhiteRed(50));
#'   title("gamma = 1")
#'   displayColors(blueWhiteRed(50, 3));
#'   title("gamma = 3")
#'   displayColors(blueWhiteRed(50, 0.5));
#'   title("gamma = 0.5")
#'
blueWhiteRed <- function(n, gamma = 1, endSaturation = 1) {
    if (endSaturation  > 1  | endSaturation < 0)
        stop("'endSaturation' must be between 0 and 1.")
    es = 1 - endSaturation
    blueEnd = c(0.05 + es * 0.45, 0.55 + es * 0.25, 1.00)
    redEnd = c(1.0, 0.2 + es * 0.6, 0.6 * es)
    middle = c(1, 1, 1)

    half = as.integer(n / 2)
    if (n %% 2 == 0) {
        index1 = c(1:half)
        index2 = c(1:half) + half
        frac1 = ((index1 - 1) / (half - 1)) ^ (1 / gamma)
        frac2 = rev(frac1)
    } else {
        index1 = c(1:(half + 1))
        index2 = c(1:half) + half + 1
        frac1 = (c(0:half) / half) ^ (1 / gamma)
        frac2 = rev((c(1:half) / half) ^ (1 / gamma))
    }
    cols = matrix(0, n, 3)
    for (c in 1:3) {
        cols[index1, c] = blueEnd[c] + (middle[c] - blueEnd[c]) * frac1
        cols[index2, c] = redEnd[c] + (middle[c] - redEnd[c]) * frac2
    }

    rgb(cols[, 1], cols[, 2], cols[, 3], maxColorValue = 1)
}

# A set of global variables and functions that should help handling color names
# for some 400 + modules. A vector called .GlobalStandardColors is defined that
# holds color names with first few entries being the well - known and  - loved
# colors. The rest is randomly chosen from the color names of R, excluding grey
# colors.

#
# GlobalStandardColors  ####
#
# This code forms a vector of color names in which the first entries are given
# by BaseColors and the rest is "randomly" chosen from the rest of R color names
# that do not contain "grey" nor "gray".

BaseColors <-c("turquoise", "blue", "brown", "yellow", "green", "red", "black",
               "pink", "magenta", "purple", "greenyellow", "tan", "salmon",
               "cyan", "midnightblue", "lightcyan", "grey60", "lightgreen",
               "lightyellow", "royalblue", "darkred", "darkgreen",
               "darkturquoise", "darkgrey", "orange", "darkorange", "white",
               "skyblue", "saddlebrown", "steelblue", "paleturquoise", "violet",
               "darkolivegreen", "darkmagenta")

RColors <- colors()[-grep("grey", colors())]
RColors <- RColors[-grep("gray", RColors)]
InBase <- match(BaseColors, RColors)
ExtraColors <- RColors[-c(InBase[!is.na(InBase)])]
nExtras <- length(ExtraColors)

# Here is the vector of colors that should be used by all functions:

.GlobalStandardColors = c(BaseColors,
                          ExtraColors[rank(
                              sin(13 * c(1:nExtras) + sin(13 * c(1:nExtras))))])

rm(BaseColors, RColors, ExtraColors, nExtras, InBase)
