
# initProgInd ####
#' Progress indicators
#'
#' Describes how much/how long does it takes
#' @param leadStr Leading characters.
#' @param trailStr Last characters.
#' @param quiet Logical: Run it in interactive mode or not?.
#' @export
initProgInd <-function(leadStr = "..", trailStr = "", quiet = !interactive()) {
    oldStr = " "
    cat(oldStr)
    progInd = list(oldStr = oldStr,
                   leadStr = leadStr,
                   trailStr = trailStr)
    class(progInd) = "progressIndicator"
    updateProgInd(0, progInd, quiet)
}

# updateProgInd  ####
#' Update the progress as it goes
#'
#' Should work on all OS to indicate the evolution of indicator
#' @rdname initProgInd
#' @param newFrac Unkown parameter.
#' @param progInd Object of class progressIndicator obtained with initProgInd.
#' @export
updateProgInd <- function(newFrac, progInd, quiet = !interactive()) {
    if (class(progInd) != "progressIndicator") {
        stop( "Parameter progInd is not of class 'progressIndicator'.\n",
              "Use initProgInd() to initialize",
              "it prior to use.")
    }

    newStr = paste0(progInd$leadStr,
                    as.integer(newFrac * 100),
                    "% ",
                    progInd$trailStr)
    if (newStr != progInd$oldStr) {
        if (quiet) {
            progInd$oldStr = newStr
        } else {
            cat(paste(rep("\b", nchar(
                progInd$oldStr
            )), collapse = ""))
            cat(newStr)
            if (exists("flush.console"))
                flush.console()
            progInd$oldStr = newStr
        }
    }
    progInd
}
