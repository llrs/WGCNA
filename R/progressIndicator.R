
# initProgInd ####
#' Inline display of progress
#'
#' These functions provide an inline display of pregress.
#'
#' A progress indicator is a simple inline display of progress intended to
#' satisfy impatient users during lengthy operations. The function
#' \code{initProgInd} initializes a progress indicator (at zero)
#' \code{updateProgInd} updates it to a specified fraction.
#'
#' Note that excessive use of \code{updateProgInd} may lead to a performance
#' penalty (see examples).
#'
#' @aliases initProgInd updateProgInd
#' @param leadStr character string that will be printed before the actual
#' progress number.
#' @param trailStr character string that will be printed after the actual
#' progress number.
#' @param quiet can be used to silence the indicator for non-interactive
#' sessions whose output is typically redirected to a file.
#' @param newFrac new fraction of progress to be displayed.
#' @param progInd an object of class \code{progressIndicator} that encodes
#' previously printed message.
#' @return Both functions return an object of class \code{progressIndicator}
#' that holds information on the last printed value and should be used for
#' subsequent updates of the indicator.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' max = 10
#' prog = initProgInd("Counting: ", "done")
#' for (c in 1:max) {
#'   Sys.sleep(0.10)
#'   prog = updateProgInd(c/max, prog)
#' }
#' printFlush("")
#'
#' printFlush("Example 2:")
#' prog = initProgInd()
#' for (c in 1:max) {
#'   Sys.sleep(0.10)
#'   prog = updateProgInd(c/max, prog)
#' }
#' printFlush("")
#'
#' ## Example of a significant slowdown:
#'
#' ## Without progress indicator:
#'
#' system.time( {a = 0; for (i in 1:10000) a = a+i; } )
#'
#' ## With progress indicator, some 50 times slower:
#'
#' system.time({
#'     prog = initProgInd("Counting: ", "done")
#'     a = 0
#'     for (i in 1:10000) {
#'       a = a+i
#'       prog = updateProgInd(i/10000, prog)
#'     }
#'   }
#' )
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
#' @rdname initProgInd
#' @export
updateProgInd <- function(newFrac, progInd, quiet = !interactive()) {
    if (class(progInd) != "progressIndicator") {
        stop( "Parameter progInd is not of class 'progressIndicator'.\n",
              "Use initProgInd() to initialize",
              "it prior to use.")
    }

    newStr = paste0(progInd$leadStr,
                    as.integer(newFrac * 100),
                    "% ", progInd$trailStr)
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
