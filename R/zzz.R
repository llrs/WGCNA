# first and last lib functions
#
#' @useDynLib WGCNA
#' @import foreach
#' @import doParallel
#' @import dynamicTreeCut
#' @import fastcluster
#' @import GO.db
#' @import matrixStats
#' @import survival
#' @import parallel
#' @import AnnotationDbi
#' @importFrom Hmisc rcorr.cens errbar
#' @import impute
#' @import splines
#' @import preprocessCore
#' @import grDevices
#' @import graphics
#' @importFrom stats anova as.dendrogram as.dist as.hclust coef cov cutree dist
#' fisher.test glm heatmap kruskal.test lm median model.matrix na.exclude
#' order.dendrogram pchisq phyper pnorm predict pt qnorm quantile reorder
#' residuals rexp rnorm runif sd smooth.spline t.test var weighted.mean
#' @importFrom utils packageVersion compareVersion data flush.console read.csv
#' write.csv write.table
#' @exportPattern "^[^\\.]"
#'
.onAttach = function(libname, pkgname) {
    ourVer = try(gsub("[^0-9_.-] ", packageVersion("WGCNA"),
                      packageVersion("WGCNA"), fixed = FALSE))

    if (inherits(ourVer, "try-error")) {
        ourVer = ""
    }

    packageStartupMessage("=====================================================",
                          "=====================\n*\n*  Package WGCNA ", ourVer,
                          " loaded.\n*")

    if (.useNThreads()==1 && .nProcessorsOnline() > 1) {
        packageStartupMessage(
            "*    Important note: It appears that your system supports ",
            "multi-threading,\n*    but it is not enabled within WGCNA in R.",
            "\n*    To allow multi-threading within WGCNA with all available ",
            "cores, use \n*\n*          allowWGCNAThreads()\n*\n*    within R.",
            " Use disableWGCNAThreads() to disable threading if necessary.\n",
            "*    Alternatively, set the following environment variable on ",
            "your system:\n*\n*          ", .threadAllowVar,
            "=<number_of_processors>\n*\n*    for example \n*\n*          ",
            .threadAllowVar, "=", .nProcessorsOnline(), "\n*\n*    To set the ",
            "environment variable in linux bash shell, type \n*\n*           ",
            "export ", .threadAllowVar, "=", .nProcessorsOnline(),
            "\n*\n*     before running R. Other operating systems or shells ",
            "will\n*     have a similar command to achieve the same aim.\n*")
    }
    packageStartupMessage("===================================================",
                          "=======================")


    imputeVer = try(gsub("[^0-9_.-]", "", packageVersion("impute"),
                         fixed = FALSE))

    if (!inherits(imputeVer, "try-error")) {
        if (compareVersion(imputeVer, "1.12")< 0) {
            packageStartupMessage(
                "*!*!*!*!*!*!* Caution: installed package 'impute' is too ",
                "old.\nOld versions of this package can occasionally crash ",
                "the code or the entire R session.\nIf you already have the ",
                "newest version available from CRAN, \nand you still see this ",
                "warning, please download the impute package \nfrom ",
                "Bioconductor at \nhttp://www.bioconductor.org/packages/",
                "release/bioc/html/impute.html . \n",
                "If the above link is dead, search for package 'impute' \n",
                "in the Downloads -> Software section of",
                " http://www.bioconductor.org .\n")
        }
    }
}

.bucketOrder = function(data, min = min(data), max = max(data), nIntervals, exact = FALSE) {

  if (any(!is.finite(data))) stop("'data' cannot contain missing or infinite values.")
  .Call("bucketOrder_R", as.numeric(data), as.double(min), as.double(max), as.integer(nIntervals),
        as.integer(exact))

}
