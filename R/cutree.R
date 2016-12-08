
# cutreeStatic ####
#' Constant-height tree cut
#'
#' Module detection in hierarchical dendrograms using a constant-height tree
#' cut. Only branches whose size is at least \code{minSize} are retained.
#'
#' This function performs a straightforward constant-height cut as implemented
#' by \code{\link{cutree}}, then calculates the number of objects on each
#' branch and only keeps branches that have at least \code{minSize} objects on
#' them.
#'
#' @param dendro a hierarchical clustering dendrogram such as returned by
#' \code{\link{hclust}}.
#' @param cutHeight height at which branches are to be cut.
#' @param minSize minimum number of object on a branch to be considered a
#' cluster.
#' @return A numeric vector giving labels of objects, with 0 meaning
#' unassigned. The largest cluster is conventionally labeled 1, the next
#' largest 2, etc.
#' @author Peter Langfelder
#' @seealso
#' \code{\link[stats]{hclust}} for hierarchical clustering,
#' \code{\link[stats]{cutree}} for other constant-height branch cuts,
#' \code{\link{standardColors}} to convert the
#' retuned numerical lables into colors for easier visualization or
#' \code{\link{cutreeStaticColor}}.
#' @keywords misc
#' @examples
#' datExpr <- matrix(rnorm(150), ncol = 5)
#' hc <- hclust(dist(datExpr))
#' cutreeStatic(hc, minSize = 2)
cutreeStatic <- function(dendro,
                         cutHeight = 0.9,
                         minSize = 50) {
    normalizeLabels(moduleNumber(dendro, cutHeight, minSize))
}

# cutreeStaticColor ####
#' Constant height tree cut using color labels
#'
#' Cluster detection by a constant height cut of a hierarchical clustering
#' dendrogram.
#'
#' This function performs a straightforward constant-height cut as implemented
#' by \code{\link{cutree}}, then calculates the number of objects on each
#' branch and only keeps branches that have at least \code{minSize} objects on
#' them.
#'
#' @inheritParams cutreeStatic
#' @return A character vector giving color labels of objects, with "grey"
#' meaning unassigned. The largest cluster is conventionally labeled
#' "turquoise", next "blue" etc. Run \code{standardColors()} to see the
#' sequence of standard color labels.
#' @author Peter Langfelder
#' @seealso
#' for hierarchical clustering: \code{\link[stats]{hclust}} ,
#' for for other constant-height branch cuts\code{\link[stats]{cutree}} and
#' \code{\link{cutreeStatic}}.
#' @examples
#' datExpr <- matrix(rnorm(150), ncol = 5)
#' hc <- hclust(dist(datExpr))
#' cutreeStaticColor(hc, minSize = 2)
#' cutreeStatic(hc, minSize = 2)
cutreeStaticColor <- function(dendro, cutHeight = 0.9, minSize = 50) {
    labels2colors(cutreeStatic(dendro, cutHeight, minSize))
}
