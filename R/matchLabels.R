# Relabel the labels in source such that modules with high overlap with those
# in reference will have the same labels
# overlapTable ####


#' Calculate overlap of modules
#' 
#' The function calculates overlap counts and Fisher exact test p-values for
#' the given two sets of module assignments.
#' 
#' 
#' @param labels1 a vector containing module labels.
#' @param labels2 a vector containing module labels to be compared to
#' \code{labels1}.
#' @param na.rm logical: should entries missing in either \code{labels1} or
#' \code{labels2} be removed?
#' @param ignore an optional vector giving label levels that are to be ignored.
#' @return A list with the following components: \item{countTable}{a matrix
#' whose rows correspond to modules (unique labels) in \code{labels1} and whose
#' columns correspond to modules (unique labels) in \code{labels2}, giving the
#' number of objects in the intersection of the two respective modules. }
#' 
#' \item{pTable}{a matrix whose rows correspond to modules (unique labels) in
#' \code{labels1} and whose columns correspond to modules (unique labels) in
#' \code{labels2}, giving Fisher's exact test significance p-values for the
#' overlap of the two respective modules. }
#' @author Peter Langfelder
#' @seealso \code{\link{fisher.test}}, \code{\link{matchLabels}}
#' @keywords misc
#' @examples
#' 
#' set.seed(1)
#' nModules <- 10
#' nGenes <- 1000
#' labels1 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' labels2 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' overlapTable(labels1, labels2)
#' 
overlapTable <- function(labels1, labels2, na.rm = TRUE, ignore = NULL) {
    labels1 <- as.factor(labels1)
    labels2 <- as.factor(labels2)
    if (na.rm) {
        keep <- !is.na(labels1) & !is.na(labels2)
        labels1 <- labels1[keep]
        labels2 <- labels2[keep]
    }

    levels1 <- levels(labels1)
    levels2 <- levels(labels2)
    n1 <- length(levels1)
    n2 <- length(levels2)
    ignore <- as.character(ignore)
    countMat <- table(labels1, labels2, exclude = c(NA, NaN, ignore),
                      dnn = NULL)

    pMat <- matrix(0, n1, n2)

    for (m1 in 1:n1) {
        if (levels1[m1] %in% ignore) {
            next
        }
        for (m2 in 1:n2) {
            if (levels2[m2] %in% ignore) {
                next
            }
            m1Members <- (labels1 == levels1[m1])
            m2Members <- (labels2 == levels2[m2])
            tab <- table(m1Members, m2Members)
            pMat[m1, m2] <- fisher.test(tab, alternative = "greater")$p.value
        }
    }

    dimnames(pMat) <- list(levels1, levels2)
    pMat <- pMat[!rownames(pMat) %in% ignore, !colnames(pMat) %in% ignore]
    pMat[is.na(pMat)] <- 1

    list(countTable = countMat, pTable = pMat)
}

# matchLabels ####


#' Relabel module labels to best match the given reference labels
#' 
#' Given a \code{source} and \code{reference} vectors of module labels, the
#' function produces a module labeling that is equivalent to \code{source}, but
#' individual modules are re-labeled so that modules with significant overlap
#' in \code{source} and \code{reference} have the same labels.
#' 
#' Each column of \code{source} is treated separately. Unlike in previous
#' version of this function, source and reference labels can be any labels, not
#' necessarily of the same type.
#' 
#' The function calculates the overlap of the \code{source} and
#' \code{reference} modules using Fisher's exact test. It then attempts to
#' relabel \code{source} modules such that each \code{source} module gets the
#' label of the \code{reference} module that it overlaps most with, subject to
#' not renaming two \code{source} modules to the same \code{reference} module.
#' (If two \code{source} modules point to the same \code{reference} module, the
#' one with the more significant overlap is chosen.)
#' 
#' Those \code{source} modules that cannot be matched to a \code{reference}
#' module are labeled using those labels from \code{extraLabels} that do not
#' occur in either of \code{source}, \code{reference} or \code{ignoreLabels}.
#' 
#' @aliases matchLabels match.Labels
#' @param source a vector or a matrix of reference labels. The labels may be
#' numeric or character.
#' @param reference a vector of reference labels.
#' @param pThreshold threshold of Fisher's exact test for considering modules
#' to have a significant overlap.
#' @param na.rm logical: should missing values in either \code{source} or
#' \code{reference} be removed? If not, missing values may be treated as a
#' standard label or the function may throw an error (exact behaviour depends
#' on whether the input labels are numeric or not).
#' @param ignoreLabels labels in \code{source} and \code{reference} to be
#' considered unmatchable. These labels are excluded from the re-labeling
#' procedure.
#' @param extraLabels a vector of labels for modules in \code{source} that
#' cannot be matched to any modules in \code{reference}. The user should ensure
#' that this vector contains enough labels since the function automatically
#' removes a values that occur in either \code{source}, \code{reference} or
#' \code{ignoreLabels}, to avoid possible confusion.
#' @return A vector (if the input \code{source} labels are a vector) or a
#' matrix (if the input \code{source} labels are a matrix) of the new labels.
#' @author Peter Langfelder
#' @seealso \code{\link{overlapTable}} for calculation of overlap counts and
#' p-values \code{\link{standardColors}} for standard non-numeric WGCNA labels.
#' @keywords misc
#' @examples
#' 
#' set.seed(1)
#' nModules <- 10
#' nGenes <- 1000
#' labels1 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' labels2 <- sample(c(1:nModules), nGenes, replace = TRUE)
#' sample(matchLabels(labels1, labels2), 20)
#' 
#' @export matchLabels
matchLabels <- function(source, reference, pThreshold = 5e-2, na.rm = TRUE,
                        ignoreLabels = ifelse(is.numeric(reference), 0, "grey"),
                        extraLabels = ifelse(is.numeric(reference), c(1:1000),
                                             standardColors())) {

    source <- as.matrix(source)
    if (nrow(source) != length(reference)) {
        stop("Number of rows of 'source' must equal the length of 'reference'.")
    }

    result <- array(NA, dim = dim(source))
    #refMods <- as.numeric(sort(unique(reference)))
    #refMods <- refMods[!refMods %in% ignoreLabels]
    for (col in 1:ncol(source)) {
        src <- source[, col]
        tab <- overlapTable(src, reference, na.rm = na.rm,
                            ignore = ignoreLabels)
        pTab <- tab$pTable
        pOrder <- apply(pTab, 2, order)
        bestOrder <- order(apply(pTab, 2, min))

        refMods <- colnames(pTab)
        if (is.numeric(reference)) {
            refMods <- as.numeric(refMods)
        }
        sourceMods <- rownames(pTab)
        newLabels <- rep(NA, length(sourceMods))
        names(newLabels) <- sourceMods
        for (rm in 1:length(bestOrder)) {
            bestInd <- 1
            done <- FALSE
            #printFlush(paste("Looking for best match for reference module ",
            #refMods[bestOrder[rm]]))
            while (!done && bestInd < length(sourceMods)) {
                bm <- pOrder[bestInd, bestOrder[rm]]
                bp <- pTab[bm, bestOrder[rm]]
                if (bp > pThreshold) {
                    done <- TRUE
                } else if (is.na(newLabels[bm])) {
                    #newLabels[bm] <- as.numeric(refMods[bestOrder[rm]])
                    #printFlush(paste("Labeling old module ", sourceMods[bm],
                    #                 "as new module",
                    #                 refMods[bestOrder[rm]], "with p=", bp))
                    newLabels[bm] <- refMods[bestOrder[rm]]
                    done <- TRUE
                }
                bestInd <- bestInd + 1
            }
        }
        if (length(ignoreLabels) > 0) {
            newLabels.ignore <- ignoreLabels
            names(newLabels.ignore) <- ignoreLabels
            newLabels <- c(newLabels.ignore, newLabels)
        }

        unassigned <- src %in% names(newLabels)[is.na(newLabels)]
        if (any(unassigned)) {
            unassdSrcTab <- table(src[!src %in% names(newLabels)])
            unassdRank <- rank(-unassdSrcTab, ties.method = "first")

            nExtra <- sum(is.na(newLabels))
            keep <- !extraLabels %in% c(refMods, ignoreLabels, names(newLabels))
            newLabels[is.na(newLabels)] <- extraLabels[keep][1:nExtra]
        }

        result[, col] <- newLabels[match(src, names(newLabels))]
    }

    result
}
