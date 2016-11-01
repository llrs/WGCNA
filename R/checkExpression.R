# Checking the quality of the expression matrix or expression multi set

# goodGenes ####
#' Filter genes with too many missing entries
#'
#' This function checks data for missing entries and returns a list of genes
#' that have non-zero variance and pass two criteria on maximum number of
#' missing values: the fraction of missing values must be below a given
#' threshold and the total number of missing samples must be below a given
#' threshold.
#'
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.  For most data sets, the fraction of missing samples criterion
#' will be much more stringent than the absolute number of missing samples
#' criterion.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples.
#' @param useSamples optional specifications of which samples to use for the
#' check. Should be a logical vector; samples whose entries are \code{FALSE}
#' will be ignored for the missing value counts. Defaults to using all samples.
#' @param useGenes optional specifications of genes for which to perform the
#' check. Should be a logical vector; genes whose entries are \code{FALSE} will
#' be ignored. Defaults to using all genes.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of non-missing samples for a gene to be
#' considered good.
#' @param minNGenes minimum number of good genes for the data set to be
#' considered fit for analysis. If the actual number of good genes falls below
#' this threshold, an error will be issued.
#' @param tol an optional 'small' number to compare the variance against.
#' Defaults to the square of \code{1e-10 * max(abs(datExpr), na.rm = TRUE)}.
#' The reason of comparing the variance to this number, rather than zero, is
#' that the fast way of computing variance used by this function sometimes
#' causes small numerical overflow errors which make variance of constant
#' vectors slightly non-zero; comparing the variance to \code{tol} rather than
#' zero prevents the retaining of such genes as 'good genes'.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A logical vector with one entry per gene that is \code{TRUE} if the
#' gene is considered good and \code{FALSE} otherwise. Note that all genes
#' excluded by \code{useGenes} are automatically assigned \code{FALSE}.
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{goodSamples}}, \code{\link{goodSamplesGenes}}
#' @keywords misc
goodGenes <- function(datExpr,
                      useSamples = NULL,
                      useGenes = NULL,
                      minFraction = 1 / 2,
                      minNSamples = ..minNSamples,
                      minNGenes = ..minNGenes,
                      tol = NULL,
                      verbose = 1,
                      indent = 0) {
    datExpr = as.matrix(datExpr)
    if (is.atomic(datExpr) && (mode(datExpr) != 'numeric'))
        stop("datExpr must contain numeric data.")

    if (is.null(tol))
        tol = 1e-10 * max(abs(datExpr), na.rm = TRUE)
    if (is.null(useGenes))
        useGenes = rep(TRUE, ncol(datExpr))
    if (is.null(useSamples))
        useSamples = rep(TRUE, nrow(datExpr))

    if (length(useGenes) != ncol(datExpr))
        stop("Length of nGenes is not compatible with number of columns in
             datExpr.")
    if (length(useSamples) != nrow(datExpr))
        stop("Length of nSamples is not compatible with number of rows in
             datExpr.")

    nSamples = sum(useSamples)
    nGenes = sum(useGenes)
    nPresent = colSums(!is.na(datExpr[useSamples, useGenes]))
    gg = useGenes
    gg[useGenes][nPresent < minNSamples] = FALSE
    var = colVars(datExpr[useSamples, gg, drop = FALSE], na.rm = TRUE)
    var[is.na(var)] = 0
    nNAsGenes = colSums(is.na(datExpr[useSamples, gg]))
    gg[gg] = (
        nNAsGenes < (1 - minFraction) * nSamples & var > tol ^ 2 & (nSamples - nNAsGenes >= minNSamples)
    )
    if (sum(gg) < minNGenes)
        stop("Too few genes with valid expression levels in the required
             number of samples.")

    if (verbose > 0 & (nGenes - sum(gg) > 0))
        printFlush(
            paste(
                "  ..Excluding",
                nGenes - sum(gg),
                "genes from the calculation due to too many missing
                samples or zero variance."
            )
            )

    gg
}

# goodSamples ####
#' Filter samples with too many missing entries
#'
#' This function checks data for missing entries and returns a list of samples
#' that pass two criteria on maximum number of missing values: the fraction of
#' missing values must be below a given threshold and the total number of
#' missing genes must be below a given threshold.
#'
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.  For most data sets, the fraction of missing samples criterion
#' will be much more stringent than the absolute number of missing samples
#' criterion.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples.
#' @param useSamples optional specifications of which samples to use for the
#' check. Should be a logical vector; samples whose entries are \code{FALSE}
#' will be ignored for the missing value counts. Defaults to using all samples.
#' @param useGenes optional specifications of genes for which to perform the
#' check. Should be a logical vector; genes whose entries are \code{FALSE} will
#' be ignored. Defaults to using all genes.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of good samples for the data set to be
#' considered fit for analysis. If the actual number of good samples falls
#' below this threshold, an error will be issued.
#' @param minNGenes minimum number of non-missing samples for a sample to be
#' considered good.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A logical vector with one entry per sample that is \code{TRUE} if
#' the sample is considered good and \code{FALSE} otherwise. Note that all
#' samples excluded by \code{useSamples} are automatically assigned
#' \code{FALSE}.
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{goodSamples}}, \code{\link{goodSamplesGenes}}
#' @keywords misc
goodSamples <- function(datExpr,
                        useSamples = NULL,
                        useGenes = NULL,
                        minFraction = 1 / 2,
                        minNSamples = ..minNSamples,
                        minNGenes = ..minNGenes,
                        verbose = 1,
                        indent = 0) {
    if (is.null(useGenes))
        useGenes = rep(TRUE, ncol(datExpr))
    if (is.null(useSamples))
        useSamples = rep(TRUE, nrow(datExpr))

    if (length(useGenes) != ncol(datExpr))
        stop("Length of nGenes is not compatible with number of columns
             in datExpr.")
    if (length(useSamples) != nrow(datExpr))
        stop("Length of nSamples is not compatible with number of rows
             in datExpr.")

    nSamples = sum(useSamples)
    nGenes = sum(useGenes)
    nNAsSamples = rowSums(is.na(datExpr[useSamples, useGenes, drop = FALSE]))
    goodSamples = useSamples
    goodSamples[useSamples] = ((nNAsSamples < (1 - minFraction) * nGenes) &
                                   (nGenes - nNAsSamples >= minNGenes))
    if (sum(goodSamples) < minNSamples)
        stop("Too few samples with valid expression levels for the required
             number of genes.")

    if (verbose > 0 & (nSamples - sum(goodSamples) > 0))
        printFlush(
            paste(
                "  ..Excluding",
                nSamples - sum(goodSamples),
                "samples from the calculation due to too many
                missing genes."
            )
            )

    goodSamples
}

# goodGenesMS ####
#' Filter genes with too many missing entries across multiple sets
#'
#' This function checks data for missing entries and returns a list of genes
#' that have non-zero variance in all sets and pass two criteria on maximum
#' number of missing values in each given set: the fraction of missing values
#' must be below a given threshold and the total number of missing samples must
#' be below a given threshold
#'
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.  For most data sets, the fraction of missing samples criterion
#' will be much more stringent than the absolute number of missing samples
#' criterion.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param useSamples optional specifications of which samples to use for the
#' check. Should be a logical vector; samples whose entries are \code{FALSE}
#' will be ignored for the missing value counts. Defaults to using all samples.
#' @param useGenes optional specifications of genes for which to perform the
#' check. Should be a logical vector; genes whose entries are \code{FALSE} will
#' be ignored. Defaults to using all genes.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of non-missing samples for a gene to be
#' considered good.
#' @param minNGenes minimum number of good genes for the data set to be
#' considered fit for analysis. If the actual number of good genes falls below
#' this threshold, an error will be issued.
#' @param tol an optional 'small' number to compare the variance against. For
#' each set in \code{multiExpr}, the default value is \code{1e-10 *
#' max(abs(multiExpr[[set]]$data), na.rm = TRUE)}.  The reason of comparing the
#' variance to this number, rather than zero, is that the fast way of computing
#' variance used by this function sometimes causes small numerical overflow
#' errors which make variance of constant vectors slightly non-zero; comparing
#' the variance to \code{tol} rather than zero prevents the retaining of such
#' genes as 'good genes'.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A logical vector with one entry per gene that is \code{TRUE} if the
#' gene is considered good and \code{FALSE} otherwise. Note that all genes
#' excluded by \code{useGenes} are automatically assigned \code{FALSE}.
#' @author Peter Langfelder
#' @seealso \code{\link{goodGenes}}, \code{\link{goodSamples}},
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately
#'
#' \code{\link{goodSamplesMS}}, \code{\link{goodSamplesGenesMS}} for additional
#' cleaning of multiple data sets together.
#' @keywords misc
goodGenesMS <- function(multiExpr,
                        useSamples = NULL,
                        useGenes = NULL,
                        minFraction = 1 / 2,
                        minNSamples = ..minNSamples,
                        minNGenes = ..minNGenes,
                        tol = NULL,
                        verbose = 1,
                        indent = 0) {
    dataSize = checkSets(multiExpr)
    nSets = dataSize$nSets
    if (is.null(useGenes))
        useGenes = rep(TRUE, dataSize$nGenes)
    if (is.null(useSamples))
    {
        useSamples = list()
        for (set in 1:nSets)
            useSamples[[set]] = rep(TRUE, dataSize$nSamples[set])
    }

    if (length(useGenes) != dataSize$nGenes)
        stop("Length of nGenes is not compatible with number of genes in
             multiExpr.")
    if (length(useSamples) != nSets)
        stop("Length of nSamples is not compatible with number of sets in
             multiExpr.")

    for (set in 1:nSets)
        if (length(useSamples[[set]]) != dataSize$nSamples[set])
            stop(
                paste(
                    "Number of samples in useSamples[[",
                    set,
                    "]] incompatible\n   ",
                    "with number of samples in the corresponding set of multiExpr."
                )
            )

    nSamples = sapply(useSamples, sum)
    nGenes = sum(useGenes)

    goodGenes = useGenes
    for (set in 1:nSets) {
        if (is.null(tol)) {
            tol1 = 1e-10 * max(abs(multiExpr[[set]]$data),
                               na.rm = TRUE)
        } else {
            tol1 = tol
        }
        if (sum(goodGenes) == 0)
            break
        if (sum(useSamples[[set]]) == 0)
            next
        expr1 = multiExpr[[set]]$data[useSamples[[set]], goodGenes, drop = FALSE]
        if (mode(expr1) == "list")
            expr1 = as.matrix(expr1)
        nPresent = colSums(!is.na(expr1))
        goodGenes[goodGenes] = (nPresent >= minNGenes)
        expr1 = expr1[, nPresent >= minNGenes, drop = FALSE]
        if (any(goodGenes)) {
            var = colVars(expr1, na.rm = TRUE)
            nNAsGenes = colSums(is.na(expr1))
            goodGenes[goodGenes][nNAsGenes > (1 - minFraction) * nSamples[set] |
                                     var <= tol1 ^ 2 |
                                     (nSamples[set] - nNAsGenes < minNSamples)] = FALSE
        }
    }
    if (sum(goodGenes) < minNGenes)
        stop(
            "Too few genes with valid expression levels in the required number
            of samples in all sets."
        )

    if (verbose > 0 & (nGenes - sum(goodGenes) > 0))
        printFlush(
            paste(
                "  ..Excluding",
                nGenes - sum(goodGenes),
                "genes from the calculation due to too many missing
                samples or zero variance."
            )
            )
    goodGenes
}

# goodSamplesMS ####
#' Filter samples with too many missing entries across multiple data sets
#'
#' This function checks data for missing entries and returns a list of samples
#' that pass two criteria on maximum number of missing values: the fraction of
#' missing values must be below a given threshold and the total number of
#' missing genes must be below a given threshold.
#'
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.  For most data sets, the fraction of missing samples criterion
#' will be much more stringent than the absolute number of missing samples
#' criterion.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param useSamples optional specifications of which samples to use for the
#' check. Should be a logical vector; samples whose entries are \code{FALSE}
#' will be ignored for the missing value counts. Defaults to using all samples.
#' @param useGenes optional specifications of genes for which to perform the
#' check. Should be a logical vector; genes whose entries are \code{FALSE} will
#' be ignored. Defaults to using all genes.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of good samples for the data set to be
#' considered fit for analysis. If the actual number of good samples falls
#' below this threshold, an error will be issued.
#' @param minNGenes minimum number of non-missing samples for a sample to be
#' considered good.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with one component per input set. Each component is a logical
#' vector with one entry per sample in the corresponding set, indicating
#' whether the sample passed the missing value criteria.
#' @author Peter Langfelder and Steve Horvath
#' @seealso \code{\link{goodGenes}}, \code{\link{goodSamples}},
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately
#'
#' \code{\link{goodGenesMS}}, \code{\link{goodSamplesGenesMS}} for additional
#' cleaning of multiple data sets together.
#' @keywords misc
goodSamplesMS <- function(multiExpr,
                          useSamples = NULL,
                          useGenes = NULL,
                          minFraction = 1 / 2,
                          minNSamples = ..minNSamples,
                          minNGenes = ..minNGenes,
                          verbose = 1,
                          indent = 0) {
    dataSize = checkSets(multiExpr)
    nSets = dataSize$nSets
    if (is.null(useGenes))
        useGenes = rep(TRUE, dataSize$nGenes)
    if (is.null(useSamples))
    {
        useSamples = list()
        for (set in 1:nSets)
            useSamples[[set]] = rep(TRUE, dataSize$nSamples[set])
    }

    if (length(useGenes) != dataSize$nGenes)
        stop("Length of nGenes is not compatible with number of genes in
             multiExpr.")
    if (length(useSamples) != dataSize$nSets)
        stop("Length of nSamples is not compatible with number of sets
             in multiExpr.")

    for (set in 1:nSets)
        if (length(useSamples[[set]]) != dataSize$nSamples[set])
            stop(
                paste(
                    "Number of samples in useSamples[[",
                    set,
                    "]] incompatible\n   ",
                    "with number of samples in the corresponding set of multiExpr."
                )
            )

    nSamples = sapply(useSamples, sum)
    nGenes = sum(useGenes)

    goodSamples = useSamples
    for (set in 1:nSets) {
        if (sum(useGenes) == 0)
            break
        if (sum(goodSamples[[set]]) == 0)
            next
        nNAsSamples = rowSums(is.na(multiExpr[[set]]$data[useSamples[[set]],
                                                          useGenes, drop = FALSE]))
        goodSamples[[set]][useSamples[[set]]]  =
            ((nNAsSamples < (1 - minFraction) * nGenes) & (nGenes - nNAsSamples >= minNGenes))
        if (sum(goodSamples[[set]]) < minNSamples)
            stop(
                "Too few samples with valid expression levels for the required n
                umber of genes in set",
                set
            )
        if (verbose > 0 &
            (nSamples[set] - sum(goodSamples[[set]]) > 0))
            printFlush(
                paste(
                    "  ..Set",
                    set,
                    ": Excluding",
                    nSamples[set] -
                        sum(goodSamples[[set]]),
                    "samples from the calculation due to too many missing genes."
                )
            )
    }
    goodSamples
}

# goodSamplesGenes ####
#' Iterative filtering of samples and genes with too many missing entries
#'
#' This function checks data for missing entries and zero-variance genes, and
#' returns a list of samples and genes that pass criteria maximum number of
#' missing values. If necessary, the filtering is iterated.
#'
#' This function iteratively identifies samples and genes with too many missing
#' entries and genes with zero variance. Iterations may be required since
#' excluding samples effectively changes criteria on genes and vice versa. The
#' process is repeated until the lists of good samples and genes are stable.
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of non-missing samples for a gene to be
#' considered good.
#' @param minNGenes minimum number of good genes for the data set to be
#' considered fit for analysis. If the actual number of good genes falls below
#' this threshold, an error will be issued.
#' @param tol an optional 'small' number to compare the variance against.
#' Defaults to the square of \code{1e-10 * max(abs(datExpr), na.rm = TRUE)}.
#' The reason of comparing the variance to this number, rather than zero, is
#' that the fast way of computing variance used by this function sometimes
#' causes small numerical overflow errors which make variance of constant
#' vectors slightly non-zero; comparing the variance to \code{tol} rather than
#' zero prevents the retaining of such genes as 'good genes'.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return
#'
#' A list with the foolowing components: \item{goodSamples}{ A logical vector
#' with one entry per sample that is \code{TRUE} if the sample is considered
#' good and \code{FALSE} otherwise.  }
#'
#' \item{goodGenes}{ A logical vector with one entry per gene that is
#' \code{TRUE} if the gene is considered good and \code{FALSE} otherwise.  }
#' @author Peter Langfelder
#' @seealso \code{\link{goodSamples}}, \code{\link{goodGenes}}
#' @keywords misc
goodSamplesGenes <- function(datExpr,
                             minFraction = 1 / 2,
                             minNSamples = ..minNSamples,
                             minNGenes = ..minNGenes,
                             tol = NULL,
                             verbose = 1,
                             indent = 0) {
    spaces = indentSpaces(indent)
    goodGenes = NULL
    goodSamples = NULL
    nBadGenes = 0
    nBadSamples = 0
    changed = TRUE
    iter = 1
    if (verbose > 0)
        printFlush(paste(
            spaces,
            "Flagging genes and samples with too many
            missing values..."
        ))
    while (changed) {
        if (verbose > 0)
            printFlush(paste(spaces, " ..step", iter))
        goodGenes = goodGenes(
            datExpr,
            goodSamples,
            goodGenes,
            minFraction = minFraction,
            minNSamples = minNSamples,
            minNGenes = minNGenes,
            tol = tol,
            verbose = verbose - 1,
            indent = indent + 1
        )
        goodSamples = goodSamples(
            datExpr,
            goodSamples,
            goodGenes,
            minFraction = minFraction,
            minNSamples = minNSamples,
            minNGenes = minNGenes,
            verbose = verbose - 1,
            indent = indent + 1
        )
        changed = ((sum(!goodGenes) > nBadGenes) |
                       (sum(!goodSamples) > nBadSamples))
        nBadGenes = sum(!goodGenes)
        nBadSamples = sum(!goodSamples)
        iter = iter + 1
    }
    allOK = (sum(c(nBadGenes, nBadSamples)) == 0)
    list(goodGenes = goodGenes,
         goodSamples = goodSamples,
         allOK = allOK)
}

# goodSamplesGenesMS ####
#' Iterative filtering of samples and genes with too many missing entries
#' across multiple data sets
#'
#' This function checks data for missing entries and zero variance across
#' multiple data sets and returns a list of samples and genes that pass
#' criteria maximum number of missing values. If necessary, the filtering is
#' iterated.
#'
#' This function iteratively identifies samples and genes with too many missing
#' entries, and genes with zero variance. Iterations may be required since
#' excluding samples effectively changes criteria on genes and vice versa. The
#' process is repeated until the lists of good samples and genes are stable.
#' The constants \code{..minNSamples} and \code{..minNGenes} are both set to
#' the value 4.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param minFraction minimum fraction of non-missing samples for a gene to be
#' considered good.
#' @param minNSamples minimum number of non-missing samples for a gene to be
#' considered good.
#' @param minNGenes minimum number of good genes for the data set to be
#' considered fit for analysis. If the actual number of good genes falls below
#' this threshold, an error will be issued.
#' @param tol an optional 'small' number to compare the variance against. For
#' each set in \code{multiExpr}, the default value is \code{1e-10 *
#' max(abs(multiExpr[[set]]$data), na.rm = TRUE)}.  The reason of comparing the
#' variance to this number, rather than zero, is that the fast way of computing
#' variance used by this function sometimes causes small numerical overflow
#' errors which make variance of constant vectors slightly non-zero; comparing
#' the variance to \code{tol} rather than zero prevents the retaining of such
#' genes as 'good genes'.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the foolowing components: \item{goodSamples}{ A list
#' with one component per given set. Each component is a logical vector with
#' one entry per sample in the corresponding set that is \code{TRUE} if the
#' sample is considered good and \code{FALSE} otherwise.  }
#'
#' \item{goodGenes}{ A logical vector with one entry per gene that is
#' \code{TRUE} if the gene is considered good and \code{FALSE} otherwise.  }
#' @author Peter Langfelder
#' @seealso \code{\link{goodGenes}}, \code{\link{goodSamples}},
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately
#'
#' \code{\link{goodSamplesMS}}, \code{\link{goodGenesMS}} for additional
#' cleaning of multiple data sets together.
#' @keywords misc
goodSamplesGenesMS <- function(multiExpr,
                               minFraction = 1 / 2,
                               minNSamples = ..minNSamples,
                               minNGenes = ..minNGenes,
                               tol = NULL,
                               verbose = 2,
                               indent = 0) {
    spaces = indentSpaces(indent)
    size = checkSets(multiExpr)
    nSets = size$nSets
    goodGenes = NULL
    goodSamples = NULL
    nBadGenes = 0
    nBadSamples = rep(0, nSets)
    changed = TRUE
    iter = 1
    if (verbose > 0)
        printFlush(paste(
            spaces,
            "Flagging genes and samples with too many missing
            values..."
        ))
    while (changed) {
        if (verbose > 0)
            printFlush(paste(spaces, " ..step", iter))
        goodGenes = goodGenesMS(
            multiExpr,
            goodSamples,
            goodGenes,
            minFraction = minFraction,
            minNSamples = minNSamples,
            minNGenes = minNGenes,
            tol = tol,
            verbose = verbose - 1,
            indent = indent + 1
        )
        goodSamples = goodSamplesMS(
            multiExpr,
            goodSamples,
            goodGenes,
            minFraction = minFraction,
            minNSamples = minNSamples,
            minNGenes = minNGenes,
            verbose = verbose - 1,
            indent = indent + 1
        )
        changed = FALSE
        for (set in 1:nSets)
            changed = (changed |
                           (sum(!goodGenes) > nBadGenes) |
                           (sum(!goodSamples[[set]]) > nBadSamples[set]))
        nBadGenes = sum(!goodGenes)
        for (set in 1:nSets)
            nBadSamples[set] = sum(!goodSamples[[set]])
        iter = iter + 1
        if (verbose > 2)
            printFlush(paste(
                spaces,
                "   ..bad gene count: ",
                nBadGenes,
                ", bad sample counts: ",
                paste0(nBadSamples, collapse = ", ")
            ))
    }
    allOK = (sum(c(nBadGenes, nBadSamples)) == 0)
    list(goodGenes = goodGenes,
         goodSamples = goodSamples,
         allOK = allOK)
}

# modified heatmap plot: allow specifying the hang parameter for both side and
# top dendrograms
.heatmap <- function(x,
                     Rowv = NULL,
                     Colv = if (symm)
                         "Rowv"
                     else
                         NULL,
                     distfun = dist,
                     hclustfun = fastcluster::hclust,
                     reorderfun = function(d, w) {
                         reorder(d, w)
                     },
                     add.expr,
                     symm = FALSE,
                     revC = identical(Colv, "Rowv"),
                     scale = c("row", "column", "none"),
                     na.rm = TRUE,
                     margins = c(1.2, 1.2),
                     ColSideColors,
                     RowSideColors,
                     cexRow = 0.2  +
                         1 / log10(nr),
                     cexCol = 0.2 + 1 / log10(nc),
                     labRow = NULL,
                     labCol = NULL,
                     main = NULL,
                     xlab = NULL,
                     ylab = NULL,
                     keep.dendro = FALSE,
                     verbose = getOption("verbose"),
                     setLayout = TRUE,
                     hang = 0.04,
                     ...) {
    scale <- if (symm && missing(scale))
        "none"
    else
        match.arg(scale)
    if (length(di <- dim(x))  != 2 || !is.numeric(x))
        stop("'x' must be a numeric matrix")
    nr <- di[1L]
    nc <- di[2L]
    if (nr <= 1 || nc <= 1)
        stop("'x' must have at least 2 rows and 2 columns")
    if (!is.numeric(margins) || length(margins)  != 2L)
        stop("'margins' must be a numeric vector of length 2")

    doRdend <- !identical(Rowv, NA)
    doCdend <- !identical(Colv, NA)
    if (!doRdend && identical(Colv, "Rowv"))
        doCdend <- FALSE
    ## by default order by row/col means
    if (is.null(Rowv))
        Rowv <- rowMeans(x, na.rm = na.rm)
    if (is.null(Colv))
        Colv <- colMeans(x, na.rm = na.rm)

    ## get the dendrograms and reordering indices
    if (doRdend) {
        if (inherits(Rowv, "dendrogram"))
            ddr <- Rowv
        else {
            hcr <- hclustfun(distfun(x))
            if (class(hcr) == 'hclust')
            {
                hcr$height = hcr$height - min(hcr$height) + hang * (max(hcr$height) - min(hcr$height))
            }
            ddr <- as.dendrogram(hcr, hang = hang)
            if (!is.logical(Rowv) || Rowv)
                ddr <- reorderfun(ddr, Rowv)
        }
        if (nr  != length(rowInd <- order.dendrogram(ddr)))
            stop("row dendrogram ordering gave index of wrong length")
    }
    else
        rowInd <- 1:nr
    if (doCdend) {
        if (inherits(Colv, "dendrogram"))
            ddc <- Colv
        else if (identical(Colv, "Rowv")) {
            if (nr  != nc)
                stop("Colv = \"Rowv\" but nrow(x)  != ncol(x)")
            ddc <- ddr
        }
        else {
            hcc <- hclustfun(distfun(if (symm)
                x
                else
                    t(x)))
            if (class(hcr) == 'hclust')
            {
                hcc$height = hcc$height - min(hcc$height) + hang * (max(hcc$height) - min(hcc$height))
            }
            ddc <- as.dendrogram(hcc, hang = hang)
            if (!is.logical(Colv) ||
                Colv)
                ddc <- reorderfun(ddc, Colv)
        }
        if (nc  != length(colInd <- order.dendrogram(ddc)))
            stop("column dendrogram ordering gave index of wrong length")
    }
    else
        colInd <- 1:nc

    ## reorder x
    x <- x[rowInd, colInd]

    labRow <- if (is.null(labRow))
        if (is.null(rownames(x)))
            (1:nr)[rowInd]
    else
        rownames(x)
    else
        labRow[rowInd]
    labCol <- if (is.null(labCol))
        if (is.null(colnames(x)))
            (1:nc)[colInd]
    else
        colnames(x)
    else
        labCol[colInd]
    if (scale == "row") {
        x <- sweep(x, 1, rowMeans(x, na.rm = na.rm))
        sx <- apply(x, 1, sd, na.rm = na.rm)
        x <- sweep(x, 1, sx, "/")
    }
    else if (scale == "column") {
        x <- sweep(x, 2, colMeans(x, na.rm = na.rm))
        sx <- apply(x, 2, sd, na.rm = na.rm)
        x <- sweep(x, 2, sx, "/")
    }

    ## Calculate the plot layout
    lmat <- rbind(c(NA, 3), 2:1)
    lwid <- c(if (doRdend)
        1
        else
            0.05, 4)
    lhei <-
        c((if (doCdend)
            1
           else
               0.05) + if (!is.null(main))
                   0.5
          else
              0, 4)
    if (!missing(ColSideColors)) {
        if (!is.character(ColSideColors) || length(ColSideColors)  != nc)
            stop("'ColSideColors' must be a character vector of length ncol(x)")
        lmat <- rbind(lmat[1,] + 1, c(NA, 1), lmat[2,] + 1)
        lhei <- c(lhei[1], 0.2, lhei[2])
    }
    if (!missing(RowSideColors)) {
        if (!is.character(RowSideColors) || length(RowSideColors)  != nr)
            stop("'RowSideColors' must be a character vector of length nrow(x)")
        lmat <- cbind(lmat[, 1] + 1, c(rep(NA, nrow(lmat) - 1), 1),
                      lmat[, 2] + 1)
        lwid <- c(lwid[1], 0.2, lwid[2])
    }
    lmat[is.na(lmat)] <- 0
    if (verbose) {
        cat("layout: widths = ", lwid, ", heights = ", lhei, "; lmat = \n")
        print(lmat)
    }
    if (!symm || scale  != "none")
        x <- t(x)
    op <- par(no.readonly = TRUE)
    if (revC) {
        iy <- nc:1
        ddr <- rev(ddr)
        rowInd.colors = rev(rowInd)
        x <- x[, iy]
    } else
        iy <- 1:nr
    #on.exit(par(op))
    # print(paste("main:", main))
    if (setLayout)
        layout(lmat,
               widths = lwid,
               heights = lhei,
               respect = TRUE)
    if (!missing(RowSideColors)) {
        par(mar = c(margins[1], 0, 0, 0.5))
        image(rbind(1:nr), col = RowSideColors[rowInd.colors], axes = FALSE)
    }
    if (!missing(ColSideColors)) {
        par(mar = c(0.5, 0, 0, margins[2]))
        image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    }
    par(mar = c(margins[1], 0, 0, margins[2]))
    image(
        1:nc,
        1:nr,
        x,
        xlim = 0.5 + c(0, nc),
        ylim = 0.5 + c(0, nr),
        axes = FALSE,
        xlab = "",
        ylab = "",
        ...
    )
    axis(
        1,
        1:nc,
        labels = labCol,
        las = 2,
        line = -0.5,
        tick = 0,
        cex.axis = cexCol
    )
    if (!is.null(xlab))
        mtext(xlab, side = 1, line = margins[1] - 1.25)
    axis(
        4,
        iy,
        labels = labRow,
        las = 2,
        line = -0.5,
        tick = 0,
        cex.axis = cexRow
    )
    if (!is.null(ylab))
        mtext(ylab, side = 4, line = margins[2] - 1.25)
    if (!missing(add.expr))
        eval.parent(substitute(add.expr))
    par(mar = c(margins[1], 0, 0, 0))
    if (doRdend) {
        .plotDendrogram(
            as.hclust(ddr),
            horiz = TRUE,
            labels = FALSE,
            axes = FALSE
        )
        #    plot(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    }
    else
        frame()
    par(mar = c(0, 0, if (!is.null(main))
        1.8
        else
            0, margins[2]))
    if (doCdend)
    {
        .plotDendrogram(
            as.hclust(ddc),
            horiz = FALSE,
            labels = FALSE,
            axes = FALSE
        )
        #    plot(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    }
    else if (!is.null(main))
        frame()
    if (!is.null(main))
        title(main, cex.main = 1.2 * op[["cex.main"]])
    invisible(list(
        rowInd = rowInd,
        colInd = colInd,
        Rowv = if (keep.dendro && doRdend)
            ddr,
        Colv = if (keep.dendro && doCdend)
            ddc
    ))
}
