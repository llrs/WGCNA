#===============================================================================
# The function PickHardThreshold can help one to estimate the cut - off value
# when using the signum (step) function.
# The first column lists the threshold ("cut"), the second column lists the
# corresponding p-value based on the Fisher Transform of the correlation.
# The third column reports the resulting scale free topology fitting index R^2.
# The fourth column reports the slope of the fitting line, it shoud be negative
# for biologically meaningul networks.
# The fifth column reports the fitting index for the truncated exponential model
# usually we ignore it.
# The remaining columns list the mean, median and maximum resulting connectivity.
# To pick a hard threshold (cut) with the scale free topology criterion:
# aim for high scale free R^2 (column 3), high connectivity (col 6) and negative
# slope (around  - 1, col 4).
# The output is a list with 2 components. The first component lists a sugggested
#  cut-off while the second component contains the whole table.
# The removeFirst option removes the first point (k = 0, P(k = 0)) from the
# regression fit.
# nBreaks specifies how many intervals used to estimate the frequency p(k) i.e.
# the no. of points in the scale free topology plot.



#' Analysis of scale free topology for hard-thresholding.
#'
#' Analysis of scale free topology for multiple hard thresholds. The aim is to
#' help the user pick an appropriate threshold for network construction.
#'
#' The function calculates unsigned networks by thresholding the correlation
#' matrix using thresholds given in \code{cutVector}. For each power the scale
#' free topology fit index is calculated and returned along with other
#' information on connectivity.
#'
#' @aliases pickHardThreshold pickHardThreshold.fromSimilarity
#' @param data expression data in a matrix or data frame. Rows correspond to
#' samples and columns to genes.
#' @param dataIsExpr logical: should the data be interpreted as expression (or
#' other numeric) data, or as a similarity matrix of network nodes?
#' @param RsquaredCut desired minimum scale free topology fitting index
#' \eqn{R^2}.
#' @param cutVector a vector of hard threshold cuts for which the scale free
#' topology fit indices are to be calculated.
#' @param moreNetworkConcepts logical: should additional network concepts be
#' calculated? If \code{TRUE}, the function will calculate how the network
#' density, the network heterogeneity, and the network centralization depend on
#' the power. For the definition of these additional network concepts, see
#' Horvath and Dong (2008).  PloS Comp Biol.
#' @param removeFirst should the first bin be removed from the connectivity
#' histogram?
#' @param nBreaks number of bins in connectivity histograms
#' @param corFnc a character string giving the correlation function to be used
#' in adjacency calculation.
#' @param corOptions further options to the correlation function specified in
#' \code{corFnc}.
#' @return A list with the following components:
#'
#' \item{cutEstimate}{ estimate of an appropriate hard-thresholding cut: the
#' lowest cut for which the scale free topology fit \eqn{R^2} exceeds
#' \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} for all cuts,
#' \code{NA} is returned. }
#'
#' \item{fitIndices}{ a data frame containing the fit indices for scale free
#' topology. The columns contain the hard threshold, Student p-value for the
#' correlation threshold, adjusted \eqn{R^2} for the linear fit, the linear
#' coefficient, adjusted \eqn{R^2} for a more complicated fit models, mean
#' connectivity, median connectivity and maximum connectivity.  If input
#' \code{moreNetworkConcepts} is \code{TRUE}, 3 additional columns containing
#' network density, centralization, and heterogeneity.}
#' @author Steve Horvath
#' @seealso \code{\link{signumAdjacency}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @keywords misc
pickHardThreshold <- function(data, dataIsExpr = TRUE, RsquaredCut = 0.85,
                              cutVector = seq(0.1, 0.9,
                                              by = 0.05), moreNetworkConcepts = FALSE, removeFirst = FALSE, nBreaks = 10,
                              corFnc = "cor",
                              corOptions = "use = 'p'") {
    nGenes = dim(data)[[2]]
    colname1 = c("Cut", "p-value", "SFT.R.sq", "slope = ",
                 "truncated R^2", "mean(k)", "median(k)", "max(k)")
    if (moreNetworkConcepts) {
        colname1 = c(colname1, "Density", "Centralization", "Heterogeneity")
    }
    if (!dataIsExpr)
    {
        checkAdjMat(data)
        if (any(diag(data) != 1)) diag(data) = 1
    } else
        nSamples = dim(data)[[1]]

    datout = data.frame(matrix(nrow = length(cutVector),
                               ncol = length(colname1)))
    names(datout) = colname1
    datout[, 1] = cutVector
    if (dataIsExpr)
    {
        for (i in 1:length(cutVector))
        {
            cut1 = cutVector[i]
            datout[i, 2] = 2 * (1 - pt(sqrt(nSamples - 1) * cut1/sqrt(1  -
                                                                          cut1^2), nSamples - 1))
        }
    } else
        datout[, 2] = NA

    fun1 <- function(x, dataIsExpr) {
        if (dataIsExpr)
        {
            corExpr = parse(text = paste(corFnc, "(x, data",
                                         prepComma(corOptions), ")"))
            corx = abs(eval(corExpr))
        } else
            corx = x
        out1 = rep(NA, length(cutVector))
        for (j in c(1:length(cutVector))) {
            out1[j] = sum(corx > cutVector[j], na.rm = TRUE)
        }
        out1
    }
    datk = t(apply(data, 2, fun1, dataIsExpr))
    for (i in c(1:length(cutVector))) {
        khelp = datk[, i] - 1
        SFT1 = scaleFreeFitIndex(k = khelp, nBreaks = nBreaks,
                                 removeFirst = removeFirst)
        datout[i, 3] = SFT1$Rsquared.SFT
        datout[i, 4] = SFT1$slope.SFT
        datout[i, 5] = SFT1$truncatedExponentialAdjRsquared
        datout[i, 6] = mean(khelp, na.rm = TRUE)
        datout[i, 7] = median(khelp, na.rm = TRUE)
        datout[i, 8] = max(khelp, na.rm = TRUE)
        if (moreNetworkConcepts) {
            Density = sum(khelp)/(nGenes * (nGenes - 1))
            datout[i, 9]  = Density
            Centralization = nGenes * (max(khelp) - mean(khelp))/(
                (nGenes - 1) * (nGenes - 2))
            datout[i, 10] = Centralization
            Heterogeneity = sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 1)
            datout[i, 11] = Heterogeneity
        }
    }
    print(signif (data.frame(datout), 3))
    ind1 = datout[, 3] > RsquaredCut
    indcut = NA
    indcut = if (sum(ind1) > 0) min(c(1:length(ind1))[ind1]) else indcut
    cutEstimate = cutVector[indcut][[1]]
    list(cutEstimate = cutEstimate, fitIndices = data.frame(datout))
} # end of function pickHardThreshold


#===============================================================================
#
# pickSoftThreshold
#
#===============================================================================
# The function pickSoftThreshold allows one to estimate the power parameter when
# using a soft thresholding approach with the use of the power function
# AF(s) = s^Power
# The removeFirst option removes the first point (k = 1, P(k = 1)) from the
# regression fit.
# PL: a rewrite that splits the data into a few blocks.
# SH: more netowkr concepts added.
# PL: re - written for parallel processing



#' Analysis of scale free topology for soft-thresholding
#'
#' Analysis of scale free topology for multiple soft thresholding powers. The
#' aim is to help the user pick an appropriate soft-thresholding power for
#' network construction.
#'
#' The function calculates weighted networks either by interpreting \code{data}
#' directly as similarity, or first transforming it to similarity of the type
#' specified by \code{networkType}. The weighted networks are obtained by
#' raising the similarity to the powers given in \code{powerVector}.  For each
#' power the scale free topology fit index is calculated and returned along
#' with other information on connectivity.
#'
#' On systems with multiple cores or processors, the function pickSoftThreshold
#' takes advantage of parallel processing if the function
#' \code{\link{enableWGCNAThreads}} has been called to allow parallel
#' processing and set up the parallel calculation back-end.
#'
#' @aliases pickSoftThreshold pickSoftThreshold.fromSimilarity
#' @rdname pickSoftThreshold
#' @param data expression data in a matrix or data frame. Rows correspond to
#' samples and columns to genes.
#' @param dataIsExpr logical: should the data be interpreted as expression (or
#' other numeric) data, or as a similarity matrix of network nodes?
#' @param RsquaredCut desired minimum scale free topology fitting index
#' \eqn{R^2}.
#' @param powerVector a vector of soft thresholding powers for which the scale
#' free topology fit indices are to be calculated.
#' @param removeFirst should the first bin be removed from the connectivity
#' histogram?
#' @param nBreaks number of bins in connectivity histograms
#' @param blockSize block size into which the calculation of connectivity
#' should be broken up. If not given, a suitable value will be calculated using
#' function \code{blockSize} and printed if \code{verbose>0}. If R runs into
#' memory problems, decrease this value.
#' @param corFnc the correlation function to be used in adjacency calculation.
#' @param corOptions a list giving further options to the correlation function
#' specified in \code{corFnc}.
#' @param networkType network type. Allowed values are (unique abbreviations
#' of) \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}. See
#' \code{\link{adjacency}}.
#' @param moreNetworkConcepts logical: should additional network concepts be
#' calculated? If \code{TRUE}, the function will calculate how the network
#' density, the network heterogeneity, and the network centralization depend on
#' the power. For the definition of these additional network concepts, see
#' Horvath and Dong (2008). PloS Comp Biol.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{powerEstimate}{ estimate of an appropriate soft-thresholding power:
#' the lowest power for which the scale free topology fit \eqn{R^2} exceeds
#' \code{RsquaredCut}. If \eqn{R^2} is below \code{RsquaredCut} for all powers,
#' \code{NA} is returned. }
#'
#' \item{fitIndices}{ a data frame containing the fit indices for scale free
#' topology. The columns contain the soft-thresholding power, adjusted
#' \eqn{R^2} for the linear fit, the linear coefficient, adjusted \eqn{R^2} for
#' a more complicated fit models, mean connectivity, median connectivity and
#' maximum connectivity. If input \code{moreNetworkConcepts} is \code{TRUE}, 3
#' additional columns containing network density, centralization, and
#' heterogeneity.}
#' @author Steve Horvath and Peter Langfelder
#' @seealso \code{\link{adjacency}}, \code{\link{softConnectivity}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @keywords misc
pickSoftThreshold <- function(data, dataIsExpr = TRUE, RsquaredCut = 0.85,
                              powerVector = c(seq(1, 10, by = 1),
                                              seq(12, 20, by = 2)),
                              removeFirst = FALSE, nBreaks = 10,
                              blockSize = NULL, corFnc = cor,
                              corOptions = list(use = 'p'),
                              networkType = "unsigned",
                              moreNetworkConcepts = FALSE, verbose = 0,
                              indent = 0) {
    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType)) {
        stop("Unrecognized 'networkType'. Recognized values are",
             paste(.networkTypes, collapse = ", "))
    }
    nGenes <- ncol(data)
    if (nGenes < 3) {
        stop("The input data data contain fewer than 3 rows (nodes).",
             "\nThis would result in a trivial correlation network.")
    }
    if (!dataIsExpr) {
        checkSimilarity(data)
        if (any(diag(data) != 1)) {
            diag(data) <- 1
        }
    }

    if (is.null(blockSize)) {
        blockSize <- blockSize(nGenes, rectangularBlocks = TRUE,
                               maxMemoryAllocation = 2^30)
        if (verbose > 0) {
            printFlush(paste0("pickSoftThreshold: will use block size ",
                              blockSize, "."))
        }
    }

    colname1 <- c("Power", "SFT.R.sq", "slope", "truncated R.sq",
                  "mean(k)", "median(k)", "max(k)")
    if (moreNetworkConcepts) {
        colname1 <- c(colname1, "Density", "Centralization", "Heterogeneity")
    }
    datout <- data.frame(matrix(666, nrow = length(powerVector),
                                ncol = length(colname1)))
    names(datout) <- colname1
    datout$Power <- powerVector
    spaces <- indentSpaces(indent)
    if (verbose > 0) {
        cat(paste(spaces,
                  "pickSoftThreshold: calculating ",
                  "connectivity for given powers..."))
        if (verbose == 1) {
            pind = initProgInd()
        } else {
            cat("\n")
        }
    }

    # if we're using one of WGNCA's own correlation functions, set the number of
    #  threads to 1.
    corFnc <- match.fun(corFnc)
    corFormals <- formals(corFnc)
    if ("nThreads" %in% names(corFormals)) {
        corOptions$nThreads <- 1
    }

    # Resulting connectivities
    datk <- matrix(0, nrow = nGenes, ncol = length(powerVector))

    # Number of threads. In this case I need this explicitly.
    nThreads <- WGCNAnThreads()

    nPowers <- length(powerVector)

    # Main loop
    startG <- 1
    while (startG <= nGenes) {
        endG <- min (startG + blockSize - 1, nGenes)

        if (verbose > 1) {
            printFlush(paste(spaces, "  ..working on genes", startG, "through",
                             endG, "of", nGenes))
        }

        nBlockGenes <- endG - startG + 1
        jobs <- allocateJobs(nBlockGenes, nThreads)
        # This assumes that the non - zero length allocations
        # precede the zero - length ones
        actualThreads <- which(sapply(jobs, length) > 0)

        datk[c(startG:endG), ] <- foreach(
            t <- actualThreads, .combine = rbind) %dopar% {
                useGenes <- c(startG:endG)[jobs[[t]]]
                nGenes1 <- length(useGenes)
                if (dataIsExpr) {
                    corOptions$x <- data
                    corOptions$y <- data[, useGenes]
                    corx <- do.call(corFnc, corOptions)
                    if (intType == 1) {
                        corx <- abs(corx)
                    } else if (intType == 2) {
                        corx <- (1 + corx)/2
                    } else if (intType == 3) {
                        corx[corx < 0] <- 0
                    }
                    if (sum(is.na(corx))  != 0)
                        warning(paste("Some correlations are NA in block",
                                      startG, ":", endG, "."))
                } else {
                    corx <- data[, useGenes]
                }
                datk.local <- matrix(nGenes1, nPowers)
                for (j in 1:nPowers) {
                    datk.local[, j] <- colSums(corx^powerVector[j],
                                               na.rm = TRUE) - 1
                }
                datk.local
            } # End of %dopar% evaluation
        # Move to the next block of genes.
        startG <- endG + 1
        if (verbose == 1) {
            pind <- updateProgInd(endG/nGenes, pind)
        }
    }
    if (verbose == 1) {
        printFlush("")
    }

    for (i in c(1:length(powerVector))) {
        khelp <- datk[, i]
        SFT1 <- scaleFreeFitIndex(k = khelp, nBreaks = nBreaks,
                                  removeFirst = removeFirst)
        datout[i, "SFT.R.sq"] <- SFT1$Rsquared.SFT
        datout[i, "slope"] <- SFT1$slope.SFT
        datout[i, "truncated R.sq"] <- SFT1$truncatedExponentialAdjRsquared
        datout[i, "mean(k)"] <- mean(khelp, na.rm = TRUE)
        datout[i, "median(k)"] <- median(khelp, na.rm = TRUE)
        datout[i, "max(k)"] <- max(khelp, na.rm = TRUE)
        if (moreNetworkConcepts) {
            Density <- sum(khelp)/(nGenes * (nGenes - 1))
            datout[i, "Density"]  <- Density
            Centralization <- nGenes * (max(khelp) - mean(khelp))/(
                (nGenes - 1) * (nGenes - 2))
            datout[i, "Centralization"] <- Centralization
            Heterogeneity <- sqrt(nGenes * sum(khelp^2)/sum(khelp)^2 - 1)
            datout[i, "Heterogenity"] <- Heterogeneity
        }
    }
    print(signif (data.frame(datout), 3))
    ind1 <- datout[, "SFT.R.sq"] > RsquaredCut
    indcut <- NA
    indcut = ifelse(sum(ind1) > 0, min(c(1:length(ind1))[ind1]), indcut)
    powerEstimate <- powerVector[indcut][[1]]
    list(powerEstimate = powerEstimate, fitIndices = data.frame(datout))
}
