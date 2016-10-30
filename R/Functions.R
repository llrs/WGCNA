# A collection of functions

# getMEprefix ####
#' Get the prefix used to label module eigengenes.
#'
#' Returns the currently used prefix used to label module eigengenes.  When
#' returning module eigengenes in a dataframe, names of the corresponding
#' columns will start with the given prefix.
#'
#' Returns the prefix used to label module eigengenes. When returning module
#' eigengenes in a dataframe, names of the corresponding columns will consist
#' of the corresponfing color label preceded by the given prefix. For example,
#' if the prefix is "PC" and the module is turquoise, the corresponding module
#' eigengene will be labeled "PCturquoise". Most of old code assumes "PC", but
#' "ME" is more instructive and used in some newer analyses.
#'
#' @return A character string.
#' @note Currently the standard prefix is \code{"ME"} and there is no way to
#' change it.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{moduleEigengenes}}
#' @keywords misc
getMEprefix <- function() {
    .moduleColorOptions$MEprefix
}

# moduleEigengenes ####
#' Calculate module eigengenes.
#'
#' Calculates module eigengenes (1st principal component) of modules in a given
#' single dataset.
#'
#' Module eigengene is defined as the first principal component of the
#' expression matrix of the corresponding module. The calculation may fail if
#' the expression data has too many missing entries. Handling of such errors is
#' controlled by the arguments \code{subHubs} and \code{trapErrors}.  If
#' \code{subHubs==TRUE}, errors in principal component calculation will be
#' trapped and a substitute calculation of hubgenes will be attempted. If this
#' fails as well, behaviour depends on \code{trapErrors}: if \code{TRUE}, the
#' offending module will be ignored and the return value will allow the user to
#' remove the module from further analysis; if \code{FALSE}, the function will
#' stop.
#'
#' From the user's point of view, setting \code{trapErrors=FALSE} ensures that
#' if the function returns normally, there will be a valid eigengene (principal
#' component or hubgene) for each of the input colors. If the user sets
#' \code{trapErrors=TRUE}, all calculational (but not input) errors will be
#' trapped, but the user should check the output (see below) to make sure all
#' modules have a valid returned eigengene.
#'
#' While the principal component calculation can fail even on relatively sound
#' data (it does not take all that many "well-placed" \code{NA} to torpedo the
#' calculation), it takes many more irregularities in the data for the hubgene
#' calculation to fail. In fact such a failure signals there likely is
#' something seriously wrong with the data.
#'
#' @param datExpr Expression data for a single set in the form of a data frame
#' where rows are samples and columns are genes (probes).
#' @param colors A vector of the same length as the number of probes in
#' \code{expr}, giving module color for all probes (genes). Color \code{"grey"}
#' is reserved for unassigned genes.
#' @param impute If \code{TRUE}, expression data will be checked for the
#' presence of \code{NA} entries and if the latter are present, numerical data
#' will be imputed, using function \code{impute.knn} and probes from the same
#' module as the missing datum. The function \code{impute.knn} uses a fixed
#' random seed giving repeatable results.
#' @param nPC Number of principal components and variance explained entries to
#' be calculated. Note that only the first principal component is returned; the
#' rest are used only for the calculation of proportion of variance explained.
#' The number of returned variance explained entries is currently
#' \code{min(nPC, 10)}. If given \code{nPC} is greater than 10, a warning is
#' issued.
#' @param align Controls whether eigengenes, whose orientation is undetermined,
#' should be aligned with average expression (\code{align = "along average"},
#' the default) or left as they are (\code{align = ""}). Any other value will
#' trigger an error.
#' @param excludeGrey Should the improper module consisting of 'grey' genes be
#' excluded from the eigengenes?
#' @param grey Value of \code{colors} designating the improper module. Note
#' that if \code{colors} is a factor of numbers, the default value will be
#' incorrect.
#' @param subHubs Controls whether hub genes should be substituted for missing
#' eigengenes. If \code{TRUE}, each missing eigengene (i.e., eigengene whose
#' calculation failed and the error was trapped) will be replaced by a weighted
#' average of the most connected hub genes in the corresponding module. If this
#' calculation fails, or if \code{subHubs==FALSE}, the value of
#' \code{trapErrors} will determine whether the offending module will be
#' removed or whether the function will issue an error and stop.
#' @param trapErrors Controls handling of errors from that may arise when there
#' are too many \code{NA} entries in expression data. If \code{TRUE}, errors
#' from calling these functions will be trapped without abnormal exit.  If
#' \code{FALSE}, errors will cause the function to stop. Note, however, that
#' \code{subHubs} takes precedence in the sense that if \code{subHubs==TRUE}
#' and \code{trapErrors==FALSE}, an error will be issued only if both the
#' principal component and the hubgene calculations have failed.
#' @param returnValidOnly logical; controls whether the returned data frame of
#' module eigengenes contains columns corresponding only to modules whose
#' eigengenes or hub genes could be calculated correctly (\code{TRUE}), or
#' whether the data frame should have columns for each of the input color
#' labels (\code{FALSE}).
#' @param softPower The power used in soft-thresholding the adjacency matrix.
#' Only used when the hubgene approximation is necessary because the principal
#' component calculation failed. It must be non-negative. The default value
#' should only be changed if there is a clear indication that it leads to
#' incorrect results.
#' @param scale logical; can be used to turn off scaling of the expression data
#' before calculating the singular value decomposition. The scaling should only
#' be turned off if the data has been scaled previously, in which case the
#' function can run a bit faster. Note however that the function first imputes,
#' then scales the expression data in each module. If the expression contain
#' missing data, scaling outside of the function and letting the function
#' impute missing data may lead to slightly different results than if the data
#' is scaled within the function.
#' @param verbose Controls verbosity of printed progress messages. 0 means
#' silent, up to (about) 5 the verbosity gradually increases.
#' @param indent A single non-negative integer controlling indentation of
#' printed messages. 0 means no indentation, each unit above that adds two
#' spaces.
#' @return A list with the following components: \item{eigengenes}{Module
#' eigengenes in a dataframe, with each column corresponding to one eigengene.
#' The columns are named by the corresponding color with an \code{"ME"}
#' prepended, e.g., \code{MEturquoise} etc. If \code{returnValidOnly==FALSE},
#' module eigengenes whose calculation failed have all components set to
#' \code{NA}.} \item{averageExpr}{If \code{align == "along average"}, a
#' dataframe containing average normalized expression in each module. The
#' columns are named by the corresponding color with an \code{"AE"} prepended,
#' e.g., \code{AEturquoise} etc.} \item{varExplained}{A dataframe in which each
#' column corresponds to a module, with the component \code{varExplained[PC,
#' module]} giving the variance of module \code{module} explained by the
#' principal component no. \code{PC}. The calculation is exact irrespective of
#' the number of computed principal components. At most 10 variance explained
#' values are recorded in this dataframe.} \item{nPC}{A copy of the input
#' \code{nPC}.} \item{validMEs}{A boolean vector. Each component (corresponding
#' to the columns in \code{data}) is \code{TRUE} if the corresponding eigengene
#' is valid, and \code{FALSE} if it is invalid. Valid eigengenes include both
#' principal components and their hubgene approximations. When
#' \code{returnValidOnly==FALSE}, by definition all returned eigengenes are
#' valid and the entries of \code{validMEs} are all \code{TRUE}. }
#' \item{validColors}{A copy of the input colors with entries corresponding to
#' invalid modules set to \code{grey} if given, otherwise 0 if \code{colors} is
#' numeric and "grey" otherwise.} \item{allOK}{Boolean flag signalling whether
#' all eigengenes have been calculated correctly, either as principal
#' components or as the hubgene average approximation.} \item{allPC}{Boolean
#' flag signalling whether all returned eigengenes are principal components.}
#' \item{isPC}{Boolean vector. Each component (corresponding to the columns in
#' \code{eigengenes}) is \code{TRUE} if the corresponding eigengene is the
#' first principal component and \code{FALSE} if it is the hubgene
#' approximation or is invalid.} \item{isHub}{Boolean vector. Each component
#' (corresponding to the columns in \code{eigengenes}) is \code{TRUE} if the
#' corresponding eigengene is the hubgene approximation and \code{FALSE} if it
#' is the first principal component or is invalid.} \item{validAEs}{Boolean
#' vector. Each component (corresponding to the columns in \code{eigengenes})
#' is \code{TRUE} if the corresponding module average expression is valid.}
#' \item{allAEOK}{Boolean flag signalling whether all returned module average
#' expressions contain valid data. Note that \code{returnValidOnly==TRUE} does
#' not imply \code{allAEOK==TRUE}: some invalid average expressions may be
#' returned if their corresponding eigengenes have been calculated correctly.}
#' @author Steve Horvath \email{SHorvath@@mednet.ucla.edu}, Peter Langfelder
#' \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{svd}}, \code{impute.knn}
#' @references Zhang, B. and Horvath, S. (2005), "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
moduleEigengenes <- function(datExpr, colors, impute = TRUE, nPC = 1,
                             align = "along average", excludeGrey = FALSE,
                             grey = ifelse(is.numeric(colors), 0, "grey"),
                             subHubs = TRUE, trapErrors = FALSE,
                             returnValidOnly = trapErrors, softPower = 6,
                             scale = TRUE, verbose = 0, indent = 0) {
    spaces = indentSpaces(indent)

    if (verbose == 1) {
        printFlush(spaces, " moduleEigengenes: Calculating ",
                   nlevels(as.factor(colors)),
                   " module eigengenes in given set.")
    }
    if (is.null(datExpr)) {
        stop("moduleEigengenes: Error: datExpr is NULL.")
    }
    if (is.null(colors)) {
        stop("moduleEigengenes: Error: colors is NULL.")
    }
    if (is.null(dim(datExpr)) || length(dim(datExpr)) != 2)
        stop("moduleEigengenes: Error: datExpr must be two - dimensional.")

    if (dim(datExpr)[2] != length(colors)) {
        stop("moduleEigengenes: Error: ncol(datExpr) and length(colors)",
             "must be equal (one color per gene).")
    }

    if (is.factor(colors)) {
        nl = nlevels(colors)
        nlDrop = nlevels(colors[, drop = TRUE])
        if (nl > nlDrop) {
            stop(paste("Argument 'colors' contains unused levels ",
                       "(empty modules).",
                       "Use colors[, drop = TRUE] to get rid of them."))
        }
    }

    if (softPower < 0) {
        stop("softPower must be non - negative")
    }

    alignRecognizedValues = c("", "along average")
    if (!is.element(align, alignRecognizedValues)) {
        stop(paste("ModulePrincipalComponents: Error:",
                   "parameter align has an unrecognised value:",
                   align, "; Recognized values are ",
                   alignRecognizedValues)
        )
    }

    maxVarExplained = 10
    if (nPC > maxVarExplained) {
        warning("Given nPC is too large. Will use value", maxVarExplained)
    }

    nVarExplained <- min(nPC, maxVarExplained)
    modlevels <- levels(factor(colors))
    if (excludeGrey) {
        if (sum(as.character(modlevels) != as.character(grey)) > 0) {
            modlevels <-
                modlevels[as.character(modlevels) != as.character(grey)]
        } else {
            stop("Color levels are empty. Possible reason: ",
                 "the only color is grey and grey module is excluded ",
                 "from the calculation.")
        }
    }
    lml <- length(modlevels)
    PrinComps <- data.frame(matrix(nrow = dim(datExpr)[[1]], ncol = lml))
    averExpr <- data.frame(matrix(nrow = dim(datExpr)[[1]], ncol = lml))
    varExpl <- data.frame(matrix(nrow = nVarExplained, ncol = lml))
    validMEs <- rep(TRUE, lml)
    validAEs <- rep(FALSE, lml)
    isPC <- validMEs
    isHub <- validAEs
    validColors <- colors
    names(PrinComps) <- paste0(getMEprefix(), modlevels)
    names(averExpr) <- paste0("AE", modlevels)
    for (i in c(1:lml)) {
        if (verbose > 1) {
            printFlush(paste(spaces,
                             "moduleEigengenes : Working on ME for module",
                             modlevels[i]))
        }
        modulename <- modlevels[i]
        restrict1 <- as.character(colors) == as.character(modulename)
        if (verbose > 2) {
            printFlush(paste(spaces, " ...", sum(restrict1), "genes"))
        }
        datModule <- as.matrix(t(datExpr[, restrict1]))
        n <- nrow(datModule)
        p <- ncol(datModule)
        pc <- try({
            if (nrow(datModule) > 1 && impute) {
                seedSaved <- FALSE
                if (exists(".Random.seed")) {
                    saved.seed <- .Random.seed
                    seedSaved <- TRUE
                }
                if (any(is.na(datModule))) {
                    if (verbose > 5) {
                        printFlush(paste(spaces,
                                         " ...imputing missing data"))
                    }
                    datModule <- impute.knn(datModule,
                                            k = min(10,
                                                    nrow(datModule) - 1))
                    # some versions of impute.knn return a list and we
                    # need the data component:
                    try({
                        if (!is.null(datModule$data)) {
                            datModule <- datModule$data
                        }
                    }, silent = TRUE)
                }
                # The <<- in the next line is extremely important.
                # Using = or <- will create a local variable of
                # the name .Random.seed and will leave the important global
                # .Random.seed untouched.
                if (seedSaved) {
                    .Random.seed <<- saved.seed
                }
            }

            if (verbose > 5) {
                printFlush(paste(spaces, " ...scaling"))
            }
            if (scale) {
                datModule <- t(scale(t(datModule)))
            }
            if (verbose > 5) {
                printFlush(paste(spaces, " ...calculating SVD"))
            }
            svd1 <- svd(datModule,
                        nu = min(n, p, nPC),
                        nv = min(n, p, nPC))
            # varExpl[, i] = (svd1$d[1:min(n, p,
            # nVarExplained)])^2/sum(svd1$d^2)
            if (verbose > 5) {
                printFlush(paste(spaces, " ...calculating PVE"))
            }
            veMat <- cor(svd1$v[, c(1:min(n, p, nVarExplained))],
                         t(datModule), use = "p")
            varExpl[c(1:min(n, p, nVarExplained)),
                    i] <-  rowMeans(veMat ^ 2, na.rm = TRUE)
            # this is the first principal component
            svd1$v[, 1]
        }, silent = TRUE)
        if (class(pc) == 'try - error') {
            if ((!subHubs) && (!trapErrors)) {
                stop(pc)
            }
            if (subHubs) {
                if (verbose > 0) {
                    printFlush(
                        paste(spaces,
                              " ..principal component calculation for module",
                              modulename, "failed with the following error:"))
                    printFlush(paste(spaces, "     ", pc, spaces,
                                     "..hub genes will be used instead of ",
                                     "principal components."))
                }
                isPC[i] <- FALSE
                pc <- try({
                    scaledExpr <- scale(t(datModule))
                    covEx <- cov(scaledExpr, use = "p")
                    covEx[!is.finite(covEx)] <- 0
                    modAdj <- abs(covEx) ^ softPower
                    kIM <- (rowMeans(modAdj, na.rm = TRUE)) ^ 3
                    if (max(kIM, na.rm = TRUE) > 1) {
                        kIM <- kIM - 1
                    }
                    kIM[is.na(kIM)] <- 0
                    hub <- which.max(kIM)
                    alignSign <- sign(covEx[, hub])
                    alignSign[is.na(alignSign)] <- 0
                    isHub[i] <- TRUE
                    pcxMat <- scaledExpr  *
                        matrix(
                            kIM * alignSign,
                            nrow = nrow(scaledExpr),
                            ncol = ncol(scaledExpr),
                            byrow = TRUE
                        ) / sum(kIM)
                    pcx <- rowMeans(pcxMat, na.rm = TRUE)
                    varExpl[1, i] <- mean(cor(pcx, t(datModule),
                                              use = "p") ^ 2, na.rm = TRUE)
                    pcx
                }, silent = TRUE)
            }
        }

        if (class(pc) == 'try - error') {
            if (!trapErrors) {
                stop(pc)
            }
            if (verbose > 0) {
                printFlush(paste(spaces, " ..ME calculation of module",
                                 modulename, "failed with the following error:"))
                printFlush(paste(spaces, "     ", pc, spaces,
                                 " ..the offending module has been removed."))
            }
            warning("Eigengene calculation of module ", modulename,
                    " failed with the following error \n     ", pc,
                    " The offending module has been removed.\n")
            validMEs[i] <- FALSE
            isPC[i] <- FALSE
            isHub[i] <- FALSE
            validColors[restrict1] <- grey
        } else {
            PrinComps[, i] <- pc
            ae <- try({if (isPC[i]) {
                scaledExpr <- scale(t(datModule))
            }
                averExpr[, i] <- rowMeans(scaledExpr, na.rm = TRUE)
                if (align == "along average") {
                    if (verbose > 4) {
                        printFlush(paste(spaces,
                                         " .. aligning module eigengene",
                                         "with average expression."))
                    }
                    corAve <- cor(averExpr[, i], PrinComps[, i], use = "p")
                    if (!is.finite(corAve)) {
                        corAve <- 0
                    }
                    if (corAve < 0) {
                        PrinComps[, i] <- -PrinComps[, i]
                    }
                }
                0
            }, silent = TRUE)
            if (class(ae) == 'try - error') {
                if (!trapErrors) {
                    stop(ae)
                }
                if (verbose > 0) {
                    printFlush(paste(spaces,
                                     " ..Average expression calculation of module",
                                     modulename, "failed with the following error:"))
                    printFlush(paste(spaces, "     ", ae, spaces,
                                     " ..the returned average expression",
                                     "vector will be invalid."))
                }
                warning(paste("Average expression calculation of module",
                              modulename,
                              "failed with the following error \n     ", ae,
                              "The returned average expression vector will",
                              "be invalid.\n"))
            }
            validAEs[i] <- !(class(ae) == 'try - error')
        }
    }
    allOK <- sum(!validMEs) == 0
    if (returnValidOnly && sum(!validMEs) > 0) {
        PrinComps <- PrinComps[, validMEs]
        averExpr <- averExpr[, validMEs]
        varExpl <- varExpl[, validMEs]
        validMEs <- rep(TRUE, times = ncol(PrinComps))
        isPC <- isPC[validMEs]
        isHub <- isHub[validMEs]
        validAEs <- validAEs[validMEs]
    }
    allPC <- (sum(!isPC) == 0)
    allAEOK <- (sum(!validAEs) == 0)
    list(eigengenes = PrinComps,
         averageExpr = averExpr,
         varExplained = varExpl,
         nPC = nPC,
         validMEs = validMEs,
         validColors = validColors,
         allOK = allOK,
         allPC = allPC,
         isPC = isPC,
         isHub = isHub,
         validAEs = validAEs,
         allAEOK = allAEOK)
}

# removeGrey ####
#' Removes the grey eigengene from a given collection of eigengenes.
#'
#' Given module eigengenes either in a single data frame or in a multi-set
#' format, removes the grey eigengenes from each set. If the grey eigengenes
#' are not found, a warning is issued.
#'
#'
#' @param MEs Module eigengenes, either in a single data frame (typicaly for a
#' single set), or in a multi-set format. See \code{\link{checkSets}} for a
#' description of the multi-set format.
#' @param greyMEName Name of the module eigengene (in each corresponding data
#' frame) that corresponds to the grey color. This will typically be "PCgrey"
#' or "MEgrey". If the module eigengenes were calculated using standard
#' functions in this library, the default should work.
#' @return Module eigengenes in the same format as input (either a single data
#' frame or a vector of lists) with the grey eigengene removed.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
removeGreyME <- function(MEs, greyMEName = paste0(getMEprefix(), "grey")) {
    newMEs = MEs
    if (is.vector(MEs) & mode(MEs) == "list") {
        warned <- 0
        newMEs <- vector(mode = "list", length = length(MEs))
        for (set in 1:length(MEs)) {
            if (!is.data.frame(MEs[[set]]$data)) {
                stop("MEs is a vector list but the list structure ",
                     "is missing the correct 'data' component.")
            }
            newMEs[[set]] <- MEs[[set]]
            if (greyMEName %in% names(MEs[[set]]$data)) {
                newMEs[[set]]$data <- MEs[[set]]$data[,
                                                      names(MEs[[set]]$data) != greyMEName]
            } else {
                if (warned == 0) {
                    warning("removeGreyME: The given grey ME name was ",
                            "not found ",
                            "among the names of given MEs.")
                    warned <- 1
                }
            }
        }
    } else {
        if (length(dim(MEs)) != 2) {
            stop("Argument 'MEs' has incorrect dimensions.")
        }
        MEs = as.data.frame(MEs)
        if (greyMEName %in% names(MEs)) {
            newMEs = MEs[, names(MEs) != greyMEName]
        } else {
            warning(
                "removeGreyME: The given grey ME name was not ",
                "found among the names of given MEs."
            )
        }
    }
    newMEs
}

# collectGarbage ####
#' Iterative garbage collection.
#'
#' Performs garbage collection until free memory idicators show no change.
#'
#'
#' @return None.
#' @author Steve Horvath
#' @keywords utilities
collectGarbage <- function() {
    while (gc()[2, 4]  != gc()[2, 4] | gc()[1, 4]  != gc()[1, 4]) {
    }
}

# orderMEs ####
#' Put close eigenvectors next to each other
#'
#' Reorder given (eigen-)vectors such that similar ones (as measured by
#' correlation) are next to each other.
#'
#' Ordering module eigengenes is useful for plotting purposes. For this
#' function the order can be specified explicitly, or a set can be given in
#' which the correlations of the eigengenes will determine the order. For the
#' latter, a hierarchical dendrogram is calculated and the order given by the
#' dendrogram is used for the eigengenes in all other sets.
#'
#' @param MEs Module eigengenes in a multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, with each list corresponding to
#' one dataset and the module eigengenes in the component \code{data}, that is
#' \code{MEs[[set]]$data[sample, module]} is the expression of the eigengene of
#' module \code{module} in sample \code{sample} in dataset \code{set}. The
#' number of samples can be different between the sets, but the modules must be
#' the same.
#' @param greyLast Normally the color grey is reserved for unassigned genes;
#' hence the grey module is not a proper module and it is conventional to put
#' it last. If this is not desired, set the parameter to \code{FALSE}.
#' @param greyName Name of the grey module eigengene.
#' @param orderBy Specifies the set by which the eigengenes are to be ordered
#' (in all other sets as well). Defaults to the first set in \code{useSets} (or
#' the first set, if \code{useSets} is not given).
#' @param order Allows the user to specify a custom ordering.
#' @param useSets Allows the user to specify for which sets the eigengene
#' ordering is to be performed.
#' @param verbose Controls verbostity of printed progress messages. 0 means
#' silent, nonzero verbose.
#' @param indent A single non-negative integer controling indentation of
#' printed messages. 0 means no indentation, each unit above zero adds two
#' spaces.
#' @return A vector of lists of the same type as \code{MEs} containing the
#' re-ordered eigengenes.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso
#' \code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}},
#' \code{\link{consensusOrderMEs}}
#' @keywords misc
orderMEs <- function(MEs,
                     greyLast = TRUE,
                     greyName = paste0(getMEprefix(), "grey"),
                     orderBy = 1,
                     order = NULL,
                     useSets = NULL,
                     verbose = 0,
                     indent = 0) {
    spaces <- indentSpaces(indent)

    if ("eigengenes" %in% names(MEs)) {
        if (is.null(order)) {
            if (verbose > 0) {
                printFlush(paste(spaces,
                                 "orderMEs: order not given, calculating ",
                                 "using given set", orderBy))
            }
            corPC <- cor(MEs$eigengenes, use = "p")
            disPC <- 1 - corPC
            order <- .clustOrder(disPC, greyLast = greyLast,
                                 greyName = greyName)
        }

        if (length(order) != dim(MEs$eigengenes)[2]) {
            stop("orderMEs: given MEs and order have incompatible dimensions.")
        }
        orderedMEs <- MEs
        orderedMEs$eigengenes <- as.data.frame(MEs$eigengenes[, order])
        colnames(orderedMEs$eigengenes) <- colnames(MEs$eigengenes)[order]
        if (!is.null(MEs$averageExpr)) {
            orderedMEs$averageExpr <- as.data.frame(MEs$averageExpr[, order])
            colnames(orderedMEs$averageExpr) <- colnames(MEs$data)[order]
        }
        if (!is.null(MEs$varExplained)) {
            orderedMEs$varExplained <- as.data.frame(MEs$varExplained[, order])
            colnames(orderedMEs$varExplained) <- colnames(MEs$data)[order]
        }
        return(orderedMEs)
    } else {
        check = checkSets(MEs, checkStructure = TRUE, useSets = useSets)
        if (check$structureOK) {
            multiSet <- TRUE
        } else {
            multiSet <- FALSE
            MEs <- fixDataStructure(MEs)
            useSets <- NULL
            orderBy <- 1
        }

        if (!is.null(useSets)) {
            if (is.na(match(orderBy, useSets))) {
                orderBy = useSets[1]
            }

            if (is.null(order)) {
                if (verbose > 0) {
                    printFlush(paste(spaces,
                            "orderMEs: order not given, calculating using given set",
                            orderBy))
                }
                corPC = cor(MEs[[orderBy]]$data, use = "p")
                disPC = 1 - corPC
                order = .clustOrder(disPC,
                                    greyLast = greyLast,
                                    greyName = greyName)
            }

            if (length(order) != dim(MEs[[orderBy]]$data)[2])
                stop("orderMEs: given MEs and order have incompatible dimensions.")

            nSets <- length(MEs)
            orderedMEs <- MEs
            if (is.null(useSets)) {
                useSets = c(1:nSets)
            }
            for (set in useSets) {
                orderedMEs[[set]]$data = as.data.frame(MEs[[set]]$data[, order])
                colnames(orderedMEs[[set]]$data) = colnames(MEs[[set]]$data)[order]
                if (!is.null(MEs[[set]]$averageExpr)) {
                    orderedMEs[[set]]$averageExpr = as.data.frame(MEs[[set]]$averageExpr[, order])
                    colnames(orderedMEs[[set]]$averageExpr) = colnames(MEs[[set]]$data)[order]
                }
                if (!is.null(MEs[[set]]$varExplained)) {
                    orderedMEs[[set]]$varExplained = as.data.frame(MEs[[set]]$varExplained[, order])
                    colnames(orderedMEs[[set]]$varExplained) = colnames(MEs[[set]]$data)[order]
                }
            }
            if (multiSet) {
                return(orderedMEs)
            } else {
                return(orderedMEs[[1]]$data)
            }
        }
    }
}

# normalizeLabels ####
#' Transform numerical labels into normal order.
#'
#' Transforms numerical labels into normal order, that is the largest group
#' will be labeled 1, next largest 2 etc. Label 0 is optionally preserved.
#'
#'
#' @param labels Numerical labels.
#' @param keepZero If \code{TRUE} (the default), labels 0 are preserved.
#' @return A vector of the same length as input, containing the normalized
#' labels.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
normalizeLabels <- function(labels, keepZero = TRUE) {
    if (keepZero) {
        NonZero = (labels != 0)
    } else {
        NonZero = rep(TRUE, length(labels))
    }
    f = as.numeric(factor(labels[NonZero]))
    t = table(labels[NonZero])
    # print(t)
    r = rank(-as.vector(t), ties.method = "first")
    norm_labs = rep(0, times = length(labels))
    norm_labs[NonZero] = r[f]
    norm_labs
}

# moduleNumber ####
#' Fixed-height cut of a dendrogram.
#'
#' Detects branches of on the input dendrogram by performing a fixed-height
#' cut.
#'
#' All contiguous branches below the height \code{cutHeight} that contain at
#' least \code{minSize} objects are assigned unique positive numerical labels;
#' all unassigned objects are assigned label 0.
#'
#' @param dendro a hierarchical clustering dendorgram such as one returned by
#' \code{hclust}.
#' @param cutHeight Maximum joining heights that will be considered.
#' @param minSize Minimum cluster size.
#' @return A vector of numerical labels giving the assigment of each object.
#' @note The numerical labels may not be sequential. See
#' \code{\link{normalizeLabels}} for a way to put the labels into a standard
#' order.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso
#' \code{\link{hclust}}, \code{\link{cutree}},
#' \code{\link{normalizeLabels}}
#' @keywords cluster
moduleNumber <- function(dendro,
                         cutHeight = 0.9,
                         minSize = 50) {
    Branches = cutree(dendro, h = cutHeight)
    NOnBranches = table(Branches)
    TrueBranch = NOnBranches >= minSize
    Branches[!TrueBranch[Branches]] = 0

    Branches
}

# fixDataStructure ####
#' Put single-set data into a form useful for multiset calculations.
#'
#' Encapsulates single-set data in a wrapper that makes the data suitable for
#' functions working on multiset data collections.
#'
#' For multiset calculations, many quantities (such as expression data, traits,
#' module eigengenes etc) are presented by a common structure, a vector of
#' lists (one list for each set) where each list has a component \code{data}
#' that contains the actual (expression, trait, eigengene) data for the
#' corresponding set in the form of a dataframe. This funtion creates a vector
#' of lists of length 1 and fills the component \code{data} with the content of
#' parameter \code{data}.
#'
#' @param data A dataframe, matrix or array with two dimensions to be
#' encapsulated.
#' @param verbose Controls verbosity. 0 is silent.
#' @param indent Controls indentation of printed progress messages. 0 means no
#' indentation, every unit adds two spaces.
#' @return As described above, input data in a format suitable for functions
#' operating on multiset data collections.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{checkSets}}
#' @keywords misc
#' @examples
#'
#'
#' singleSetData = matrix(rnorm(100), 10,10);
#' encapsData = fixDataStructure(singleSetData);
#' length(encapsData)
#' names(encapsData[[1]])
#' dim(encapsData[[1]]$data)
#' all.equal(encapsData[[1]]$data, singleSetData);
#'
#' @export
fixDataStructure <- function(data, verbose = 0, indent = 0) {
    spaces = indentSpaces(indent)
    if ((class(data) != "list") || (class(data[[1]]) != "list")) {
        if (verbose > 0)
            printFlush(paste(spaces,
                             "fixDataStructure: data is not a vector of lists:",
                             "converting it into one."))
        x = data
        data = vector(mode = "list", length = 1)
        data[[1]] = list(data = x)
        rm(x)
    }
    data
}

# checkSets ####
#' Check structure and retrieve sizes of a group of datasets.
#'
#' Checks whether given sets have the correct format and retrieves dimensions.
#'
#' For multiset calculations, many quantities (such as expression data, traits,
#' module eigengenes etc) are presented by a common structure, a vector of
#' lists (one list for each set) where each list has a component \code{data}
#' that contains the actual (expression, trait, eigengene) data for the
#' corresponding set in the form of a dataframe. This funtion checks whether
#' \code{data} conforms to this convention and retrieves some basic dimension
#' information (see output).
#'
#' @param data A vector of lists; in each list there must be a component named
#' \code{data} whose content is a matrix or dataframe or array of dimension 2.
#' @param checkStructure If \code{FALSE}, incorrect structure of \code{data}
#' will trigger an error. If \code{TRUE}, an appropriate flag (see output) will
#' be set to indicate whether \code{data} has correct structure.
#' @param useSets Optional specification of entries of the vector \code{data}
#' that are to be checked. Defaults to all components. This may be useful when
#' \code{data} only contains information for some of the sets.
#' @return A list with components \item{nSets}{Number of sets (length of the
#' vector \code{data}).} \item{nGenes}{Number of columns in the \code{data}
#' components in the lists. This number must be the same for all sets.}
#' \item{nSamples}{A vector of length \code{nSets} giving the number of rows in
#' the \code{data} components.} \item{structureOK}{Only set if the argument
#' \code{checkStructure} equals \code{TRUE}.  The value is \code{TRUE} if the
#' paramter \code{data} passes a few tests of its structure, and \code{FALSE}
#' otherwise. The tests are not exhaustive and are meant to catch obvious user
#' errors rather than be bulletproof.}
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
checkSets <- function(data, checkStructure = FALSE, useSets = NULL) {
        nSets = length(data)
        if (is.null(useSets)) {
            useSets = c(1:nSets)
        }
        if (nSets <= 0) {
            stop("No data given.")
        }
        structureOK = TRUE
        if ((class(data) != "list") || (class(data[[useSets[1]]]) != "list")) {
            if (checkStructure) {
                structureOK = FALSE
                nGenes = 0
                nSamples = 0
            } else {
                stop("data does not appear to have the correct format. ",
                     "Consider using fixDataStructure or setting ",
                     "checkStructure = TRUE when calling this function.")
            }
        } else {
            nSamples = vector(length = nSets)
            nGenes = dim(data[[useSets[1]]]$data)[2]
            for (set in useSets) {
                if (nGenes != dim(data[[set]]$data)[2]) {
                    if (checkStructure) {
                        structureOK = FALSE
                    } else {
                        stop("Incompatible number of genes in set 1 and ", set)
                    }
                }
                nSamples[set] = dim(data[[set]]$data)[1]
            }
        }

        list(nSets = nSets, nGenes = nGenes, nSamples = nSamples,
             structureOK = structureOK)
    }

# multiSetMEs ####
#' Calculate module eigengenes.
#'
#' Calculates module eigengenes for several sets.
#'
#' This function calls \code{\link{moduleEigengenes}} for each set in
#' \code{exprData}.
#'
#' Module eigengene is defined as the first principal component of the
#' expression matrix of the corresponding module. The calculation may fail if
#' the expression data has too many missing entries. Handling of such errors is
#' controlled by the arguments \code{subHubs} and \code{trapErrors}. If
#' \code{subHubs==TRUE}, errors in principal component calculation will be
#' trapped and a substitute calculation of hubgenes will be attempted. If this
#' fails as well, behaviour depends on \code{trapErrors}: if \code{TRUE}, the
#' offending module will be ignored and the return value will allow the user to
#' remove the module from further analysis; if \code{FALSE}, the function will
#' stop. If \code{universalColors} is given, any offending module will be
#' removed from all sets (see \code{validMEs} in return value below).
#'
#' From the user's point of view, setting \code{trapErrors=FALSE} ensures that
#' if the function returns normally, there will be a valid eigengene (principal
#' component or hubgene) for each of the input colors. If the user sets
#' \code{trapErrors=TRUE}, all calculational (but not input) errors will be
#' trapped, but the user should check the output (see below) to make sure all
#' modules have a valid returned eigengene.
#'
#' While the principal component calculation can fail even on relatively sound
#' data (it does not take all that many "well-placed" \code{NA} to torpedo the
#' calculation), it takes many more irregularities in the data for the hubgene
#' calculation to fail. In fact such a failure signals there likely is
#' something seriously wrong with the data.
#'
#' @param exprData Expression data in a multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, with each list corresponding to
#' one microarray dataset and expression data in the component \code{data},
#' that is \code{expr[[set]]$data[sample, probe]} is the expression of probe
#' \code{probe} in sample \code{sample} in dataset \code{set}. The number of
#' samples can be different between the sets, but the probes must be the same.
#' @param colors A matrix of dimensions (number of probes, number of sets)
#' giving the module assignment of each gene in each set. The color "grey" is
#' interpreted as unassigned.
#' @param universalColors Alternative specification of module assignment. A
#' single vector of length (number of probes) giving the module assignment of
#' each gene in all sets (that is the modules are common to all sets). If
#' given, takes precedence over \code{color}.
#' @param useSets If calculations are requested in (a) selected set(s) only,
#' the set(s) can be specified here. Defaults to all sets.
#' @param useGenes Can be used to restrict calculation to a subset of genes
#' (the same subset in all sets). If given, \code{validColors} in the returned
#' list will only contain colors for the genes specified in \code{useGenes}.
#' @param impute Logical. If \code{TRUE}, expression data will be checked for
#' the presence of \code{NA} entries and if the latter are present, numerical
#' data will be imputed, using function \code{impute.knn} and probes from the
#' same module as the missing datum. The function \code{impute.knn} uses a
#' fixed random seed giving repeatable results.
#' @param nPC Number of principal components to be calculated. If only
#' eigengenes are needed, it is best to set it to 1 (default). If variance
#' explained is needed as well, use value \code{NULL}. This will cause all
#' principal components to be computed, which is slower.
#' @param align Controls whether eigengenes, whose orientation is undetermined,
#' should be aligned with average expression (\code{align = "along average"},
#' the default) or left as they are (\code{align = ""}). Any other value will
#' trigger an error.
#' @param excludeGrey Should the improper module consisting of 'grey' genes be
#' excluded from the eigengenes?
#' @param grey Value of \code{colors} or \code{universalColors} (whichever
#' applies) designating the improper module. Note that if the appropriate
#' colors argument is a factor of numbers, the default value will be incorrect.
#' @param subHubs Controls whether hub genes should be substituted for missing
#' eigengenes. If \code{TRUE}, each missing eigengene (i.e., eigengene whose
#' calculation failed and the error was trapped) will be replaced by a weighted
#' average of the most connected hub genes in the corresponding module. If this
#' calculation fails, or if \code{subHubs==FALSE}, the value of
#' \code{trapErrors} will determine whether the offending module will be
#' removed or whether the function will issue an error and stop.
#' @param trapErrors Controls handling of errors from that may arise when there
#' are too many \code{NA} entries in expression data. If \code{TRUE}, errors
#' from calling these functions will be trapped without abnormal exit. If
#' \code{FALSE}, errors will cause the function to stop. Note, however, that
#' \code{subHubs} takes precedence in the sense that if \code{subHubs==TRUE}
#' and \code{trapErrors==FALSE}, an error will be issued only if both the
#' principal component and the hubgene calculations have failed.
#' @param returnValidOnly Boolean. Controls whether the returned data frames of
#' module eigengenes contain columns corresponding only to modules whose
#' eigengenes or hub genes could be calculated correctly in every set
#' (\code{TRUE}), or whether the data frame should have columns for each of the
#' input color labels (\code{FALSE}).
#' @param softPower The power used in soft-thresholding the adjacency matrix.
#' Only used when the hubgene approximation is necessary because the principal
#' component calculation failed. It must be non-negative. The default value
#' should only be changed if there is a clear indication that it leads to
#' incorrect results.
#' @param verbose Controls verbosity of printed progress messages. 0 means
#' silent, up to (about) 5 the verbosity gradually increases.
#' @param indent A single non-negative integer controlling indentation of
#' printed messages. 0 means no indentation, each unit above that adds two
#' spaces.
#' @return A vector of lists similar in spirit to the input \code{exprData}.
#' For each set there is a list with the following components:
#' \item{data}{Module eigengenes in a data frame, with each column
#' corresponding to one eigengene. The columns are named by the corresponding
#' color with an \code{"ME"} prepended, e.g., \code{MEturquoise} etc. Note
#' that, when \code{trapErrors == TRUE} and \code{returnValidOnly==FALSE}, this
#' data frame also contains entries corresponding to removed modules, if any.
#' (\code{validMEs} below indicates which eigengenes are valid and \code{allOK}
#' whether all module eigengens were successfully calculated.) }
#' \item{averageExpr}{If \code{align == "along average"}, a dataframe
#' containing average normalized expression in each module. The columns are
#' named by the corresponding color with an \code{"AE"} prepended, e.g.,
#' \code{AEturquoise} etc.} \item{varExplained}{A dataframe in which each
#' column corresponds to a module, with the component \code{varExplained[PC,
#' module]} giving the variance of module \code{module} explained by the
#' principal component no. \code{PC}. This is only accurate if all principal
#' components have been computed (input \code{nPC = NULL}). At most 5 principal
#' components are recorded in this dataframe.} \item{nPC}{A copy of the input
#' \code{nPC}.} \item{validMEs}{A boolean vector. Each component (corresponding
#' to the columns in \code{data}) is \code{TRUE} if the corresponding eigengene
#' is valid, and \code{FALSE} if it is invalid. Valid eigengenes include both
#' principal components and their hubgene approximations. When
#' \code{returnValidOnly==FALSE}, by definition all returned eigengenes are
#' valid and the entries of \code{validMEs} are all \code{TRUE}. }
#' \item{validColors}{A copy of the input colors (\code{universalColors} if
#' set, otherwise \code{colors[, set]}) with entries corresponding to invalid
#' modules set to \code{grey} if given, otherwise 0 if the appropriate input
#' colors are numeric and "grey" otherwise.} \item{allOK}{Boolean flag
#' signalling whether all eigengenes have been calculated correctly, either as
#' principal components or as the hubgene approximation. If
#' \code{universalColors} is set, this flag signals whether all eigengenes are
#' valid in all sets.} \item{allPC}{Boolean flag signalling whether all
#' returned eigengenes are principal components. This flag (as well as the
#' subsequent ones) is set independently for each set.} \item{isPC}{Boolean
#' vector. Each component (corresponding to the columns in \code{eigengenes})
#' is \code{TRUE} if the corresponding eigengene is the first principal
#' component and \code{FALSE} if it is the hubgene approximation or is invalid.
#' } \item{isHub}{Boolean vector. Each component (corresponding to the columns
#' in \code{eigengenes}) is \code{TRUE} if the corresponding eigengene is the
#' hubgene approximation and \code{FALSE} if it is the first principal
#' component or is invalid.} \item{validAEs}{Boolean vector. Each component
#' (corresponding to the columns in \code{eigengenes}) is \code{TRUE} if the
#' corresponding module average expression is valid.} \item{allAEOK}{Boolean
#' flag signalling whether all returned module average expressions contain
#' valid data. Note that \code{returnValidOnly==TRUE} does not imply
#' \code{allAEOK==TRUE}: some invalid average expressions may be returned if
#' their corresponding eigengenes have been calculated correctly.}
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{moduleEigengenes}}
#' @keywords misc
multiSetMEs <- function(exprData,
                        colors,
                        universalColors = NULL,
                        useSets = NULL,
                        useGenes = NULL,
                        impute = TRUE,
                        nPC = 1,
                        align = "along average",
                        excludeGrey = FALSE,
                        grey = if (is.null(universalColors)) {
                            if (is.numeric(colors)) {0} else{"grey"}
                        } else {
                            if (is.numeric(universalColors)) {0} else {"grey"}},
                        subHubs = TRUE,
                        trapErrors = FALSE,
                        returnValidOnly = trapErrors,
                        softPower = 6,
                        verbose = 1,
                        indent = 0) {
    spaces = indentSpaces(indent)
    nSets = length(exprData)
    setsize = checkSets(exprData, useSets = useSets)
    nGenes = setsize$nGenes
    nSamples = setsize$nSamples
    if (verbose > 0) {
        printFlush(paste(spaces, "multiSetMEs: Calculating module MEs."))
    }
    MEs = vector(mode = "list", length = nSets)
    consValidMEs = NULL
    if (!is.null(universalColors)) {
        consValidColors = universalColors}
    if (is.null(useSets)) {
        useSets = c(1:nSets)}
    if (is.null(useGenes)) {
        for (set in useSets) {
            if (verbose > 0)
                printFlush(paste(spaces, "  Working on set", as.character(set), "..."))
            if (is.null(universalColors)) {
                setColors = colors[, set]
            } else {
                setColors = universalColors
            }
            setMEs = moduleEigengenes(
                datExpr = exprData[[set]]$data,
                colors = setColors,
                impute = impute,
                nPC = nPC,
                align = align,
                excludeGrey = excludeGrey,
                grey = grey,
                trapErrors = trapErrors,
                subHubs = subHubs,
                returnValidOnly = FALSE,
                softPower = softPower,
                verbose = verbose - 1,
                indent = indent + 1
            )
            if (!is.null(universalColors) && (!setMEs$allOK)) {
                if (is.null(consValidMEs)) {
                    consValidMEs = setMEs$validMEs
                } else {
                    consValidMEs = consValidMEs * setMEs$validMEs
                }
                consValidColors[setMEs$validColors != universalColors]  =
                    setMEs$validColors[setMEs$validColors != universalColors]
            }
            MEs[[set]] = setMEs
            names(MEs[[set]])[names(setMEs) == 'eigengenes'] = 'data'
            # Here's what moduleEigengenes returns:
            #
            #  list(eigengenes = PrinComps, averageExpr = averExpr,
            #  varExplained = varExpl, nPC = nPC,
            #       validMEs = validMEs, validColors = validColors, allOK = allOK,
            #       allPC = allPC, isPC = isPC,
            #       isHub = isHub, validAEs = validAEs, allAEOK = allAEOK)
        }
    } else {
        for (set in useSets) {
            if (verbose > 0) {
                printFlush(spaces, "  Working on set", as.character(set),
                           "...")
            }
            if (is.null(universalColors)) {
                setColors = colors[useGenes, set]
            } else {
                setColors = universalColors[useGenes]
            }
            setMEs = moduleEigengenes(
                datExpr = exprData[[set]]$data[, useGenes],
                colors = setColors,
                impute = impute,
                nPC = nPC,
                align = align,
                excludeGrey = excludeGrey,
                grey = grey,
                trapErrors = trapErrors,
                subHubs = subHubs,
                returnValidOnly = FALSE,
                softPower = softPower,
                verbose = verbose - 1,
                indent = indent + 1
            )
            if (!is.null(universalColors) && (!setMEs$allOK)) {
                if (is.null(consValidMEs)) {
                    consValidMEs = setMEs$validMEs
                } else {
                    consValidMEs = consValidMEs * setMEs$validMEs
                }
                consValidColors[
                    setMEs$validColors != universalColors[useGenes]]  =
                    setMEs$validColors[
                        setMEs$validColors != universalColors[useGenes]]
            }
            MEs[[set]] = setMEs
            names(MEs[[set]])[names(setMEs) == 'eigengenes'] = 'data'
        }
    }
    if (!is.null(universalColors)) {
        for (set in 1:nSets) {
            if (!is.null(consValidMEs))
                MEs[[set]]$validMEs = consValidMEs
            MEs[[set]]$validColors = consValidColors
        }
    }
    for (set in 1:nSets) {
        MEs[[set]]$allOK = (sum(!MEs[[set]]$validMEs) == 0)
        if (returnValidOnly) {
            valid = (MEs[[set]]$validMEs > 0)
            MEs[[set]]$data = MEs[[set]]$data[, valid]
            MEs[[set]]$averageExpr = MEs[[set]]$averageExpr[, valid]
            MEs[[set]]$varExplained = MEs[[set]]$varExplained[, valid]
            MEs[[set]]$isPC = MEs[[set]]$isPC[valid]
            MEs[[set]]$allPC = (sum(!MEs[[set]]$isPC) == 0)
            MEs[[set]]$isHub = MEs[[set]]$isHub[valid]
            MEs[[set]]$validAEs = MEs[[set]]$validAEs[valid]
            MEs[[set]]$allAEOK = (sum(!MEs[[set]]$validAEs) == 0)
            MEs[[set]]$validMEs = rep(TRUE, times = ncol(MEs[[set]]$data))
        }
    }
    MEs
}

# MergeCloseModules ####
#' Merge close modules in gene expression data
#'
#' Merges modules in gene expression networks that are too close as measured by
#' the correlation of their eigengenes.
#'
#' This function returns the color labels for modules that are obtained from
#' the input modules by merging ones that are closely related. The
#' relationships are quantified by correlations of module eigengenes; a
#' ``consensus'' measure is defined as the ``consensus quantile'' over the
#' corresponding relationship in each set. Once the (dis-)similarity is
#' calculated, average linkage hierarchical clustering of the module eigengenes
#' is performed, the dendrogram is cut at the height \code{cutHeight} and
#' modules on each branch are merged. The process is (optionally) repeated
#' until no more modules are merged.
#'
#' If, for a particular module, the module eigengene calculation fails, a
#' hubgene approximation will be used.
#'
#' The user should be aware that if a computational error occurs and
#' \code{trapErrors==TRUE}, the returned list (see below) will not contain all
#' of the components returned upon normal execution.
#'
#' @param exprData Expression data, either a single data frame with rows
#' corresponding to samples and columns to genes, or in a multi-set format (see
#' \code{\link{checkSets}}). See \code{checkDataStructure} below.
#' @param colors A vector (numeric, character or a factor) giving module colors
#' for genes.  The method only makes sense when genes have the same color label
#' in all sets, hence a single vector.
#' @param MEs If module eigengenes have been calculated before, the user can
#' save some computational time by inputting them. \code{MEs} should have the
#' same format as \code{exprData}.  If they are not given, they will be
#' calculated.
#' @param useSets A vector of scalar allowing the user to specify which sets
#' will be used to calculate the consensus dissimilarity of module eigengenes.
#' Defaults to all given sets.
#' @param impute Should missing values be imputed in eigengene calculation? If
#' imputation is disabled, the presence of \code{NA} entries will cause the
#' eigengene calculation to fail and eigengenes will be replaced by their
#' hubgene approximation. See \code{\link{moduleEigengenes}} for more details.
#' @param checkDataFormat If TRUE, the function will check \code{exprData} and
#' \code{MEs} for correct multi-set structure. If single set data is given, it
#' will be converted into a format usable for the function. If FALSE, incorrect
#' structure of input data will trigger an error.
#' @param unassdColor Specifies the string that labels unassigned genes. Module
#' of this color will not enter the module eigengene clustering and will not be
#' merged with other modules.
#' @param corFnc Correlation function to be used to calculate correlation of
#' module eigengenes.
#' @param corOptions Can be used to specify options to the correlation
#' function, in addition to argument \code{x} which is used to pass the actual
#' data to calculate the correlation of.
#' @param useAbs Specifies whether absolute value of correlation or plain
#' correlation (of module eigengenes) should be used in calculating module
#' dissimilarity.
#' @param equalizeQuantiles Logical: should quantiles of the eigengene
#' dissimilarity matrix be equalized ("quantile normalized")? The default is
#' \code{FALSE} for reproducibility of old code, but better results will
#' probably be achieved if quantile equalizatio is used.
#' @param quantileSummary One of \code{"mean"} or \code{"median"}. Controls how
#' a reference dissimilarity is computed from the input ones (using mean or
#' median, respectively).
#' @param consensusQuantile A number giving the desired quantile to use in the
#' consensus similarity calculation (see details).
#' @param cutHeight Maximum dissimilarity (i.e., 1-correlation) that qualifies
#' modules for merging.
#' @param iterate Controls whether the merging procedure should be repeated
#' until there is no change. If FALSE, only one iteration will be executed.
#' @param relabel Controls whether, after merging, color labels should be
#' ordered by module size.
#' @param colorSeq Color labels to be used for relabeling. Defaults to the
#' standard color order used in this package if \code{colors} are not numeric,
#' and to integers starting from 1 if \code{colors} is numeric.
#' @param getNewMEs Controls whether module eigengenes of merged modules should
#' be calculated and returned.
#' @param getNewUnassdME When doing module eigengene manipulations, the
#' function does not normally calculate the eigengene of the 'module' of
#' unassigned ('grey') genes. Setting this option to \code{TRUE} will force the
#' calculation of the unassigned eigengene in the returned newMEs, but not in
#' the returned oldMEs.
#' @param trapErrors Controls whether computational errors in calculating
#' module eigengenes, their dissimilarity, and merging trees should be trapped.
#' If \code{TRUE}, errors will be trapped and the function will return the
#' input colors. If \code{FALSE}, errors will cause the function to stop.
#' @param verbose Controls verbosity of printed progress messages. 0 means
#' silent, up to (about) 5 the verbosity gradually increases.
#' @param indent A single non-negative integer controlling indentation of
#' printed messages. 0 means no indentation, each unit above that adds two
#' spaces.
#' @return If no errors occurred, a list with components \item{colors}{Color
#' labels for the genes corresponding to merged modules. The function attempts
#' to mimic the mode of the input \code{colors}: if the input \code{colors} is
#' numeric, character and factor, respectively, so is the output. Note,
#' however, that if the fnction performs relabeling, a standard sequence of
#' labels will be used: integers starting at 1 if the input \code{colors} is
#' numeric, and a sequence of color labels otherwise (see \code{colorSeq}
#' above).}
#'
#' \item{dendro}{Hierarchical clustering dendrogram (average linkage) of the
#' eigengenes of the most recently computed tree. If \code{iterate} was set
#' TRUE, this will be the dendrogram of the merged modules, otherwise it will
#' be the dendrogram of the original modules.}
#'
#' \item{oldDendro}{Hierarchical clustering dendrogram (average linkage) of the
#' eigengenes of the original modules.}
#'
#' \item{cutHeight}{The input cutHeight.}
#'
#' \item{oldMEs}{Module eigengenes of the original modules in the sets given by
#' \code{useSets}.}
#'
#' \item{newMEs}{Module eigengenes of the merged modules in the sets given by
#' \code{useSets}.}
#'
#' \item{allOK}{A boolean set to \code{TRUE}.}
#'
#' If an error occurred and \code{trapErrors==TRUE}, the list only contains
#' these components:
#'
#' \item{colors}{A copy of the input colors.}
#'
#' \item{allOK}{a boolean set to \code{FALSE}.}
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
mergeCloseModules <- function(# input data
    exprData,
    colors,

    # Optional starting eigengenes
    MEs = NULL,

    # Optional restriction to a subset of all sets
    useSets = NULL,

    # If missing data are present, impute them?
    impute = TRUE,

    # Input handling options
    checkDataFormat = TRUE,
    unassdColor = if (is.numeric(colors))
        0
    else
        "grey",

    # Options for eigengene network construction
    corFnc = cor,
    corOptions = list(use = 'p'),
    useAbs = FALSE,

    # Options for constructing the consensus
    equalizeQuantiles = FALSE,
    quantileSummary = "mean",
    consensusQuantile = 0,

    # Merging options
    cutHeight = 0.2,
    iterate = TRUE,

    # Output options
    relabel = FALSE,
    colorSeq = NULL,
    getNewMEs = TRUE,
    getNewUnassdME = TRUE,

    # Options controlling behaviour of the function
    trapErrors = FALSE,
    verbose = 1,
    indent = 0) {
    MEsInSingleFrame = FALSE
    spaces = indentSpaces(indent)

    #numCols = is.numeric(colors)
    #facCols = is.factor(colors)
    #charCols = is.character(colors)

    origColors = colors

    colors = colors[, drop = TRUE]

    greyMEname = paste0(getMEprefix(), unassdColor)

    if (verbose > 0) {
        printFlush(
            paste(
                spaces,
                "mergeCloseModules: Merging modules whose distance is less than",
                cutHeight
            )
        )
    }

    if (verbose > 3)
        printFlush(paste(spaces,
                         "  .. will look for grey label", greyMEname))

    if (!checkSets(exprData, checkStructure = TRUE,
                   useSets = useSets)$structureOK)
    {
        if (checkDataFormat)
        {
            exprData = fixDataStructure(exprData)
            MEsInSingleFrame = TRUE
        } else {
            stop("Given exprData appear to be misformatted.")
        }
    }

    setsize = checkSets(exprData, useSets = useSets)
    nSets = setsize$nSets

    if (!is.null(MEs))
    {
        checkMEs = checkSets(MEs, checkStructure = TRUE, useSets = useSets)
        if (checkMEs$structureOK)
        {
            if (nSets != checkMEs$nSets)
                stop("Input error: numbers of sets in exprData and MEs differ.")
            for (set in 1:nSets)
            {
                if (checkMEs$nSamples[set] != setsize$nSamples[set])
                    stop(
                        paste(
                            "Number of samples in MEs is incompatible with subset
                            length for set",
                            set
                        )
                    )
            }
        } else {
            if (MEsInSingleFrame)
            {
                MEs = fixDataStructure(MEs)
                checkMEs = checkSets(MEs)
            } else {
                stop("MEs do not have the appropriate structure (same as exprData).")
            }
        }
    }

    if (setsize$nGenes != length(colors))
        stop(
            "Number of genes in exprData is different from the length of original
            colors. They must equal."
        )

    if ((cutHeight < 0) | (cutHeight > (1 + as.integer(useAbs))))
        stop(paste(
            "Given cutHeight is out of sensible range between 0 and",
            1 + as.integer(useAbs)
        ))

    done = FALSE
    iteration = 1

    MergedColors = colors
    ok = try({
        while (!done)
        {
            if (is.null(MEs))
            {
                MEs = multiSetMEs(
                    exprData,
                    colors = NULL,
                    universalColors = colors,
                    useSets = useSets,
                    impute = impute,
                    subHubs = TRUE,
                    trapErrors = FALSE,
                    excludeGrey = TRUE,
                    grey = unassdColor,
                    verbose = verbose - 1,
                    indent = indent + 1
                )
                MEs = consensusOrderMEs(
                    MEs,
                    useAbs = useAbs,
                    useSets = useSets,
                    greyLast = FALSE
                )
                collectGarbage()
            } else if (nlevels(as.factor(colors)) != checkMEs$nGenes)
            {
                if ((iteration == 1) & (verbose > 0)) {
                    printFlush(
                        paste(
                            spaces,
                            "Number of given module colors",
                            "does not match number of given MEs  = > recalculating the MEs."
                        )
                    )
                }
                MEs = multiSetMEs(
                    exprData,
                    colors = NULL,
                    universalColors = colors,
                    useSets = useSets,
                    impute = impute,
                    subHubs = TRUE,
                    trapErrors = FALSE,
                    excludeGrey = TRUE,
                    grey = unassdColor,
                    verbose = verbose - 1,
                    indent = indent + 1
                )
                MEs = consensusOrderMEs(
                    MEs,
                    useAbs = useAbs,
                    useSets = useSets,
                    greyLast = FALSE
                )
                collectGarbage()
            }
            if (iteration == 1)
                oldMEs = MEs

            # Check colors for number of distinct colors that are not grey

            colLevs = as.character(levels(as.factor(colors)))
            if (length(colLevs[colLevs != as.character(unassdColor)]) <
                2)
            {
                printFlush(paste(
                    spaces,
                    "mergeCloseModules: less than two proper modules."
                ))
                printFlush(paste(
                    spaces,
                    " ..color levels are",
                    paste(colLevs, collapse = ", ")
                ))
                printFlush(paste(spaces, " ..there is nothing to merge."))
                MergedNewColors = colors
                MergedColors = colors
                nOldMods = 1
                nNewMods = 1
                oldTree = NULL
                Tree = NULL
                break
            }

            # Cluster the found module eigengenes and merge ones that are too close
            # according to the specified quantile.

            nOldMods = nlevels(as.factor(colors))

            ConsDiss = .consensusMEDissimilarity(MEs,
                                                 equalizeQuantiles = equalizeQuantiles,
                                                 quantileSummary = quantileSummary,
                                                 consensusQuantile = consensusQuantile,
                                                 useAbs = useAbs,
                                                 corFnc = corFnc,
                                                 corOptions = corOptions,
                                                 useSets = useSets,
                                                 greyMEname = greyMEname)

            Tree = fastcluster::hclust(as.dist(ConsDiss), method = "average")
            if (iteration == 1) {
                oldTree = Tree
            }
            TreeBranches = as.factor(moduleNumber(dendro = Tree,
                                                  cutHeight = cutHeight,
                                                  minSize = 1))
            UniqueBranches = levels(TreeBranches)
            nBranches = nlevels(TreeBranches)
            NumberOnBranch = table(TreeBranches)
            MergedColors = colors

            # Merge modules on the same branch

            for (branch in 1:nBranches)
                if (NumberOnBranch[branch] > 1) {
                    ModulesOnThisBranch = names(TreeBranches)[TreeBranches == UniqueBranches[branch]]
                    ColorsOnThisBranch = substring(ModulesOnThisBranch, 3)
                    if (is.numeric(origColors)) {
                        ColorsOnThisBranch = as.numeric(ColorsOnThisBranch)
                    }
                    if (verbose > 3)
                        printFlush(paste(
                            spaces,
                            "  Merging original colors",
                            paste(ColorsOnThisBranch, collapse = ", ")
                        ))
                    for (color in 2:length(ColorsOnThisBranch)) {
                        MergedColors[MergedColors == ColorsOnThisBranch[color]] =
                            ColorsOnThisBranch[1]
                    }
                }

            MergedColors = MergedColors[, drop = TRUE]

            nNewMods = nlevels(as.factor(MergedColors))

            if (nNewMods < nOldMods & iterate)
            {
                colors = MergedColors
                MEs = NULL
            } else {
                done = TRUE
            }
            iteration = iteration + 1
        }
        if (relabel)
        {
            RawModuleColors = levels(as.factor(MergedColors))
            # relabel the merged colors to the usual order based on the number of
            # genes in each module
            if (is.null(colorSeq))
            {
                if (is.numeric(origColors)) {
                    colorSeq = c(1:length(table(origColors)))
                } else {
                    nNewColors = length(RawModuleColors)
                    colorSeq = labels2colors(c(1:nNewColors))
                }
            }

            # nGenesInModule = rep(0, nNewMods)
            # for (mod in 1:nNewMods) nGenesInModule[mod] =
            # sum(MergedColors == RawModuleColors[mod])
            nGenesInModule = table(MergedColors)

            SortedRawModuleColors = RawModuleColors[order(-nGenesInModule)]

            # Change the color names to the standard sequence, but leave grey grey
            # (that's why rank in general does not equal color)
            MergedNewColors = MergedColors
            if (is.factor(MergedNewColors)) {
                MergedNewColors = as.character(MergedNewColors)
            }
            if (verbose > 3)
                printFlush(paste(spaces, "   Changing original colors:"))
            rank = 0
            for (color in 1:length(SortedRawModuleColors)) {
                if (SortedRawModuleColors[color] != unassdColor) {
                    rank = rank + 1
                    if (verbose > 3) {
                        printFlush(paste(
                            spaces,
                            "      ",
                            SortedRawModuleColors[color],
                            "to ",
                            colorSeq[rank]
                        ))
                    }
                    MergedNewColors[MergedColors == SortedRawModuleColors[color]] =
                        colorSeq[rank]
                }

            }
            if (is.factor(MergedColors))
                MergedNewColors = as.factor(MergedNewColors)
        } else {
            MergedNewColors = MergedColors
        }
        MergedNewColors = MergedNewColors[, drop = TRUE]

        if (getNewMEs)
        {
            if (nNewMods < nOldMods | relabel | getNewUnassdME)
            {
                if (verbose > 0)
                    printFlush(paste(spaces, "  Calculating new MEs..."))
                NewMEs = multiSetMEs(
                    exprData,
                    colors = NULL,
                    universalColors = MergedNewColors,
                    useSets = useSets,
                    impute = impute,
                    subHubs = TRUE,
                    trapErrors = FALSE,
                    excludeGrey = !getNewUnassdME,
                    grey = unassdColor,
                    verbose = verbose - 1,
                    indent = indent + 1
                )
                newMEs = consensusOrderMEs(
                    NewMEs,
                    useAbs = useAbs,
                    useSets = useSets,
                    greyLast = TRUE,
                    greyName = greyMEname
                )

                ConsDiss = .consensusMEDissimilarity(
                    newMEs,
                    equalizeQuantiles = equalizeQuantiles,
                    quantileSummary = quantileSummary,
                    consensusQuantile = consensusQuantile,
                    useAbs = useAbs,
                    corFnc = corFnc,
                    corOptions = corOptions,
                    useSets = useSets,
                    greyMEname = greyMEname
                )
                if (length(ConsDiss) > 1)
                {
                    Tree = fastcluster::hclust(as.dist(ConsDiss), method = "average")
                } else
                    Tree = NULL
            } else {
                newMEs = MEs
            }
        } else {
            newMEs = NULL
        }
        if (MEsInSingleFrame)
        {
            newMEs = newMEs[[1]]$data
            oldMEs = oldMEs[[1]]$data
        }
    }, silent = TRUE)

    if (class(ok) == 'try - error')
    {
        if (!trapErrors)
            stop(ok)
        if (verbose > 0)
        {
            printFlush(
                paste(
                    spaces,
                    "Warning: merging of modules failed with the following error:"
                )
            )
            printFlush(paste('   ', spaces, ok))
            printFlush(paste(
                spaces,
                " - - > returning unmerged modules and  * no *  eigengenes."
            ))
        }
        warning(
            paste(
                "mergeCloseModules: merging of modules failed with the following error:\n",
                "    ",
                ok,
                " - - > return(code)ing unmerged modules and  * no *  eigengenes.\n"
            )
        )
        list(colors = origColors, allOK = FALSE)
    } else {
        list(
            colors = MergedNewColors,
            dendro = Tree,
            oldDendro = oldTree,
            cutHeight = cutHeight,
            oldMEs = oldMEs,
            newMEs = newMEs,
            allOK = TRUE
        )
    }
}

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

# checkAdjMat ####
#' Check adjacency matrix
#'
#' Checks a given matrix for properties that an adjacency matrix must satisfy.
#'
#' The function checks whether the given matrix really is a 2-dimensional
#' numeric matrix, whether it is square, symmetric, and all finite entries are
#' between \code{min} and \code{max}.  If any of the conditions is not met, the
#' function issues an error.
#'
#' @aliases checkAdjMat checkSimilarity
#' @param adjMat matrix to be checked
#' @param similarity matrix to be checked
#' @param min minimum allowed value for entries of the input
#' @param max maximum allowed value for entries of the input
#' @return None. The function returns normally if all conditions are met.
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}}
#' @keywords misc
checkAdjMat <- function(adjMat, min = 0, max = 1) {
    dim = dim(adjMat)
    if (is.null(dim) || length(dim) != 2)
        stop("adjacency is not two - dimensional")
    if (!is.numeric(adjMat))
        stop("adjacency is not numeric")
    if (dim[1] != dim[2])
        stop("adjacency is not square")
    if (max(abs(adjMat - t(adjMat)), na.rm = TRUE) > 1e-12)
        stop("adjacency is not symmetric")
    if (min(adjMat, na.rm = TRUE) < min ||
        max(adjMat, na.rm = TRUE) > max)
        stop("some entries are not between", min, "and", max)
}

# signedKME ####
#' Signed eigengene-based connectivity
#'
#' Calculation of (signed) eigengene-based connectivity, also known as module
#' membership.
#'
#' Signed eigengene-based connectivity of a gene in a module is defined as the
#' correlation of the gene with the corresponding module eigengene.  The
#' samples in \code{datExpr} and \code{datME} must be the same.
#'
#' @param datExpr a data frame containing the gene expression data. Rows
#' correspond to samples and columns to genes. Missing values are allowed and
#' will be ignored.
#' @param datME a data frame containing module eigengenes. Rows correspond to
#' samples and columns to module eigengenes.
#' @param outputColumnName a character string specifying the prefix of column
#' names of the output.
#' @param corFnc character string specifying the function to be used to
#' calculate co-expression similarity. Defaults to Pearson correlation. Any
#' function returning values between -1 and 1 can be used.
#' @param corOptions character string specifying additional arguments to be
#' passed to the function given by \code{corFnc}. Use \code{"use = 'p', method
#' = 'spearman'"} to obtain Spearman correlation.
#' @return A data frame in which rows correspond to input genes and columns to
#' module eigengenes, giving the signed eigengene-based connectivity of each
#' gene with respect to each eigengene.
#' @author Steve Horvath
#' @references
#'
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#' @keywords misc
signedKME <- function(datExpr,
                      datME,
                      outputColumnName = "kME",
                      corFnc = "cor",
                      corOptions = "use = 'p'") {
    datExpr = data.frame(datExpr)
    datME = data.frame(datME)
    output = list()
    if (dim(as.matrix(datME))[[1]]  != dim(as.matrix(datExpr))[[1]]) {
        stop("Number of samples (rows) in 'datExpr' and 'datME' must be the ",
             "same.")
    }
    varianceZeroIndicatordatExpr = as.vector(apply(as.matrix(datExpr), 2, var,
                                                   na.rm = TRUE)) == 0
    varianceZeroIndicatordatME  = as.vector(apply(as.matrix(datME), 2, var,
                                                  na.rm = TRUE)) == 0
    if (sum(varianceZeroIndicatordatExpr, na.rm = TRUE) > 0) {
        warning("Some genes are constant. Hint: consider removing constant ",
                "columns from datExpr.")
    }
    if (sum(varianceZeroIndicatordatME, na.rm = TRUE) > 0) {
        warning("Some module eigengenes are constant, which is suspicious.\n",
                "    Hint: consider removing constant columns from datME.")
    }
    no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)), 2, sum))
    if (min(no.presentdatExpr) < ..minNSamples) {
        warning("Some gene expressions have fewer than 4 observations.\n",
                "    Hint: consider removing genes with too many missing values or
                collect more arrays.")
    }

    #output = data.frame(cor(datExpr, datME, use = "p"))
    corExpr = parse(text = paste("data.frame(", corFnc, "(datExpr, datME ",
                                 prepComma(corOptions), "))"))
    output = eval(corExpr)

    output[no.presentdatExpr < ..minNSamples,] = NA
    names(output) = paste0(outputColumnName, substring(names(datME), first = 3,
                                                       last = 100))
    dimnames(output)[[1]] = names(datExpr)
    output
} # end of function signedKME

# clusterCoef ####
#' Clustering coefficient calculation
#'
#' This function calculates the clustering coefficients for all nodes in the
#' network given by the input adjacency matrix.
#'
#'
#' @param adjMat adjacency matrix
#' @return A vector of clustering coefficients for each node.
#' @author Steve Horvath
#' @keywords misc
clusterCoef <- function(adjMat) {
    checkAdjMat(adjMat)
    diag(adjMat) = 0
    nNodes = dim(adjMat)[[1]]
    computeLinksInNeighbors <- function(x, imatrix) {
        x %*% imatrix %*% x
    }
    nolinksNeighbors <- c(rep(-666, nNodes))
    total.edge <- c(rep(-666, nNodes))
    maxh1 = max(as.dist(adjMat))
    minh1 = min(as.dist(adjMat))
    if (maxh1 > 1 | minh1 < 0) {
        stop("The adjacency matrix contains entries that are larger than 1 or
                smaller than 0: max  = ", maxh1, ", min  = ", minh1)
    }
    nolinksNeighbors <- apply(adjMat, 1, computeLinksInNeighbors,
                              imatrix = adjMat)
    plainsum  <- apply(adjMat, 1, sum)
    squaresum <- apply(adjMat ^ 2, 1, sum)
    total.edge = plainsum ^ 2 - squaresum
    CChelp = rep(-666, nNodes)
    CChelp = ifelse(total.edge == 0, 0, nolinksNeighbors / total.edge)
    CChelp
} # end of function

# addErrorBars ####
#' Add error bars to a barplot.
#'
#' This function adds error bars to an existing barplot.
#'
#'
#' @param means vector of means plotted in the barplot
#' @param errors vector of standard errors (signle positive values) to be
#' plotted.
#' @param two.side should the error bars be two-sided?
#' @return None.
#' @author Steve Horvath and Peter Langfelder
#' @keywords hplot
addErrorBars <- function(means, errors, two.side = FALSE) {
    if (!is.numeric(means)) {
        stop("All arguments must be numeric")
    }

    if (is.null(dim(means)) || length(dim(means)) == 1) {
        xval <- (cumsum(c(0.7, rep(
            1.2, length(means) - 1
        ))))
    } else{
        if (length(dim(means)) == 2) {
            xval <- cumsum(array(c(1, rep(
                0, dim(means)[1] - 1
            )),
            dim = c(1, length(means)))) + 0:(length(means) - 1) + .5
        } else{
            stop("First argument must either be a vector or a matrix")
        }
    }
    MW <- 0.25 * (max(xval) / length(xval))
    ERR1 <- means + errors
    ERR2 <- means - errors
    for (i in 1:length(means)) {
        segments(xval[i], means[i], xval[i], ERR1[i])
        segments(xval[i] - MW, ERR1[i], xval[i] + MW, ERR1[i])
        if (two.side) {
            segments(xval[i], means[i], xval[i], ERR2[i])
            segments(xval[i] - MW, ERR2[i], xval[i] + MW, ERR2[i])
        }
    }
}

# stdErr ####
#' Standard error of the mean of a given vector.
#'
#' Returns the standard error of the mean of a given vector. Missing values are
#' ignored.
#'
#'
#' @param x a numeric vector
#' @return Standard error of the mean of x.
#' @author Steve Horvath
#' @keywords misc
stdErr <- function(x) {
    sqrt(var(x, na.rm = TRUE) / sum(!is.na(x)))
}

#===============================================================================
# The following two functions are for displaying the pair - wise correlation in
# a panel when using the command "pairs()"
# Typically, we use "pairs(DATA, upper.panel = panel.smooth,
# lower.panel = .panel.cor, diag.panel = panel.hist)" to
# put the correlation coefficients on the lower panel.

.panel.hist <- function(x, ...) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(usr[1:2], 0, 1.5))
    h <- hist(x, plot = FALSE)
    breaks <- h$breaks
    nB <- length(breaks)
    y <- h$counts
    y <- y / max(y)
    rect(breaks[-nB], 0, breaks[-1], y, col = "cyan", ...)
}

#===============================================================================
# This function is used in "pairs()" function. The problem of the original
# panel.cor is that when the correlation coefficient is very small, the lower
# panel will have a large font instead of a mini-font in a saved .ps file.
# This new function uses a format for corr = 0.2 when corr<0.2, but it still
# reports the original value of corr, with a minimum format.

.panel.cor <- function(x,
                       y,
                       digits = 2,
                       prefix = "",
                       cex.cor) {
    usr <- par("usr")
    on.exit(par(usr))
    par(usr = c(0, 1, 0, 1))
    r <- abs(cor(x, y))
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    txt1 = txt
    r1 = r
    if (r < 0.2) {
        r1 = 0.2
        txt1 <- format(c(r1, 0.123456789), digits = digits)[1]
        txt1 <- paste0(prefix, txt1)
    }
    if (missing(cex.cor))
        cex <- 0.8 / strwidth(txt1)
    cex = cex * r1
    r <- round(r, digits)
    txt <- format(c(r, 0.123456789), digits = digits)[1]
    txt <- paste0(prefix, txt)
    text(0.5, 0.5, txt, cex = cex)
}

# dynamicMergeCut ####
#' Threshold for module merging
#'
#' Calculate a suitable threshold for module merging based on the number of
#' samples and a desired Z quantile.
#'
#' This function calculates the threshold for module merging. The threshold is
#' calculated as the lower boundary of the interval around the theoretical
#' correlation \code{mergeCor} whose width is given by the Z value
#' \code{Zquantile}.
#'
#' @param n number of samples
#' @param mergeCor theoretical correlation threshold for module merging
#' @param Zquantile Z quantile for module merging
#' @return The correlation threshold for module merging; a single number.
#' @author Steve Horvath
#' @seealso \code{\link{moduleEigengenes}}, \code{\link{mergeCloseModules}}
#' @keywords misc
#' @examples
#' dynamicMergeCut(20)
#' dynamicMergeCut(50)
#' dynamicMergeCut(100)
dynamicMergeCut <- function(n, mergeCor = .9, Zquantile = 2.35) {
    if (mergeCor > 1 | mergeCor < 0) {
        stop("'mergeCor' must be between 0 and 1.")
    }
    if (mergeCor == 1) {
        printFlush("dynamicMergeCut: given mergeCor = 1 will be set to .999.")
        mergeCor = .999
    }
    if (n < 4) {
        printFlush(
            paste(
                "Warning in function dynamicMergeCut: too few observations
                for the dynamic",
                "assignment of the merge threshold.\n
                Will set the threshold to .35"
            )
        )
        mergethreshold = .35
    } else {
        # Fisher transform of the true merge correlation
        FishermergeCor = .5 * log((1 + mergeCor) / (1 - mergeCor))
        E = exp(2 * (FishermergeCor  - Zquantile / sqrt(n - 3)))
        LowerBoundCIcor = (E - 1) / (E + 1)
        mergethreshold = 1 -  LowerBoundCIcor
    }
    if (mergethreshold > 1) {
        1
    }
    else {
        mergethreshold
    }
}# end of function dynamicMergeCut

# propVarExplained ####
#' Proportion of variance explained by eigengenes.
#'
#' This function calculates the proportion of variance of genes in each module
#' explained by the respective module eigengene.
#'
#' For compatibility with other functions, entries in \code{color} are matched
#' to a substring of \code{names(MEs)} starting at position 3. For example, the
#' entry \code{"turquoise"} in \code{colors} will be matched to the eigengene
#' named \code{"MEturquoise"}. The first two characters of the eigengene name
#' are ignored and can be arbitrary.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. NAs are allowed and will be ignored.
#' @param colors a vector giving module assignment for genes given in
#' \code{datExpr}. Unique values should correspond to the names of the
#' eigengenes in \code{MEs}.
#' @param MEs a data frame of module eigengenes in which each column is an
#' eigengene and each row corresponds to a sample.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @return A vector with one entry per eigengene containing the proportion of
#' variance of the module explained by the eigengene.
#' @author Peter Langfelder
#' @seealso \code{\link{moduleEigengenes}}
#' @keywords misc
propVarExplained <- function(datExpr, colors, MEs, corFnc = "cor",
                             corOptions = "use = 'p'") {
    fc = as.factor(colors)
    mods = levels(fc)
    nMods = nlevels(fc)
    nGenes = ncol(datExpr)
    if (nMods != ncol(MEs)) {
        stop("Input error: number of distinct 'colors' differs from\n",
                   "     the number of module eigengenes given in ME.")
    }

    if (ncol(datExpr) != length(colors)) {
        stop("Input error: number of probes (columns) in 'datExpr' differs ",
             "from the length of goven 'colors'.")
    }

    if (nrow(datExpr) != nrow(MEs))
        stop("Input error: number of observations (rows) in 'datExpr' and 'MEs'
             differ.")

    PVE = rep(0, nMods)

    col2MEs = match(mods, substring(names(MEs), 3))

    if (sum(is.na(col2MEs)) > 0) {
        stop("Input error: not all given colors could be matched to names of",
             " module eigengenes.")
    }
    for (mod in 1:nMods) {
        modGenes = c(1:nGenes)[as.character(colors) == mods[mod]]
        corExpr = parse(text = paste(
            corFnc,
            "(datExpr[, modGenes], MEs[, col2MEs[mod]]",
            prepComma(corOptions),
            ")"
        ))
        PVE[mod] = mean(as.vector(eval(corExpr) ^ 2))
    }

    names(PVE) = paste0("PVE", mods)
    PVE
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
#'   addGrid();
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

# alignExpr ####
# If y is supplied, it multiplies columns of datExpr by +/- 1 to make all
# correlations with y positive.
# If y is not supplied, the first column of datExpr is used as the reference
# direction.
#' Align expression data with given vector
#'
#' Multiplies genes (columns) in given expression data such that their
#' correlation with given reference vector is non-negative.
#'
#' The function basically multiplies each column in \code{datExpr} by the sign
#' of its correlation with \code{y}. If \code{y} is not given, the first column
#' in \code{datExpr} will be used as the reference vector.
#'
#' @param datExpr expression data to be aligned. A data frame with columns
#' corresponding to genes and rows to samples.
#' @param y reference vector of length equal the number of samples (rows) in
#' \code{datExpr}
#' @return A data frame containing the aligned expression data, of the same
#' dimensions as the input data frame.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
alignExpr <- function(datExpr, y = NULL) {
    if (!is.null(y) & dim(as.matrix(datExpr))[[1]]  != length(y))
        stop("Incompatible number of samples in 'datExpr' and 'y'.")
    if (is.null(y))
        y <- as.numeric(datExpr[, 1])
    sign1 <- sign(as.numeric(cor(y, datExpr, use = "p")))
    as.data.frame(scale(t(t(datExpr) * sign1)))
} # end of function alignExpr


# this function can be used to rank the values in x. Ties are broken by the
# method first.
# This function does not appear to be used anywhere in these functions.
# rank1 <- function(x){
#    rank(x, ties.method = "first")
#}

# hubGeneSignificance ####
# Input a data frame with possibly signed module membership measures (also known
# as module eigengene based connectivity kME.
# The input to this function can include the sign of the correlation.
#' Hubgene significance
#'
#' Calculate approximate hub gene significance for all modules in network.
#'
#' In \code{datKME} rows correspond to genes and columns to modules.
#'
#' @param datKME a data frame (or a matrix-like object) containing
#' eigengene-based connectivities of all genes in the network.
#' @param GS a vector with one entry for every gene containing its gene
#' significance.
#' @return A vector whose entries are the hub gene significances for each
#' module.
#' @author Steve Horvath
#' @references Dong J, Horvath S (2007) Understanding Network Concepts in
#' Modules, BMC Systems Biology 2007, 1:24
#' @keywords misc
hubGeneSignificance <- function(datKME, GS) {
    nMEs = dim(as.matrix(datKME))[[2]]
    nGenes = dim(as.matrix(datKME))[[1]]
    if (length(GS)  != nGenes)
        stop("Numbers of genes in 'datKME' and 'GS' are not compatible.")
    Kmax = as.numeric(apply(as.matrix(abs(datKME)), 2, max, na.rm = TRUE))
    Kmax[Kmax == 0] = 1
    datKME = scale(datKME, center = FALSE, scale = Kmax)
    sumKsq = as.numeric(apply(as.matrix(datKME ^ 2), 2, sum, na.rm = TRUE))
    sumKsq[sumKsq == 0] = 1
    HGS = as.numeric(apply(I(GS) * datKME, 2, sum, na.rm = TRUE)) / sumKsq
    as.numeric(HGS)
} #end of function hubGeneSignificance

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

# addTraitToPCs ####
#' Add trait information to multi-set module eigengene structure
#'
#' Adds trait information to multi-set module eigengene structure.
#'
#' The function simply \code{cbind}'s the module eigengenes and traits for each
#' set. The number of sets and numbers of samples in each set must be
#' consistent between \code{multiMEs} and \code{multiTraits}.
#'
#' @param multiME Module eigengenes in multi-set format. A vector of lists, one
#' list per set. Each list must contain an element named \code{data} that is a
#' data frame with module eigengenes.
#' @param multiTraits Microarray sample trait(s) in multi-set format. A vector
#' of lists, one list per set. Each list must contain an element named
#' \code{data} that is a data frame in which each column corresponds to a
#' trait, and each row to an individual sample.
#' @return A multi-set structure analogous to the input: a vector of lists, one
#' list per set. Each list will contain a component \code{data} with the merged
#' eigengenes and traits for the corresponding set.
#' @author Peter Langfelder
#' @seealso \code{\link{checkSets}}, \code{\link{moduleEigengenes}}
#' @keywords misc
addTraitToMEs <- function(multiME, multiTraits) {
    nSets = length(multiTraits)
    setsize = checkSets(multiTraits)
    nTraits = setsize$nGenes
    nSamples = setsize$nSamples

    if (length(multiME) != nSets) {
        stop("Numbers of sets in multiME and multiTraits parameters differ -
             must be the same.")
    }

    multiMETs = vector(mode = "list", length = nSets)
    for (set in 1:nSets) {
        trait.subs = multiTraits[[set]]$data
        multiMET = as.data.frame(cbind(multiME[[set]]$data, trait.subs))
        colnames(multiMET) = c(colnames(multiME[[set]]$data),
                               colnames(trait.subs))
        if (!is.null(multiME[[set]]$AET)) {
            AET = as.data.frame(cbind(multiME[[set]]$averageExpr, trait.subs))
            colnames(AET) = c(colnames(multiME[[set]]$averageExpr),
                              colnames(trait.subs))
        }
        multiMETs[[set]] = list(data = multiMET)
    }
    multiMETs
}

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
#' heatmaps and barplots;
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
    par(mfrow = c(nPlotRows, nPlotCols))
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

# numbers2colors ####
#' Color representation for a numeric variable
#'
#' The function creates a color represenation for the given numeric input.
#'
#' Each column of \code{x} is processed individually, meaning that the color
#' palette is adjusted individually for each column of \code{x}.
#'
#' @param x a vector or matrix of numbers. Missing values are allowed and will
#' be assigned the color given in \code{naColor}. If a matrix, each column of
#' the matrix is processed separately and the return value will be a matrix of
#' colors.
#' @param signed logical: should \code{x} be considered signed? If \code{TRUE},
#' the default setting is to use to use a palette that starts with green for
#' the most negative values, continues with white for values around zero and
#' turns red for positive values. If \code{FALSE}, the default palette ranges
#' from white for minimum values to red for maximum values. If not given, the
#' behaviour is controlled by values in \code{x}: if there are both positive
#' and negative values, \code{signed} will be considered \code{TRUE}, otherwise
#' \code{FALSE}.
#' @param centered logical. If \code{TRUE} and \code{signed==TRUE}, numeric
#' value zero will correspond to the middle of the color palette. If
#' \code{FALSE} or \code{signed==FALSE}, the middle of the color palette will
#' correspond to the average of the minimum and maximum value. If neither
#' \code{signed} nor \code{centered} are given, \code{centered} will follow
#' \code{signed} (see above).
#' @param lim optional specification of limits, that is numeric values that
#' should correspond to the first and last entry of \code{colors}.
#' @param commonLim logical: should limits be calculated separately for each
#' column of x, or should the limits be the same for all columns? Only applies
#' if \code{lim} is \code{NULL}.
#' @param colors color palette to represent the given numbers.
#' @param naColor color to represent missing values in \code{x}.
#' @return A vector or matrix (of the same dimensions as \code{x}) of colors.
#' @author Peter Langfelder
#' @seealso \code{\link{labels2colors}} for color coding of ordinal labels.
#' @keywords misc
numbers2colors <- function(x,
                           signed = NULL,
                           centered = signed,
                           lim = NULL,
                           commonLim = FALSE,
                           colors = if (signed) {
                               blueWhiteRed(100)
                           } else {
                               blueWhiteRed(100)[51:100]
                           },
                           naColor = "grey") {
    x = as.matrix(x)
    if (!is.numeric(x))
        stop("'x' must be numeric. For a factor, please use as.numeric(x) in
             the call.")
    if (is.null(signed)) {
        if (any(x < 0, na.rm = TRUE) & any(x > 0, na.rm = TRUE))
        {
            signed = TRUE
        } else
            signed = FALSE
    }
    if (is.null(centered))
        centered = signed

    if (is.null(lim)) {
        if (signed & centered) {
            max = apply(abs(x), 2, max, na.rm = TRUE)
            lim = as.matrix(cbind(-max, max))
        } else {
            lim = as.matrix(cbind(
                apply(x, 2, min, na.rm = TRUE),
                apply(x, 2, max, na.rm = TRUE)
            ))
        }
        if (commonLim)
            lim = c(min(lim[, 1], na.rm = TRUE), max(lim[, 2], na.rm = TRUE))
    }
    if (is.null(dim(lim))) {
        if (length(lim) != 2)
            stop("'lim' must be a vector of length 2 or a matrix with 2 columns.")
        if (!is.numeric(lim))
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != 2)
            stop("'lim' must be finite.")
        lim = t(as.matrix(lim))
    } else {
        if (ncol(x) != nrow(lim))
            stop("Incompatible numbers of columns in 'x' and rows in 'lim'.")
        if (!is.numeric(lim))
            stop("'lim' must be numeric")
        if (sum(is.finite(lim)) != length(lim))
            stop("'lim' must be finite.")
    }

    xMin = matrix(lim[, 1],
                  nrow = nrow(x),
                  ncol = ncol(x),
                  byrow = TRUE)
    xMax = matrix(lim[, 2],
                  nrow = nrow(x),
                  ncol = ncol(x),
                  byrow = TRUE)

    if (sum(xMin == xMax) > 0)
        warning("(some columns in) 'x' are constant. Their color will be the
                color of NA.")

    xx = x
    xx[is.na(xx)] = ((xMin + xMax)[is.na(xx)]) / 2
    if (sum(x < xMin, na.rm = TRUE) > 0)
    {
        warning("Some values of 'x' are below given minimum and will be
                truncated to the minimum.")
        x[xx < xMin] = xMin[xx < xMin]
    }

    if (sum(x > xMax, na.rm = TRUE) > 0)
    {
        warning("Some values of 'x' are above given maximum and will be
                truncated to the maximum.")
        x[xx > xMax] = xMax[xx > xMax]
    }

    mmEq = xMin == xMax

    nColors = length(colors)

    xCol = array(naColor, dim = dim(x))

    xInd = (x - xMin) / (xMax - xMin)
    xInd[xInd == 1] = 1 - 0.5 / nColors
    xCol[!mmEq] = colors[as.integer(xInd[!mmEq] * nColors) + 1]
    xCol[is.na(xCol)] = naColor

    xCol
}

# randIndex ####
#  Rand index calculation
# this function is used for computing the Rand index below...
#
.choosenew <- function(n, k) {
    n <- c(n)
    out1 <- rep(0, length(n))
    for (i in c(1:length(n))) {
        if (n[i] < k) {
            out1[i] <- 0
        }
        else {
            out1[i] <- choose(n[i], k)
        }
    }
    out1
}
#' Rand index of two partitions
#'
#' Computes the Rand index, a measure of the similarity between two
#' clusterings.
#'
#'
#' @param tab a matrix giving the cross-tabulation table of two clusterings.
#' @param adjust logical: should the "adjusted" version be computed?
#' @return the Rand index of the input table.
#' @author Steve Horvath
#' @references W. M. Rand (1971). "Objective criteria for the evaluation of
#' clustering methods". Journal of the American Statistical Association 66:
#' 846-850
#' @keywords misc
randIndex <- function(tab, adjust = TRUE) {
    a <- 0
    b <- 0
    c <- 0
    d <- 0
    nn <- 0
    m <- nrow(tab)
    n <- ncol(tab)
    for (i in 1:m) {
        c <- 0
        for (j in 1:n) {
            a <- a + .choosenew(tab[i, j], 2)
            nj <- sum(tab[, j])
            c <- c + .choosenew(nj, 2)
        }
        ni <- sum(tab[i,])
        b <- b + .choosenew(ni, 2)
        nn <- nn + ni
    }
    if (adjust) {
        d <- .choosenew(nn, 2)
        adrand <- (a - (b * c) / d) / (0.5 * (b + c) - (b * c) / d)
        adrand
    } else {
        b <- b - a
        c <- c - a
        d <- .choosenew(nn, 2) - a - b - c
        rand <- (a + d) / (a + b + c + d)
        rand
    }
}

# vectorizeMatrix ####
#' Turn a matrix into a vector of non-redundant components
#'
#' A convenient function to turn a matrix into a vector of non-redundant
#' components. If the matrix is non-symmetric, returns a vector containing all
#' entries of the matrix. If the matrix is symmetric, only returns the upper
#' triangle and optionally the diagonal.
#'
#'
#' @param M the matrix or data frame to be vectorized.
#' @param diag logical: should the diagonal be included in the output?
#' @return A vector containing the non-redundant entries of the input matrix.
#' @author Steve Horvath
#' @keywords misc
vectorizeMatrix <- function(M, diag = FALSE) {
    if (is.null(dim(M)))
        stop("The input of the vectorize function is not a matrix or data frame.")
    if (length(dim(M)) != 2)
        stop("The input of the vectorize function is not a matrix or data frame.")
    # now we check whether the matrix is symmetrical
    if (dim(M)[[1]] == dim(M)[[2]]) {
        M = as.matrix(M)
        Mtranspose = t(M)
        abs.difference = max(abs(M - Mtranspose), na.rm = TRUE)
        if (abs.difference < 10 ^ (-14)) {
            out = M[upper.tri(M, diag)]
        }
        else
            out = as.vector(M)
    } else
        out = as.vector(M)
    out
} # end

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

# metaZfunction ####
#' Meta-analysis Z statistic
#'
#' The function calculates a meta analysis Z statistic based on an input data
#' frame of Z statistics.
#'
#' For example, if datZ has 3 columns whose columns are labelled Z1,Z2,Z3 then
#' ZMeta= (Z1+Z2+Z3)/sqrt(3). Under the null hypothesis (where all Z statistics
#' follow a standard normal distribution and the Z statistics are independent),
#' ZMeta also follows a standard normal distribution.  To calculate a 2 sided
#' p-value, one an use the following code pvalue=2*pnorm(-abs(ZMeta) )
#'
#' @param datZ Matrix or data frame of Z statistics (assuming standard normal
#' distribution under the null hypothesis). Rows correspond to genes, columns
#' to independent data sets.
#' @param columnweights optional vector of non-negative numbers for weighing
#' the columns of datZ.
#' @return Vector of meta analysis Z statistic. Under the null hypothesis this
#' should follow a standard normal distribution.
#' @author Steve Horvath
#' @keywords misc
metaZfunction <- function(datZ, columnweights = NULL) {
    if (!is.null(columnweights))  {
        datZ = t(t(datZ) *  columnweights)
    }
    datZpresent = !is.na(datZ) + 0.0
    if (!is.null(columnweights))  {
        datZpresent = t(t(datZpresent) *  columnweights)
    }
    sumZ = as.numeric(rowSums(datZ, na.rm = TRUE))
    variance = as.numeric(rowSums(datZpresent ^ 2))
    sumZ / sqrt(variance)
}

# prepComma ####
#' Prepend a comma to a non-empty string
#'
#' Utility function that prepends a comma before the input string if the string
#' is non-empty.
#'
#'
#' @param s Character string.
#' @return If \code{s} is non-empty, returns \code{paste(",", s)}, otherwise
#' returns s.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' prepComma("abc")
#' prepComma("")
#'
prepComma <- function(s) {
    if (s == "")
        return (s)
    paste(", ", s)
}

# multiUnion ####
#' Union and intersection of multiple sets
#'
#' Union and intersection of multiple sets. These function generalize the
#' standard functions \code{\link{union}} and \code{\link{intersect}}.
#'
#'
#' @aliases multiUnion multiIntersect
#' @param setList A list containing the sets to be performed upon.
#' @return The union or intersection of the given sets.
#' @author Peter Langfelder
#' @seealso
#' The "standard" functions \code{\link{union}} and \code{\link{intersect}}.
#' @keywords misc
multiUnion <- function(setList) {
    len = length(setList)
    if (len == 0) {
        return(NULL)
    }
    if (len == 1) {
        return(setList[[1]])
    }

    out = setList[[1]]
    for (elem in 2:len) {
        out <- union(out, setList[[elem]])
    }
    out
}

#' @export
#' @rdname multiUnion
multiIntersect <- function(setList) {
    len = length(setList)
    if (len == 0) {
        return(NULL)
    }
    if (len == 1) {
        return(setList[[1]])
    }

    out = setList[[1]]
    for (elem in 2:len) {
        out <- intersect(out, setList[[elem]])
    }
    out
}

# prependZeros ####
#' Pad numbers with leading zeros to specified total width
#'
#' This function pads the specified numbers with zeros to a specified total
#' width.
#'
#'
#' @param x Vector of numbers to be padded.
#' @param width Width to pad the numbers to.
#' @return Character vector with the 0-padded numbers.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' prependZeros(1:10)
#' prependZeros(1:10, 4)
#'
prependZeros <- function(x, width = max(nchar(x))) {
    lengths = nchar(x)
    if (width < max(lengths))
        stop("Some entries of 'x' are too long.")
    out = as.character(x)
    n = length(x)
    for (i in 1:n) {
        if (lengths[i] < width) {
            out[i] <- paste0(paste(rep("0", width - lengths[i]), collapse = ""),
                            x[i])
        }
    }
    out
}

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
