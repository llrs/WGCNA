# Categories of functions:
# . network construction (including connectivity calculation)
# . module detection
# . gene screening
# . data simulation
# . general statistical functions
# . visualization



#-------------------------------------------------------------------------------
#
# Overall options and settings for the package
#
#-------------------------------------------------------------------------------

.moduleColorOptions = list(MEprefix = "ME")


# moduleColor.getMEprefix ####
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
moduleColor.getMEprefix <- function() {
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
            printFlush(
                paste(
                    spaces,
                    "moduleEigengenes: Calculating",
                    nlevels(as.factor(colors)),
                    "module eigengenes in given set."
                )
            )
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
            stop(
                "moduleEigengenes: Error: ncol(datExpr) and length(colors)",
                "must be equal (one color per gene)."
            )
        }

        if (is.factor(colors)) {
            nl = nlevels(colors)
            nlDrop = nlevels(colors[, drop = TRUE])
            if (nl > nlDrop)
                stop(
                    paste(
                        "Argument 'colors' contains unused levels ",
                        "(empty modules).",
                        "Use colors[, drop = TRUE] to get rid of them."
                    )
                )
        }

        if (softPower < 0) {
            stop("softPower must be non - negative")
        }

        alignRecognizedValues = c("", "along average")
        if (!is.element(align, alignRecognizedValues)) {
            stop(
                paste(
                    "ModulePrincipalComponents: Error:",
                    "parameter align has an unrecognised value:",
                    align,
                    "; Recognized values are ",
                    alignRecognizedValues
                )
            )
        }

        maxVarExplained = 10
        if (nPC > maxVarExplained) {
            warning(paste("Given nPC is too large. Will use value",
                          maxVarExplained))
        }

        nVarExplained <- min(nPC, maxVarExplained)
        modlevels <- levels(factor(colors))
        if (excludeGrey) {
            if (sum(as.character(modlevels) != as.character(grey)) > 0) {
                modlevels <-
                    modlevels[as.character(modlevels) != as.character(grey)]
            } else {
                stop(
                    "Color levels are empty. Possible reason: ",
                    "the only color is grey and grey module is excluded ",
                    "from the calculation."
                )
            }
        }
        lml <- length(modlevels)
        PrinComps <-
            data.frame(matrix(nrow = dim(datExpr)[[1]], ncol = lml))
        averExpr <-
            data.frame(matrix(nrow = dim(datExpr)[[1]], ncol = lml))
        varExpl <- data.frame(matrix(nrow = nVarExplained, ncol = lml))
        validMEs <- rep(TRUE, lml)
        validAEs <- rep(FALSE, lml)
        isPC <- validMEs
        isHub <- validAEs
        validColors <- colors
        names(PrinComps) <-
            paste0(moduleColor.getMEprefix(), modlevels)
        names(averExpr) <- paste0("AE", modlevels)
        for (i in c(1:lml)) {
            if (verbose > 1) {
                printFlush(paste(
                    spaces,
                    "moduleEigengenes : Working on ME for module",
                    modlevels[i]
                ))
            }
            modulename <- modlevels[i]
            restrict1 <-
                as.character(colors) == as.character(modulename)
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
                        }, silent <- TRUE)
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
                svd1 <-
                    svd(datModule,
                        nu = min(n, p, nPC),
                        nv = min(n, p, nPC))
                # varExpl[, i] = (svd1$d[1:min(n, p,
                # nVarExplained)])^2/sum(svd1$d^2)
                if (verbose > 5) {
                    printFlush(paste(spaces, " ...calculating PVE"))
                }
                veMat <-
                    cor(svd1$v[, c(1:min(n, p, nVarExplained))],
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
                            paste(
                                spaces,
                                " ..principal component calculation for module",
                                modulename,
                                "failed with the following error:"
                            )
                        )
                        printFlush(
                            paste(
                                spaces,
                                "     ",
                                pc,
                                spaces,
                                "..hub genes will be used instead of ",
                                "principal components."
                            )
                        )
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
                    printFlush(
                        paste(
                            spaces,
                            " ..ME calculation of module",
                            modulename,
                            "failed with the following error:"
                        )
                    )
                    printFlush(
                        paste(
                            spaces,
                            "     ",
                            pc,
                            spaces,
                            " ..the offending module has been removed."
                        )
                    )
                }
                warning(
                    paste(
                        "Eigengene calculation of module",
                        modulename,
                        "failed with the following error \n     ",
                        pc,
                        "The offending module has been removed.\n"
                    )
                )
                validMEs[i] <- FALSE
                isPC[i] <- FALSE
                isHub[i] <- FALSE
                validColors[restrict1] <- grey
            } else {
                PrinComps[, i] <- pc
                ae <- try({
                    if (isPC[i]) {
                        scaledExpr <- scale(t(datModule))
                    }
                    averExpr[, i] <-
                        rowMeans(scaledExpr, na.rm = TRUE)
                    if (align == "along average") {
                        if (verbose > 4) {
                            printFlush(
                                paste(
                                    spaces,
                                    " .. aligning module eigengene",
                                    "with average expression."
                                )
                            )
                        }
                        corAve <-
                            cor(averExpr[, i], PrinComps[, i], use = "p")
                        if (!is.finite(corAve)) {
                            corAve <- 0
                        }
                        if (corAve < 0) {
                            PrinComps[, i] <- -PrinComps[, i]
                        }
                    }
                    0
                }, silent = TRUE)
                if (class(ae) == 'try - error')
                {
                    if (!trapErrors) {
                        stop(ae)
                    }
                    if (verbose > 0) {
                        printFlush(
                            paste(
                                spaces,
                                " ..Average expression calculation of module",
                                modulename,
                                "failed with the following error:"
                            )
                        )
                        printFlush(
                            paste(
                                spaces,
                                "     ",
                                ae,
                                spaces,
                                " ..the returned average expression vector will be invalid."
                            )
                        )
                    }
                    warning(
                        paste(
                            "Average expression calculation of module",
                            modulename,
                            "failed with the following error \n     ",
                            ae,
                            "The returned average expression vector will be invalid.\n"
                        )
                    )
                }
                validAEs[i] <- !(class(ae) == 'try - error')
            }
        }
        allOK <- (sum(!validMEs) == 0)
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
        list(
            eigengenes = PrinComps,
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
            allAEOK = allAEOK
        )
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
removeGreyME <- function(MEs, greyMEName = paste0(moduleColor.getMEprefix(),
                                             "grey")) {
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
#' @seealso \code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}},
#' \code{\link{consensusOrderMEs}}
#' @keywords misc
orderMEs <- function(MEs,
                     greyLast = TRUE,
                     greyName = paste0(moduleColor.getMEprefix(), "grey"),
                     orderBy = 1,
                     order = NULL,
                     useSets = NULL,
                     verbose = 0,
                     indent = 0) {
    spaces <- indentSpaces(indent)

    if ("eigengenes" %in% names(MEs)) {
        if (is.null(order)) {
            if (verbose > 0) {
                printFlush(
                    paste(
                        spaces,
                        "orderMEs: order not given, calculating ",
                        "using given set",
                        orderBy
                    )
                )
            }
            corPC <- cor(MEs$eigengenes, use = "p")
            disPC <- 1 - corPC
            order <-
                .clustOrder(disPC, greyLast = greyLast, greyName = greyName)
        }

        if (length(order) != dim(MEs$eigengenes)[2]) {
            stop("orderMEs: given MEs and order have incompatible dimensions.")
        }
        orderedMEs <- MEs
        orderedMEs$eigengenes <-
            as.data.frame(MEs$eigengenes[, order])
        colnames(orderedMEs$eigengenes) <-
            colnames(MEs$eigengenes)[order]
        if (!is.null(MEs$averageExpr)) {
            orderedMEs$averageExpr <- as.data.frame(MEs$averageExpr[, order])
            colnames(orderedMEs$averageExpr) <-
                colnames(MEs$data)[order]
        }
        if (!is.null(MEs$varExplained)) {
            orderedMEs$varExplained <- as.data.frame(MEs$varExplained[, order])
            colnames(orderedMEs$varExplained) <-
                colnames(MEs$data)[order]
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
                    printFlush(
                        paste(
                            spaces,
                            "orderMEs: order not given, calculating using given set",
                            orderBy
                        )
                    )
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

.clustOrder <- function(distM,
                        greyLast = TRUE,
                        greyName = paste0(moduleColor.getMEprefix(), "grey")) {
    distM = as.matrix(distM)
    distNames = dimnames(distM)[[1]]
    greyInd = match(greyName, distNames)
    if (greyLast && !is.na(greyInd)) {
        clusterMEs = (greyName != distNames)
        if (sum(clusterMEs) > 1) {
            h = fastcluster::hclust(as.dist(distM[clusterMEs, clusterMEs]),
                                    method = "average")
            order = h$order
            if (sum(order >= greyInd) > 0) {
                order[order >= greyInd] = order[order >= greyInd] + 1
            }
            order = c(order, greyInd)
        } else if (ncol(distM) > 1) {
            if (greyInd == 1) {
                order = c(2, 1)
            } else {
                order = c(1, 2)
            }
        } else {
            order = 1
        }
    } else {
        if (length(distM) > 1) {
            h = fastcluster::hclust(as.dist(distM), method = "average")
            order = h$order
        } else {
            order = 1
        }
    }
    order

    # print(paste("names:", names(distM), collapse = ", "))
    # print(paste("order:", order, collapse = ", "))
}

# consensusOrderMEs ####
#' Put close eigenvectors next to each other in several sets.
#'
#' Reorder given (eigen-)vectors such that similar ones (as measured by
#' correlation) are next to each other. This is a multi-set version of
#' \code{\link{orderMEs}}; the dissimilarity used can be of consensus type (for
#' each pair of eigenvectors the consensus dissimilarity is the maximum of
#' individual set dissimilarities over all sets) or of majority type (for each
#' pair of eigenvectors the consensus dissimilarity is the average of
#' individual set dissimilarities over all sets).
#'
#' Ordering module eigengenes is useful for plotting purposes. This function
#' calculates the consensus or majority dissimilarity of given eigengenes over
#' the sets specified by \code{useSets} (defaults to all sets). A hierarchical
#' dendrogram is calculated using the dissimilarity and the order given by the
#' dendrogram is used for the eigengenes in all other sets.
#'
#' @param MEs Module eigengenes of several sets in a multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, with each list corresponding to
#' one dataset and the module eigengenes in the component \code{data}, that is
#' \code{MEs[[set]]$data[sample, module]} is the expression of the eigengene of
#' module \code{module} in sample \code{sample} in dataset \code{set}. The
#' number of samples can be different between the sets, but the modules must be
#' the same.
#' @param useAbs Controls whether vector similarity should be given by absolute
#' value of correlation or plain correlation.
#' @param useSets Allows the user to specify for which sets the eigengene
#' ordering is to be performed.
#' @param greyLast Normally the color grey is reserved for unassigned genes;
#' hence the grey module is not a proper module and it is conventional to put
#' it last. If this is not desired, set the parameter to \code{FALSE}.
#' @param greyName Name of the grey module eigengene.
#' @param method A character string giving the method to be used calculating
#' the consensus dissimilarity. Allowed values are (abbreviations of)
#' \code{"consensus"} and \code{"majority"}. The consensus dissimilarity is
#' calculated as the maximum of given set dissimilarities for
#' \code{"consensus"} and as the average for \code{"majority"}.
#' @return A vector of lists of the same type as \code{MEs} containing the
#' re-ordered eigengenes.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{moduleEigengenes}}, \code{\link{multiSetMEs}},
#' \code{\link{orderMEs}}
#' @keywords misc
consensusOrderMEs <- function(MEs,
                              useAbs = FALSE,
                              useSets = NULL,
                              greyLast = TRUE,
                              greyName = paste0(moduleColor.getMEprefix(),
                                                "grey"),
                              method = "consensus") {
    # Debugging code:
    #printFlush("consensusOrderMEs:")
    #size = checkSets(MEs)
    #print(size)
    # end debuging code
    Diss = consensusMEDissimilarity(MEs,
                                    useAbs = useAbs,
                                    useSets = useSets,
                                    method = method)
    order = .clustOrder(Diss, greyLast, greyName)
    #print(order)
    orderMEs(
        MEs,
        greyLast = greyLast,
        greyName = greyName,
        order = order,
        useSets = useSets
    )
}

# consensusMEDissimilarity ####
#' Consensus dissimilarity of module eigengenes.
#'
#' Calculates consensus dissimilarity \code{(1-cor)} of given module eigengenes
#' relaized in several sets.
#'
#' This function calculates the individual set dissimilarities of the given
#' eigengenes in each set, then takes the (parallel) maximum or average over
#' all sets. For details on the structure of imput data, see
#' \code{\link{checkSets}}.
#'
#' @param MEs Module eigengenes of the same modules in several sets.
#' @param useAbs Controls whether absolute value of correlation should be used
#' instead of correlation in the calculation of dissimilarity.
#' @param useSets If the consensus is to include only a selection of the given
#' sets, this vector (or scalar in the case of a single set) can be used to
#' specify the selection. If \code{NULL}, all sets will be used.
#' @param method A character string giving the method to use. Allowed values
#' are (abbreviations of) \code{"consensus"} and \code{"majority"}. The
#' consensus dissimilarity is calculated as the minimum of given set
#' dissimilarities for \code{"consensus"} and as the average for
#' \code{"majority"}.
#' @return A dataframe containing the matrix of dissimilarities, with
#' \code{names} and \code{rownames} set appropriately.
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link{checkSets}}
#' @keywords misc
consensusMEDissimilarity <- function(MEs, useAbs = FALSE, useSets = NULL,
                                     method = "consensus") {
        methods = c("consensus", "majority")
        m = charmatch(method, methods)
        if (is.na(m))
            stop("Unrecognized method given. Recognized values are",
                 paste(methods, collapse  = ", "))

        nSets = length(MEs)
        MEDiss = vector(mode = "list", length = nSets)
        if (is.null(useSets))
            useSets = c(1:nSets)
        for (set in useSets)
        {
            if (useAbs)
            {
                diss = 1 - abs(cor(MEs[[set]]$data, use = "p"))
            } else
            {
                diss = 1 - cor(MEs[[set]]$data, use = "p")
            }
            MEDiss[[set]] = list(Diss = diss)
        }

        for (set in useSets)
            if (set == useSets[1])
            {
                ConsDiss = MEDiss[[set]]$Diss
            } else {
                if (m == 1) {
                    ConsDiss = pmax(ConsDiss, MEDiss[[set]]$Diss)
                } else {
                    ConsDiss = ConsDiss + MEDiss[[set]]$Diss
                }
            }

        if (m == 2)
            ConsDiss = ConsDiss / nSets

        ConsDiss = as.data.frame(ConsDiss)
        names(ConsDiss) = names(MEs[[useSets[1]]]$data)
        rownames(ConsDiss) = names(MEs[[useSets[1]]]$data)

        ConsDiss
    }

# Quantile normalization
# normalize each column such that (column) quantiles are the same
# The final value for each quantile is the 'summaryType' of the corresponding
# quantiles across the columns

.equalizeQuantiles <- function(data, summaryType = c("median", "mean")) {
        summaryType = match.arg(summaryType)
        data.sorted = apply(data, 2, sort)

        if (summaryType == "median") {
            refSample = rowMedians(data.sorted, na.rm = TRUE)
        } else if (summaryType == "mean") {
            refSample = rowMeans(data.sorted, na.rm = TRUE)
        }
        ranks = round(colRanks(data, ties.method = "average", preserveShape = TRUE))
        out = refSample [ranks]
        dim(out) = dim(data)
        dimnames(out) = dimnames(data)

        out
    }

.turnVectorIntoDist <- function(x, size, Diag, Upper) {
    attr(x, "Size") = size
    attr(x, "Diag") = FALSE
    attr(x, "Upper") = FALSE
    class(x) = c("dist", class(x))
    x
}

.turnDistVectorIntoMatrix <- function(x, size, Diag, Upper, diagValue) {
        mat = as.matrix(.turnVectorIntoDist(x, size, Diag, Upper))
        if (!Diag) {
            diag(mat) = diagValue
        }
        mat
}

# This function calculates consensus dissimilarity of module eigengenes

.consensusMEDissimilarity <- function(multiMEs,
                                      useSets = NULL,
                                      corFnc = cor,
                                      corOptions = list(use = 'p'),
                                      equalizeQuantiles = FALSE,
                                      quantileSummary = "mean",
                                      consensusQuantile = 0,
                                      useAbs = FALSE,
                                      greyMEname = "ME0") {
    nSets <- checkSets(multiMEs)$nSets
    init <- multiMEs[[1]]$data
    useMEs <- c(1:ncol(init))[names(init) != greyMEname]
    useNames <- names(init)[useMEs]
    nUseMEs <- length(useMEs)
    #  if (nUseMEs<2)
    #    stop("Something is wrong: there are two or more proper modules,
    #    but less than two proper",
    #         "eigengenes. Please check that the grey color label and module
    #         eigengene label",
    #         "are correct.")

    if (is.null(useSets)) {
        useSets <- c(1:nSets)
    }
    nUseSets <- length(useSets)
    MEDiss <- array(NA, dim = c(nUseMEs, nUseMEs, nUseSets))
    for (set in useSets) {
        corOptions$x <- multiMEs[[set]]$data[, useMEs]
        if (useAbs) {
            diss <- 1 - abs(do.call(corFnc, corOptions))
        } else {
            diss <- 1 - do.call(corFnc, corOptions)
        }
        MEDiss[, , set] <- diss
    }

    if (equalizeQuantiles) {
        distMat <- apply(MEDiss, 3, function(x) {
            as.numeric(as.dist(x))
        })
        dim(distMat) <- c(nUseMEs * (nUseMEs - 1) / 2, nUseSets)
        normalized <- .equalizeQuantiles(distMat, summaryType = quantileSummary)
        MEDiss <- apply(normalized, 2, .turnDistVectorIntoMatrix,
                        size = nUseMEs, Diag = FALSE, Upper = FALSE,
                        diagValue = 0)
    }

    ConsDiss <- apply(MEDiss, c(1:2), quantile, probs = 1 - consensusQuantile,
                      names = FALSE, na.rm = TRUE)
    colnames(ConsDiss) <- useNames
    rownames(ConsDiss) <- useNames
    ConsDiss
}

#===============================================================================
# ColorHandler.R
#===============================================================================

# A set of global variables and functions that should help handling color names
# for some 400 + modules. A vector called .GlobalStandardColors is defined that
# holds color names with first few entries being the well - known and  - loved
# colors. The rest is randomly chosen from the color names of R, excluding grey
# colors.

#-------------------------------------------------------------------------------
#
# .GlobalStandardColors
#
#-------------------------------------------------------------------------------
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

.GlobalStandardColors = c(BaseColors, ExtraColors[rank(sin(13 * c(1:nExtras) + sin(13 * c(1:nExtras))))])

rm(BaseColors, RColors, ExtraColors, nExtras, InBase)

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
#' @seealso \code{\link{hclust}}, \code{\link{cutree}},
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
            printFlush(
                paste(
                    spaces,
                    "fixDataStructure: data is not a vector of lists: converting it into one."
                )
            )
        x = data
        data = vector(mode = "list", length = 1)
        data[[1]] = list(data = x)
        rm(x)
    }
    data
}

# checkSets ####
.permissiveDim <- function(x) {
    d = dim(x)
    if (is.null(d))
        return(c(length(x), 1))
    return(d)
}



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
            nGenes = .permissiveDim(data[[useSets[1]]]$data)[2]
            for (set in useSets) {
                if (nGenes != .permissiveDim(data[[set]]$data)[2]) {
                    if (checkStructure) {
                        structureOK = FALSE
                    } else {
                        stop("Incompatible number of genes in set 1 and ", set))
                    }
                }
                nSamples[set] = .permissiveDim(data[[set]]$data)[1]
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
                            if (is.numeric(colors))
                                0
                            else
                                "grey"
                        } else
                            if (is.numeric(universalColors)) {
                                0
                            } else {
                                "grey"
                            },
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
                           "..."))
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

    greyMEname = paste0(moduleColor.getMEprefix(), unassdColor)

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

            ConsDiss = .consensusMEDissimilarity(
                MEs,
                equalizeQuantiles = equalizeQuantiles,
                quantileSummary = quantileSummary,
                consensusQuantile = consensusQuantile,
                useAbs = useAbs,
                corFnc = corFnc,
                corOptions = corOptions,
                useSets = useSets,
                greyMEname = greyMEname
            )

            Tree = fastcluster::hclust(as.dist(ConsDiss), method = "average")
            if (iteration == 1)
                oldTree = Tree
            TreeBranches = as.factor(moduleNumber(
                dendro = Tree,
                cutHeight = cutHeight,
                minSize = 1
            ))
            UniqueBranches = levels(TreeBranches)
            nBranches = nlevels(TreeBranches)
            NumberOnBranch = table(TreeBranches)
            MergedColors = colors

            # Merge modules on the same branch

            for (branch in 1:nBranches)
                if (NumberOnBranch[branch] > 1)
                {
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
#' @seealso \code{\link{softConnectivity}} for connectivity calculation in
#' weigheted networks.
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
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
# GTOMdist ####
# This function computes a TOMk dissimilarity
# which generalizes the topological overlap matrix (Ravasz et al)
# Input: an Adjacency matrix with entries in [0, 1]
# WARNING:  ONLY FOR UNWEIGHTED NETWORKS, i.e. the adjacency matrix contains
# binary entries...
# This function is explained in Yip and Horvath (2005)
# http://www.genetics.ucla.edu/labs/horvath/GTOM/
#' Generalized Topological Overlap Measure
#'
#' Generalized Topological Overlap Measure, taking into account interactions of
#' higher degree.
#'
#' @param adjMat adjacency matrix. See details below.
#' @param degree integer specifying the maximum degree to be calculated.
#' @return Matrix of the same dimension as the input \code{adjMat}.
#' @author Steve Horvath and Andy Yip
#' @references Yip A, Horvath S (2007) Gene network interconnectedness and the
#' generalized topological overlap measure. BMC Bioinformatics 2007, 8:22
#' @keywords misc
GTOMdist <- function(adjMat, degree = 1) {
    maxh1 = max(as.dist(adjMat))
    minh1 = min(as.dist(adjMat))
    if (degree != round(abs(degree)))
        stop("'degree' must be a positive integer.")
    if (maxh1 > 1 | minh1 < 0)
        stop(
            paste(
                "Entries of the adjacency matrix are not between 0 and 1: max  = ",
                maxh1,
                ", min  = ",
                minh1
            )
        )

    if (max(c(as.dist(abs(
        adjMat - t(adjMat)
    )))) > 0)
        stop("Given adjacency matrix is not symmetric.")

    B <- adjMat
    if (degree >= 2)
        for (i in 2:degree) {
            diag(B) <- diag(B) + 1
            # Calculates the number of paths with length at most degree connecting
            #  a pair
            B = B %*% adjMat
        }
    # this gives the degree - step reachability from a node to another
    B <- (B > 0)
    diag(B) <- 0 # exclude each node being its own neighbor
    # this gives the number of common degree-step-neighbor that a pair of nodes
    # share
    B <- B %*% B

    Nk <- diag(B)
    B <- B + adjMat # numerator
    diag(B) <- 1
    denomTOM = outer(Nk, Nk, FUN = "pmin") + 1 - adjMat
    diag(denomTOM) <- 1
    1 - B / denomTOM   # this turns the TOM matrix into a dissimilarity
}

# vectorTOM ####
# vectorTOM: calculate TOM of a vector (or a 'small' matrix) with expression
# data. If the number of columns in vect is small (or 1), number of columns in
# datExpr can be large.
#' Topological overlap for a subset of the whole set of genes
#'
#' This function calculates topological overlap of a small set of vectors with
#' respect to a whole data set.
#'
#' Topological overlap can be viewed as the normalized count of shared
#' neighbors encoded in an adjacency matrix. In this case, the adjacency matrix
#' is calculated between the columns of \code{vect} and \code{datExpr} and the
#' topological overlap of vectors in \code{vect} measures the number of shared
#' neighbors in \code{datExpr} that vectors of \code{vect} share.
#'
#' @param datExpr a data frame containing the expression data of the whole set,
#' with rows corresponding to samples and columns to genes.
#' @param vect a single vector or a matrix-like object containing vectors whose
#' topological overlap is to be calculated.
#' @param subtract1 logical: should calculation be corrected for
#' self-correlation? Set this to \code{TRUE} if \code{vect} contains a subset
#' of \code{datExpr}.
#' @param blockSize maximum block size for correlation calculations. Only
#' important if \code{vect} contains a large number of columns.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function.
#' @param networkType character string giving network type. Allowed values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, \code{"signed
#' hybrid"}. See \code{\link{adjacency}}.
#' @param power soft-thresholding power for network construction.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of dimensions \code{n*n}, where \code{n} is the number of
#' columns in \code{vect}.
#' @author Peter Langfelder
#' @seealso \code{\link{TOMsimilarity}} for standard calculation of topological
#' overlap.
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
vectorTOM <-function(datExpr, vect, subtract1 = FALSE, blockSize = 2000,
                     corFnc = "cor", corOptions = "use = 'p'",
                     networkType = "unsigned", power = 6, verbose = 1,
                     indent = 0) {
        spaces = indentSpaces(indent)

        intType = charmatch(networkType, .networkTypes)
        if (is.na(intType))
            stop(paste(
                "Unrecognized 'networkType'. Recognized values are",
                paste(.networkTypes, collapse = ", ")
            ))

        if (is.null(dim(vect)))
        {
            vect = as.matrix(vect)
            vectIsVector = TRUE
        } else
            vectIsVector = FALSE

        if (nrow(vect) != nrow(datExpr))
            stop("Input error: numbers of samples in 'vect' and
                 'datExpr' must be the same.")

        if (ncol(vect) > blockSize)
            stop(
                paste(
                    "Input error: number of columns in 'vect' is too large.",
                    "If you are certain you want to try anyway, increase 'blockSize'
                    to at least the number of columns in 'vect'."
                )
            )

        corEval = parse(text = paste(corFnc,
                                     "(datExpr, vect ", prepComma(corOptions), ")"))
        corVE = eval(corEval)
        if (intType == 1)
        {
            corVE = abs(corVE)
        } else if (intType == 2)
        {
            corVE = (1 + corVE) / 2
        } else if (intType == 3)
        {
            corVE[corVE < 0] = 0
        } else
            stop(
                "Unrecognized networkType argument. Recognized values are 'unsigned',
                'signed', and 'signed hybrid'."
            )

        corVE = corVE ^ power

        subtract1 = as.numeric(subtract1)

        nVect = ncol(vect)
        nGenes = ncol(datExpr)
        TOM = matrix(nrow = nGenes, ncol = nVect)

        if (verbose > 0) {
            if (verbose > 1) {
                cat(paste(
                    spaces,
                    "Calculating TOM of a set of vectors with genes"
                ))
            }
            pind = initProgInd()
        }
        start = 1
        denomArr = array(0, dim = c(2, blockSize, nVect))
        while (start <= nGenes)
        {
            end = min(start + blockSize - 1, nGenes)
            blockInd = c(start:end)
            corEval = parse(text = paste(
                corFnc,
                "(datExpr[, blockInd], datExpr ",
                prepComma(corOptions),
                ")"
            ))
            corEE = eval(corEval)
            if (intType == 1) {
                corEE = abs(corEE)
            } else if (intType == 2) {
                corEE = (1 + corEE) / 2
            } else if (intType == 3) {
                corEE[corEE < 0] = 0
            }
            corEE = corEE ^ power
            num = corEE %*% corVE  - subtract1 * corVE[blockInd,]
            kV = apply(corVE, 2, sum, na.rm = TRUE) - subtract1
            kE = apply(corEE, 1, sum, na.rm = TRUE) - 1
            denomArr[1, 1:(end - start + 1),] = matrix(kV,
                                                       nrow = end - start + 1,
                                                       ncol = nVect,
                                                       byrow = TRUE)
            denomArr[2, 1:(end - start + 1),] = matrix(kE, nrow = end - start + 1,
                                                       ncol = nVect)
            denom = apply(denomArr[, 1:(end - start + 1),], c(2, 3), min) +
                1 - corVE[blockInd,]
            TOM[blockInd,] = num / denom
            if (verbose > 0)
                pind = updateProgInd(end / nGenes, pind)
            start = end + 1
            collectGarbage()
        }
        if (verbose > 0)
            printFlush(" ")

        TOM
    }

# subsetTOM  ####
# subsetTOM: calculate TOM of a subset of vectors with respect to a full set of
# vectors.
#' Topological overlap for a subset of a whole set of genes
#'
#' This function calculates topological overlap of a subset of vectors with
#' respect to a whole data set.
#'
#' This function is designed to calculated topological overlaps of small
#' subsets of large expression data sets, for example in individual modules.
#'
#' @param datExpr a data frame containing the expression data of the whole set,
#' with rows corresponding to samples and columns to genes.
#' @param subset a single logical or numeric vector giving the indices of the
#' nodes for which the TOM is to be calculated.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function.
#' @param networkType character string giving network type. Allowed values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, \code{"signed
#' hybrid"}. See \code{\link{adjacency}}.
#' @param power soft-thresholding power for network construction.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of dimensions \code{n*n}, where \code{n} is the number of
#' entries selected by \code{block}.
#' @author Peter Langfelder
#' @seealso \code{\link{TOMsimilarity}} for standard calculation of topological
#' overlap.
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
subsetTOM <- function(datExpr,
                      subset,
                      corFnc = "cor",
                      corOptions = "use = 'p'",
                      networkType = "unsigned",
                      power = 6,
                      verbose = 1,
                      indent = 0) {
    spaces = indentSpaces(indent)

    if (!is.null(dim(subset)))
        stop("'subset' must be a dimensionless vector.")

    if (is.null(dim(datExpr)))
        stop("'datExpr' must be a matrix or data frame.")
    if (length(dim(datExpr)) != 2)
        stop("'datExpr' must be two - dimensional.")

    nGenes = ncol(datExpr)

    if (is.logical(subset))
        subset = c(1:nGenes)[subset]

    nBlock = length(subset)

    if (any(!is.finite(subset)))
        stop("Entries of 'subset' must all be finite.")

    if (min(subset) < 1 | max(subset) > nGenes)
        stop(
            paste(
                "Some entries of 'subset' are out of range.",
                "\nNote: 'subset' must contain indices of the subset
                for which the TOM is calculated."
            )
        )

    intType = charmatch(networkType, .networkTypes)
    if (is.na(intType))
        stop(paste(
            "Unrecognized 'networkType'. Recognized values are",
            paste(.networkTypes, collapse = ", ")
        ))

    adj = adjacency(
        datExpr,
        subset,
        power = power,
        type = networkType,
        corFnc = corFnc,
        corOptions = corOptions
    )

    adj[is.na(adj)] = 0
    num = t(adj) %*% adj - adj[subset,]

    k = apply(adj, 2, sum)

    kMat = matrix(k, nBlock, nBlock)

    denom = pmin(kMat, t(kMat)) - adj[subset,]

    TOM = num / denom
    diag(TOM) = 1

    TOM
}

#-------------------------------------------------------------------------------
#
# adjacency
#
#-------------------------------------------------------------------------------
# Computes the adjacency from the expression data: takes cor, transforms it as
# appropriate and possibly adds a sign if requested. No subselection on datExpr
# is performed.
# A slighly reworked version that assumes one wants the adjacency matrix of data
# with itself or a subset. The data are given only once, and an additional
# selection index for columns is given.
# Caution: no checking of selectCols validity is performed.

.compiledAdjacency <- function(expr,
                               corType = "pearson",
                               networkType = "unsigned",
                               power = 6,
                               maxPOutliers = 1,
                               quickCor = 0,
                               pearsonFallback = "individual",
                               cosineCorrelation = FALSE,
                               nThreads = 0,
                               verbose = 1,
                               indent = 0) {
    corTypeC = as.integer(pmatch(corType, .corTypes) - 1)
    if (is.na(corTypeC))
        stop(paste(
            "Invalid 'corType'. Recognized values are",
            paste(.corTypes, collapse = ", ")
        ))

    if ((maxPOutliers < 0) | (maxPOutliers > 1))
        stop("maxPOutliers must be between 0 and 1.")
    if (quickCor < 0)
        stop("quickCor must be positive.")
    if ((maxPOutliers < 0) | (maxPOutliers > 1))
        stop("maxPOutliers must be between 0 and 1.")

    fallback = as.integer(pmatch(pearsonFallback, .pearsonFallbacks))
    if (is.na(fallback))
        stop(
            paste(
                "Unrecognized 'pearsonFallback'. Recognized values are
                (unique abbreviations of)\n",
                paste(.pearsonFallbacks, collapse = ", ")
            )
        )

    if (nThreads < 0)
        stop("nThreads must be positive.")
    if (is.null(nThreads) ||
        (nThreads == 0))
        nThreads = .useNThreads()

    if ((power < 1) |
        (power > 30))
        stop("power must be between 1 and 30.")

    networkTypeC = as.integer(charmatch(networkType, .networkTypes) - 1)
    if (is.na(networkTypeC))
        stop(
            paste(
                "Unrecognized networkType argument.",
                "Recognized values are (unique abbreviations of)",
                paste(.networkTypes, collapse = ", ")
            )
        )
    dimEx = dim(expr)
    if (length(dimEx) != 2)
        stop("expr has incorrect dimensions.")
    nGenes = dimEx[2]
    nSamples = dimEx[1]
    warn = as.integer(0)

    expr = as.matrix(expr)

    adj = matrix(0, nGenes, nGenes)
    err = as.integer(0)

    res = .C(
        "testAdjacency",
        as.double(expr),
        as.integer(nSamples),
        as.integer(nGenes),
        as.integer(corTypeC),
        as.integer(networkTypeC),
        as.double(power),
        as.double(maxPOutliers),
        as.double(quickCor),
        as.integer(fallback),
        as.integer(cosineCorrelation),
        adj = as.double(adj),
        as.integer(err),
        as.integer(warn),
        as.integer(nThreads),
        NAOK = TRUE
    )

    adj = res$adj
    dim(adj) = c(nGenes, nGenes)
    adj
}




# A presumably faster and less memory - intensive version, only for "unsigned"
# networks.





################################################################################
################################################################################
# C) Defining gene modules using clustering procedures
################################################################################
################################################################################




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
#' @seealso \code{\link{hclust}} for hierarchical clustering,
#' \code{\link{cutree}} and \code{\link{cutreeStatic}} for other
#' constant-height branch cuts, \code{\link{standardColors}} to convert the
#' retuned numerical lables into colors for easier visualization.
#' @keywords misc
cutreeStatic <- function(dendro,
                         cutHeight = 0.9,
                         minSize = 50) {
    normalizeLabels(moduleNumber(dendro, cutHeight, minSize))
}



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
#' @param dendro a hierarchical clustering dendrogram such as returned by
#' \code{\link{hclust}}.
#' @param cutHeight height at which branches are to be cut.
#' @param minSize minimum number of object on a branch to be considered a
#' cluster.
#' @return A character vector giving color labels of objects, with "grey"
#' meaning unassigned. The largest cluster is conventionally labeled
#' "turquoise", next "blue" etc. Run \code{standardColors()} to see the
#' sequence of standard color labels.
#' @author Peter Langfelder
#' @seealso \code{\link{hclust}} for hierarchical clustering,
#' \code{\link{cutree}} and \code{\link{cutreeStatic}} for other
#' constant-height branch cuts, \code{\link{standardColors}} to see the
#' sequence of color labels that can be assigned.
#' @keywords misc
cutreeStaticColor <-
    function(dendro,
             cutHeight = 0.9,
             minSize = 50) {
        labels2colors(normalizeLabels(moduleNumber(dendro, cutHeight, minSize)))
    }




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
#' @author Steve Horvath \email{SHorvath@@mednet.ucla.edu} and Peter Langfelder
#' \email{Peter.Langfelder@@gmail.com}
#' @seealso \code{\link[dynamicTreeCut]{cutreeDynamic}} for module detection in
#' a dendrogram;
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
    if (!is.null(rowText))
    {
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

#===============================================================================
# This function can be used to create an average linkage hierarchical
# clustering tree
# or the microarray samples. The rows of datExpr correspond to the samples and
# the columns to the genes
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
#' @seealso \code{\link[stats]{dist}}, \code{\link[stats]{hclust}},
#' \code{\link{plotDendroAndColors}}
#' @keywords hplot misc
plotClusterTreeSamples <-
    function(datExpr,
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
            plot(
                dendro,
                main = main,
                sub = "",
                xlab = "",
                labels = dendroLabels,
                cex = cex.dendroLabels
            )
            if (saveMar)
                par(oldMar)
        } else {
            if (is.null(traitLabels))
                traitLabels = names(as.data.frame(y))
            y = as.matrix(y)
            if (!is.numeric(y)) {
                warning(paste(
                    "The microarray sample trait y will be transformed to
                    numeric."
                ))
                dimy = dim(y)
                y = as.numeric(y)
                dim(y) = dimy
            } # end of if (!is.numeric(y))
            if (nrow(as.matrix(datExpr))  != nrow(y))
                stop(
                    paste(
                        "Input Error: dim(as.matrix(datExpr))[[1]]  != length(y)\n",
                        "  In plain English: The number of microarray sample arrays does
                        not match the number",
                        "of samples for the trait.\n",
                        "   Hint: Make sure rows of 'datExpr' (and 'y', if it is a
                        matrix) correspond to samples."
                    )
                )

            if (is.integer(y)) {
                y = y - min(0, min(y, na.rm = TRUE)) + 1
            } else {
                y = (y >= median(y, na.rm = TRUE)) + 1
            }
            plotDendroAndColors(
                dendro,
                colors = y,
                groupLabels = traitLabels,
                rowText = yLabels,
                setLayout = setLayout,
                autoColorHeight = autoColorHeight,
                colorHeight = colorHeight,
                addGuide = addGuide,
                guideAll = guideAll,
                guideCount = guideCount,
                guideHang = guideHang,
                cex.colorLabels = cex.traitLabels,
                cex.dendroLabels = cex.dendroLabels,
                marAll = marAll,
                saveMar = saveMar,
                abHeight = abHeight,
                abCol = abCol,
                main = main,
                ...
            )
        }
    }

#===============================================================================
# The function TOMplot creates a TOM plot
# Inputs:  distance measure, hierarchical (hclust) object, color label = colors



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
TOMplot <-
    function(dissim,
             dendro,
             Colors = NULL,
             ColorsLeft = Colors,
             terrainColors = FALSE,
             setLayout = TRUE,
             ...) {
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
plotNetworkHeatmap <-
    function(datExpr,
             plotGenes,
             useTOM = TRUE,
             power = 6 ,
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

################################################################################
################################################################################
# E) Relating a measure of gene significance to the modules
################################################################################
################################################################################

#===============================================================================
# The function ModuleEnrichment1 creates a bar plot that shows whether modules
# are enriched with significant genes.
# More specifically, it reports the mean gene significance for each module.
# The gene significance can be a binary variable or a quantitative variable.
# It also plots the 95% confidence interval of the mean
# (CI = mean +/-  1.96 * standard error).
# It also reports a Kruskal Wallis p-value.



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
plotModuleSignificance <-
    function(geneSignificance,
             colors,
             boxplot = FALSE,
             main = "Gene significance across modules, ",
             ylab = "Gene Significance",
             ...)
    {
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

################################################################################
################################################################################
# G) Miscellaneous other functions, e.g. for computing the cluster coefficient.
################################################################################
################################################################################


# The function signedKME computes the module eigengene based connectivity.
# Input: datExpr = a possibly very large gene expression data set where the rows
# correspond to samples and the columns represent genes
# datME = data frame of module eigengenes (columns correspond to module
# eigengenes or MEs)
# A module eigengene based connectivity KME value will be computed if the gene
# has a non missing expression value in at least minNSamples arrays.
# Output a data frame where columns are the KME values
# corresponding to different modules.
# By splitting the expression data into different blocks, the function can deal
# with tens of thousands of gene expression data.
# If there are many eigengenes (say hundreds) consider decreasing the block size



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
    if (dim(as.matrix(datME))[[1]]  != dim(as.matrix(datExpr))[[1]])
        stop("Number of samples (rows) in 'datExpr' and 'datME' must be the same.")
    varianceZeroIndicatordatExpr = as.vector(apply(as.matrix(datExpr), 2, var,
                                                   na.rm = TRUE)) == 0
    varianceZeroIndicatordatME  = as.vector(apply(as.matrix(datME), 2, var,
                                                  na.rm = TRUE)) == 0
    if (sum(varianceZeroIndicatordatExpr, na.rm = TRUE) > 0)
        warning("Some genes are constant. Hint: consider removing constant columns
                from datExpr.")
    if (sum(varianceZeroIndicatordatME, na.rm = TRUE) > 0)
        warning(
            paste(
                "Some module eigengenes are constant, which is suspicious.\n",
                "    Hint: consider removing constant columns from datME."
            )
        )
    no.presentdatExpr = as.vector(apply(!is.na(as.matrix(datExpr)), 2, sum))
    if (min(no.presentdatExpr) < ..minNSamples)
        warning(
            paste(
                "Some gene expressions have fewer than 4 observations.\n",
                "    Hint: consider removing genes with too many missing values or
                collect more arrays."
            )
        )

    #output = data.frame(cor(datExpr, datME, use = "p"))
    corExpr = parse(text = paste(
        "data.frame(",
        corFnc,
        "(datExpr, datME ",
        prepComma(corOptions),
        "))"
    ))
    output = eval(corExpr)

    output[no.presentdatExpr < ..minNSamples,] = NA
    names(output) = paste0(outputColumnName, substring(names(datME), first = 3,
                                                       last = 100))
    dimnames(output)[[1]] = names(datExpr)
    output
} # end of function signedKME




#===============================================================================
# The function clusterCoef computes the cluster coefficients.
# Input is an adjacency matrix



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
    computeLinksInNeighbors <-
        function(x, imatrix) {
            x %*% imatrix %*% x
        }
    nolinksNeighbors <- c(rep(-666, nNodes))
    total.edge <- c(rep(-666, nNodes))
    maxh1 = max(as.dist(adjMat))
    minh1 = min(as.dist(adjMat))
    if (maxh1 > 1 | minh1 < 0)
        stop(
            paste(
                "The adjacency matrix contains entries that are larger than 1 or
                smaller than 0: max  = ",
                maxh1,
                ", min  = ",
                minh1
            )
        )
    nolinksNeighbors <-
        apply(adjMat, 1, computeLinksInNeighbors, imatrix = adjMat)
    plainsum  <- apply(adjMat, 1, sum)
    squaresum <- apply(adjMat ^ 2, 1, sum)
    total.edge = plainsum ^ 2 - squaresum
    CChelp = rep(-666, nNodes)
    CChelp = ifelse(total.edge == 0, 0, nolinksNeighbors / total.edge)
    CChelp
} # end of function



#===============================================================================
# The function addErrorBars  is used to create error bars in a barplot
# usage: addErrorBars(as.vector(means), as.vector(stderrs), two.side = FALSE)


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

#===============================================================================
# this function computes the standard error


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
dynamicMergeCut <- function(n,
                            mergeCor = .9,
                            Zquantile = 2.35) {
    if (mergeCor > 1 |
        mergeCor < 0)
        stop("'mergeCor' must be between 0 and 1.")
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
    if (mergethreshold > 1)
        1
    else
        mergethreshold
}# end of function dynamicMergeCut



#===============================================================================
#
# print.flush
#
#===============================================================================

#print.flush <- function(...)
#{
#   printFlush(...)
#}



#===============================================================================
#
# Correlation p-value for multiple correlation values
#
#===============================================================================



#' Fisher's asymptotic p-value for correlation
#'
#' Calculates Fisher's asymptotic p-value for given correlations.
#'
#'
#' @param cor A vector of correlation values whose corresponding p-values are
#' to be calculated
#' @param nSamples Number of samples from which the correlations were
#' calculated
#' @param twoSided logical: should the calculated p-values be two sided?
#' @return A vector of p-values of the same length as the input correlations.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
corPvalueFisher <- function(cor, nSamples, twoSided = TRUE) {
    if (sum(abs(cor) > 1, na.rm = TRUE) > 0)
        stop("Some entries in 'cor' are out of normal range  - 1 to 1.")
    if (twoSided) {
        z = abs(0.5 * log((1 + cor) / (1 - cor)) * sqrt(nSamples - 3))
        2 * pnorm(-z)
    } else {
        # return a small p-value for positive correlations
        z = -0.5 * log((1 + cor) / (1 - cor)) * sqrt(nSamples - 3)
        pnorm(-z)
    }
}

# this function compute an asymptotic p-value for a given correlation (r) and
# sample size (n). Needs a new name before we commit it to the package.



#' Student asymptotic p-value for correlation
#'
#' Calculates Student asymptotic p-value for given correlations.
#'
#'
#' @param cor A vector of correlation values whose corresponding p-values are
#' to be calculated
#' @param nSamples Number of samples from which the correlations were
#' calculated
#' @return A vector of p-values of the same length as the input correlations.
#' @author Steve Horvath and Peter Langfelder
#' @keywords misc
corPvalueStudent <- function(cor, nSamples) {
    T = sqrt(nSamples - 2) * cor / sqrt(1 - cor ^ 2)
    2 * pt(abs(T), nSamples - 2, lower.tail = FALSE)
}


################################################################################



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
propVarExplained <- function(datExpr,
                             colors,
                             MEs,
                             corFnc = "cor",
                             corOptions = "use = 'p'") {
    fc = as.factor(colors)
    mods = levels(fc)
    nMods = nlevels(fc)
    nGenes = ncol(datExpr)
    if (nMods != ncol(MEs))
        stop(
            paste(
                "Input error: number of distinct 'colors' differs from\n",
                "     the number of module eigengenes given in ME."
            )
        )

    if (ncol(datExpr) != length(colors))
        stop(
            "Input error: number of probes (columns) in 'datExpr' differs from the
            length of goven 'colors'."
        )

    if (nrow(datExpr) != nrow(MEs))
        stop("Input error: number of observations (rows) in 'datExpr' and 'MEs'
             differ.")

    PVE = rep(0, nMods)

    col2MEs = match(mods, substring(names(MEs), 3))

    if (sum(is.na(col2MEs)) > 0)
        stop("Input error: not all given colors could be matched to names of module
             eigengenes.")

    for (mod in 1:nMods)
    {
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


#===============================================================================
#
# addGrid
#
#===============================================================================
# This function adds a horizontal grid to a plot



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
addGrid <- function(linesPerTick = NULL,
                    horiz = TRUE,
                    vert = FALSE,
                    col = "grey30",
                    lty = 3) {
    box = par("usr")
    if (horiz) {
        ticks = par("yaxp")
        nTicks = ticks[3]
        if (is.null(linesPerTick)) {
            if (nTicks < 6)
                linesPerTick = 5
            else
                linesPerTick = 2
        }
        spacing = (ticks[2] - ticks[1]) / (linesPerTick * nTicks)
        first = ceiling((box[3] - ticks[1]) / spacing)
        last = floor((box[4] - ticks[1]) / spacing)
        #print(paste("addGrid: first = ", first, ", last  = ", last, "box = ",
        #paste(signif (box, 2), collapse = ", "),
        #"ticks = ", paste(signif (ticks, 2), collapse = ", "),
        #"spacing  = ", spacing))
        for (k in first:last) {
            lines(
                x = box[c(1, 2)],
                y = rep(ticks[1] + spacing * k, 2),
                col = col,
                lty = lty
            )
        }
    }
    if (vert) {
        ticks = par("xaxp")
        nTicks = ticks[3]
        if (is.null(linesPerTick)) {
            if (nTicks < 6)
                linesPerTick = 5
            else
                linesPerTick = 2
        }
        spacing = (ticks[2] - ticks[1]) / (linesPerTick * ticks[3])
        first = ceiling((box[1] - ticks[1]) / spacing)
        last = floor((box[2] - ticks[1]) / spacing)
        #print(paste("addGrid: first = ", first, ", last  = ", last, "box = ",
        #paste(signif (box, 2), collapse = ", "),
        #            "ticks = ", paste(signif (ticks, 2), collapse = ", "),
        #            "spacing  = ", spacing))
        for (l in first:last) {
            lines(
                x = rep(ticks[1] + spacing * l, 2),
                y = box[c(3, 4)],
                col = col,
                lty = lty
            )
        }
    }
}

#-------------------------------------------------------------------------------
#
# Add vertical "guide" lines to a dendrogram to facilitate identification of
# clusters with color bars
#
#-------------------------------------------------------------------------------



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
addGuideLines <-
    function(dendro,
             all = FALSE,
             count = 50,
             positions = NULL,
             col = "grey30",
             lty = 3,
             hang = 0) {
        if (all) {
            positions = 1:(length(dendro$height) + 1)
        } else {
            if (is.null(positions)) {
                lineSpacing = (length(dendro$height) + 1) / count
                positions = (1:count) *  lineSpacing
            }
        }
        objHeights = rep(0, length(dendro$height + 1))
        objHeights[-dendro$merge[dendro$merge[, 1] < 0, 1]] = dendro$height[dendro$merge[, 1] <
                                                                                0]
        objHeights[-dendro$merge[dendro$merge[, 2] <  0, 2]] = dendro$height[dendro$merge[, 2] <
                                                                                 0]
        box = par("usr")
        ymin = box[3]
        ymax = box[4]
        objHeights = objHeights - hang * (ymax - ymin)
        objHeights[objHeights < ymin] = ymin
        posHeights = pmin(objHeights[dendro$order][floor(positions)],
                          objHeights[dendro$order][ceiling(positions)])
        for (line in 1:length(positions)) {
            # The last guide line is superfluous
            lines(
                x = rep(positions[line], 2),
                y = c(ymin, posHeights[line]),
                lty = 3,
                col = col
            )
        }
    }

#-------------------------------------------------------------------------------
#
# nearestNeighborConnectivity
#
#-------------------------------------------------------------------------------
# This function takes expression data (rows = samples, colummns = genes)
# and the power exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of
# dimensions nGenes * nSets containing the connectivities of each gene in each
# subset.



#' Connectivity to a constant number of nearest neighbors
#'
#' Given expression data and basic network parameters, the function calculates
#' connectivity of each gene to a given number of nearest neighbors.
#'
#' Connectivity of gene \code{i} is the sum of adjacency strengths between gene
#' \code{i} and other genes; in this case we take the \code{nNeighbors} nodes
#' with the highest connection strength to gene \code{i}. The adjacency
#' strengths are calculated by correlating the given expression data using the
#' function supplied in \code{corFNC} and transforming them into adjacency
#' according to the given network \code{type} and \code{power}.
#'
#' @param datExpr a data frame containing expression data, with rows
#' corresponding to samples and columns to genes. Missing values are allowed
#' and will be ignored.
#' @param nNeighbors number of nearest neighbors to use.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param type a character string encoding network type. Recognized values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored upon
#' exit.
#' @param verbose integer controlling the level of verbosity. 0 means silent.
#' @param indent integer controlling indentation of output. Each unit above 0
#' adds two spaces.
#' @return A vector with one component for each gene containing the nearest
#' neighbor connectivity.
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}}, \code{\link{softConnectivity}}
#' @keywords misc
nearestNeighborConnectivity <-
    function(datExpr,
             nNeighbors = 50,
             power = 6,
             type = "unsigned",
             corFnc = "cor",
             corOptions = "use = 'p'",
             blockSize = 1000,
             sampleLinks = NULL,
             nLinks = 5000,
             setSeed = 38457,
             verbose = 1,
             indent = 0)
    {
        spaces = indentSpaces(indent)
        nGenes = dim(datExpr)[2]
        nSamples = dim(datExpr)[1]

        if (is.null(sampleLinks)) {
            sampleLinks = (nGenes > nLinks)
        }

        if (sampleLinks)
            nLinks = min(nLinks, nGenes)
        else
            nLinks = nGenes

        #printFlush(paste("blockSize  = ", blockSize))
        #printFlush(paste("nGenes  = ", nGenes))
        #printFlush(paste(".largestBlockSize  = ", .largestBlockSize))

        if (blockSize * nLinks > .largestBlockSize) {
            blockSize = as.integer(.largestBlockSize / nLinks)
        }
        intNetworkType = charmatch(type, .networkTypes)
        if (is.na(intNetworkType))
            stop(
                paste(
                    "Unrecognized networkType argument. Recognized values are
                    (unique abbreviations of)",
                    paste(.networkTypes, collapse = ", ")
                )
            )

        subtract = rep(1, nGenes)
        if (sampleLinks) {
            if (verbose > 0) {
                printFlush(
                    paste(
                        spaces,
                        "nearestNeighborConnectivity: selecting sample pool
                        of size",
                        nLinks,
                        ".."
                    )
                )
            }
            sd = apply(datExpr, 2, sd, na.rm = TRUE)
            order = order(-sd)
            saved = FALSE
            if (exists(".Random.seed")) {
                saved = TRUE
                savedSeed = .Random.seed
                if (is.numeric(setSeed))
                    set.seed(setSeed)
            }
            samplePool = order[sample(x = nGenes, size = nLinks)]
            if (saved) {
                .Random.seed <<- savedSeed
            }
            poolExpr = datExpr[, samplePool]
            subtract[-samplePool] = 0
        }

        if (verbose > 0) {
            printFlush(
                paste(
                    spaces,
                    "nearestNeighborConnectivity: received",
                    "dataset with nGenes  = ",
                    as.character(nGenes)
                )
            )
            cat(
                paste(
                    spaces,
                    "..using nNeighbors  = ",
                    nNeighbors,
                    "and blockSize  = ",
                    blockSize,
                    "  "
                )
            )
            pind = initProgInd(trailStr = " done")
        }

        nearestNeighborConn = rep(0, nGenes)

        nBlocks = as.integer((nGenes - 1) / blockSize)
        SetRestrConn = NULL
        start = 1
        if (sampleLinks) {
            corEval = parse(text = paste(
                corFnc,
                "(poolExpr, datExpr[, blockIndex] ",
                prepComma(corOptions),
                ")"
            ))
        } else {
            corEval = parse(text = paste(
                corFnc,
                "(datExpr, datExpr[, blockIndex] ",
                prepComma(corOptions),
                ")"
            ))
        }

        while (start <= nGenes) {
            end = start + blockSize - 1
            if (end > nGenes)
                end = nGenes
            blockIndex = c(start:end)
            #if (verbose > 1) printFlush(paste(spaces, "..working on genes", start,
            #"through", end, "of", nGenes))
            c = eval(corEval)
            if (intNetworkType == 1) {
                c = abs(c)
            } else if (intNetworkType == 2) {
                c = (1 + c) / 2
            } else if (intNetworkType == 3) {
                c[c < 0] = 0
            } else {
                stop(
                    "Internal error: intNetworkType has wrong value:",
                    intNetworkType,
                    ". Sorry!"
                )
            }
            adj_mat = as.matrix(c ^ power)
            adj_mat[is.na(adj_mat)] = 0
            sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors + 1),])
            nearestNeighborConn[blockIndex] = apply(sortedAdj, 2, sum) - subtract[blockIndex]
            start = end + 1
            if (verbose > 0)
                pind = updateProgInd(end / nGenes, pind)
            collectGarbage()
        }
        if (verbose > 0)
            printFlush(" ")
        nearestNeighborConn
    }


#Try to merge this with the single - set function.
#-------------------------------------------------------------------------------
#
# nearestNeighborConnectivityMS
#
#-------------------------------------------------------------------------------
# This function takes expression data (rows = samples, colummns = genes) in the
# multi-set format and the power exponent used in weighting the
# correlations to get the network adjacency matrix, and returns an array of
# dimensions nGenes * nSets containing the connectivities of each gene in each
# subset.



#' Connectivity to a constant number of nearest neighbors across multiple data
#' sets
#'
#' Given expression data from several sets and basic network parameters, the
#' function calculates connectivity of each gene to a given number of nearest
#' neighbors in each set.
#'
#' Connectivity of gene \code{i} is the sum of adjacency strengths between gene
#' \code{i} and other genes; in this case we take the \code{nNeighbors} nodes
#' with the highest connection strength to gene \code{i}. The adjacency
#' strengths are calculated by correlating the given expression data using the
#' function supplied in \code{corFNC} and transforming them into adjacency
#' according to the given network \code{type} and \code{power}.
#'
#' @param multiExpr expression data in multi-set format. A vector of lists, one
#' list per set. In each list there must be a component named \code{data} whose
#' content is a matrix or dataframe or array of dimension 2 containing the
#' expression data. Rows correspond to samples and columns to genes (probes).
#' @param nNeighbors number of nearest neighbors to use.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param type a character string encoding network type. Recognized values are
#' (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored after.
#' @param verbose integer controlling the level of verbosity. 0 means silent.
#' @param indent integer controlling indentation of output. Each unit above 0
#' adds two spaces.
#' @return A matrix in which columns correspond to sets and rows to genes; each
#' entry contains the nearest neighbor connectivity of the corresponding gene.
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}}, \code{\link{softConnectivity}},
#' \code{\link{nearestNeighborConnectivity}}
#' @keywords misc
nearestNeighborConnectivityMS <-
    function(multiExpr,
             nNeighbors = 50,
             power = 6,
             type = "unsigned",
             corFnc = "cor",
             corOptions = "use = 'p'",
             blockSize = 1000,
             sampleLinks = NULL,
             nLinks = 5000,
             setSeed = 36492,
             verbose = 1,
             indent = 0) {
        spaces = indentSpaces(indent)
        setsize = checkSets(multiExpr)
        nGenes = setsize$nGenes
        nSamples = setsize$nSamples
        nSets = setsize$nSets

        if (is.null(sampleLinks)) {
            sampleLinks = (nGenes > nLinks)
        }

        if (sampleLinks)
            nLinks = min(nLinks, nGenes)
        else
            nLinks = nGenes

        #printFlush(paste("blockSize  = ", blockSize))
        #printFlush(paste("nGenes  = ", nGenes))
        #printFlush(paste(".largestBlockSize  = ", .largestBlockSize))

        if (blockSize * nLinks > .largestBlockSize)
            blockSize = as.integer(.largestBlockSize / nLinks)

        if (length(power) == 1)
        {
            power = rep(power, nSets)
        } else if (length(power) != nSets)
            stop("Invalid arguments: length of 'power' must equal number sets in
                 'multiExpr'")

        intNetworkType = charmatch(type, .networkTypes)
        if (is.na(intNetworkType))
            stop(
                paste(
                    "Unrecognized networkType argument. Recognized values are
                    (unique abbreviations of)",
                    paste(.networkTypes, collapse = ", ")
                )
            )

        subtract = rep(1, nGenes)
        if (sampleLinks) {
            if (verbose > 0)
                printFlush(
                    paste(
                        spaces,
                        "nearestNeighborConnectivityMS: selecting
                        sample pool of size",
                        nLinks,
                        ".."
                    )
                )
            sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE)
            order = order(-sd)
            saved = FALSE
            if (exists(".Random.seed"))
            {
                saved = TRUE
                savedSeed = .Random.seed
                if (is.numeric(setSeed))
                    set.seed(setSeed)
            }
            samplePool = order[sample(x = nGenes, size = nLinks)]
            if (saved) {
                .Random.seed <<- savedSeed
            }
            subtract[-samplePool] = 0
        }

        if (verbose > 0)
            printFlush(
                paste(
                    spaces,
                    "nearestNeighborConnectivityMS: received",
                    nSets,
                    "datasets with nGenes  = ",
                    as.character(nGenes)
                )
            )
        if (verbose > 0)
            printFlush(paste(spaces, "  Using nNeighbors  = ",
                             nNeighbors))

        nearestNeighborConn = matrix(nrow = nGenes, ncol = nSets)

        if (sampleLinks) {
            corEval = parse(
                text = paste(
                    corFnc,
                    "(multiExpr[[set]]$data[, samplePool],
                    multiExpr[[set]]$data[, blockIndex] ",
                    prepComma(corOptions),
                    ")"
                )
            )
        } else {
            corEval = parse(
                text = paste(
                    corFnc,
                    "(multiExpr[[set]]$data,
                    multiExpr[[set]]$data[, blockIndex] ",
                    prepComma(corOptions),
                    ")"
                )
            )
        }


        for (set in 1:nSets) {
            if (verbose > 0) {
                cat(paste(spaces, "  Working on set", set))
                pind = initProgInd(trailStr = " done")
            }
            nBlocks = as.integer((nGenes - 1) / blockSize)
            SetRestrConn = NULL
            start = 1
            while (start <= nGenes) {
                end = start + blockSize - 1
                if (end > nGenes)
                    end = nGenes
                blockIndex = c(start:end)
                #if (verbose > 1) printFlush(paste(spaces, " .. working on genes",
                #start, "through", end, "of", nGenes))
                c = eval(corEval)
                if (intNetworkType == 1) {
                    c = abs(c)
                } else if (intNetworkType == 2) {
                    c = (1 + c) / 2
                } else if (intNetworkType == 3) {
                    c[c < 0] = 0
                } else {
                    stop(
                        "Internal error: intNetworkType has wrong value:",
                        intNetworkType,
                        ". Sorry!"
                    )
                }
                adj_mat = as.matrix(c ^ power[set])
                adj_mat[is.na(adj_mat)] = 0
                sortedAdj = as.matrix(apply(adj_mat, 2, sort, decreasing = TRUE)[1:(nNeighbors + 1),])
                nearestNeighborConn[blockIndex, set] = apply(sortedAdj, 2,
                                                             sum) -
                    subtract[blockIndex]
                collectGarbage()
                start = end + 1
                if (verbose > 0)
                    pind = updateProgInd(end / nGenes, pind)
                collectGarbage()
            }
            if (verbose > 0)
                printFlush(" ")
        }
        nearestNeighborConn
    }

#===============================================================================
#
# Nifty display of progress.
#
#===============================================================================
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

#' Update the progress as it goes
#'
#' Should work on all OS to indicate the evolution of indicator
#' @rdname initProgInd
#' @param newFrac Unkown parameter.
#' @param progInd Object of class progressIndicator obtained with initProgInd.
#' @export
updateProgInd <- function(newFrac, progInd, quiet = !interactive())
{
    if (class(progInd) != "progressIndicator")
        stop(
            "Parameter progInd is not of class 'progressIndicator'.
            Use initProgInd() to initialize",
            "it prior to use."
        )

    newStr = paste0(progInd$leadStr,
                    as.integer(newFrac * 100),
                    "% ",
                    progInd$trailStr)
    if (newStr != progInd$oldStr)
    {
        if (quiet)
        {
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

#===============================================================================
#
# Plot a dendrogram and a set of labels underneath
#
#===============================================================================
#



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

################################################################################
#
#  Functions included from NetworkScreeningFunctions
#
################################################################################

# this function creates pairwise scatter plots between module eigengenes (above
# the diagonal) Below the diagonal are the absolute values of the Pearson
# correlation coefficients.
# The diagonal contains histograms of the module eigengene expressions.



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


#-------------------------------------------------------------------------------
#
# corPredictionSuccess
#
#-------------------------------------------------------------------------------

# The function corPredictionSuccess can be used to determine which method is
# best for predicting correlations in a new test set. corTestSet should be a
# vector of correlations in the test set.
# The parameter topNumber specifies that the top number most positive and the
# top most negative predicted correlations
#  TopNumber is a vector of integers.
# corPrediction should be a data frame of predictions for the correlations.
# Output a list with the following components:
# meancorTestSetPositive = mean test set correlation among the topNumber of
# genes which are predicted to have positive correlations.
# meancorTestSetNegative = mean test set correlation among the topNumber of
# genes which are predicted to have negative correlations.
# meancorTestSetOverall = (meancorTestSetPositive - meancorTestSetNegative)/2



#' Qunatification of success of gene screening
#'
#' This function calculates the success of gene screening.
#'
#' For each column in \code{corPrediction}, the function evaluates the mean
#' \code{corTestSet} for the number of top genes (ranked by the column in
#' \code{corPrediction}) given in \code{topNumber}. The higher the mean
#' \code{corTestSet} (for positive \code{corPrediction}) or negative (for
#' negative \code{corPrediction}), the more successful the prediction.
#'
#' @param corPrediction a vector or a matrix of prediction statistics
#' @param corTestSet correlation or other statistics on test set
#' @param topNumber a vector of the number of top genes to consider
#' @return \item{meancorTestSetOverall }{ difference of
#' \code{meancorTestSetPositive} and \code{meancorTestSetNegative} below }
#' \item{meancorTestSetPositive}{ mean \code{corTestSet} on top genes with
#' positive \code{corPrediction} } \item{meancorTestSetNegative}{ mean
#' \code{corTestSet} on top genes with negative \code{corPrediction} } ...
#' @author Steve Horvath
#' @seealso \code{\link{relativeCorPredictionSuccess}}
#' @keywords misc
corPredictionSuccess <-
    function(corPrediction, corTestSet, topNumber = 100)
    {
        nPredictors = dim(as.matrix(corPrediction))[[2]]
        nGenes = dim(as.matrix(corPrediction))[[1]]
        if (length(as.numeric(corTestSet)) != nGenes)
            stop("non - compatible dimensions of 'corPrediction' and 'corTestSet'")
        out1 = rep(NA, nPredictors)
        meancorTestSetPositive = matrix(ncol = nPredictors,
                                        nrow = length(topNumber))
        meancorTestSetNegative = matrix(ncol = nPredictors,
                                        nrow = length(topNumber))
        for (i in c(1:nPredictors)) {
            rankpositive = rank(-as.matrix(corPrediction)[, i],
                                ties.method = "first")
            ranknegative = rank(as.matrix(corPrediction)[, i],
                                ties.method = "first")
            for (j in c(1:length(topNumber))) {
                meancorTestSetPositive[j, i] = mean(corTestSet[rankpositive <= topNumber[j]], na.rm = TRUE)
                meancorTestSetNegative[j, i] = mean(corTestSet[ranknegative <= topNumber[j]], na.rm = TRUE)
            } # end of j loop over topNumber
        } # end of i loop over predictors
        meancorTestSetOverall = data.frame((meancorTestSetPositive - meancorTestSetNegative) /
                                               2)
        dimnames(meancorTestSetOverall)[[2]] = names(data.frame(corPrediction))
        meancorTestSetOverall = data.frame(topNumber = topNumber,
                                           meancorTestSetOverall)
        meancorTestSetPositive = data.frame(meancorTestSetPositive)
        dimnames(meancorTestSetPositive)[[2]] = names(data.frame(corPrediction))
        meancorTestSetPositive = data.frame(topNumber = topNumber,
                                            meancorTestSetPositive)
        meancorTestSetNegative = data.frame(meancorTestSetNegative)
        dimnames(meancorTestSetNegative)[[2]] = names(data.frame(corPrediction))
        meancorTestSetNegative = data.frame(topNumber = topNumber,
                                            meancorTestSetNegative)
        datout = list(
            meancorTestSetOverall = meancorTestSetOverall,
            meancorTestSetPositive = meancorTestSetPositive,
            meancorTestSetNegative  = meancorTestSetNegative
        )
        datout
    } # end of function corPredictionSuccess



#-------------------------------------------------------------------------------
#
# relativeCorPredictionSuccess
#
#-------------------------------------------------------------------------------

# The function relativeCorPredictionSuccess can be used to test whether a gene
# screening method is significantly better than a standard method.
# For each gene screening method (column of corPredictionNew) it provides a
# Kruskal Wallis test p-value for comparison with the vector
# corPredictionStandard, TopNumber is a vector of integers.
# corTestSet should be a vector of correlations in the test set.
# corPredictionNew should be a data frame of predictions for the correlations.
# corPredictionStandard should be the standard prediction (correlation in the
# training data).
# The function outputs a p-value for the Kruskal test that the new correlation
# prediction methods outperform the standard correlation prediction method.



#' Compare prediction success
#'
#' Compare prediction success of several gene screening methods.
#'
#'
#' @param corPredictionNew Matrix of predictor statistics
#' @param corPredictionStandard Reference presdictor statistics
#' @param corTestSet Correlations of predictor variables with trait in test set
#' @param topNumber A vector giving the numbers of top genes to consider
#' @return Data frame with components \item{topNumber}{copy of the input
#' \code{topNumber}} \item{kruskalp}{Kruskal-Wallis p-values}
#' @author Steve Horvath
#' @seealso \code{\link{corPredictionSuccess}}
#' @keywords misc
relativeCorPredictionSuccess <- function(corPredictionNew,
                                         corPredictionStandard,
                                         corTestSet,
                                         topNumber = 100) {
    nPredictors = dim(as.matrix(corPredictionNew))[[2]]
    nGenes = dim(as.matrix(corPredictionNew))[[1]]
    if (length(as.numeric(corTestSet)) != nGenes)
        stop("non - compatible dimensions of 'corPrediction' and 'corTestSet'.")
    if (length(as.numeric(corTestSet)) != length(corPredictionStandard))
        stop("non - compatible dimensions of 'corTestSet' and
             'corPredictionStandard'.")
    kruskalp = matrix(nrow = length(topNumber), ncol = nPredictors)
    for (i in c(1:nPredictors)) {
        rankhighNew = rank(-as.matrix(corPredictionNew)[, i],
                           ties.method = "first")
        ranklowNew = rank(as.matrix(corPredictionNew)[, i],
                          ties.method = "first")
        for (j in c(1:length(topNumber))) {
            highCorNew = as.numeric(corTestSet[rankhighNew <= topNumber[j]])
            lowCorNew = as.numeric(corTestSet[ranklowNew <= topNumber[j]])
            highCorStandard = as.numeric(corTestSet[rank(-as.numeric(corPredictionStandard),
                                                         ties.method = "first") <= topNumber[j]])
            lowCorStandard = as.numeric(corTestSet[rank(as.numeric(corPredictionStandard), ties.method = "first")
                                                   <= topNumber[j]])
            signedCorNew = c(highCorNew,-lowCorNew)
            signedCorStandard = c(highCorStandard,-lowCorStandard)
            x1 = c(signedCorNew, signedCorStandard)
            Grouping = rep(c(2, 1), c(length(signedCorNew),
                                      length(signedCorStandard)))
            sign1 = sign(cor(Grouping, x1, use = "p"))
            if (sign1 == 0)
                sign1 = 1
            kruskalp[j, i] = kruskal.test(x = x1, g = Grouping)$p.value * sign1
            #print(names(data.frame(corPredictionNew))[[i]])
            #print(paste("This correlation is positive if the new method is
            #better than the old method" ,
            # signif (cor(Grouping, x1, use = "p"), 3)))
        } # end of j loop
    } # end of i loop
    kruskalp[kruskalp < 0] = 1
    kruskalp = data.frame(kruskalp)
    dimnames(kruskalp)[[2]] = paste0(names(data.frame(corPredictionNew)),
                                     ".kruskalP")
    kruskalp = data.frame(topNumber = topNumber, kruskalp)
    kruskalp
} # end of function relativeCorPredictionSuccess

#-------------------------------------------------------------------------------
#
# alignExpr
#
#-------------------------------------------------------------------------------

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
        y = as.numeric(datExpr[, 1])
    sign1 = sign(as.numeric(cor(y, datExpr, use = "p")))
    as.data.frame(scale(t(t(datExpr) * sign1)))
} # end of function alignExpr


# this function can be used to rank the values in x. Ties are broken by the
# method first.
# This function does not appear to be used anywhere in these functions.
# rank1 <- function(x){
#    rank(x, ties.method = "first")
#}

################################################################################
#
# Gene expression simulations (functions by P.L.)
#
################################################################################

#-------------------------------------------------------------------------------
#
# .causalChildren
#
#-------------------------------------------------------------------------------
# Note: The returned vector may contain multiple occurences of the same child.

.causalChildren <- function(parents, causeMat) {
    nNodes = dim(causeMat)[[1]]

    # print(paste("Length of parents: ", length(parents)))
    if (length(parents) == 0)
        return(NULL)

    Child_ind = apply(as.matrix(abs(causeMat[, parents])), 1, sum) > 0
    if (sum(Child_ind) > 0)
    {
        children = c(1:nNodes)[Child_ind]
    } else {
        children = NULL
    }
    children
}


#-------------------------------------------------------------------------------
#
# simulateEigengeneNetwork
#
#-------------------------------------------------------------------------------
#
# Given a set of causal anchors, this function creates a network of vectors that
# should satisfy the causal relations encoded in the causal matrix causeMat,
# i.e. causeMat[j, i] is the causal effect of vector i on vector j.

# The function starts by initializing all vectors to noise given in the noise
# specification. (The noise can be specified for each vector separately.) Then
# it runs the standard causal network signal propagation and returns the
# resulting vectors.



#' Simulate eigengene network from a causal model
#'
#' Simulates a set of eigengenes (vectors) from a given set of causal anchors
#' and a causal matrix.
#'
#' The algorithm starts with the anchor vectors and iteratively generates the
#' rest from the path coefficients given in the matrix \code{causeMat}.
#'
#' @param causeMat causal matrix. The entry \code{[i,j]} is the influence (path
#' coefficient) of vector \code{j} on vector \code{i}.
#' @param anchorIndex specifies the indices of the anchor vectors.
#' @param anchorVectors a matrix giving the actual anchor vectors as columns.
#' Their number must equal the length of \code{anchorIndex}.
#' @param noise standard deviation of the noise added to each simulated vector.
#' @param verbose level of verbosity. 0 means silent.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation; each unit adds two spaces.
#' @return A list with the following components: \item{eigengenes }{ generated
#' eigengenes. } \item{causeMat }{ a copy of the input causal matrix}
#' \item{levels}{ useful for debugging. A vector with one entry for each
#' eigengene giving the number of generations of parents of the eigengene.
#' Anchors have level 0, their direct causal children have level 1 etc.}
#' \item{anchorIndex}{a copy of the input \code{anchorIndex}. }
#' @author Peter Langfelder
#' @keywords misc
simulateEigengeneNetwork <-
    function(causeMat,
             anchorIndex,
             anchorVectors,
             noise = 1,
             verbose = 0,
             indent = 0) {
        spaces = indentSpaces(indent)

        if (verbose > 0)
            printFlush(paste(spaces, "Creating seed vectors..."))
        nNodes = dim(causeMat)[[1]]
        nSamples = dim(anchorVectors)[[1]]


        if (length(anchorIndex) != dim(anchorVectors)[[2]])
            stop(paste(
                "Length of anchorIndex must equal the number of vectors in
                anchorVectors."
            ))

        if (length(noise) == 1)
            noise = rep(noise, nNodes)
        if (length(noise) != nNodes)
            stop(
                paste(
                    "Length of noise must equal",
                    "the number of nodes as given by the dimension of the
                    causeMat matrix."
                )
            )

        # Initialize all node vectors to noise with given standard deviation

        NodeVectors = matrix(0, nrow = nSamples, ncol = nNodes)
        for (i in 1:nNodes)
            NodeVectors[, i] = rnorm(n = nSamples, mean = 0,
                                     sd = noise[i])

        Levels = rep(0, times = nNodes)

        # Calculate levels for all nodes: start from anchors and go through each
        # successive level of children

        level = 0
        parents = anchorIndex
        Children = .causalChildren(parents = parents, causeMat = causeMat)
        if (verbose > 1)
            printFlush(paste(spaces, "..Determining level structure..."))
        while (!is.null(Children))
        {
            # print(paste("level:", level))
            # print(paste("   parents:", parents))
            # print(paste("   Children:", Children))
            level = level + 1
            if ((verbose > 1) & (level / 10 == as.integer(level / 10)))
                printFlush(paste(spaces, "  ..Detected level", level))
            #printFlush(paste("Detected level", level))
            Levels[Children] = level
            parents = Children
            Children = .causalChildren(parents = parents, causeMat = causeMat)
        }

        HighestLevel = level

        # Generate the whole network

        if (verbose > 1)
            printFlush(paste(spaces, "..Calculating network..."))
        NodeVectors[, anchorIndex] = NodeVectors[, anchorIndex] + anchorVectors
        for (level in (1:HighestLevel))
        {
            if ((verbose > 1) & (level / 10 == as.integer(level / 10)))
                printFlush(paste(spaces, " .Working on level", level))
            #printFlush(paste("Working on level", level))
            LevelChildren = c(1:nNodes)[Levels == level]
            for (child in LevelChildren)
            {
                LevelParents = c(1:nNodes)[causeMat[child,] != 0]
                for (parent in LevelParents)
                    NodeVectors[, child] = scale(NodeVectors[, child] + causeMat[child, parent]  *
                                                     NodeVectors[, parent])
            }
        }

        Nodes = list(
            eigengenes = NodeVectors,
            causeMat = causeMat,
            levels = Levels,
            anchorIndex = anchorIndex
        )
        Nodes
    }

#-------------------------------------------------------------------------------
#
# simulateModule
#
#-------------------------------------------------------------------------------
# The resulting data is normalized.
# Attributes contain the component trueKME giving simulated correlation with
# module eigengene for both module genes and near - module genes.
# corPower controls how fast the correlation drops with index i in the module;
# the curve is roughly x^{1/corPower} with x<1 and x~0 near the "center", so the
#  higher the power, the faster the curve rises.



#' Simulate a gene co-expression module
#'
#' Simulation of a single gene co-expression module.
#'
#' Module genes are simulated around the eigengene by choosing them such that
#' their (expected) correlations with the seed eigengene decrease progressively
#' from (just below) \code{maxCor} to \code{minCor}. The genes are otherwise
#' independent from one another. The variable \code{corPower} determines how
#' fast the correlation drops towards \code{minCor}. Higher powers lead to a
#' faster frop-off; \code{corPower} must be above zero but need not be integer.
#'
#' If \code{signed} is \code{FALSE}, the genes are simulated so as to be part
#' of an unsigned network module, that is some genes will be simulated with a
#' negative correlation with the seed eigengene (but of the same absolute value
#' that a positively correlated gene would be simulated with). The proportion
#' of genes with negative correlation is controlled by \code{propNegativeCor}.
#'
#' Optionally, the function can also simulate genes that are "near" the module,
#' meaning they are simulated with a low but non-zero correlation with the seed
#' eigengene. The correlations run between \code{minCor} and zero.
#'
#' @param ME seed module eigengene.
#' @param nGenes number of genes in the module to be simulated. Must be
#' non-zero.
#' @param nNearGenes number of genes to be simulated with low correlation with
#' the seed eigengene.
#' @param minCor minimum correlation of module genes with the eigengene. See
#' details.
#' @param maxCor maximum correlation of module genes with the eigengene. See
#' details.
#' @param corPower controls the dropoff of gene-eigengene correlation. See
#' details.
#' @param signed logical: should the genes be simulated as belonging to a
#' signed network? If \code{TRUE}, all genes will be simulated to have positive
#' correlation with the eigengene. If \code{FALSE}, a proportion given by
#' \code{propNegativeCor} will be simulated with negative correlations of the
#' same absolute values.
#' @param propNegativeCor proportion of genes to be simulated with negative
#' gene-eigengene correlations. Only effective if \code{signed} is
#' \code{FALSE}.
#' @param geneMeans optional vector of length \code{nGenes} giving desired mean
#' expression for each gene. If not given, the returned expression profiles
#' will have mean zero.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix containing the expression data with rows corresponding to
#' samples and columns to genes.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{simulateEigengeneNetwork}} for a simulation of eigengenes with a
#' given causal structure;
#'
#' \code{\link{simulateDatExpr}} for simulations of whole datasets consisting
#' of multiple modules;
#'
#' \code{\link{simulateDatExpr5Modules}} for a simplified interface to
#' expression simulations;
#'
#' \code{\link{simulateMultiExpr}} for a simulation of several related data
#' sets.
#' @references A short description of the simulation method can also be found
#' in the Supplementary Material to the article
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co-expression modules. BMC Systems Biology 2007, 1:54.
#'
#' The material is posted at
#' http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf.
#' @keywords misc
simulateModule <- function(ME, nGenes, nNearGenes = 0, minCor = 0.3, maxCor = 1,
             corPower = 1, signed = FALSE, propNegativeCor = 0.3,
             geneMeans = NULL, verbose = 0, indent = 0) {
        nSamples = length(ME)

        datExpr = matrix(rnorm((nGenes + nNearGenes) * nSamples),
                         nrow = nSamples,
                         ncol = nGenes + nNearGenes)

        VarME = var(ME)

        # generate the in - module genes
        CorME = maxCor - (c(1:nGenes) / nGenes) ^ (1 / corPower) * (maxCor - minCor)
        noise = sqrt(VarME * (1 - CorME ^ 2) / CorME ^ 2)
        sign = rep(1, nGenes)
        if (!signed)
        {
            negGenes = as.integer(
                seq(
                    from = 1 / propNegativeCor,
                    by = 1 / propNegativeCor,
                    length.out = nGenes * propNegativeCor
                )
            )
            negGenes = negGenes[negGenes <= nGenes]
            sign[negGenes] = -1
        }
        for (gene in 1:nGenes)
        {
            datExpr[, gene] = sign[gene] * (ME + rnorm(nSamples, sd = noise[gene]))
        }

        trueKME = CorME
        # generate the near - module genes

        if (nNearGenes > 0) {
            CorME = c(1:nNearGenes) / nNearGenes * minCor
            noise = sqrt(VarME * (1 - CorME ^ 2) / CorME ^ 2)
            sign = rep(1, nNearGenes)
            if (!signed) {
                negGenes = as.integer(
                    seq(
                        from = 1 / propNegativeCor,
                        by = 1 / propNegativeCor,
                        length.out = nNearGenes * propNegativeCor
                    )
                )
                negGenes = negGenes[negGenes <= nNearGenes]
                sign[negGenes] = -1
            }
            for (gene in 1:nNearGenes)
                datExpr[, nGenes + gene] = ME + sign[gene] * rnorm(nSamples,
                                                                   sd = noise[gene])
            trueKME = c(trueKME, CorME)
        }

        datExpr = scale(datExpr)
        if (!is.null(geneMeans)) {
            if (any(is.na(geneMeans)))
                stop("All entries of 'geneMeans' must be finite.")
            if (length(geneMeans) != nGenes + nNearGenes)
                stop("The lenght of 'geneMeans' must equal nGenes + nNearGenes.")
            datExpr = datExpr + matrix(geneMeans, nSamples, nGenes + nNearGenes,
                                       byrow = TRUE)
        }

        attributes(datExpr)$trueKME = trueKME

        datExpr

    }

#.SimulateModule <- function(ME, size, minimumCor = .3) {
#if (size<3) print("WARNING: module size smaller than 3")
#if (minimumCor == 0) minimumCor = 0.0001
#maxnoisevariance = var(ME, na.rm = TRUE) * (1/minimumCor^2 - 1)
#SDvector = sqrt(c(1:size)/size * maxnoisevariance)
#datSignal = suppressWarnings(matrix(c(ME, ME, - ME), nrow = size,
#ncol = length(ME), byrow = TRUE))
#datNoise = SDvector *  matrix(rnorm(size * length(ME)), nrow = size,
# ncol = length(ME))
#datModule = datSignal + datNoise
#t(datModule)
#} # end of function



#-------------------------------------------------------------------------------
#
# simulateSmallLayer
#
#-------------------------------------------------------------------------------
# Simulates a bunch of small and weakly expressed modules.



#' Simulate small modules
#'
#' This function simulates a set of small modules. The primary purpose is to
#' add a submodule structure to the main module structure simulated by
#' \code{\link{simulateDatExpr}}.
#'
#' Module eigenvectors are chosen randomly and independently. Module sizes are
#' chosen randomly from an exponential distribution with mean equal
#' \code{averageModuleSize}. Two thirds of genes in each module are simulated
#' as proper module genes and one third as near-module genes (see
#' \code{\link{simulateModule}} for details).  Between each successive pairs of
#' modules a number of genes given by \code{moduleSpacing} will be left
#' unsimulated (zero expression). Module expression, that is the expected
#' standard deviation of the module expression vectors, is chosen randomly from
#' an exponential distribution with mean equal \code{averageExpr}. The
#' expression profiles are chosen such that their correlations with the
#' eigengene run from just below \code{maxCor} to \code{minCor * maxCor} (hence
#' minCor must be between 0 and 1, not including the bounds). The parameter
#' \code{corPower} can be chosen to control the behaviour of the simulated
#' correlation with the gene index; values higher than 1 will result in the
#' correlation approaching \code{minCor * maxCor} faster and lower than 1
#' slower.
#'
#' The simulated genes will be returned in the order given in \code{order}.
#'
#' @param order a vector giving the simulation order for vectors. See details.
#' @param nSamples integer giving the number of samples to be simulated.
#' @param minCor a multiple of \code{maxCor} (see below) giving the minimum
#' correlation of module genes with the corresponding eigengene. See details.
#' @param maxCor maximum correlation of module genes with the corresponding
#' eigengene. See details.
#' @param corPower controls the dropoff of gene-eigengene correlation. See
#' details.
#' @param averageModuleSize average number of genes in a module. See details.
#' @param averageExpr average strength of module expression vectors.
#' @param moduleSpacing a number giving module spacing: this multiple of the
#' module size will lie between the module and the next one.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A matrix of simulated gene expressions, with dimension
#' \code{(nSamples, length(order))}.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{simulateModule}} for simulation of individual modules;
#'
#' \code{\link{simulateDatExpr}} for the main gene expression simulation
#' function.
#' @keywords misc
simulateSmallLayer <- function(order,
                               nSamples,
                               minCor = 0.3,
                               maxCor = 0.5,
                               corPower = 1,
                               averageModuleSize,
                               averageExpr,
                               moduleSpacing,
                               verbose = 4,
                               indent = 0) {
    spaces = indentSpaces(indent)
    nGenes = length(order)
    datExpr = matrix(0, nrow = nSamples, ncol = nGenes)

    maxCorN0 = averageModuleSize

    if (verbose > 0)
        printFlush(
            paste(
                spaces,
                "simulateSmallLayer: simulating
                modules with min corr",
                minCor,
                ", average expression",
                averageExpr,
                ", average module size",
                averageModuleSize,
                ", inverse density",
                moduleSpacing
            )
        )

    index = 0
    while (index < nGenes)
    {
        ModSize = as.integer(rexp(1, 1 / averageModuleSize))
        if (ModSize < 3)
            ModSize = 3
        if (index + ModSize > nGenes)
            ModSize = nGenes - index
        if (ModSize > 2)
            # Otherwise don't bother :)
        {
            ModuleExpr = rexp(1, 1 / averageExpr)
            if (verbose > 4)
                printFlush(
                    paste(
                        spaces,
                        "  Module of size",
                        ModSize,
                        ", expression",
                        ModuleExpr,
                        ", min corr",
                        minCor,
                        "inserted at index",
                        index + 1
                    )
                )
            ME = rnorm(nSamples, sd = ModuleExpr)
            NInModule = as.integer(ModSize * 2 / 3)
            nNearModule = ModSize - NInModule
            EffMinCor = minCor * maxCor
            datExpr[, order[(index + 1):(index + ModSize)]]  =
                ModuleExpr * simulateModule(ME,
                                            NInModule,
                                            nNearModule,
                                            EffMinCor,
                                            maxCor,
                                            corPower)
        }
        index = index + ModSize * moduleSpacing
    }
    datExpr
}


#-------------------------------------------------------------------------------
#
# simulateDatExpr
#
#-------------------------------------------------------------------------------
#
# Caution: the last Mod.Props entry gives the number of "true grey" genes
# the corresponding minCor entry must be absent (i.e.
# length(minCor) = length(modProportions) - 1

# SubmoduleLayers: layers of small modules with weaker correlation, ordered in
# the same order as the genes in the big modules. Needs average number of genes
# in a module (exponential distribution), average expression strength (
# exponential density) and inverse density.

# ScatteredModuleLayers: Layers of small modules whose order is random.



#' Simulation of expression data
#'
#' Simulation of expression data with a customizable modular structure and
#' several different types of noise.
#'
#'
#' Given \code{eigengenes} can be unrelated or they can exhibit non-trivial
#' correlations. Each module is simulated separately from others. The
#' expression profiles are chosen such that their correlations with the
#' eigengene run from just below \code{maxCor} to \code{minCor} (hence minCor
#' must be between 0 and 1, not including the bounds). The parameter
#' \code{corPower} can be chosen to control the behaviour of the simulated
#' correlation with the gene index; values higher than 1 will result in the
#' correlation approaching \code{minCor} faster and lower than 1 slower.
#'
#' Numbers of genes in each module are specified (as fractions of the total
#' number of genes \code{nGenes}) by \code{modProportions}. The last entry in
#' \code{modProportions} corresponds to the genes that will be simulated as
#' unrelated to anything else ("grey" genes). The proportion must add up to 1
#' or less. If the sum is less than one, the remaining genes will be
#' partitioned into groups and simulated to be "close" to the proper modules,
#' that is with small but non-zero correlations (between \code{minCor} and 0)
#' with the module eigengene.
#'
#' If \code{signed} is set \code{FALSE}, the correlation for some of the module
#' genes is chosen negative (but the absolute values remain the same as they
#' would be for positively correlated genes). To ensure consistency for
#' simulations of multiple sets, the indices of the negatively correlated genes
#' are fixed and distributed evenly.
#'
#' In addition to the primary module structure, a secondary structure can be
#' optionally simulated. Modules in the secondary structure have sizes chosen
#' from an exponential distribution with mean equal
#' \code{averageNGenesInSubmodule}. Expression vectors simulated in the
#' secondary structure are simulated with expected standard deviation chosen
#' from an exponential distribution with mean equal
#' \code{averageExprInSubmodule}; the higher this coefficient, the more
#' pronounced will the submodules be in the main modules. The secondary
#' structure can be simulated in several layers; their number is given by
#' \code{SubmoduleLayers}. Genes in these submodules are ordered in the same
#' order as in the main modules.
#'
#' In addition to the ordered submodule structure, a scattered submodule
#' structure can be simulated as well. This structure can be viewed as noise
#' that tends to correlate random groups of genes. The size and effect
#' parameters are the same as for the ordered submodules, and the number of
#' layers added is controlled by \code{nScatteredModuleLayers}.
#'
#' @param eigengenes a data frame containing the seed eigengenes for the
#' simulated modules. Rows correspond to samples and columns to modules.
#' @param nGenes total number of genes to be simulated.
#' @param modProportions a numeric vector with length equal the number of
#' eigengenes in \code{eigengenes} plus one, containing fractions of the total
#' number of genes to be put into each of the modules and into the "grey
#' module", which means genes not related to any of the modules. See details.
#' @param minCor minimum correlation of module genes with the corresponding
#' eigengene. See details.
#' @param maxCor maximum correlation of module genes with the corresponding
#' eigengene. See details.
#' @param corPower controls the dropoff of gene-eigengene correlation. See
#' details.
#' @param signed logical: should the genes be simulated as belonging to a
#' signed network? If \code{TRUE}, all genes will be simulated to have positive
#' correlation with the eigengene. If \code{FALSE}, a proportion given by
#' \code{propNegativeCor} will be simulated with negative correlations of the
#' same absolute values.
#' @param propNegativeCor proportion of genes to be simulated with negative
#' gene-eigengene correlations. Only effective if \code{signed} is
#' \code{FALSE}.
#' @param geneMeans optional vector of length \code{nGenes} giving desired mean
#' expression for each gene. If not given, the returned expression profiles
#' will have mean zero.
#' @param backgroundNoise amount of background noise to be added to the
#' simulated expression data.
#' @param leaveOut optional specification of modules that should be left out of
#' the simulation, that is their genes will be simulated as unrelated ("grey").
#' This can be useful when simulating several sets, in some which a module is
#' present while in others it is absent.
#' @param nSubmoduleLayers number of layers of ordered submodules to be added.
#' See details.
#' @param nScatteredModuleLayers number of layers of scattered submodules to be
#' added. See details.
#' @param averageNGenesInSubmodule average number of genes in a submodule. See
#' details.
#' @param averageExprInSubmodule average strength of submodule expression
#' vectors.
#' @param submoduleSpacing a number giving submodule spacing: this multiple of
#' the submodule size will lie between the submodule and the next one.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components: \item{datExpr}{ simulated
#' expression data in a data frame whose columns correspond genes and rows to
#' samples. }
#'
#' \item{setLabels}{ simulated module assignment. Module labels are numeric,
#' starting from 1. Genes simulated to be outside of proper modules have label
#' 0.  Modules that are left out (specified in \code{leaveOut}) are indicated
#' as 0 here. }
#'
#' \item{allLabels}{ simulated module assignment. Genes that belong to leftout
#' modules (specified in \code{leaveOut}) are indicated by their would-be
#' assignment here. }
#'
#' \item{labelOrder}{ a vector specifying the order in which labels correspond
#' to the given eigengenes, that is \code{labelOrder[1]} is the label assigned
#' to module whose seed is \code{eigengenes[, 1]} etc.  }
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{simulateEigengeneNetwork}} for a simulation of eigengenes with a
#' given causal structure;
#'
#' \code{\link{simulateModule}} for simulations of individual modules;
#'
#' \code{\link{simulateDatExpr5Modules}} for a simplified interface to
#' expression simulations;
#'
#' \code{\link{simulateMultiExpr}} for a simulation of several related data
#' sets.
#' @references A short description of the simulation method can also be found
#' in the Supplementary Material to the article
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co-expression modules. BMC Systems Biology 2007, 1:54.
#'
#' The material is posted at
#' http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf.
#' @keywords misc
simulateDatExpr <- function(eigengenes,
                            nGenes,
                            modProportions,
                            minCor = 0.3,
                            maxCor = 1,
                            corPower = 1,
                            signed = FALSE,
                            propNegativeCor = 0.3,
                            geneMeans = NULL,
                            backgroundNoise = 0.1,
                            leaveOut = NULL,
                            nSubmoduleLayers = 0,
                            nScatteredModuleLayers = 0,
                            averageNGenesInSubmodule = 10,
                            averageExprInSubmodule = 0.2,
                            submoduleSpacing = 2,
                            verbose = 1,
                            indent = 0) {
    spaces = indentSpaces(indent)

    nMods = length(modProportions) - 1

    nSamples = dim(eigengenes)[[1]]

    if (length(minCor) == 1)
        minCor = rep(minCor, nMods)
    if (length(maxCor) == 1)
        maxCor = rep(maxCor, nMods)

    if (length(minCor) != nMods)
        stop(
            paste(
                "Input error: minCor is an array of different lentgh than",
                "the length - 1 of modProportions array."
            )
        )

    if (length(maxCor) != nMods)
        stop(
            paste(
                "Input error: maxCor is an array of different lentgh than",
                "the length - 1 of modProportions array."
            )
        )

    if (dim(eigengenes)[[2]] != nMods)
        stop(
            paste(
                "Input error: Number of seed vectors must equal the",
                "length of modProportions."
            )
        )

    if (is.null(geneMeans))
        geneMeans = rep(0, nGenes)
    if (length(geneMeans) != nGenes)
        stop("Length of 'geneMeans' must equal 'nGenes'.")

    if (any(is.na(geneMeans)))
        stop("All entries of 'geneMeans' must be finite.")

    grey = 0
    moduleLabels = c(1:nMods)

    if (sum(modProportions) > 1) {
        stop("Input error: the sum of Mod.Props must be less than 1")
    }
    #if (sum(modProportions[c(1:(length(modProportions) - 1))]) >= 0.5)
    #print(paste("SimulateExprData: Input warning: the sum of modProportions for
    #proper modules",
    #"should ideally be less than 0.5."))

    no.in.modules = as.integer(nGenes * modProportions)
    no.in.proper.modules = no.in.modules[c(1:(length(modProportions) - 1))]
    no.near.modules = as.integer((nGenes - sum(no.in.modules))  *
                                     no.in.proper.modules / sum(no.in.proper.modules))

    simulate.module = rep(TRUE, times = nMods)
    if (!is.null(leaveOut))
        simulate.module[leaveOut] = FALSE

    no.in.modules[nMods + 1] = nGenes - sum(no.in.proper.modules[simulate.module]) -
        sum(no.near.modules[simulate.module])

    labelOrder = moduleLabels[rank(-modProportions[-length(modProportions)],
                                   ties.method = "first")]
    labelOrder = c(labelOrder, grey)

    if (verbose > 0)
        printFlush(
            paste(
                spaces,
                "simulateDatExpr: simulating",
                nGenes,
                "genes in",
                nMods,
                "modules."
            )
        )

    if (verbose > 1) {
        #  printFlush(paste(spaces, "    Minimum correlation in a module is",
        #  minCor,
        #              " and its dropoff is characterized by power", corPower))
        printFlush(paste(
            spaces,
            "    Simulated labels:",
            paste(labelOrder[1:nMods], collapse = ", "),
            " and ",
            grey
        ))
        printFlush(paste(
            spaces,
            "    Module sizes:",
            paste(no.in.modules,
                  collapse = ", ")
        ))
        printFlush(paste(
            spaces,
            "    near module sizes:",
            paste(no.near.modules, collapse = ", ")
        ))
        printFlush(paste(
            spaces,
            "    Min correaltion:",
            paste(minCor, collapse = ", ")
        ))
        if (!is.null(leaveOut))
            printFlush(paste(
                spaces,
                "    _leaving out_ modules",
                paste(labelOrder[leaveOut],
                      collapse = ", ")
            ))
    }

    truemodule = rep(grey, nGenes)
    # These have the colors for left - out modules as well.
    allLabels = rep(grey, nGenes)

    # This matrix contains the simulated expression values (rows are genes,
    # columns samples)
    # Each simulated cluster has a distinct mean expression across the samples

    datExpr = matrix(rnorm(nGenes * nSamples), nrow = nSamples, ncol = nGenes)
    trueKME = rep(NA, nGenes)
    trueKME.whichMod = rep(0, nGenes)

    gene.index = 0 # Where to put the current gene into datExpr

    for (mod in c(1:nMods)) {
        nModGenes = no.in.modules[mod]
        nNearGenes = no.near.modules[mod]
        if (simulate.module[mod]) {
            ME = eigengenes[, mod]
            EffMaxCor = maxCor[mod]
            EffMinCor = minCor[mod]
            range = (gene.index + 1):(gene.index + nModGenes + nNearGenes)
            temp = simulateModule(
                ME,
                nModGenes,
                nNearGenes,
                minCor[mod],
                maxCor[mod],
                corPower,
                signed = signed,
                propNegativeCor = propNegativeCor,
                geneMeans = NULL,
                verbose = verbose - 2,
                indent = indent + 2
            )
            datExpr[, range] = temp
            truemodule[(gene.index + 1):(gene.index + nModGenes)] = labelOrder[mod]
            trueKME[range] = attributes(temp)$trueKME
            trueKME.whichMod[range] = mod
        }
        allLabels[(gene.index + 1):(gene.index + nModGenes)] = labelOrder[mod]
        gene.index = gene.index + nModGenes + nNearGenes
    }

    if (nSubmoduleLayers > 0) {
        OrderVector = c(1:nGenes)
        for (layer in 1:nSubmoduleLayers) {
            if (verbose > 1)
                printFlush(paste(spaces,
                                 "Simulating ordereded extra layer",
                                 layer))
            datExpr = datExpr + simulateSmallLayer(
                OrderVector,
                nSamples,
                minCor[1],
                maxCor[1],
                corPower,
                averageNGenesInSubmodule,
                averageExprInSubmodule,
                submoduleSpacing,
                verbose - 1,
                indent + 1
            )
        }
    }
    if (nScatteredModuleLayers > 0)
        for (layer in 1:nScatteredModuleLayers) {
            if (verbose > 1)
                printFlush(paste(spaces,
                                 "Simulating unordereded extra layer",
                                 layer))
            OrderVector = sample(nGenes)
            datExpr = datExpr + simulateSmallLayer(
                OrderVector,
                nSamples,
                minCor[1],
                maxCor[1],
                corPower,
                averageNGenesInSubmodule,
                averageExprInSubmodule,
                submoduleSpacing,
                verbose = verbose - 1,
                indent = indent + 1
            )
        }
    collectGarbage()
    if (verbose > 1)
        printFlush(paste(
            spaces,
            "  Adding background noise with amplitude",
            backgroundNoise
        ))
    datExpr = datExpr + rnorm(n = nGenes * nSamples, sd = backgroundNoise)
    means = colMeans(datExpr)

    datExpr = datExpr + matrix(geneMeans - means, nSamples, nGenes,
                               byrow = TRUE)

    colnames(datExpr) = paste0("Gene.", c(1:nGenes))
    rownames(datExpr) = paste0("Sample.", c(1:nSamples))

    list(
        datExpr = datExpr,
        setLabels = truemodule,
        allLabels = allLabels,
        labelOrder = labelOrder,
        trueKME = trueKME,
        trueKME.whichMod = trueKME.whichMod
    )
} # end of function

#-------------------------------------------------------------------------------
#
# simulateMultiExpr
#
#-------------------------------------------------------------------------------
# simulate several sets with some of the modules left out.
# eigengenes are specified in a standard multi - set data format.
# leaveOut must be a matrix of No.Modules x No.Sets of TRUE/FALSE values
# minCor must be a single number here modProportions are a single vector,
# since the proportions should be the same for all sets.
# nSamples is a vector specifying the number of samples in each set this must
# be compatible with the dimensions of the eigengenes.



#' Simulate multi-set expression data
#'
#' Simulation of expression data in several sets with relate module structure.
#'
#' For details of simulation of individual data sets and the meaning of
#' individual set simulation arguments, see \code{\link{simulateDatExpr}}. This
#' function simulates several data sets at a time and puts the result in a
#' multi-set format. The number of genes is the same for all data sets. Module
#' memberships are also the same, but modules can optionally be ``dissolved'',
#' that is their genes will be simulated as unassigned. Such ``dissolved'', or
#' left out, modules can be specified in the matrix \code{leaveOut}.
#'
#' @param eigengenes the seed eigengenes for the simulated modules in a
#' multi-set format. A list with one component per set. Each component is again
#' a list that must contain a component \code{data}. This is a data frame of
#' seed eigengenes for the corresponding data set. Columns correspond to
#' modules, rows to samples. Number of samples in the simulated data is
#' determined from the number of samples of the eigengenes.
#' @param nGenes integer specifyin the number of simulated genes.
#' @param modProportions a numeric vector with length equal the number of
#' eigengenes in \code{eigengenes} plus one, containing fractions of the total
#' number of genes to be put into each of the modules and into the "grey
#' module", which means genes not related to any of the modules. See details.
#' @param minCor minimum correlation of module genes with the corresponding
#' eigengene. See details.
#' @param maxCor maximum correlation of module genes with the corresponding
#' eigengene. See details.
#' @param corPower controls the dropoff of gene-eigengene correlation. See
#' details.
#' @param backgroundNoise amount of background noise to be added to the
#' simulated expression data.
#' @param leaveOut optional specification of modules that should be left out of
#' the simulation, that is their genes will be simulated as unrelated ("grey").
#' A logical matrix in which columns correspond to sets and rows to modules.
#' Wherever \code{TRUE}, the corresponding module in the corresponding data set
#' will not be simulated, that is its genes will be simulated independently of
#' the eigengene.
#' @param signed logical: should the genes be simulated as belonging to a
#' signed network? If \code{TRUE}, all genes will be simulated to have positive
#' correlation with the eigengene. If \code{FALSE}, a proportion given by
#' \code{propNegativeCor} will be simulated with negative correlations of the
#' same absolute values.
#' @param propNegativeCor proportion of genes to be simulated with negative
#' gene-eigengene correlations. Only effective if \code{signed} is
#' \code{FALSE}.
#' @param geneMeans optional vector of length \code{nGenes} giving desired mean
#' expression for each gene. If not given, the returned expression profiles
#' will have mean zero.
#' @param nSubmoduleLayers number of layers of ordered submodules to be added.
#' See details.
#' @param nScatteredModuleLayers number of layers of scattered submodules to be
#' added. See details.
#' @param averageNGenesInSubmodule average number of genes in a submodule. See
#' details.
#' @param averageExprInSubmodule average strength of submodule expression
#' vectors.
#' @param submoduleSpacing a number giving submodule spacing: this multiple of
#' the submodule size will lie between the submodule and the next one.
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{multiExpr }{simulated expression data in multi-set format analogous to
#' that of the input \code{eigengenes}.  A list with one component per set.
#' Each component is again a list that must contains a component \code{data}.
#' This is a data frame of expression data for the corresponding data set.
#' Columns correspond to genes, rows to samples.}
#'
#' \item{setLabels}{a matrix of dimensions (number of genes) times (number of
#' sets) that contains module labels for each genes in each simulated data set.
#' }
#'
#' \item{allLabels}{a matrix of dimensions (number of genes) times (number of
#' sets) that contains the module labels that would be simulated if no module
#' were left out using \code{leaveOut}. This means that all columns of the
#' matrix are equal; the columns are repeated for convenience so
#' \code{allLabels} has the same dimensions as \code{setLabels}. }
#'
#' \item{labelOrder}{a matrix of dimensions (number of modules) times (number
#' of sets) that contains the order in which module labels were assigned to
#' genes in each set. The first label is assigned to genes 1...(module size of
#' module labeled by first label), the second label to the following batch of
#' genes etc.}
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{simulateEigengeneNetwork}} for a simulation of eigengenes with a
#' given causal structure;
#'
#' \code{\link{simulateDatExpr}} for simulation of individual data sets;
#'
#' \code{\link{simulateDatExpr5Modules}} for a simple simulation of a data set
#' consisting of 5 modules;
#'
#' \code{\link{simulateModule}} for simulations of individual modules;
#' @references A short description of the simulation method can also be found
#' in the Supplementary Material to the article
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co-expression modules. BMC Systems Biology 2007, 1:54.
#'
#' The material is posted at
#' http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf.
#' @keywords misc
simulateMultiExpr <- function(eigengenes,
                              nGenes,
                              modProportions,
                              minCor = 0.5,
                              maxCor = 1,
                              corPower = 1,
                              backgroundNoise = 0.1,
                              leaveOut = NULL,
                              signed = FALSE,
                              propNegativeCor = 0.3,
                              geneMeans = NULL,
                              nSubmoduleLayers = 0,
                              nScatteredModuleLayers = 0,
                              averageNGenesInSubmodule = 10,
                              averageExprInSubmodule = 0.2,
                              submoduleSpacing = 2,
                              verbose = 1,
                              indent = 0) {
    MEsize = checkSets(eigengenes)
    nSets = MEsize$nSets
    nMods = MEsize$nGenes
    nSamples = MEsize$nSamples

    nAllSamples = sum(nSamples)

    if (is.null(geneMeans)) {
        geneMeans = matrix(0, nGenes, nSets)
    } else {
        geneMeans = as.matrix(geneMeans)
        if (nrow(geneMeans) != nGenes) {
            stop("Number of rows (or entries) in 'geneMeans' must equal
                 'nGenes'.")
        } else if (ncol(geneMeans) == 1) {
            geneMeans = matrix(geneMeans, nGenes, nSets)
        } else if (ncol(geneMeans) != nSets)
            stop("Number of columns in geneMeans must either equal the number
                 of sets or be 1.")
    }

    if (any(is.na(geneMeans)))
        stop("All entries of 'geneMeans' must be finite.")

    d2 = length(modProportions) - 1
    if (d2  != nMods)
        stop(paste(
            "Incompatible numbers of modules in
            'eigengenes' and 'modProportions'"
        ))
    if (is.null(leaveOut)) {
        leaveOut = matrix(FALSE, nMods, nSets)
    } else {
        d3 = dim(leaveOut)
        if ((d3[1]  != nMods) | (d3[2]  != nSets))
            stop(paste(
                "Incompatible dimensions of 'leaveOut' and set
                eigengenes."
            ))
    }

    multiExpr = vector(mode = "list", length = nSets)
    setLabels = NULL
    allLabels = NULL
    labelOrder = NULL

    for (set in 1:nSets) {
        SetEigengenes = scale(eigengenes[[set]]$data)
        setLeaveOut = leaveOut[, set]
        # Convert setLeaveOut from boolean to a list of indices where it's TRUE
        # SetMinCor = rep(minCor, nMods)
        # SetMaxCor = rep(maxCor, nMods)
        SetLO = c(1:nMods)[setLeaveOut]
        setData = simulateDatExpr(
            SetEigengenes,
            nGenes,
            modProportions,
            minCor = minCor,
            maxCor = maxCor,
            corPower = corPower,
            signed = signed,
            propNegativeCor = propNegativeCor,
            backgroundNoise = backgroundNoise,
            leaveOut = SetLO,
            nSubmoduleLayers = nSubmoduleLayers,
            nScatteredModuleLayers  = nScatteredModuleLayers ,
            averageNGenesInSubmodule = averageNGenesInSubmodule,
            averageExprInSubmodule = averageExprInSubmodule,
            submoduleSpacing = submoduleSpacing,
            verbose = verbose - 1,
            indent = indent + 1
        )
        multiExpr[[set]] = list(data = setData$datExpr)
        setLabels = cbind(setLabels, setData$setLabels)
        allLabels = cbind(allLabels, setData$allLabels)
        labelOrder = cbind(labelOrder, setData$labelOrder)
    }
    list(
        multiExpr = multiExpr,
        setLabels = setLabels,
        allLabels = allLabels,
        labelOrder = labelOrder
    )
}

#-------------------------------------------------------------------------------
#
# simulateDatExpr5Modules
#
#-------------------------------------------------------------------------------



#' Simplified simulation of expression data
#'
#' This function provides a simplified interface to the expression data
#' simulation, at the cost of considerably less flexibility.
#'
#' Roughly one-third of the genes are simulated with a negative correlation to
#' their seed eigengene. See the functions \code{\link{simulateModule}} and
#' \code{\link{simulateDatExpr}} for more details.
#'
#' @param nGenes total number of genes to be simulated.
#' @param colorLabels labels for simulated modules.
#' @param simulateProportions a vector of length 5 giving proportions of the
#' total number of genes to be placed in each individual module. The entries
#' must be positive and sum to at most 1. If the sum is less than 1, the
#' leftover genes will be simulated outside of modules.
#' @param MEturquoise seed module eigengene for the first module.
#' @param MEblue seed module eigengene for the second module.
#' @param MEbrown seed module eigengene for the third module.
#' @param MEyellow seed module eigengene for the fourth module.
#' @param MEgreen seed module eigengene for the fifth module.
#' @param SDnoise level of noise to be added to the simulated expressions.
#' @param backgroundCor backgrond correlation. If non-zero, a component will be
#' added to all genes such that the average correlation of otherwise unrelated
#' genes will be \code{backgroundCor}.
#' @return
#'
#' A list with the following components: \item{datExpr}{ the simulated
#' expression data in a data frame, with rows corresponding to samples and
#' columns to genes. }
#'
#' \item{truemodule}{ a vector with one entry per gene containing the simulated
#' module membership. }
#'
#' \item{datME}{a data frame containing a copy of the input module eigengenes.
#' }
#' @author Steve Horvath and Peter Langfelder
#' @seealso
#'
#' \code{\link{simulateModule}} for simulation of individual modules;
#'
#' \code{\link{simulateDatExpr}} for a more comprehensive data simulation
#' interface.
#' @keywords misc
simulateDatExpr5Modules <- function(nGenes = 2000,
                                    colorLabels = c("turquoise", "blue", "brown", "yellow", "green"),
                                    simulateProportions = c(0.10, 0.08, 0.06, 0.04, 0.02),
                                    MEturquoise,
                                    MEblue,
                                    MEbrown,
                                    MEyellow,
                                    MEgreen,
                                    SDnoise = 1,
                                    backgroundCor = 0.3) {
    nSamples = length(MEturquoise)
    if (length(MEturquoise)  != length(MEblue) |
        length(MEturquoise)  != length(MEbrown) |
        length(MEturquoise)  != length(MEyellow) |
        length(MEturquoise)  != length(MEgreen))
        stop("Numbers of samples in module eigengenes (MEs) are not consistent")
    if (sum(simulateProportions) > 1) {
        stop(
            "Sum of module proportions is larger than 1. Please ensure
            sum(simulateProportions) <= 1."
        )
        # simulateProportions = rep(1/10, 5)
    }
    modulesizes = round(nGenes * c(simulateProportions, 1 - sum(simulateProportions)))
    truemodule = rep(c(as.character(colorLabels), "grey"), modulesizes)
    ModuleEigengenes = data.frame(MEturquoise, MEblue, MEbrown, MEyellow, MEgreen)
    no.MEs = dim(ModuleEigengenes)[[2]]
    # This matrix contains the simulated expression values
    #(rows are samples, columns genes)
    # it contains some background noise
    datExpr = matrix(rnorm(nSamples * nGenes, mean = 0, sd = SDnoise),
                     nrow = nSamples,
                     ncol = nGenes)

    if (is.logical(backgroundCor))
        backgroundCor = 0.3 * as.numeric(backgroundCor)

    if (as.numeric(backgroundCor) > 0) {
        MEbackground = MEturquoise
        datSignal = matrix(
            MEbackground,
            nrow = length(MEturquoise),
            ncol = nGenes,
            byrow = FALSE
        )
        datExpr = datExpr + as.numeric(backgroundCor) * datSignal
    }# end of if backgroundCor

    for (i in c(1:no.MEs)) {
        restrict1 = truemodule == colorLabels[i]
        datModule = simulateModule(ModuleEigengenes[, i],
                                   nGenes = modulesizes[i],
                                   corPower = 2.5)
        datExpr[, restrict1] = datModule
    } # end of for loop
    # this is the output of the function
    list(datExpr  = datExpr,
         truemodule  = truemodule,
         datME = ModuleEigengenes)
} # end of simulation function


#-------------------------------------------------------------------------------
#
# automaticNetworkScreening
#
#-------------------------------------------------------------------------------




#' One-step automatic network gene screening
#'
#' This function performs gene screening based on a given trait and gene
#' network properties
#'
#' Network screening is a method for identifying genes that have a high gene
#' significance and are members of important modules at the same time.  If
#' \code{datME} is given, the function calls \code{\link{networkScreening}}
#' with the default parameters. If \code{datME} is not given, module eigengenes
#' are first calculated using network analysis based on supplied parameters.
#'
#' @param datExpr data frame containing the expression data, columns
#' corresponding to genes and rows to samples
#' @param y vector containing trait values for all samples in \code{datExpr}
#' @param power soft thresholding power used in network construction
#' @param networkType character string specifying network type. Allowed values
#' are (unique abbreviations of) \code{"unsigned"}, \code{"signed"},
#' \code{"hybrid"}.
#' @param detectCutHeight cut height of the gene hierarchical clustering
#' dendrogram. See \code{cutreeDynamic} for details.
#' @param minModuleSize minimum module size to be used in module detection
#' procedure.
#' @param datME optional specification of module eigengenes. A data frame whose
#' columns are the module eigengenes. If given, module analysis will not be
#' performed.
#' @param getQValues logical: should q-values (local FDR) be calculated?
#' @param ... other arguments to the module identification function
#' \code{\link{blockwiseModules}}
#' @return A list with the following components: \item{networkScreening}{a data
#' frame containing results of the network screening procedure. See
#' \code{\link{networkScreening}} for more details.} \item{datME}{ calculated
#' module eigengenes (or a copy of the input \code{datME}, if given).}
#' \item{hubGeneSignificance}{ hub gene significance for all calculated
#' modules. See \code{\link{hubGeneSignificance}}. }
#' @author Steve Horvath
#' @seealso \code{\link{networkScreening}}, \code{\link{hubGeneSignificance}},
#' \code{\link{networkScreening}}, \code{\link[dynamicTreeCut]{cutreeDynamic}}
#' @keywords misc
automaticNetworkScreening <- function(datExpr,
                                      y,
                                      power = 6,
                                      networkType = "unsigned",
                                      detectCutHeight = 0.995,
                                      minModuleSize = min(20, ncol(as.matrix(datExpr)) / 2),
                                      datME = NULL,
                                      getQValues = TRUE,
                                      ...) {
    y = as.numeric(as.character(y))
    if (length(y)  != dim(as.matrix(datExpr))[[1]])
        stop(
            "Number of samples in 'y' and 'datExpr' disagree: length(y)  !=
            dim(as.matrix(datExpr))[[1]] "
        )

    nAvailable = apply(as.matrix(!is.na(datExpr)), 2, sum)
    ExprVariance = apply(as.matrix(datExpr), 2, var, na.rm = TRUE)
    restrictGenes = (nAvailable >= ..minNSamples) &
        (ExprVariance > 0)
    numberUsefulGenes = sum(restrictGenes, na.rm = TRUE)
    if (numberUsefulGenes < 3) {
        stop(
            paste(
                "IMPORTANT: there are not enough useful genes. \n",
                "    Your input genes have fewer than 4 observations or they
                are constant.\n",
                "    WGCNA cannot be used for these data. Hint: collect more
                arrays or input genes that vary."
            )
        )
        #warning(paste("IMPORTANT: there are not enough useful genes. \n",
        #   "    Your input genes have fewer than 4 observations or they are
        #   constant.\n",
        #   "    WGCNA cannot be used for these data. Hint: collect more arrays
        #   or input genes that vary."))
        #output = list(NetworkScreening = data.frame(NS1 = rep(NA,
        # dim(as.matrix(datExpr))[[2]])),
        #            datME = rep(NA, dim(as.matrix(datExpr))[[1]]),
        #            EigengeneSignificance = NA, AAcriterion = NA)
        #return(output)
    }

    datExprUsefulGenes = as.matrix(datExpr)[, restrictGenes &
                                                !is.na(restrictGenes)]
    if (is.null(datME)) {
        mergeCutHeight1 = dynamicMergeCut(n = dim(as.matrix(datExprUsefulGenes))[[1]])
        B = blockwiseModules(
            datExprUsefulGenes,
            mergeCutHeight = mergeCutHeight1,
            TOMType = "none",
            power = power,
            networkType = networkType,
            detectCutHeight = detectCutHeight,
            minModuleSize = minModuleSize
        )
        datME = data.frame(B$MEs)
    }

    if (dim(as.matrix(datME))[[1]]  != dim(as.matrix(datExpr))[[1]])
        stop(
            paste(
                "Numbers of samples in 'datME' and 'datExpr' are incompatible:",
                "dim(as.matrix(datME))[[1]]  != dim(as.matrix(datExpr))[[1]]"
            )
        )

    MMdata = signedKME(datExpr = datExpr,
                       datME = datME,
                       outputColumnName = "MM.")
    MMdataPvalue = as.matrix(corPvalueStudent(as.matrix(MMdata), nSamples = dim(as.matrix(datExpr))[[1]]))
    dimnames(MMdataPvalue)[[2]] = paste("Pvalue", names(MMdata), sep = ".")

    NS1 = networkScreening(
        y = y,
        datME = datME,
        datExpr = datExpr,
        getQValues = getQValues
    )
    # here we compute the eigengene significance measures
    ES = data.frame(cor(y, datME, use = "p"))

    ESvector = as.vector(as.matrix(ES))
    EScounts = tapply(abs(ESvector), cut(abs(ESvector),
                                         seq(
                                             from = 0, to = 1, by = .1
                                         )), length)
    EScounts[is.na(EScounts)] = 0

    rr = max(abs(ES), na.rm = TRUE)
    AAcriterion = sqrt(length(y) - 2) * rr / sqrt(1 - rr ^ 2)


    ESy = (1 + max(abs(ES), na.rm = TRUE)) / 2
    ES = data.frame(ES, ESy = ESy)

    # to avoid dividing by zero, we set correlation that are 1 equal to .9999
    ES.999 = as.numeric(as.vector(ES))
    ES.999[!is.na(ES) &  ES > 0.9999] = .9999
    ES.pvalue = corPvalueStudent(cor = abs(ES.999), nSamples = sum(!is.na(y)))
    ES.pvalue[length(ES.999)] = 0
    EigengeneSignificance.pvalue = data.frame(matrix(ES.pvalue, nrow = 1))
    names(EigengeneSignificance.pvalue) = names(ES)

    datME = data.frame(datME, y = y)
    names(ES) = paste0("ES", substr(names(ES), 3, 100))

    print(signif (ES, 2))

    output = list(
        networkScreening = data.frame(NS1, MMdata, MMdataPvalue),
        datME = data.frame(datME),
        eigengeneSignificance = data.frame(ES) ,
        EScounts = EScounts,
        eigengeneSignificance.pvalue = EigengeneSignificance.pvalue,
        AAcriterion = AAcriterion
    )

    output
} # end of function automaticNetworkScreening


#-------------------------------------------------------------------------------
#
# automaticNetworkScreeningGS
#
#-------------------------------------------------------------------------------



#' One-step automatic network gene screening with external gene significance
#'
#' This function performs gene screening based on external gene significance
#' and their network properties.
#'
#' Network screening is a method for identifying genes that have a high gene
#' significance and are members of important modules at the same time.  If
#' \code{datME} is given, the function calls \code{\link{networkScreeningGS}}
#' with the default parameters. If \code{datME} is not given, module eigengenes
#' are first calculated using network analysis based on supplied parameters.
#'
#' @param datExpr data frame containing the expression data, columns
#' corresponding to genes and rows to samples
#' @param GS vector containing gene significance for all genes given in
#' \code{datExpr}
#' @param power soft thresholding power used in network construction
#' @param networkType character string specifying network type. Allowed values
#' are (unique abbreviations of) \code{"unsigned"}, \code{"signed"},
#' \code{"hybrid"}.
#' @param detectCutHeight cut height of the gene hierarchical clustering
#' dendrogram. See \code{cutreeDynamic} for details.
#' @param minModuleSize minimum module size to be used in module detection
#' procedure.
#' @param datME optional specification of module eigengenes. A data frame whose
#' columns are the module eigengenes. If given, module analysis will not be
#' performed.
#' @return A list with the following components: \item{networkScreening}{a data
#' frame containing results of the network screening procedure. See
#' \code{\link{networkScreeningGS}} for more details.} \item{datME}{ calculated
#' module eigengenes (or a copy of the input \code{datME}, if given).}
#' \item{hubGeneSignificance}{ hub gene significance for all calculated
#' modules. See \code{\link{hubGeneSignificance}}. }
#' @author Steve Horvath
#' @seealso \code{\link{networkScreening}}, \code{\link{hubGeneSignificance}},
#' \code{\link{networkScreening}}, \code{\link[dynamicTreeCut]{cutreeDynamic}}
#' @keywords misc
automaticNetworkScreeningGS <- function(datExpr, GS, power = 6,
                                        networkType = "unsigned",
                                        detectCutHeight = 0.995,
                                        minModuleSize = min(20, ncol(as.matrix(datExpr)) / 2),
                                        datME = NULL) {
    if (!is.numeric(GS))
        stop("Gene significance 'GS' is not numeric.")
    if (dim(as.matrix(datExpr))[[2]]  != length(GS))
        stop(
            "length of gene significance variable GS does not equal the number
            of columns of datExpr."
        )

    mergeCutHeight1 = dynamicMergeCut(n = dim(as.matrix(datExpr))[[1]])
    nAvailable = apply(as.matrix(!is.na(datExpr)), 2, sum)
    ExprVariance = apply(as.matrix(datExpr), 2, var, na.rm = TRUE)
    restrictGenes = nAvailable >= 4 & ExprVariance > 0
    numberUsefulGenes = sum(restrictGenes, na.rm = TRUE)
    if (numberUsefulGenes < 3) {
        stop(
            paste(
                "IMPORTANT: there are not enough useful genes. \n",
                "    Your input genes have fewer than 4 observations or they
                are constant.\n",
                "    WGCNA cannot be used for these data. Hint: collect more
                arrays or input genes that vary."
            )
        )
        #output = list(NetworkScreening = data.frame(NS1 = rep(NA,
        #dim(as.matrix(datExpr))[[2]])) , datME = rep(NA,
        #dim(as.matrix(datExpr))[[1]])   , hubGeneSignificance = NA)
    } # end of if
    datExprUsefulGenes = as.matrix(datExpr)[, restrictGenes &
                                                !is.na(restrictGenes)]

    if (is.null(datME)) {
        B = blockwiseModules(
            datExprUsefulGenes,
            mergeCutHeight = mergeCutHeight1,
            TOMType = "none",
            power = power,
            networkType = networkType,
            detectCutHeight = detectCutHeight,
            minModuleSize = minModuleSize
        )
        datME = data.frame(B$MEs)
    } #end of if
    MMdata = signedKME(datExpr = datExpr,
                       datME = datME,
                       outputColumnName = "MM.")
    MMdataPvalue = as.matrix(corPvalueStudent(as.matrix(MMdata), nSamples = dim(as.matrix(datExpr))[[1]]))
    dimnames(MMdataPvalue)[[2]] = paste("Pvalue", names(MMdata), sep = ".")

    NS1 = networkScreeningGS(datExpr = datExpr,
                             datME = datME,
                             GS = GS)
    # here we compute the eigengene significance measures
    HGS1 = data.frame(as.matrix(t(hubGeneSignificance(MMdata ^ 3, GS ^ 3)), nrow = 1))
    datME = data.frame(datME)
    names(HGS1) = paste0("HGS", substr(names(MMdata), 4, 100))
    # now we compute the AA criterion
    print(signif (HGS1, 2))
    output = list(
        networkScreening = data.frame(NS1, MMdata, MMdataPvalue),
        datME = data.frame(datME),

        hubGeneSignificance = data.frame(HGS1)
    )
    output
} # end of function automaticNetworkScreeningGS


#-------------------------------------------------------------------------------
#
#  hubGeneSignificance
#
#-------------------------------------------------------------------------------

# The following function computes the hub gene significance as defined in
# in the paper Horvath and Dong. Input a data frame with possibly signed
# module membership measures (also known as module eigengene based connectivity
#kME. Further it requires a possibly signed gene significance measure.
# GS = 0 means that the gene is not significant, high positive or negative values
#  mean that it is significant.
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


#-------------------------------------------------------------------------------
#
#  networkScreeningGS
#
#-------------------------------------------------------------------------------



#' Network gene screening with an external gene significance measure
#'
#' This function blends standard and network approaches to selecting genes (or
#' variables in general) with high gene significance
#'
#' This function should be considered experimental. It takes into account both
#' the "standard" and the network measures of gene importance for the trait.
#'
#' @param datExpr data frame of expression data
#' @param datME data frame of module eigengenes
#' @param GS numeric vector of gene significances
#' @param oddPower odd integer used as a power to raise module memberships and
#' significances
#' @param blockSize block size to use for calculations with large data sets
#' @param minimumSampleSize minimum acceptable number of samples. Defaults to
#' the default minimum number of samples used throughout the WGCNA package,
#' currently 4.
#' @param addGS logical: should gene significances be added to the screening
#' statistics?
#' @return \item{GS.Weighted }{weighted gene significance } \item{GS }{copy of
#' the input gene significances (only if \code{addGS=TRUE})}
#' @author Steve Horvath
#' @seealso \code{\link{networkScreening}},
#' \code{\link{automaticNetworkScreeningGS}}
#' @keywords misc
networkScreeningGS <- function(datExpr,
                               datME,
                               GS ,
                               oddPower = 3,
                               blockSize = 1000,
                               minimumSampleSize = ..minNSamples,
                               addGS = TRUE) {
    oddPower = as.integer(oddPower)
    if (as.integer(oddPower / 2) == oddPower / 2) {
        oddPower = oddPower + 1
    }
    nMEs = dim(as.matrix(datME))[[2]]
    nGenes = dim(as.matrix(datExpr))[[2]]
    GS.Weighted = rep(0, nGenes)

    if (dim(as.matrix(datExpr))[[1]]  != dim(as.matrix(datME))[[1]])
        stop(
            paste(
                "Expression data and the module eigengenes have different\n",
                "      numbers of observations (arrays). Specifically:\n",
                "      dim(as.matrix(datExpr))[[1]]  != dim(as.matrix(datME))[[1]] "
            )
        )

    if (dim(as.matrix(datExpr))[[2]]  != length(GS))
        stop(
            paste(
                "The number of genes in the expression data does not match\n",
                "      the length of the genes significance variable.
                Specifically:\n",
                "       dim(as.matrix(datExpr))[[2]]  != length(GS)   "
            )
        )

    nAvailable = apply(as.matrix(!is.na(datExpr)), 2, sum)
    ExprVariance = apply(as.matrix(datExpr), 2, var, na.rm = TRUE)
    restrictGenes = nAvailable >= 4 & ExprVariance > 0
    numberUsefulGenes = sum(restrictGenes, na.rm = TRUE)
    if (numberUsefulGenes < 3)
    {
        stop(
            paste(
                "IMPORTANT: there are fewer than 3 useful genes. \n",
                "    Violations: either fewer than 4 observations or they are
                constant.\n",
                "    WGCNA cannot be used for these data. Hint: collect more arrays
                or input genes that vary."
            )
        )
        # datout = data.frame(GS.Weighted = rep(NA,
        # dim(as.matrix(datExpr))[[2]]), GS = GS)
    } # end of if

    nBlocks = as.integer(nMEs / blockSize)
    if (nBlocks > 0)
        for (i in 1:nBlocks)
        {
            printFlush(paste("block number = ", i))
            index1 = c(1:blockSize) + (i - 1) *  blockSize
            datMEBatch = datME[, index1]
            datKMEBatch = as.matrix(signedKME(datExpr, datMEBatch,
                                              outputColumnName = "MM."))
            ESBatch = hubGeneSignificance(datKMEBatch ^ oddPower, GS ^ oddPower)
            # the following omits the diagonal when datME = datExpr
            if (nGenes == nMEs) {
                diag(datKMEBatch[index1,]) = 0
                # missing values will not be used
                datKMEBatch[is.na(datKMEBatch)] = 0
                ESBatch[is.na(ESBatch)] = 0
            } # end of if
            GS.WeightedBatch = as.matrix(datKMEBatch) ^ oddPower %*% as.matrix(ESBatch)
            GS.Weighted = GS.Weighted + GS.WeightedBatch
        } # end of for (i in 1:nBlocks
    if (nMEs - nBlocks * blockSize > 0) {
        restindex = c((nBlocks * blockSize + 1):nMEs)
        datMEBatch = datME[, restindex]
        datKMEBatch = as.matrix(signedKME(datExpr, datMEBatch,
                                          outputColumnName = "MM."))
        ESBatch = hubGeneSignificance(datKMEBatch ^ oddPower, GS ^ oddPower)
        # the following omits the diagonal when datME = datExpr
        if (nGenes == nMEs) {
            diag(datKMEBatch[restindex,]) = 0
            # missing values will not be used
            datKMEBatch[is.na(datKMEBatch)] = 0
            ESBatch[is.na(ESBatch)] = 0
        } # end of if (nGenes == nMEs)
        GS.WeightedBatch = as.matrix(datKMEBatch) ^ oddPower %*% ESBatch
        GS.Weighted = GS.Weighted + GS.WeightedBatch
    } # end of if (nMEs - nBlocks * blockSize > 0)
    GS.Weighted = GS.Weighted / nMEs
    GS.Weighted[nAvailable < minimumSampleSize] = NA

    rankGS.Weighted = rank(-GS.Weighted, ties.method = "first")
    rankGS = rank(-GS, ties.method = "first")
    printFlush(paste("Proportion of agreement between GS.Weighted and GS:"))
    for (i in c(10, 20, 50, 100, 200, 500, 1000)) {
        printFlush(paste(
            "Top ",
            i,
            " list of genes: prop. of agreement = ",
            signif (sum(
                rankGS.Weighted <= i & rankGS <= i,
                na.rm = TRUE
            ) / i, 3)
        ))
    } # end of for loop
    if (mean(abs(GS.Weighted), na.rm = TRUE) > 0)
    {
        GS.Weighted = GS.Weighted / mean(abs(GS.Weighted), na.rm = TRUE) *
            mean(abs(GS), na.rm = TRUE)
    }
    if (addGS)
        GS.Weighted = apply(data.frame(GS.Weighted, GS), 1,
                            mean, na.rm = TRUE)
    datout = data.frame(GS.Weighted, GS)

    datout
} # end of function

#-------------------------------------------------------------------------------
#
# networkScreening
#
#-------------------------------------------------------------------------------



#' Identification of genes related to a trait
#'
#' This function blends standard and network approaches to selecting genes (or
#' variables in general) highly related to a given trait.
#'
#' This function should be considered experimental. It takes into account both
#' the "standard" and the network measures of gene importance for the trait.
#'
#' @param y clinical trait given as a numeric vector (one value per sample)
#' @param datME data frame of module eigengenes
#' @param datExpr data frame of expression data
#' @param corFnc character string specifying the function to be used to
#' calculate co-expression similarity. Defaults to Pearson correlation. Any
#' function returning values between -1 and 1 can be used.
#' @param corOptions character string specifying additional arguments to be
#' passed to the function given by \code{corFnc}. Use \code{"use = 'p', method
#' = 'spearman'"} to obtain Spearman correlation.
#' @param oddPower odd integer used as a power to raise module memberships and
#' significances
#' @param blockSize block size to use for calculations with large data sets
#' @param minimumSampleSize minimum acceptable number of samples. Defaults to
#' the default minimum number of samples used throughout the WGCNA package,
#' currently 4.
#' @param addMEy logical: should the trait be used as an additional "module
#' eigengene"?
#' @param removeDiag logical: remove the diagonal?
#' @param weightESy weight to use for the trait as an additional eigengene;
#' should be between 0 and 1
#' @param getQValues logical: should q-values be calculated?
#' @return datout = data.frame(p.Weighted, q.Weighted, Cor.Weighted,
#' Z.Weighted, p.Standard, q.Standard, Cor.Standard, Z.Standard) Data frame
#' reporting the following quantities for each given gene:
#'
#' \item{p.Weighted }{weighted p-value of association with the trait}
#'
#' \item{q.Weighted }{q-value (local FDR) calculated from \code{p.Weighted}}
#'
#' \item{cor.Weighted}{correlation of trait with gene expression weighted by a
#' network term}
#'
#' \item{Z.Weighted}{ Fisher Z score of the weighted correlation}
#'
#' \item{p.Standard}{ standard Student p-value of association of the gene with
#' the trait}
#'
#' \item{q.Standard}{ q-value (local FDR) calculated from \code{p.Standard}}
#'
#' \item{cor.Standard}{ correlation of gene with the trait}
#'
#' \item{Z.Standard}{ Fisher Z score of the standard correlation}
#' @author Steve Horvath
#' @keywords misc
networkScreening <- function(y,
                             datME,
                             datExpr,
                             corFnc = "cor",
                             corOptions = "use = 'p'",
                             oddPower = 3,
                             blockSize = 1000,
                             minimumSampleSize = ..minNSamples,
                             addMEy = TRUE,
                             removeDiag = FALSE,
                             weightESy = 0.5,
                             getQValues = TRUE)
{
    oddPower = as.integer(oddPower)
    if (as.integer(oddPower / 2) == oddPower / 2) {
        oddPower = oddPower + 1
    }
    nMEs = dim(as.matrix(datME))[[2]]
    nGenes = dim(as.matrix(datExpr))[[2]]
    # Here we add y as extra ME
    if (nGenes > nMEs & addMEy) {
        datME = data.frame(y, datME)
    }
    nMEs = dim(as.matrix(datME))[[2]]
    RawCor.Weighted = rep(0, nGenes)
    #Cor.Standard = as.numeric(cor(y, datExpr, use = "p"))
    corExpr = parse(text = paste(
        "as.numeric(",
        corFnc,
        "(y, datExpr ",
        prepComma(corOptions),
        "))"
    ))
    Cor.Standard = eval(corExpr)

    NoAvailable = apply(!is.na(datExpr), 2, sum)
    Cor.Standard[NoAvailable < minimumSampleSize] = NA
    if (nGenes == 1) {
        #RawCor.Weighted = as.numeric(cor(y, datExpr, use = "p"))
        corExpr = parse(text = paste(
            "as.numeric(",
            corFnc,
            "(y, datExpr ",
            prepComma(corOptions),
            "))"
        ))
        RawCor.Weighted = eval(corExpr)
    }
    start = 1
    i = 1
    while (start <= nMEs) {
        end = min(start + blockSize - 1, nMEs)
        if (i > 1 ||
            end < nMEs)
            printFlush(paste("block number = ", i))
        index1 = c(start:end)
        datMEBatch = datME[, index1]
        datKMEBatch = as.matrix(
            signedKME(
                datExpr,
                datMEBatch,
                outputColumnName = "MM.",
                corFnc = corFnc,
                corOptions = corOptions
            )
        )
        # ES.CorBatch = as.vector(cor( as.numeric(as.character(y)),
        # datMEBatch, use = "p"))
        corExpr = parse(
            text = paste(
                "as.vector(",
                corFnc,
                "( as.numeric(as.character(y)) , datMEBatch",
                prepComma(corOptions),
                "))"
            )
        )
        ES.CorBatch = eval(corExpr)

        #weightESy
        ES.CorBatch[ES.CorBatch > .999] = weightESy * 1 + (1 -  weightESy) *
            max(abs(ES.CorBatch[ES.CorBatch < .999]), na.rm = TRUE)
        # the following omits the diagonal when datME = datExpr
        if (nGenes == nMEs &
            removeDiag) {
            diag(datKMEBatch[index1,]) = 0
        }
        if (nGenes == nMEs) {
            # missing values will not be used
            datKMEBatch[is.na(datKMEBatch)] = 0
            ES.CorBatch[is.na(ES.CorBatch)] = 0
        } # end of if
        RawCor.WeightedBatch = as.matrix(datKMEBatch) ^ oddPower %*%  as.matrix(ES.CorBatch ^
                                                                                    oddPower)
        RawCor.Weighted = RawCor.Weighted + RawCor.WeightedBatch
        start = end + 1
    } # end of while (start <= nMEs)
    RawCor.Weighted = RawCor.Weighted / nMEs
    RawCor.Weighted[NoAvailable < minimumSampleSize] = NA
    #to avoid dividing by zero we scale it as follows
    if (max(abs(RawCor.Weighted), na.rm = TRUE) == 1) {
        RawCor.Weighted = RawCor.Weighted / 1.0000001
    }
    if (max(abs(Cor.Standard), na.rm = TRUE) == 1)  {
        Cor.Standard = Cor.Standard / 1.0000001
    }
    RawZ.Weighted = sqrt(NoAvailable  - 2) *
        RawCor.Weighted / sqrt(1 - RawCor.Weighted ^ 2)
    Z.Standard = sqrt(NoAvailable  - 2) *  Cor.Standard / sqrt(1 - Cor.Standard ^
                                                                   2)

    if (sum(abs(Z.Standard), na.rm = TRUE)  > 0) {
        Z.Weighted = RawZ.Weighted / sum(abs(RawZ.Weighted), na.rm = TRUE) * sum(abs(Z.Standard), na.rm = TRUE)
    } # end of if
    h1 = Z.Weighted / sqrt(NoAvailable - 2)
    Cor.Weighted = h1 / sqrt(1 + h1 ^ 2)
    p.Weighted = as.numeric(2 * (1 - pt(abs(Z.Weighted), NoAvailable - 2)))
    p.Standard = 2 * (1 - pt(abs(Z.Standard), NoAvailable - 2))

    if (getQValues) {
        # since the function qvalue cannot handle missing data, we set missing
        # p-values to 1.
        p.Weighted2 = p.Weighted
        p.Standard2 = p.Standard
        p.Weighted2[is.na(p.Weighted)] = 1
        p.Standard2[is.na(p.Standard)] = 1

        q.Weighted = try(qvalue(p.Weighted2)$qvalues, silent = TRUE)
        q.Standard = try(qvalue(p.Standard2)$qvalues, silent = TRUE)

        if (class(q.Weighted) == "try - error") {
            warning(
                "Calculation of weighted q-values failed; the q-values
                will be returned as NAs."
            )
            q.Weighted = rep(NA, length(p.Weighted))
        }
        if (class(q.Standard) == "try - error") {
            warning(
                "Calculation of standard q-values failed; the q-values will
                be returned as NAs."
            )
            q.Standard = rep(NA, length(p.Standard))
        }
    } else {
        q.Weighted = rep(NA, length(p.Weighted))
        q.Standard = rep(NA, length(p.Standard))
        if (getQValues)
            printFlush(
                "networkScreening: Warning: package qvalue not found.
                q-values will not be calculated."
            )
    }
    rankCor.Weighted = rank(-abs(Cor.Weighted), ties.method = "first")
    rankCor.Standard = rank(-abs(Cor.Standard), ties.method = "first")
    printFlush(
        paste(
            "Proportion of agreement between lists based on
            abs(Cor.Weighted) and abs(Cor.Standard):"
        )
    )
    for (i in c(10, 20, 50, 100, 200, 500, 1000)) {
        printFlush(paste("Top ", i, " list of genes: prop. agree = ",
                         signif (
                             sum(rankCor.Weighted <= i & rankCor.Standard <= i,
                                 na.rm = TRUE) / i,
                             3
                         )))
    } # end of for loop


    datout = data.frame(
        p.Weighted,
        q.Weighted,
        Cor.Weighted,
        Z.Weighted,
        p.Standard,
        q.Standard,
        Cor.Standard,
        Z.Standard
    )
    names(datout) = sub("Cor", corFnc, names(datout), fixed = TRUE)
    datout
} # end of function


################################################################################
#
# Functions included from NetworkFunctions - PL - 07.R
# Selected ones only
#
################################################################################


#-------------------------------------------------------------------------------
# labeledHeatmap.R
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#
# .reverseRows <- function(Matrix)
#
#-------------------------------------------------------------------------------
#


.reverseRows <- function(Matrix) {
    ind = seq(from = dim(Matrix)[1],
              to = 1,
              by = -1)
    Matrix[ind,]
    #Matrix
}

.reverseVector <- function(Vector) {
    ind = seq(from = length(Vector),
              to = 1,
              by = -1)
    Vector[ind]
    #Vector
}

.extend <- function(x, n) {
    nRep = ceiling(n / length(x))
    rep(x, nRep)[1:n]
}

#-------------------------------------------------------------------------------
#
# labeledHeatmap <- function(Matrix, xLabels, yLabels, ...) {
#
#-------------------------------------------------------------------------------
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
#' To label rows and columns by color squares, use \code{colorLabels=TRUE};
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
#' # Copy and paste the whole example into an R session with an interactive plot window;
#' # alternatively, you may replace the command sizeGrWindow below by opening
#' # another graphical device such as pdf.
#'
#' # Generate a matrix to be plotted
#'
#' nCol = 8; nRow = 7;
#' mat = matrix(runif(nCol*nRow, min = -1, max = 1), nRow, nCol);
#'
#' rowColors = standardColors(nRow);
#' colColors = standardColors(nRow + nCol)[(nRow+1):(nRow + nCol)];
#'
#' rowColors;
#' colColors;
#'
#' sizeGrWindow(9,7)
#' par(mfrow = c(2,2))
#' par(mar = c(4, 5, 4, 6));
#'
#' # Label rows and columns by text:
#'
#' labeledHeatmap(mat, xLabels = colColors, yLabels = rowColors,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Text-labeled heatmap");
#'
#' # Label rows and columns by colors:
#'
#' rowLabels = paste("ME", rowColors, sep="");
#' colLabels = paste("ME", colColors, sep="");
#'
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Color-labeled heatmap");
#'
#' # Mix text and color labels:
#'
#' rowLabels[3] = "Row 3";
#' colLabels[1] = "Column 1";
#'
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Mix-labeled heatmap");
#'
#' # Color labels and additional text labels
#'
#' rowLabels = paste("ME", rowColors, sep="");
#' colLabels = paste("ME", colColors, sep="");
#'
#' extraRowLabels = paste("Row", c(1:nRow));
#' extraColLabels = paste("Column", c(1:nCol));
#'
#' # Extend margins to fit all labels
#' par(mar = c(6, 6, 4, 6));
#' labeledHeatmap(mat, xLabels = colLabels, yLabels = rowLabels,
#'                xSymbols = extraColLabels,
#'                ySymbols = extraRowLabels,
#'                colorLabels = TRUE,
#'                colors = greenWhiteRed(50),
#'                setStdMargins = FALSE,
#'                textMatrix = signif(mat, 2),
#'                main = "Text- + color-labeled heatmap");
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

#===============================================================================
#
# multi - page labeled heatmap
#
#===============================================================================



#' Labeled heatmap divided into several separate plots.
#'
#' This function produces labaled heatmaps divided into several plots. This is
#' useful for large heatmaps where labels on individual columns and rows may
#' become unreadably small (or overlap).
#'
#' The function \code{\link{labeledHeatmap}} is used to produce each plot/page;
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
#' To label rows and columns by color squares, use \code{colorLabels=TRUE};
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
#' @seealso The workhorse function \code{\link{labeledHeatmap}} for the actual
#' heatmap plot;
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

    ...)
{
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




#-------------------------------------------------------------------------------
#
# labeledBarplot <- function(Matrix, labels, ...) {
#
#-------------------------------------------------------------------------------
#
# Plots a barplot of the Matrix and writes the labels underneath such that they
# are readable.



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
labeledBarplot <-
    function(Matrix,
             labels,
             colorLabels = FALSE,
             colored = TRUE,
             setStdMargins = TRUE,
             stdErrors = NULL,
             cex.lab = NULL,
             xLabelsAngle = 45,
             ...)
    {
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

#-------------------------------------------------------------------------------
#
# sizeGrWindow
#
#-------------------------------------------------------------------------------
# if the current device isn't of the required dimensions, close it and open a
# new one.



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

#===============================================================================
# GreenToRed.R
#===============================================================================



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
#'\donotrun{
#'   par(mfrow = c(3, 1))
#'   displayColors(greenWhiteRed(50));
#'   title("gamma = 1")
#'   displayColors(greenWhiteRed(50, 3));
#'   title("gamma = 3")
#'   displayColors(greenWhiteRed(50, 0.5));
#'   title("gamma = 0.5")}
#'
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

#===============================================================================
#
# Color pallettes that are more friendly to people with common color blindness
#
#===============================================================================



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
blueWhiteRed <- function(n,
                         gamma = 1,
                         endSaturation = 1) {
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

#===============================================================================
#
# KeepCommonProbes
#
#-------------------------------------------------------------------------------
# Filters out probes that are not common to all datasets, and puts probes into the same order in each
# set. Works by creating dataframes of probe names and their indices and merging them all.



#' Keep probes that are shared among given data sets
#'
#' This function strips out probes that are not shared by all given data sets,
#' and orders the remaining common probes using the same order in all sets.
#'
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param orderBy index of the set by which probes are to be ordered.
#' @return Expression data in the same format as the input data, containing
#' only common probes.
#' @author Peter Langfelder
#' @seealso \code{\link{checkSets}}
#' @keywords misc
keepCommonProbes <- function(multiExpr, orderBy = 1) {
    size = checkSets(multiExpr)
    nSets = size$nSets
    if (nSets <= 0)
        stop("No expression data given!")

    Names = data.frame(Names = names(multiExpr[[orderBy]]$data))

    if (nSets > 1)
        for (set in (1:nSets)) {
            SetNames = data.frame(Names = names(multiExpr[[set]]$data),
                                  index = c(1:dim(multiExpr[[set]]$data)[2]))
            Names = merge(
                Names,
                SetNames,
                by.x = "Names",
                by.y = "Names",
                all = FALSE,
                sort = FALSE
            )
        }

    for (set in 1:nSets)
        multiExpr[[set]]$data = multiExpr[[set]]$data[, Names[, set + 1]]

    multiExpr
}

#-------------------------------------------------------------------------------
#
# addTraitToPCs
#
#-------------------------------------------------------------------------------

# Adds a trait vector to a set of eigenvectors.
# Caution: multiTraits is assumed to be a vector of lists with each list having
# an entry data which is a nSamples x nTraits data frame with an appropriate
# column name, not a vector.



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

    if (length(multiME) != nSets)
        stop("Numbers of sets in multiME and multiTraits parameters differ -
             must be the same.")

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


#-------------------------------------------------------------------------------
#
# CorrelationPreservation
#
#-------------------------------------------------------------------------------
#
# Given a set of multiME (or OrderedMEs), calculate the preservation values for
# each module in each pair
# of datasets and return them as a matrix



#' Preservation of eigengene correlations
#'
#' Calculates a summary measure of preservation of eigengene correlations
#' across data sets
#'
#' The function calculates the preservation of correlation of each eigengene
#' with all other eigengenes (optionally except the 'grey' eigengene) in all
#' pairs of sets.
#'
#' @param multiME consensus module eigengenes in a multi-set format. A vector
#' of lists with one list corresponding to each set. Each list must contain a
#' component \code{data} that is a data frame whose columns are consensus
#' module eigengenes.
#' @param setLabels names to be used for the sets represented in
#' \code{multiME}.
#' @param excludeGrey logical: exclude the 'grey' eigengene from preservation
#' measure?
#' @param greyLabel module label corresponding to the 'grey' module. Usually
#' this will be the character string \code{"grey"} if the labels are colors,
#' and the number 0 if the labels are numeric.
#' @return A data frame whose rows correspond to consensus module eigengenes
#' given in the input \code{multiME}, and columns correspond to all possible
#' set comparisons. The two sets compared in each column are indicated in the
#' column name.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSetMEs}} and module\code{\link{checkSets}} in
#' package moduleColor for more on eigengenes and the multi-set format
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
correlationPreservation <-
    function(multiME,
             setLabels,
             excludeGrey = TRUE,
             greyLabel = "grey") {
        nSets = length(multiME)
        if (nSets != length(setLabels))
            stop("The lengths of multiME and setLabels must equal.")
        if (nSets <= 1)
            stop("Something is wrong with argument multiME: its length is 0 or 1")
        Names = names(multiME[[1]]$data)
        if (excludeGrey) {
            Use = substring(Names, 3) != greyLabel
        } else {
            Use = rep(TRUE, times = length(Names))
        }
        No.Mods = ncol(multiME[[1]]$data[, Use])
        CP = matrix(0, nrow = No.Mods, ncol = nSets * (nSets - 1) / 2)
        diag(CP) = 1
        CPInd = 1
        CPNames = NULL
        for (i in 1:(nSets - 1))
            for (j in (i + 1):nSets) {
                corME1 = cor(multiME[[i]]$data[, Use], use = "p")
                corME2 = cor(multiME[[j]]$data[, Use], use = "p")
                d = 1 - abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2)) ^
                                     2))
                CP[, CPInd] = apply(d, 1, sum) - 1
                CPNames = c(CPNames,
                            paste(setLabels[i], "::", setLabels[j], collapse = ""))
                CPInd = CPInd + 1
            }
        CPx = as.data.frame(CP)
        names(CPx) = CPNames
        rownames(CPx) = Names[Use]
        CPx
    }


#-------------------------------------------------------------------------------
#
# setCorrelationPreservation
#
#-------------------------------------------------------------------------------
#
# Given a set of multiME (or OrderedMEs), calculate the preservation values for
# each each pair of datasets and return them as a matrix.



#' Summary correlation preservation measure
#'
#' Given consensus eigengenes, the function calculates the average correlation
#' preservation pair-wise for all pairs of sets.
#'
#' For each pair of sets, the function calculates the average preservation of
#' correlation among the eigengenes. Two preservation measures are available,
#' the abosolute preservation (high if the two correlations are similar and low
#' if they are different), and the hyperbolically scaled preservation, which
#' de-emphasizes preservation of low correlation values.
#'
#' @param multiME consensus module eigengenes in a multi-set format. A vector
#' of lists with one list corresponding to each set. Each list must contain a
#' component \code{data} that is a data frame whose columns are consensus
#' module eigengenes.
#' @param setLabels names to be used for the sets represented in
#' \code{multiME}.
#' @param excludeGrey logical: exclude the 'grey' eigengene from preservation
#' measure?
#' @param greyLabel module label corresponding to the 'grey' module. Usually
#' this will be the character string \code{"grey"} if the labels are colors,
#' and the number 0 if the labels are numeric.
#' @param method character string giving the correlation preservation measure
#' to use. Recognized values are (unique abbreviations of) \code{"absolute"},
#' \code{"hyperbolic"}.
#' @return A data frame with each row and column corresponding to a set given
#' in \code{multiME}, containing the pairwise average correlation preservation
#' values. Names and rownames are set to entries of \code{setLabels}.
#' @author Peter Langfelder
#' @seealso
#'
#' \code{\link{multiSetMEs}} for module eigengene calculation;
#'
#' \code{\link{plotEigengeneNetworks}} for eigengene network visualization.
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
setCorrelationPreservation <-
    function(multiME,
             setLabels,
             excludeGrey = TRUE,
             greyLabel = "grey",
             method = "absolute") {
        m = charmatch(method, c("absolute", "hyperbolic"))
        if (is.na(m)) {
            stop("Unrecognized method given. Recognized methods are absolute,
                 hyperbolic.")
        }
        nSets = length(multiME)
        if (nSets != length(setLabels))
            stop("The lengths of multiME and setLabels must equal.")
        if (nSets <= 1)
            stop("Something is wrong with argument multiME: its length is 0 or 1")
        Names = names(multiME[[1]]$data)
        if (excludeGrey) {
            Use = substring(Names, 3) != greyLabel
        } else {
            Use = rep(TRUE, times = length(Names))
        }
        No.Mods = ncol(multiME[[1]]$data[, Use])
        SCP = matrix(0, nrow = nSets, ncol = nSets)
        diag(SCP) = 0
        for (i in 1:(nSets - 1))
            for (j in (i + 1):nSets) {
                corME1 = cor(multiME[[i]]$data[, Use], use = "p")
                corME2 = cor(multiME[[j]]$data[, Use], use = "p")
                if (m == 1) {
                    d = 1 - abs(corME1 - corME2) / 2
                } else {
                    d = 1 - abs(tanh((corME1 - corME2) / (abs(corME1) + abs(corME2)) ^ 2))
                }
                SCP[i, j] = sum(d[upper.tri(d)]) / sum(upper.tri(d))
                SCP[j, i] = SCP[i, j]
            }
        SCPx = as.data.frame(SCP)
        names(SCPx) = setLabels
        rownames(SCPx) = setLabels
        SCPx
    }

#-------------------------------------------------------------------------------
#
# preservationNetworkDensity
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
#
# preservationNetworkConnectivity
#
#-------------------------------------------------------------------------------

# This function returns connectivities of nodes in preservation networks



#' Network preservation calculations
#'
#' This function calculates several measures of gene network preservation.
#' Given gene expression data in several individual data sets, it calculates
#' the individual adjacency matrices, forms the preservation network and
#' finally forms several summary measures of adjacency preservation for each
#' node (gene) in the network.
#'
#' The preservation network is formed from adjacencies of compared sets. For
#' 'complete' preservations, all given sets are compared at once; for
#' 'pairwise' preservations, the sets are compared in pairs. Unweighted
#' preservations are simple mean preservations for each node; their weighted
#' counterparts are weighted averages in which a preservation of adjacencies
#' \eqn{A^{(1)}_{ij}}{A[i,j; 1]} and \eqn{A^{(2)}_{ij}}{A[i,j; 2]} of nodes
#' \eqn{i,j} between sets 1 and 2 is weighted by \eqn{[ (A^{(1)}_{ij} +
#' A^{(2)}_{ij} )/2]^weightPower}{ ( (A[i,j; 1]+A[i,j; 2])/2)^weightPower}. The
#' hyperbolic preservation is based on \eqn{tanh[( max -
#' min)/(max+min)^2]}{tanh[( max - min)/(max+min)^2]}, where \eqn{max}{max} and
#' \eqn{min}{min} are the componentwise maximum and minimum of the compared
#' adjacencies, respectively.
#'
#' @param multiExpr expression data in the multi-set format (see
#' \code{\link{checkSets}}). A vector of lists, one per set. Each set must
#' contain a component \code{data} that contains the expression data, with rows
#' corresponding to samples and columns to genes or probes.
#' @param useSets optional specification of sets to be used for the
#' preservation calculation. Defaults to using all sets.
#' @param useGenes optional specification of genes to be used for the
#' preservation calculation. Defaults to all genes.
#' @param corFnc character string containing the name of the function to
#' calculate correlation. Suggested functions include \code{"cor"} and
#' \code{"bicor"}.
#' @param corOptions further argument to the correlation function.
#' @param networkType a character string encoding network type. Recognized
#' values are (unique abbreviations of) \code{"unsigned"}, \code{"signed"}, and
#' \code{"signed hybrid"}.
#' @param power soft thresholding power for network construction. Should be a
#' number greater than 1.
#' @param sampleLinks logical: should network connections be sampled
#' (\code{TRUE}) or should all connections be used systematically
#' (\code{FALSE})?
#' @param nLinks number of links to be sampled. Should be set such that
#' \code{nLinks * nNeighbors} be several times larger than the number of genes.
#' @param blockSize correlation calculations will be split into square blocks
#' of this size, to prevent running out of memory for large gene sets.
#' @param setSeed seed to be used for sampling, for repeatability. If a seed
#' already exists, it is saved before the sampling starts and restored upon
#' exit.
#' @param weightPower power with which higher adjacencies will be weighted in
#' weighted means
#' @param verbose integer level of verbosity. Zero means silent, higher values
#' make the output progressively more and more verbose.
#' @param indent indentation for diagnostic messages. Zero means no
#' indentation, each unit adds two spaces.
#' @return A list with the following components:
#'
#' \item{pairwise}{ a matrix with rows corresponding to genes and columns to
#' unique pairs of given sets, giving the pairwise preservation of the
#' adjacencies connecting the gene to all other genes.}
#'
#' \item{complete}{ a vector with one entry for each input gene containing the
#' complete mean preservation of the adjacencies connecting the gene to all
#' other genes.}
#'
#' \item{pairwiseWeighted}{ a matrix with rows corresponding to genes and
#' columns to unique pairs of given sets, giving the pairwise weighted
#' preservation of the adjacencies connecting the gene to all other genes.}
#'
#' \item{completeWeighted}{ a vector with one entry for each input gene
#' containing the complete weighted mean preservation of the adjacencies
#' connecting the gene to all other genes.}
#'
#' \item{pairwiseHyperbolic}{ a matrix with rows corresponding to genes and
#' columns to unique pairs of given sets, giving the pairwise hyperbolic
#' preservation of the adjacencies connecting the gene to all other genes.}
#'
#' \item{completeHyperbolic}{ a vector with one entry for each input gene
#' containing the complete mean hyperbolic preservation of the adjacencies
#' connecting the gene to all other genes.}
#'
#' \item{pairwiseWeightedHyperbolic}{ a matrix with rows corresponding to genes
#' and columns to unique pairs of given sets, giving the pairwise weighted
#' hyperbolic preservation of the adjacencies connecting the gene to all other
#' genes.}
#'
#' \item{completeWeightedHyperbolic}{ a vector with one entry for each input
#' gene containing the complete weighted hyperbolic mean preservation of the
#' adjacencies connecting the gene to all other genes.}
#' @author Peter Langfelder
#' @seealso \code{\link{adjacency}} for calculation of adjacency;
#' @references Langfelder P, Horvath S (2007) Eigengene networks for studying
#' the relationships between co-expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords misc
preservationNetworkConnectivity <- function(multiExpr,
                                            useSets = NULL,
                                            useGenes = NULL,
                                            corFnc = "cor",
                                            corOptions = "use = 'p'",
                                            networkType = "unsigned",
                                            power = 6,
                                            sampleLinks = NULL,
                                            nLinks = 5000,
                                            blockSize = 1000,
                                            setSeed = 12345,
                                            weightPower = 2,
                                            verbose = 2,
                                            indent = 0) {
    spaces = indentSpaces(indent)

    size = checkSets(multiExpr)
    nGenes = size$nGenes
    nSets = size$nSets
    if (!is.null(useSets) || !is.null(useGenes)) {
        if (is.null(useSets))
            useSets = c(1:nSets)
        if (is.null(useGenes))
            useGenes = c(1:nGenes)
        useExpr = vector(mode = "list", length = length(useSets))
        for (set in 1:length(useSets))
            useExpr[[set]] = list(data = multiExpr[[useSets[set]]]$data[, useGenes])
        multiExpr = useExpr
        rm(useExpr)
        collectGarbage()
    }
    size = checkSets(multiExpr)
    nGenes = size$nGenes
    nSets = size$nSets

    if (is.null(sampleLinks)) {
        sampleLinks = (nGenes > nLinks)
    }

    if (sampleLinks)
        nLinks = min(nLinks, nGenes)
    else
        nLinks = nGenes

    if (blockSize * nLinks > .largestBlockSize)
        blockSize = as.integer(.largestBlockSize / nLinks)

    intNetworkType = charmatch(networkType, .networkTypes)
    if (is.na(intNetworkType))
        stop(
            paste(
                "Unrecognized networkType argument. Recognized values are
                (unique abbreviations of)",
                paste(.networkTypes, collapse = ", ")
            )
        )

    subtract = rep(1, nGenes)
    if (sampleLinks) {
        if (verbose > 0)
            printFlush(
                paste(
                    spaces,
                    "preservationNetworkConnectivity:
                    selecting sample pool of size",
                    nLinks,
                    ".."
                )
            )
        sd = apply(multiExpr[[1]]$data, 2, sd, na.rm = TRUE)
        order = order(-sd)
        saved = FALSE
        if (exists(".Random.seed")) {
            saved = TRUE
            savedSeed = .Random.seed
            if (is.numeric(setSeed))
                set.seed(setSeed)
        }
        samplePool = order[sample(x = nGenes, size = nLinks)]
        if (saved) {
            .Random.seed <<- savedSeed
        }
        subtract[-samplePool] = 0
    }

    nPairComps = nSets * (nSets  - 1) / 2

    allPres = rep(NA, nGenes)
    allPresW = rep(NA, nGenes)
    allPresH = rep(NA, nGenes)
    allPresWH = rep(NA, nGenes)

    pairPres = matrix(nGenes, nPairComps)
    pairPresW = matrix(nGenes, nPairComps)
    pairPresH = matrix(nGenes, nPairComps)
    pairPresWH = matrix(nGenes, nPairComps)

    compNames = NULL
    for (set1 in 1:(nSets - 1))
        for (set2 in (set1 + 1):nSets)
            compNames = c(compNames, paste(set1, "vs", set2))

    dimnames(pairPres) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresW) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresH) = list(names(multiExpr[[1]]$data), compNames)
    dimnames(pairPresWH) = list(names(multiExpr[[1]]$data), compNames)

    if (verbose > 0)
    {
        pind = initProgInd(trailStr = " done")
    }

    nBlocks = as.integer((nGenes - 1) / blockSize)
    SetRestrConn = NULL
    start = 1
    if (sampleLinks)
    {
        corEval = parse(
            text = paste(
                corFnc,
                "(multiExpr[[set]]$data[, samplePool],
                multiExpr[[set]]$data[, blockIndex] ",
                prepComma(corOptions),
                ")"
            )
        )
    } else {
        corEval = parse(
            text = paste(
                corFnc,
                "(multiExpr[[set]]$data,
                multiExpr[[set]]$data[, blockIndex] ",
                prepComma(corOptions),
                ")"
            )
        )
    }

    while (start <= nGenes) {
        end = start + blockSize - 1
        if (end > nGenes)
            end = nGenes
        blockIndex = c(start:end)
        nBlockGenes = end - start + 1
        blockAdj = array(0, dim = c(nSets, nLinks, nBlockGenes))
        #if (verbose > 1) printFlush(paste(spaces, "..working on genes", start,
        #"through", end, "of", nGenes))
        for (set in 1:nSets) {
            c = eval(corEval)
            if (intNetworkType == 1) {
                c = abs(c)
            } else if (intNetworkType == 2) {
                c = (1 + c) / 2
            } else if (intNetworkType == 3) {
                c[c < 0] = 0
            } else {
                stop(
                    "Internal error: intNetworkType has wrong value:",
                    intNetworkType,
                    ". Sorry!"
                )
            }
            adj_mat = as.matrix(c ^ power)
            if (sum(is.na(adj_mat)) > 0)
                stop(
                    "NA values present in adjacency - this function cannot
                    handle them yet. Sorry!"
                )
            adj_mat[is.na(adj_mat)] = 0
            blockAdj[set, ,] = adj_mat
        }
        min = matrix(0, nLinks, nBlockGenes)
        which = matrix(0, nLinks, nBlockGenes)
        res = .C(
            "minWhichMin",
            as.double(blockAdj),
            as.integer(nSets),
            as.integer(nLinks * nBlockGenes),
            min = as.double(min),
            as.double(which)
        )
        min[,] = res$min
        max = matrix(0, nLinks, nBlockGenes)
        res = .C(
            "minWhichMin",
            as.double(-blockAdj),
            as.integer(nSets),
            as.integer(nLinks * nBlockGenes),
            min = as.double(min),
            as.double(which)
        )
        max[,] = -res$min
        rm(res)
        diff = max - min
        allPres[blockIndex] = (apply(1 - diff, 2, sum) - subtract[blockIndex]) /
            (nLinks - subtract[blockIndex])
        weight = ((max + min) / 2) ^ weightPower
        allPresW[blockIndex] = (apply((1 - diff) * weight, 2, sum) -
                                    subtract[blockIndex]) /
            (apply(weight, 2, sum) - subtract[blockIndex])
        hyp = 1 - tanh(diff / (max + min) ^ 2)
        allPresH[blockIndex] = (apply(hyp, 2, sum) - subtract[blockIndex]) /
            (nLinks - subtract[blockIndex])
        allPresWH[blockIndex] = (apply(hyp * weight, 2, sum) -
                                     subtract[blockIndex]) /
            (apply(weight, 2, sum) - subtract[blockIndex])

        compNames = NULL
        compInd = 1
        for (set1 in 1:(nSets - 1))
            for (set2 in (set1 + 1):nSets) {
                diff = abs(blockAdj[set1, ,] - blockAdj[set2, ,])
                compNames = c(compNames, paste(set1, "vs", set2))
                pairPres[blockIndex, compInd] = (apply(1 - diff, 2, sum) - subtract[blockIndex]) /
                    (nLinks - subtract[blockIndex])
                weight = ((blockAdj[set1, ,] + blockAdj[set2, ,]) / 2) ^
                    weightPower
                pairPresW[blockIndex, compInd] = (apply((1 - diff) * weight, 2, sum) - subtract[blockIndex]) /
                    (apply(weight, 2, sum) - subtract[blockIndex])
                hyp = 1 - tanh(diff / (blockAdj[set1, ,] + blockAdj[set2, ,]) ^
                                   2)
                pairPresH[blockIndex, compInd] = (apply(hyp, 2, sum) - subtract[blockIndex]) /
                    (nLinks - subtract[blockIndex])
                pairPresWH[blockIndex, compInd] = (apply(hyp * weight, 2, sum) - subtract[blockIndex]) /
                    (apply(weight, 2, sum) - subtract[blockIndex])
                compInd = compInd + 1
            }

        start = end + 1
        if (verbose > 0)
            pind = updateProgInd(end / nGenes, pind)
        collectGarbage()
    }
    if (verbose > 0)
        printFlush(" ")
    list(
        pairwise = pairPres,
        complete = allPres,
        pairwiseWeighted = pairPresW,
        completeWeighted = allPresW,
        pairwiseHyperbolic = pairPresH,
        completeHyperbolic = allPresH,
        pairwiseWeightedHyperbolic = pairPresWH,
        completeWeightedHyperbolic = allPresWH
    )
}

#-------------------------------------------------------------------------------
#
# plotEigengeneNetworks
#
#-------------------------------------------------------------------------------
# Plots a matrix plot of the ME(T)s. On the diagonal the heatmaps show
# correlation of MEs in the particular subset; off - diagonal are differences in
# the correlation matrix.
# setLabels is a vector of titles for the diagonal diagrams; the off - diagonal
# will have no title for now.



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
        Letters = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

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

#===============================================================================
#
# numbers2colors: convert a vector of numbers to colors
#
#===============================================================================

# Turn a numerical variable into a color indicator. x can be a matrix or a vector.
# For discrete variables, consider also labels2colors.



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
                           naColor = "grey")
{
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

#===============================================================================
#
#  Rand index calculation
#
#===============================================================================

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


# the following function computes the Rand index between 2 clusterings
# assesses how similar two clusterings are


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
randIndex <- function(tab, adjust = TRUE)
{
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


#===============================================================================
#
# Check expression data: mark genes and samples with too many missing entries
#
#===============================================================================
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
                      indent = 0)
{
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
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately;
#'
#' \code{\link{goodSamplesMS}}, \code{\link{goodSamplesGenesMS}} for additional
#' cleaning of multiple data sets together.
#' @keywords misc
goodGenesMS <-
    function(multiExpr,
             useSamples = NULL,
             useGenes = NULL,
             minFraction = 1 / 2,
             minNSamples = ..minNSamples,
             minNGenes = ..minNGenes,
             tol = NULL,
             verbose = 1,
             indent = 0)
    {
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
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately;
#'
#' \code{\link{goodGenesMS}}, \code{\link{goodSamplesGenesMS}} for additional
#' cleaning of multiple data sets together.
#' @keywords misc
goodSamplesMS <-
    function(multiExpr,
             useSamples = NULL,
             useGenes = NULL,
             minFraction = 1 / 2,
             minNSamples = ..minNSamples,
             minNGenes = ..minNGenes,
             verbose = 1,
             indent = 0)
    {
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
                             indent = 0)
{
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
#' \code{\link{goodSamplesGenes}} for cleaning individual sets separately;
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
                               indent = 0)
{
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

#===============================================================================
#
# modified heatmap plot: allow specifying the hang parameter for both side and
# top dendrograms
#
#===============================================================================
.heatmap <-
    function(x,
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
             ...)
    {
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


#===============================================================================
# The vectorize functions turns a matrix or data frame into a vector. If the
# matrix is not #symmetric the number of entries of the vector equals the number
# of rows times the #number of columns of the matrix.
# But if the matrix is symmetrical then it only uses the #entries in the upper
# triangular matrix.
# If the option diag  = TRUE, it also includes the diagonal elements of the
# symmetric matrix. By default it  excludes the diagonal elements of a symmetric
#  matrix.



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
vectorizeMatrix <- function(M, diag = FALSE)
{
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

#===============================================================================



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
scaleFreeFitIndex <-
    function(k,
             nBreaks = 10,
             removeFirst = FALSE) {
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

#===============================================================================



#' Standard Screening with regard to a Censored Time Variable
#'
#' The function standardScreeningCensoredTime computes association measures
#' between the columns of the input data datE and a censored time variable
#' (e.g. survival time). The censored time is specified using two input
#' variables "time" and "event". The event variable is binary where 1 indicates
#' that the event took place (e.g. the person died) and 0 indicates censored
#' (i.e. lost to follow up).  The function fits univariate Cox regression
#' models (one for each column of datE) and outputs a Wald test p-value, a
#' logrank p-value, corresponding local false discovery rates (known as
#' q-values, Storey et al 2004), hazard ratios. Further it reports the
#' concordance index (also know as area under the ROC curve) and optionally
#' results from dichotomizing the columns of datE.
#'
#'
#' If input option fastCalculation=TRUE, then the function outputs correlation
#' test p-values (and q-values) for correlating the columns of datE with the
#' expected hazard (if no covariate is fit). Specifically, the expected hazard
#' is defined as the deviance residual of an intercept only Cox regression
#' model. The results are very similar to those resulting from a univariate Cox
#' model where the censored time is regressed on the columns of dat.
#' Specifically, this computational speed up is facilitated by the insight that
#' the p-values resulting from a univariate Cox regression
#' coxph(Surv(time,event)~datE[,i]) are very similar to those from
#' corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples)
#'
#' @param time numeric variable showing time to event or time to last follow
#' up.
#' @param event Input variable \code{time} specifies the time to event or time
#' to last follow up. Input variable \code{event} indicates whether the event
#' happend (=1) or whether there was censoring (=0).
#' @param datExpr a data frame or matrix whose columns will be related to the
#' censored time.
#' @param percentiles numeric vector which is only used when
#' dichotomizationResults=T. Each value should lie between 0 and 1. For each
#' value specified in the vector percentiles, a binary vector will be defined
#' by dichotomizing the column value according to the corresponding quantile.
#' Next a corresponding p-value will be calculated.
#' @param dichotomizationResults logical. If this option is set to TRUE then
#' the values of the columns of datE will be dichotomized and corresponding Cox
#' regression p-values will be calculated.
#' @param qValues logical. If this option is set to TRUE (default) then
#' q-values will be calculated for the Cox regression p-values.
#' @param fastCalculation logical. If set to TRUE, the function outputs
#' correlation test p-values (and q-values) for correlating the columns of datE
#' with the expected hazard (if no covariate is fit). Specifically, the
#' expected hazard is defined as the deviance residual of an intercept only Cox
#' regression model. The results are very similar to those resulting from a
#' univariate Cox model where the censored time is regressed on the columns of
#' dat. Specifically, this computational speed up is facilitated by the insight
#' that the p-values resulting from a univariate Cox regression
#' coxph(Surv(time,event)~datE[,i]) are very similar to those from
#' corPvalueFisher(cor(devianceResidual,datE[,i]), nSamples).
#' @return If \code{fastCalculation} is \code{FALSE}, the function outputs a
#' data frame whose rows correspond to the columns of datE and whose columns
#' report \item{ID}{column names of the input data datExpr.}
#' \item{pvalueWald}{Wald test p-value from fitting a univariate Cox regression
#' model where the censored time is regressed on each column of datExpr.}
#' \item{qValueWald}{local false discovery rate (q-value) corresponding to the
#' Wald test p-value. } \item{pvalueLogrank}{Logrank p-value resulting from the
#' Cox regression model. Also known as score test p-value. For large sample
#' sizes this sould be similar to the Wald test p-value. }
#' \item{qValueLogrank}{local false discovery rate (q-value) corresponding to
#' the Logrank test p-value. } \item{HazardRatio}{hazard ratio resulting from
#' the Cox model. If the value is larger than 1, then high values of the column
#' are associated with shorter time, e.g. increased hazard of death. A hazard
#' ratio equal to 1 means no relationship between the column and time. HR<1
#' means that high values are associated with longer time, i.e. lower hazard.}
#' \item{CI.LowerLimitHR}{Lower bound of the 95 percent confidence interval of
#' the hazard ratio. } \item{CI.UpperLimitHR}{Upper bound of the 95 percent
#' confidence interval of the hazard ratio. } \item{C.index}{concordance index,
#' also known as C-index or area under the ROC curve. Calculated with the
#' rcorr.cens option outx=TRUE (ties are ignored).}
#' \item{MinimumDichotPvalue}{This is the smallest p-value from the
#' dichotomization results. To see which dichotomized variable (and percentile)
#' corresponds to the minimum, study the following columns. }
#' \item{pValueDichot0.1}{This columns report the p-value when the column is
#' dichotomized according to the specified percentile (here 0.1). The
#' percentiles are specified in the input option percentiles. }
#' \item{pvalueDeviance}{The p-value resulting from using a correlation test to
#' relate the expected hazard (deviance residual) with each (undichotomized)
#' column of datE. Specifically, the Fisher transformation is used to calculate
#' the p-value for the Pearson correlation. The resulting p-value should be
#' very similar to that of a univariate Cox regression model.}
#' \item{qvalueDeviance}{Local false discovery rate (q-value) corresponding to
#' pvalueDeviance.} \item{corDeviance}{Pearson correlation between the expected
#' hazard (deviance residual) with each (undichotomized) column of datExpr.}
#' @author Steve Horvath
#' @keywords misc
standardScreeningCensoredTime <- function(time,
                                          event,
                                          datExpr,
                                          percentiles = seq(from = 0.1, to = 0.9, by = 0.2),
                                          dichotomizationResults = FALSE,
                                          qValues = TRUE,
                                          fastCalculation = TRUE)
{
    datExpr = data.frame(datExpr)
    no.Columns = dim(as.matrix(datExpr))[[2]]
    m = dim(as.matrix(datExpr))[[1]]
    if (length(time)  != m)
        stop(
            "The length of the time variable does not equal the number of ",
            "rows of datExpr.\nConsider transposing datExpr."
        )
    if (length(event)  != m)
        stop(
            "The length of the event variable does not equal the number of ",
            "rows of datExpr.\nConsider transposing datExpr."
        )
    if (fastCalculation) {
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0,
                           2 * CumHazard,
                           -2 * log(CumHazard)  +
                               2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr,
                                     use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance,
                                         nSamples = no.nonMissing)
        qvalueDeviance = rep(NA, length(pvalueDeviance))
        rest1 = !is.na(pvalueDeviance)
        qvalueDeviance [rest1] = qvalue(pvalueDeviance [rest1])$qvalues

        datout = data.frame(ID = dimnames(datExpr)[[2]],
                            pvalueDeviance,
                            qvalueDeviance,
                            corDeviance)
    }
    if (!fastCalculation) {
        pvalueWald = rep(NA, no.Columns)
        HazardRatio = rep(NA, no.Columns)
        CI.UpperLimitHR = rep(NA, no.Columns)
        CI.LowerLimitHR = rep(NA, no.Columns)
        C.index = rep(NA, no.Columns)
        pvalueLogrank = rep(NA, no.Columns)
        pValuesDichotomized = data.frame(matrix(nrow = no.Columns,
                                                ncol = length(percentiles)))
        names(pValuesDichotomized) = paste0("pValueDichotPercentile",
                                            as.character(percentiles))
        fittemp = summary(coxph(Surv(time, event) ~ 1, na.action = na.exclude))
        CumHazard = predict(fittemp, type = "expected")
        martingale1 = event - CumHazard
        deviance0 = ifelse(event == 0,
                           2 * CumHazard,
                           -2 * log(CumHazard)  +
                               2 * CumHazard - 2)
        devianceresidual = sign(martingale1) * sqrt(deviance0)
        corDeviance = as.numeric(cor(devianceresidual, datExpr,
                                     use = "p"))
        no.nonMissing = sum(!is.na(time))
        pvalueDeviance = corPvalueFisher(cor = corDeviance,
                                         nSamples = no.nonMissing)


        for (i in 1:no.Columns) {
            Column = as.numeric(as.matrix(datExpr[, i]))
            var1 = var(Column, na.rm = TRUE)
            if (var1 == 0 | is.na(var1)) {
                pvalueWald[i] = NA
                pvalueLogrank[i] = NA
                HazardRatio[i] = NA
                CI.UpperLimitHR[i] = NA
                CI.LowerLimitHR[i] = NA
                C.index[i] = NA
            }  # end of              if (var1 == 0 | is.na(var1))
            if (var1  != 0 & !is.na(var1)) {
                cox1 = summary(coxph(Surv(time, event) ~ Column,
                                     na.action = na.exclude))
                pvalueWald[i] = cox1$coef[5]
                pvalueLogrank[i] = cox1$sctest[[3]]
                HazardRatio[i] = exp(cox1$coef[1])
                CI.UpperLimitHR[i] = exp(cox1$coef[1] + 1.96  *
                                             cox1$coef[3])
                CI.LowerLimitHR[i] = exp(cox1$coef[1] - 1.96  *
                                             cox1$coef[3])
                C.index[i] = rcorr.cens(Column, Surv(time, event),
                                        outx = TRUE)[[1]]
            } # end of   if (var1  != 0 & !is.na(var1))


            if (dichotomizationResults) {
                quantilesE = as.numeric(quantile(Column, prob = percentiles))
                for (j in 1:length(quantilesE)) {
                    ColumnDichot = I(Column > quantilesE[j])
                    var1 = var(ColumnDichot, na.rm = TRUE)
                    if (var1 == 0 | is.na(var1)) {
                        pValuesDichotomized[i, j] = NA
                    } # end of if
                    if (var1  != 0 & !is.na(var1)) {
                        coxh = summary(coxph(
                            Surv(time, event) ~
                                ColumnDichot,
                            na.action = na.exclude
                        ))
                        pValuesDichotomized[i, j] = coxh$coef[5]
                    } # end of if
                } # end of for (j
                MinimumDichotPvalue = apply(pValuesDichotomized,
                                            1, min, na.rm = TRUE)
            } # end of if (dichotomizationResults)



            if (!qValues) {
                datout = data.frame(
                    ID = dimnames(datExpr)[[2]],
                    pvalueWald,
                    pvalueLogrank,
                    pvalueDeviance,
                    corDeviance,
                    HazardRatio,
                    CI.LowerLimitHR,
                    CI.UpperLimitHR,
                    C.index
                )
            }      # end of      if (!qValues) {

        } # end of for (i in 1:no.Columns)


        if (qValues) {
            qvalueWald = rep(NA, length(pvalueWald))
            rest1 = !is.na(pvalueWald)
            qvalueWald [rest1] = qvalue(pvalueWald[rest1])$qvalues

            qvalueLogrank = rep(NA, length(pvalueLogrank))
            rest1 = !is.na(pvalueLogrank)
            qvalueLogrank [rest1] = qvalue(pvalueLogrank[rest1])$qvalues

            qvalueDeviance = rep(NA, length(pvalueDeviance))
            rest1 = !is.na(pvalueDeviance)
            qvalueDeviance [rest1] = qvalue(pvalueDeviance[rest1])$qvalues

            datout = data.frame(
                ID = dimnames(datExpr)[[2]],
                pvalueWald,
                qvalueWald,
                pvalueLogrank,
                qvalueLogrank,
                pvalueDeviance,
                qvalueDeviance,
                corDeviance,
                HazardRatio,
                CI.LowerLimitHR,
                CI.UpperLimitHR,
                C.index
            )
        } # end of  if (qValues)


        if (dichotomizationResults) {
            datout = data.frame(datout, MinimumDichotPvalue,
                                pValuesDichotomized)
        }
    }
    datout
} # end of function standardScreeningCensoredTime


#===============================================================================
#
# standardScreeningNumericTrait
#
#===============================================================================



#' Standard screening for numeric traits
#'
#' Standard screening for numeric traits based on Pearson correlation.
#'
#' The function calculates the correlations, associated p-values, area under
#' the ROC, and q-values
#'
#' @param datExpr data frame containing expression data (or more generally
#' variables to be screened), with rows corresponding to samples and columns to
#' genes (variables)
#' @param yNumeric a numeric vector giving the trait measurements for each
#' sample
#' @param corFnc correlation function.  Defaults to Pearson correlation but can
#' also be \code{\link{bicor}}.
#' @param corOptions list specifying additional arguments to be passed to the
#' correlation function given by \code{corFnc}.
#' @param alternative alternative hypothesis for the correlation test
#' @param qValues logical: should q-values be calculated?
#' @param areaUnderROC logical: should are under the receiver-operating curve
#' be calculated?
#' @return
#'
#' Data frame with the following components:
#'
#' \item{ID }{Gene (or variable) identifiers copied from
#' \code{colnames(datExpr)}} \item{cor}{correlations of all genes with the
#' trait} \item{Z}{Fisher Z statistics corresponding to the correlations}
#' \item{pvalueStudent }{Student p-values of the correlations}
#' \item{qvalueStudent }{(if input \code{qValues==TRUE}) q-values of the
#' correlations calculated from the p-values} \item{AreaUnderROC }{(if input
#' \code{areaUnderROC==TRUE}) area under the ROC} \item{nPresentSamples}{number
#' of samples present for the calculation of each association. }
#' @author Steve Horvath
#' @seealso \code{\link{standardScreeningBinaryTrait}},
#' \code{\link{standardScreeningCensoredTime}}
#' @keywords misc
standardScreeningNumericTrait <- function(datExpr,
                                          yNumeric,
                                          corFnc = cor,
                                          corOptions = list(use = 'p'),
                                          alternative = c("two.sided", "less", "greater"),
                                          qValues = TRUE,
                                          areaUnderROC = TRUE)
{
    datExpr = as.matrix(datExpr)
    nGenes = ncol(datExpr)
    nSamples = nrow(datExpr)
    if (length(yNumeric)  != nSamples)
        stop("the length of the sample trait y does not equal the number of ",
             "rows of datExpr")
    corPearson = rep(NA, nGenes)
    pvalueStudent = rep(NA, nGenes)
    AreaUnderROC = rep(NA, nGenes)
    nPresent = Z = rep(NA, nGenes)

    corFnc = match.fun(corFnc)
    corOptions$y = yNumeric
    corOptions$x = as.matrix(datExpr)
    cp = do.call(corFnc, corOptions)
    corPearson = as.numeric(cp)

    finMat = !is.na(datExpr)
    np = t(finMat) %*% (!is.na(as.matrix(yNumeric)))

    nPresent = as.numeric(np)

    ia = match.arg(alternative)
    T = sqrt(np - 2) * corPearson / sqrt(1 - corPearson ^ 2)
    if (ia == "two.sided") {
        p = 2 * pt(abs(T), np - 2, lower.tail = FALSE)
    }
    else if (ia == "less") {
        p = pt(T, np - 2, lower.tail = TRUE)
    }
    else if (ia == "greater") {
        p = pt(T, np - 2, lower.tail = FALSE)
    }
    pvalueStudent = as.numeric(p)

    Z = 0.5 * log((1 + corPearson) / (1 - corPearson)) * sqrt(nPresent  - 2)

    if (areaUnderROC)
        for (i in 1:dim(datExpr)[[2]])
        {
            AreaUnderROC[i] = rcorr.cens(datExpr[, i], yNumeric, outx = TRUE)[[1]]
        }

    q.Student = rep(NA, length(pvalueStudent))
    rest1 = !is.na(pvalueStudent)
    if (qValues)
    {
        x = try({
            q.Student[rest1] = qvalue(pvalueStudent[rest1])$qvalues
        },
        silent = TRUE)
        if (inherits(x, "try - error"))
            printFlush(
                paste(
                    "Warning in standardScreeningNumericTrait: function qvalue ",
                    "returned an error.\n",
                    "The returned qvalues will be invalid. The qvalue error: ",
                    x,
                    "\n"
                )
            )
    }
    if (is.null(colnames(datExpr))) {
        ID = paste0("Variable.", 1:ncol(datExpr))
    } else
        ID = colnames(datExpr)

    output = data.frame(
        ID = ID,
        cor = corPearson,
        Z = Z,
        pvalueStudent = pvalueStudent
    )
    if (qValues)
        output$qvalueStudent = q.Student
    if (areaUnderROC)
        output$AreaUnderROC = AreaUnderROC

    output$nPresentSamples = nPresent

    output
}


#===============================================================================
#
# metaZfunction
#
#===============================================================================



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

#===============================================================================
#
# rankPvalue
#
#===============================================================================



#' Estimate the p-value for ranking consistently high (or low) on multiple
#' lists
#'
#' The function rankPvalue calculates the p-value for observing that an object
#' (corresponding to a row of the input data frame \code{datS}) has a
#' consistently high ranking (or low ranking) according to multiple ordinal
#' scores (corresponding to the columns of the input data frame \code{datS}).
#'
#' The function calculates asymptotic p-values (and optionally q-values) for
#' testing the null hypothesis that the values in the columns of datS are
#' independent. This allows us to find objects (rows) with consistently high
#' (or low) values across the columns.
#'
#' Example: Imagine you have 5 vectors of Z statistics corresponding to the
#' columns of datS. Further assume that a gene has ranks 1,1,1,1,20 in the 5
#' lists. It seems very significant that the gene ranks number 1 in 4 out of
#' the 5 lists. The function rankPvalue can be used to calculate a p-value for
#' this occurrence.
#'
#' The function uses the central limit theorem to calculate asymptotic p-values
#' for two types of test statistics that measure consistently high or low
#' ordinal values. The first method (referred to as percentile rank method)
#' leads to accurate estimates of p-values if datS has at least 4 columns but
#' it can be overly conservative.  The percentile rank method replaces each
#' column datS by the ranked version rank(datS[,i]) (referred to ask low
#' ranking) and by rank(-datS[,i]) (referred to as high ranking). Low ranking
#' and high ranking allow one to find consistently small values or consistently
#' large values of datS, respectively.  All ranks are divided by the maximum
#' rank so that the result lies in the unit interval [0,1]. In the following,
#' we refer to rank/max(rank) as percentile rank. For a given object
#' (corresponding to a row of datS) the observed percentile rank follows
#' approximately a uniform distribution under the null hypothesis. The test
#' statistic is defined as the sum of the percentile ranks (across the columns
#' of datS). Under the null hypothesis that there is no relationship between
#' the rankings of the columns of datS, this (row sum) test statistic follows a
#' distribution that is given by the convolution of random uniform
#' distributions. Under the null hypothesis, the individual percentile ranks
#' are independent and one can invoke the central limit theorem to argue that
#' the row sum test statistic follows asymptotically a normal distribution.  It
#' is well-known that the speed of convergence to the normal distribution is
#' extremely fast in case of identically distributed uniform distributions.
#' Even when datS has only 4 columns, the difference between the normal
#' approximation and the exact distribution is negligible in practice (Killmann
#' et al 2001). In summary, we use the central limit theorem to argue that the
#' sum of the percentile ranks follows a normal distribution whose mean and
#' variance can be calculated using the fact that the mean value of a uniform
#' random variable (on the unit interval) equals 0.5 and its variance equals
#' 1/12.
#'
#' The second method for calculating p-values is referred to as scale method.
#' It is often more powerful but its asymptotic p-value can only be trusted if
#' either datS has a lot of columns or if the ordinal scores (columns of datS)
#' follow an approximate normal distribution.  The scale method scales (or
#' standardizes) each ordinal variable (column of datS) so that it has mean 0
#' and variance 1. Under the null hypothesis of independence, the row sum
#' follows approximately a normal distribution if the assumptions of the
#' central limit theorem are met. In practice, we find that the second approach
#' is often more powerful but it makes more distributional assumptions (if datS
#' has few columns).
#'
#' @param datS a data frame whose rows represent objects that will be ranked.
#' Each column of \code{datS} represents an ordinal variable (which can take on
#' negative values). The columns correspond to (possibly signed) object
#' significance measures, e.g., statistics (such as Z statistics), ranks, or
#' correlations.
#' @param columnweights allows the user to input a vector of non-negative
#' numbers reflecting weights for the different columns of \code{datZ}. If it
#' is set to \code{NULL} then all weights are equal.
#' @param na.last controls the treatment of missing values (NAs) in the rank
#' function. If \code{TRUE}, missing values in the data are put last (i.e. they
#' get the highest rank values). If \code{FALSE}, they are put first; if
#' \code{NA}, they are removed; if \code{"keep"} they are kept with rank NA.
#' See \code{\link{rank}} for more details.
#' @param ties.method represents the ties method used in the rank function for
#' the percentile rank method. See \code{\link{rank}} for more details.
#' @param calculateQvalue logical: should q-values be calculated? If set to
#' TRUE then the function calculates corresponding q-values (local false
#' discovery rates) using the qvalue package, see Storey JD and Tibshirani R.
#' (2003). This option assumes that qvalue package has been installed.
#' @param pValueMethod determines which method is used for calculating
#' p-values. By default it is set to "all", i.e. both methods are used. If it
#' is set to "rank" then only the percentile rank method is used. If it set to
#' "scale" then only the scale method will be used.
#' @return
#'
#' A list whose actual content depends on which p-value methods is selected,
#' and whether q0values are calculated. The following inner components are
#' calculated, organized in outer components \code{datoutrank} and
#' \code{datoutscale},:
#'
#' \item{pValueExtremeRank}{This is the minimum between pValueLowRank and
#' pValueHighRank, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowRank}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the rank method.}
#'
#' \item{pValueHighRank}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the rank method.}
#'
#' \item{pValueExtremeScale}{This is the minimum between pValueLowScale and
#' pValueHighScale, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowScale}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the Scale method.}
#'
#' \item{pValueHighScale}{Asymptotic p-value for observing a consistently low
#' value across the columns of datS based on the Scale method.}
#'
#' \item{qValueExtremeRank}{local false discovery rate (q-value) corresponding
#' to the p-value pValueExtremeRank}
#'
#' \item{qValueLowRank}{local false discovery rate (q-value) corresponding to
#' the p-value pValueLowRank}
#'
#' \item{qValueHighRank}{local false discovery rate (q-value) corresponding to
#' the p-value pValueHighRank}
#'
#' \item{qValueExtremeScale}{local false discovery rate (q-value) corresponding
#' to the p-value pValueExtremeScale}
#'
#' \item{qValueLowScale}{local false discovery rate (q-value) corresponding to
#' the p-value pValueLowScale}
#'
#' \item{qValueHighScale}{local false discovery rate (q-value) corresponding to
#' the p-value pValueHighScale}
#' @author Steve Horvath
#' @seealso \code{\link{rank}}, \code{\link{qvalue}}
#' @references Killmann F, VonCollani E (2001) A Note on the Convolution of the
#' Uniform and Related Distributions and Their Use in Quality Control. Economic
#' Quality Control Vol 16 (2001), No. 1, 17-41.ISSN 0940-5151
#'
#' Storey JD and Tibshirani R. (2003) Statistical significance for genome-wide
#' experiments. Proceedings of the National Academy of Sciences, 100:
#' 9440-9445.
#' @keywords misc
rankPvalue <- function(datS,
                       columnweights = NULL,
                       na.last = "keep",
                       ties.method = "average",
                       calculateQvalue = TRUE,
                       pValueMethod = "all")
{
    no.rows = dim(datS)[[1]]
    no.cols = dim(datS)[[2]]
    if (!is.null(columnweights) & no.cols  != length(columnweights))
        stop(
            "The number of components of the vector columnweights is unequal",
            "to the number of columns of datS. Hint: consider transposing ",
            "datS."
        )

    if (!is.null(columnweights)) {
        if (min(columnweights, na.rm = TRUE) < 0)
            stop(
                "At least one component of columnweights is negative, which ",
                "makes no sense. The entries should be positive numbers"
            )
        if (sum(is.na(columnweights)) > 0)
            stop(
                "At least one component of columnweights is missing, which ",
                "makes no sense. The entries should be positive numbers"
            )
        if (sum(columnweights) != 1) {
            # warning("The entries of columnweights do not sum to 1.
            # Therefore, they will divided by the sum. Then the resulting
            # weights sum to 1.")
            columnweights = columnweights / sum(columnweights)
        }
    }

    if (pValueMethod  != "scale") {
        percentilerank1 <- function(x) {
            R1 = rank(x, ties.method = ties.method, na.last = na.last)
            (R1 - .5) / max(R1, na.rm = TRUE)
        }

        datrankslow = apply(datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankslow = t(t(datrankslow) * columnweights)
        }
        datSpresent = !is.na(datS) + 0
        if (!is.null(columnweights)) {
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expectedsum = rowSums(datSpresent, na.rm = TRUE)  *
            0.5
        varsum = rowSums(datSpresent ^ 2, na.rm = TRUE) * 1 / 12
        observed.sumPercentileslow = as.numeric(rowSums(datrankslow,
                                                        na.rm = TRUE))
        Zstatisticlow = (observed.sumPercentileslow - expectedsum) / sqrt(varsum)
        datrankshigh = apply(-datS, 2, percentilerank1)
        if (!is.null(columnweights)) {
            datrankshigh = t(t(datrankshigh) * columnweights)
        }
        observed.sumPercentileshigh = as.numeric(rowSums(datrankshigh,
                                                         na.rm = TRUE))
        Zstatistichigh = (observed.sumPercentileshigh - expectedsum) / sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = pmin(pValueLow, pValueHigh)
        datoutrank = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme = pmin(qValueLow, qValueHigh)
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutrank = data.frame(datoutrank, datq)
            names(datoutrank) = paste0(names(datoutrank), "Rank")
        }
    }
    if (pValueMethod  != "rank") {
        datSpresent = !is.na(datS) + 0
        scaled.datS = scale(datS)
        if (!is.null(columnweights)) {
            scaled.datS = t(t(scaled.datS) * columnweights)
            datSpresent = t(t(datSpresent) * columnweights)
        }
        expected.value = rep(0, no.rows)
        varsum = rowSums(datSpresent ^ 2) * 1
        observed.sumScaleddatS = as.numeric(rowSums(scaled.datS, na.rm = TRUE))
        Zstatisticlow = (observed.sumScaleddatS - expected.value) / sqrt(varsum)
        scaled.minusdatS = scale(-datS)
        if (!is.null(columnweights)) {
            scaled.minusdatS = t(t(scaled.minusdatS) * columnweights)
        }
        observed.sumScaledminusdatS = as.numeric(rowSums(scaled.minusdatS,
                                                         na.rm = TRUE))
        Zstatistichigh = (observed.sumScaledminusdatS - expected.value) /
            sqrt(varsum)
        pValueLow = pnorm((Zstatisticlow))
        pValueHigh = pnorm((Zstatistichigh))
        pValueExtreme = 2 * pnorm(-abs(Zstatisticlow))
        datoutscale = data.frame(pValueExtreme, pValueLow, pValueHigh)
        if (calculateQvalue) {
            qValueLow = rep(NA, dim(datS)[[1]])
            qValueHigh = rep(NA, dim(datS)[[1]])
            qValueExtreme = rep(NA, dim(datS)[[1]])
            rest1 = !is.na(pValueLow)
            qValueLow[rest1] = qvalue(pValueLow[rest1])$qvalues
            rest1 = !is.na(pValueHigh)
            qValueHigh[rest1] = qvalue(pValueHigh[rest1])$qvalues
            rest1 = !is.na(pValueExtreme)
            qValueExtreme[rest1] = qvalue(pValueExtreme[rest1])$qvalues
            datq = data.frame(qValueExtreme, qValueLow, qValueHigh)
            datoutscale = data.frame(datoutscale, datq)
        }
        names(datoutscale) = paste0(names(datoutscale), "Scale")
    }
    if (pValueMethod == "rank") {
        datout = datoutrank
    }
    if (pValueMethod == "scale") {
        datout = datoutscale
    }
    if (pValueMethod  != "rank" & pValueMethod  != "scale")
        datout = data.frame(datoutrank, datoutscale)
    datout
} # End of function

#===============================================================================
#
# utility function: add a comma to string if the string is non - empty
#
#===============================================================================



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
#' prepComma("abc");
#' prepComma("");
#'
prepComma <- function(s)
{
    if (s == "")
        return (s)
    paste(", ", s)
}


#===============================================================================
#
# "restricted" q - value calculation
#
#===============================================================================



#' qvalue convenience wrapper
#'
#' This function calls \code{\link{qvalue}} on finite input p-values,
#' optionally traps errors from the q-value calculation, and returns just the q
#' values.
#'
#'
#' @param p a vector of p-values. Missing data are allowed and will be removed.
#' @param trapErrors logical: should errors generated by function
#' \code{\link{qvalue}} trapped? If \code{TRUE}, the errors will be silently
#' ignored and the returned q-values will all be \code{NA}.
#' @param \dots other arguments to function \code{\link{qvalue}}.
#' @return A vector of q-values. Entries whose corresponding p-values were not
#' finite will be \code{NA}.
#' @author Peter Langfelder
#' @seealso \code{\link{qvalue}}
#' @keywords misc
qvalue.restricted <- function(p, trapErrors = TRUE, ...) {
    fin = is.finite(p)
    qx = try(qvalue(p[fin], ...)$qvalues, silent = TRUE)
    q = rep(NA, length(p))
    if (inherits(qx, "try - error"))
    {
        if (!trapErrors)
            stop(qx)
    } else
        q[fin] = qx
    q
}


#===============================================================================
#
# consensusKME
#
#===============================================================================


.interleave <-
    function(matrices,
             nameBase = names(matrices),
             sep = ".",
             baseFirst = TRUE)
    {
        # Drop null entries in the list
        keep = sapply(matrices, function(x)
            ! is.null(x))
        nameBase = nameBase[keep]
        matrices = matrices[keep]

        nMats = length(matrices)
        nCols = ncol(matrices[[1]])

        dims = lapply(matrices, dim)

        if (baseFirst)
        {
            for (m in 1:nMats)
                colnames(matrices[[m]]) = paste0(nameBase[m],
                                                 sep, colnames(matrices[[m]]))
        } else {
            for (m in 1:nMats)
                colnames(matrices[[m]]) = paste0(colnames(matrices[[m]]),
                                                 sep, nameBase[m])
        }

        out = as.data.frame(lapply(1:nCols,
                                   function(index, matrices)
                                       as.data.frame(lapply(matrices,
                                                            function(x, i)
                                                                x[, i, drop = FALSE], index)),
                                   matrices))

        rownames(out) = rownames(matrices[[1]])
        out
    }




#' Calculate consensus kME (eigengene-based connectivities) across multiple
#' data sets.
#'
#' Calculate consensus kME (eigengene-based connectivities) across multiple
#' data sets, typically following a consensus module analysis.
#'
#' The function \code{corAndPvalueFnc} is currently is expected to accept
#' arguments \code{x} (gene expression profiles), \code{y} (eigengene
#' expression profiles), and \code{alternative} with possibilities at least
#' \code{"greater", "two.sided"}.  Any additional arguments can be passed via
#' \code{corOptions}.
#'
#' The function \code{corAndPvalueFnc} should return a list which at the least
#' contains (1) a matrix of associations of genes and eigengenes (this
#' component should have the name given by \code{corComponent}), and (2) a
#' matrix of the corresponding p-values, named "p" or "p.value". Other
#' components are optional but for full functionality should include (3)
#' \code{nObs} giving the number of observations for each association (which is
#' the number of samples less number of missing data - this can in principle
#' vary from association to association), and (4) \code{Z} giving a Z static
#' for each observation. If these are missing, \code{nObs} is calculated in the
#' main function, and calculations using the Z statistic are skipped.
#'
#' @param multiExpr Expression (or other numeric) data in a multi-set format. A
#' vector of lists; in each list there must be a component named `data' whose
#' content is a matrix or dataframe or array of dimension 2.
#' @param moduleLabels Module labels: one label for each gene in
#' \code{multiExpr}.
#' @param multiEigengenes Optional eigengenes of modules specified in
#' \code{moduleLabels}. If not given, will be calculated from \code{multiExpr}.
#' @param signed logical: should the network be considered signed? In signed
#' networks (\code{TRUE}), negative kME values are not considered significant
#' and the corresponding p-values will be one-sided. In unsigned networks
#' (\code{FALSE}), negative kME values are considered significant and the
#' corresponding p-values will be two-sided.
#' @param useModules Optional specification of module labels to which the
#' analysis should be restricted. This could be useful if there are many
#' modules, most of which are not interesting. Note that the "grey" module
#' cannot be used with \code{useModules}.
#' @param consensusQuantile Quantile for the consensus calculation. Should be a
#' number between 0 (minimum) and 1.
#' @param metaAnalysisWeights Optional specification of meta-analysis weights
#' for each input set. If given, must be a numeric vector of length equal the
#' number of input data sets (i.e., \code{length(multiExpr)}). These weights
#' will be used in addition to constant weights and weights proportional to
#' number of samples (observations) in each set.
#' @param corAndPvalueFnc Function that calculates associations between
#' expression profiles and eigengenes. See details.
#' @param corOptions List giving additional arguments to function
#' \code{corAndPvalueFnc}. See details.
#' @param corComponent Name of the component of output of
#' \code{corAndPvalueFnc} that contains the actual correlation.
#' @param getQvalues logical: should q-values (estimates of FDR) be calculated?
#' @param useRankPvalue Logical: should the \code{\link{rankPvalue}} function
#' be used to obtain alternative meta-analysis statistics?
#' @param rankPvalueOptions Additional options for function
#' \code{\link{rankPvalue}}. These include \code{na.last} (default
#' \code{"keep"}), \code{ties.method} (default \code{"average"}),
#' \code{calculateQvalue} (default copied from input \code{getQvalues}), and
#' \code{pValueMethod} (default \code{"scale"}). See the help file for
#' \code{\link{rankPvalue}} for full details.
#' @param setNames names for the input sets. If not given, will be taken from
#' \code{names(multiExpr)}. If those are \code{NULL} as well, the names will be
#' \code{"Set_1", "Set_2", ...}.
#' @param excludeGrey logical: should the grey module be excluded from the kME
#' tables? Since the grey module is typically not a real module, it makes
#' little sense to report kME values for it.
#' @param greyLabel label that labels the grey module.
#' @return Data frame with the following components (for easier readability the
#' order here is not the same as in the actual output): \item{ID}{Gene ID,
#' taken from the column names of the first input data set}
#'
#' \item{consensus.kME.1, consensus.kME.2, ...}{Consensus kME (that is, the
#' requested quantile of the kMEs in the individual data sets)in each module
#' for each gene across the input data sets. The module labels (here 1, 2,
#' etc.) correspond to those in \code{moduleLabels}.}
#'
#' \item{weightedAverage.equalWeights.kME1, weightedAverage.equalWeights.kME2,
#' ...}{ Average kME in each module for each gene across the input data sets. }
#'
#' \item{weightedAverage.RootDoFWeights.kME1,
#' weightedAverage.RootDoFWeights.kME2, ...}{ Weighted average kME in each
#' module for each gene across the input data sets. The weight of each data set
#' is proportional to the square root of the number of samples in the set. }
#'
#' \item{weightedAverage.DoFWeights.kME1, weightedAverage.DoFWeights.kME2,
#' ...}{ Weighted average kME in each module for each gene across the input
#' data sets. The weight of each data set is proportional to number of samples
#' in the set. }
#'
#' \item{weightedAverage.userWeights.kME1, weightedAverage.userWeights.kME2,
#' ...}{ (Only present if input \code{metaAnalysisWeights} is non-NULL.)
#' Weighted average kME in each module for each gene across the input data
#' sets. The weight of each data set is given in \code{metaAnalysisWeights}.}
#'
#' \item{meta.Z.equalWeights.kME1, meta.Z.equalWeights.kME2, ...}{Meta-analysis
#' Z statistic for kME in each module, obtained by weighing the Z scores in
#' each set equally. Only returned if the function \code{corAndPvalueFnc}
#' returns the Z statistics corresponding to the correlations.}
#'
#' \item{meta.Z.RootDoFWeights.kME1, meta.Z.RootDoFWeights.kME2, ...}{
#' Meta-analysis Z statistic for kME in each module, obtained by weighing the Z
#' scores in each set by the square root of the number of samples. Only
#' returned if the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the correlations.}
#'
#' \item{meta.Z.DoFWeights.kME1, meta.Z.DoFWeights.kME2, ...}{Meta-analysis Z
#' statistic for kME in each module, obtained by weighing the Z scores in each
#' set by the number of samples. Only returned if the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations.}
#'
#' \item{meta.Z.userWeights.kME1, meta.Z.userWeights.kME2, ...}{Meta-analysis Z
#' statistic for kME in each module, obtained by weighing the Z scores in each
#' set by \code{metaAnalysisWeights}.  Only returned if
#' \code{metaAnalysisWeights} is non-NULL and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations.}
#'
#' \item{meta.p.equalWeights.kME1, meta.p.equalWeights.kME2, ...}{ p-values
#' obtained from the equal-weight meta-analysis Z statistics. Only returned if
#' the function \code{corAndPvalueFnc} returns the Z statistics corresponding
#' to the correlations. }
#'
#' \item{meta.p.RootDoFWeights.kME1, meta.p.RootDoFWeights.kME2, ...}{ p-values
#' obtained from the meta-analysis Z statistics with weights proportional to
#' the square root of the number of samples. Only returned if the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations. }
#'
#' \item{meta.p.DoFWeights.kME1, meta.p.DoFWeights.kME2, ...}{ p-values
#' obtained from the degree-of-freedom weight meta-analysis Z statistics. Only
#' returned if the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the correlations. }
#'
#' \item{meta.p.userWeights.kME1, meta.p.userWeights.kME2, ...}{ p-values
#' obtained from the user-supplied weight meta-analysis Z statistics. Only
#' returned if \code{metaAnalysisWeights} is non-NULL and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the
#' correlations. }
#'
#' \item{meta.q.equalWeights.kME1, meta.q.equalWeights.kME2, ...}{ q-values
#' obtained from the equal-weight meta-analysis p-values. Only present if
#' \code{getQvalues} is \code{TRUE} and the function \code{corAndPvalueFnc}
#' returns the Z statistics corresponding to the kME values.}
#'
#' \item{meta.q.RootDoFWeights.kME1, meta.q.RootDoFWeights.kME2, ...}{ q-values
#' obtained from the meta-analysis p-values with weights proportional to the
#' square root of the number of samples. Only present if \code{getQvalues} is
#' \code{TRUE} and the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the kME values.}
#'
#' \item{meta.q.DoFWeights.kME1, meta.q.DoFWeights.kME2, ...}{ q-values
#' obtained from the degree-of-freedom weight meta-analysis p-values. Only
#' present if \code{getQvalues} is \code{TRUE} and the function
#' \code{corAndPvalueFnc} returns the Z statistics corresponding to the kME
#' values.}
#'
#' \item{meta.q.userWeights.kME1, meta.q.userWeights.kME2, ...}{ q-values
#' obtained from the user-specified weight meta-analysis p-values. Only present
#' if \code{metaAnalysisWeights} is non-NULL, \code{getQvalues} is \code{TRUE}
#' and the function \code{corAndPvalueFnc} returns the Z statistics
#' corresponding to the kME values.}
#'
#' The next set of columns contain the results of function
#' \code{\link{rankPvalue}} and are only present if input \code{useRankPvalue}
#' is \code{TRUE}. Some columns may be missing depending on the options
#' specified in \code{rankPvalueOptions}. We explicitly list columns that are
#' based on weighing each set equally; names of these columns carry the suffix
#' \code{.equalWeights}
#'
#' \item{pValueExtremeRank.ME1.equalWeights,
#' pValueExtremeRank.ME2.equalWeights, ...}{ This is the minimum between
#' pValueLowRank and pValueHighRank, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowRank.ME1.equalWeights, pValueLowRank.ME2.equalWeights, ...}{
#' Asymptotic p-value for observing a consistently low value across the columns
#' of datS based on the rank method.}
#'
#' \item{pValueHighRank.ME1.equalWeights, pValueHighRank.ME2.equalWeights,
#' ...}{ Asymptotic p-value for observing a consistently low value across the
#' columns of datS based on the rank method.}
#'
#' \item{pValueExtremeScale.ME1.equalWeights,
#' pValueExtremeScale.ME2.equalWeights, ...}{ This is the minimum between
#' pValueLowScale and pValueHighScale, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowScale.ME1.equalWeights, pValueLowScale.ME2.equalWeights,
#' ...}{ Asymptotic p-value for observing a consistently low value across the
#' columns of datS based on the Scale method.}
#'
#' \item{pValueHighScale.ME1.equalWeights, pValueHighScale.ME2.equalWeights,
#' ...}{ Asymptotic p-value for observing a consistently low value across the
#' columns of datS based on the Scale method.}
#'
#' \item{qValueExtremeRank.ME1.equalWeights,
#' qValueExtremeRank.ME2.equalWeights, ...}{ local false discovery rate
#' (q-value) corresponding to the p-value pValueExtremeRank}
#'
#' \item{qValueLowRank.ME1.equalWeights, qValueLowRank.ME2.equalWeights, ...}{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueLowRank}
#'
#' \item{qValueHighRank.ME1.equalWeights, lueHighRank.ME2.equalWeights, ...}{
#' local false discovery rate (q-value) corresponding to the p-value
#' pValueHighRank}
#'
#' \item{qValueExtremeScale.ME1.equalWeights,
#' qValueExtremeScale.ME2.equalWeights, ...}{ local false discovery rate
#' (q-value) corresponding to the p-value pValueExtremeScale}
#'
#' \item{qValueLowScale.ME1.equalWeights, qValueLowScale.ME2.equalWeights,
#' ...}{ local false discovery rate (q-value) corresponding to the p-value
#' pValueLowScale}
#'
#' \item{qValueHighScale.ME1.equalWeights,qValueHighScale.ME2.equalWeights,
#' ...}{ local false discovery rate (q-value) corresponding to the p-value
#' pValueHighScale}
#'
#' \item{...}{Analogous columns corresponding to weighing individual sets by
#' the square root of the number of samples, by number of samples, and by user
#' weights (if given). The corresponding column name suffixes are
#' \code{.RootDoFWeights}, \code{.DoFWeights}, and \code{.userWeights}.}
#'
#' The following set of columns summarize kME in individual input data sets.
#'
#' \item{kME1.Set_1, kME1.Set_2, ..., kME2.Set_1, kME2.Set_2, ...}{ kME values
#' for each gene in each module in each given data set. }
#'
#' \item{p.kME1.Set_1, p.kME1.Set_2, ..., p.kME2.Set_1, p.kME2.Set_2, ...}{
#' p-values corresponding to kME values for each gene in each module in each
#' given data set. }
#'
#' \item{q.kME1.Set_1, q.kME1.Set_2, ..., q.kME2.Set_1, q.kME2.Set_2, ...}{
#' q-values corresponding to kME values for each gene in each module in each
#' given data set. Only returned if \code{getQvalues} is \code{TRUE}. }
#'
#' \item{Z.kME1.Set_1, Z.kME1.Set_2, ..., Z.kME2.Set_1, Z.kME2.Set_2, ...}{ Z
#' statistics corresponding to kME values for each gene in each module in each
#' given data set. Only present if the function \code{corAndPvalueFnc} returns
#' the Z statistics corresponding to the kME values. }
#' @author Peter Langfelder
#' @seealso \link{signedKME} for eigengene based connectivity in a single data
#' set. \link{corAndPvalue}, \link{bicorAndPvalue} for two alternatives for
#' calculating correlations and the corresponding p-values and Z scores. Both
#' can be used with this function.
#' @references Langfelder P, Horvath S., WGCNA: an R package for weighted
#' correlation network analysis. BMC Bioinformatics. 2008 Dec 29; 9:559.
#' @keywords misc
consensusKME <-
    function(multiExpr,
             moduleLabels,
             multiEigengenes = NULL,
             consensusQuantile = 0,
             signed = TRUE,
             useModules = NULL,
             metaAnalysisWeights = NULL,
             corAndPvalueFnc = corAndPvalue,
             corOptions = list(),
             corComponent = "cor",
             getQvalues = FALSE,
             useRankPvalue = TRUE,
             rankPvalueOptions = list(calculateQvalue = getQvalues,
                                      pValueMethod = "scale"),
             setNames = NULL,
             excludeGrey = TRUE,
             greyLabel = if (is.numeric(moduleLabels))
                 0
             else
                 "grey")
    {
        corAndPvalueFnc = match.fun(corAndPvalueFnc)

        size = checkSets(multiExpr)
        nSets = size$nSets
        nGenes = size$nGenes
        nSamples = size$nSamples

        if (!is.null(metaAnalysisWeights))
            if (length(metaAnalysisWeights) != nSets)
                stop("Length of 'metaAnalysisWeights' must equal number of",
                     " input sets.")

        if (!is.null(useModules))
        {
            if (greyLabel %in% useModules)
                stop(
                    paste(
                        "Grey module (or module 0) cannot be used with ",
                        "'useModules'.\n Use 'excludeGrey = FALSE' to obtain ",
                        "results for the grey module as well."
                    )
                )
            keep = moduleLabels %in% useModules
            if (sum(keep) == 0)
                stop("Incorrectly specified 'useModules': no such module(s).")
            moduleLabels [!keep] = greyLabel
        }

        if (is.null(multiEigengenes))
            multiEigengenes = multiSetMEs(
                multiExpr,
                universalColors = moduleLabels,
                verbose = 0,
                excludeGrey = excludeGrey,
                grey = greyLabel
            )

        modLevels = substring(colnames(multiEigengenes[[1]]$data), 3)
        nModules = length(modLevels)

        kME = p = Z = nObs = array(NA, dim = c(nGenes, nModules, nSets))

        corOptions$alternative = c("two.sided", "greater")[signed + 1]

        haveZs = FALSE
        for (set in 1:nSets) {
            corOptions$x = multiExpr[[set]]$data
            corOptions$y = multiEigengenes[[set]]$data
            cp = do.call(corAndPvalueFnc, args = corOptions)
            corComp = grep(corComponent, names(cp))
            pComp = match("p", names(cp))
            if (is.na(pComp))
                pComp = match("p.value", names(cp))
            if (is.na(pComp))
                stop("Function `corAndPvalueFnc' did not return a p-value.")
            kME[, , set] = cp[[corComp]]
            p[, , set] = cp[[pComp]]
            if (!is.null(cp$Z)) {
                Z[, , set] = cp$Z
                haveZs = TRUE
            }
            if (!is.null(cp$nObs)) {
                nObs[, , set] = cp$nObs
            } else
                nObs[, , set] = t(is.na(multiExpr[[set]]$data)) %*% (!is.na(multiEigengenes[[set]]$data))
        }

        if (getQvalues)
        {
            q = apply(p, c(2:3), qvalue.restricted)
        } else
            q = NULL
        # not neccessary since weighted average also contains it
        # kME.average = rowMeans(kME, dims = 2)

        powers = c(0, 0.5, 1)
        nPowers = length(powers)
        nWeights = nPowers+!is.null(metaAnalysisWeights)
        weightNames = c("equalWeights",
                        "RootDoFWeights",
                        "DoFWeights",
                        "userWeights")[1:nWeights]
        kME.weightedAverage = array(NA, dim = c(nGenes, nWeights, nModules))
        for (m in 1:nWeights) {
            if (m <= nPowers) {
                weights = nObs ^ powers[m]
            } else
                weights = array(rep(metaAnalysisWeights, rep(nGenes * nModules,
                                                             nSets)),
                                dim = c(nGenes, nModules, nSets))
            kME.weightedAverage[, m,] = rowSums(kME * weights, na.rm = TRUE,
                                                dims = 2) /
                rowSums(weights, dims = 2, na.rm = TRUE)
        }

        dim(kME.weightedAverage) = c(nGenes * nWeights, nModules)

        if (any(is.na(kME)))
        {
            kME.consensus.1 = apply(kME,
                                    c(1, 2),
                                    quantile,
                                    prob = consensusQuantile,
                                    na.rm = TRUE)
            kME.consensus.2 = apply(kME,
                                    c(1, 2),
                                    quantile,
                                    prob = 1 - consensusQuantile,
                                    na.rm = TRUE)
            kME.median = apply(kME, c(1, 2), median, na.rm = TRUE)
        } else {
            kME.consensus.1 = matrix(colQuantileC(t(
                matrix(kME, nGenes * nModules, nSets)
            ),
            p = consensusQuantile),
            nGenes,
            nModules)
            kME.consensus.2 = matrix(colQuantileC(t(
                matrix(kME, nGenes * nModules, nSets)
            ),
            p = 1 - consensusQuantile),
            nGenes,
            nModules)
            kME.median = matrix(colQuantileC(t(
                matrix(kME, nGenes * nModules, nSets)
            ), p = 0.5),
            nGenes, nModules)
        }
        kME.consensus = ifelse(kME.median > 0, kME.consensus.1, kME.consensus.2)

        kME.consensus[kME.consensus * kME.median < 0] = 0

        # Prepare identifiers for the variables (genes)
        if (is.null(colnames(multiExpr[[1]]$data)))
        {
            ID = paste0("Variable.", 1:nGenes)
        } else
            ID = colnames(multiExpr[[1]]$data)

        # Get meta - Z, - p, - q values
        if (haveZs)
        {
            Z.kME.meta = p.kME.meta = array(0, dim = c(nGenes, nWeights, nModules))
            if (getQvalues)
                q.kME.meta = array(0, dim = c(nGenes, nWeights, nModules))
            for (m in 1:nWeights) {
                if (m <= nPowers) {
                    weights = nObs ^ powers[m]
                } else
                    weights = array(rep(metaAnalysisWeights,
                                        rep(nGenes * nModules, nSets)),
                                    dim = c(nGenes, nModules, nSets))

                Z1 = rowSums(Z * weights, na.rm = TRUE, dims = 2) / sqrt(rowSums(weights ^
                                                                                     2, na.rm = TRUE, dims = 2))
                if (signed) {
                    p1 = pnorm(Z1, lower.tail = FALSE)
                } else
                    p1 = 2 * pnorm(abs(Z1), lower.tail = FALSE)
                Z.kME.meta[, m,] = Z1
                p.kME.meta[, m,] = p1
                if (getQvalues) {
                    q1 = apply(p1, 2, qvalue.restricted)
                    q.kME.meta[, m,] = q1
                }
            }
            dim(Z.kME.meta) = dim(p.kME.meta) = c(nGenes *  nWeights, nModules)
            if (getQvalues) {
                dim(q.kME.meta) = c(nGenes * nWeights, nModules)
            } else
                q.kME.meta = NULL
        } else {
            Z.kME.meta = p.kME.meta = q.kME.meta = NULL
        }

        # Call rankPvalue

        if (useRankPvalue) {
            for (mod in 1:nModules)
                for (m in 1:nWeights)
                {
                    if (m <= nPowers) {
                        weights = nObs[, mod,] ^ powers[m]
                    } else
                        weights = matrix(metaAnalysisWeights, nGenes, nSets,
                                         byrow = TRUE)
                    # rankPvalue requires a vector of weights... so compress the weights to a vector.
                    # Output a warning if the compression loses information.
                    nDifferent = apply(weights, 2, function(x) {
                        length(unique(x))
                    })
                    if (any(nDifferent) > 1)
                        printFlush(
                            paste(
                                "Warning in consensusKME: rankPvalue requires compressed weights.\n",
                                "Some weights may not be entirely accurate."
                            )
                        )
                    cw = colMeans(weights, na.rm = TRUE)
                    rankPvalueOptions$columnweights = cw / sum(cw)

                    rankPvalueOptions$datS = kME[, mod,]
                    rp1 = do.call(rankPvalue, rankPvalueOptions)
                    colnames(rp1) = paste0(colnames(rp1),
                                           ".ME",
                                           modLevels[mod],
                                           ".",
                                           weightNames[m])
                    if (mod == 1 && m == 1) {
                        rp = rp1
                    } else
                        rp = cbind(rp, rp1)
                }
        }

        # Format the output... this will entail some rearranging of the individual set results.
        if (is.null(setNames))
            setNames = names(multiExpr)

        if (is.null(setNames))
            setNames = paste0("Set_", c(1:nSets))

        if (!haveZs)
            Z = NULL

        keep = c(TRUE, TRUE, getQvalues, haveZs)
        varNames = c("kME", "p.kME", "q.kME", "Z.kME")[keep]
        nVars = sum(keep)

        dimnames(kME) = list(mtd.colnames(multiExpr),
                             paste0("k", mtd.colnames(multiEigengenes)),
                             setNames)

        dimnames(p) = list(mtd.colnames(multiExpr),
                           paste0("p.k", mtd.colnames(multiEigengenes)),
                           setNames)

        if (getQvalues)
            dimnames(q) = list(mtd.colnames(multiExpr),
                               paste0("q.k", mtd.colnames(multiEigengenes)),
                               setNames)

        if (haveZs)
            dimnames(Z) = list(mtd.colnames(multiExpr),
                               paste0("Z.k", mtd.colnames(multiEigengenes)),
                               setNames)


        varList = list(
            kME = kME,
            p = p,
            q = if (getQvalues)
                q
            else
                NULL,
            Z = if (haveZs)
                Z
            else
                NULL
        )
        varList.interleaved = lapply(varList,
                                     function(arr) {
                                         if (!is.null(dim(arr))) {
                                             split = lapply(1:dim(arr)[3], function(i)
                                                 arr[, , i])
                                             .interleave(split, nameBase = setNames,
                                                         baseFirst = FALSE)
                                         } else
                                             NULL
                                     })

        # the following seems to choke on larger data sets, at least in R 3.2.1
        # combined = array(c (kME, p, q, Z), dim = c(nGenes, nModules, nSets, nVars))
        # recast = matrix(c(cast(melt(combined), X1~X4~X3~X2)), nGenes,
        # nSets * nModules * nVars)

        # ... so I will replace it with more cumbersome but hopefully workable code.

        recast = .interleave(varList.interleaved,
                             nameBase = rep("", 4),
                             sep = "")

        combinedMeta.0 = rbind(kME.consensus,
                               kME.weightedAverage,
                               Z.kME.meta,
                               p.kME.meta,
                               q.kME.meta)

        combinedMeta = matrix(combinedMeta.0, nGenes,
                              (1 + nWeights + (2 * haveZs + haveZs * getQvalues) *
                                   nWeights) * nModules)
        metaNames = c(
            "consensus.kME",
            paste0("weightedAverage.", weightNames, ".kME"),
            paste0("meta.Z.", weightNames, ".kME"),
            paste0("meta.p.", weightNames, ".kME"),
            paste0("meta.q.", weightNames, ".kME")
        )[c(
            rep(TRUE, nWeights + 1),
            rep(haveZs, nWeights),
            rep(haveZs, nWeights),
            rep(haveZs && getQvalues, nWeights)
        )]
        nMetaVars = length(metaNames)
        colnames(combinedMeta) = paste0 (rep(metaNames, nModules),
                                         rep(modLevels, rep(nMetaVars, nModules)))

        if (useRankPvalue) {
            out = data.frame(ID = ID, combinedMeta, rp, recast)
        } else
            out = data.frame(ID = ID, combinedMeta, recast)

        out
    }

#===============================================================================
#
# Meta - analysis
#
#===============================================================================

.isBinary <- function(multiTrait) {
    bin = TRUE
    for (set in 1:length(multiTrait))
        if (length(sort(unique(multiTrait[[set]]$data))) > 2)
            bin = FALSE
        bin
}



#' Meta-analysis of binary and continuous variables
#'
#' This is a meta-analysis complement to functions
#' \code{\link{standardScreeningBinaryTrait}} and
#' \code{\link{standardScreeningNumericTrait}}. Given expression (or other)
#' data from multiple independent data sets, and the corresponding clinical
#' traits or outcomes, the function calculates multiple screening statistics in
#' each data set, then calculates meta-analysis Z scores, p-values, and
#' optionally q-values (False Discovery Rates). Three different ways of
#' calculating the meta-analysis Z scores are provided: the Stouffer method,
#' weighted Stouffer method, and using user-specified weights.
#'
#' The Stouffer method of combines Z statistics by simply taking a mean of
#' input Z statistics and multiplying it by \code{sqrt(n)}, where \code{n} is
#' the number of input data sets. We refer to this method as
#' \code{Stouffer.equalWeights}. In general, a better (i.e., more powerful)
#' method of combining Z statistics is to weigh them by the number of degrees
#' of freedom (which approximately equals \code{n}). We refer to this method as
#' \code{weightedStouffer}. Finally, the user can also specify custom weights,
#' for example if a data set needs to be downweighted due to technical
#' concerns; however, specifying own weights by hand should be done carefully
#' to avoid possible selection biases.
#'
#' @param multiExpr Expression data (or other data) in multi-set format (see
#' \code{\link{checkSets}}). A vector of lists; in each list there must be a
#' component named \code{data} whose content is a matrix or dataframe or array
#' of dimension 2.
#' @param multiTrait Trait or ourcome data in multi-set format. Only one trait
#' is allowed; consequesntly, the \code{data} component of each component list
#' can be either a vector or a data frame (matrix, array of dimension 2).
#' @param binary Logical: is the trait binary (\code{TRUE}) or continuous
#' (\code{FALSE})? If not given, the decision will be made based on the content
#' of \code{multiTrait}.
#' @param metaAnalysisWeights Optional specification of set weights for
#' meta-analysis. If given, must be a vector of non-negative weights, one entry
#' for each set contained in \code{multiExpr}.
#' @param corFnc Correlation function to be used for screening. Should be
#' either the default \code{\link{cor}} or its robust alternative,
#' \code{\link{bicor}}.
#' @param corOptions A named list giving extra arguments to be passed to the
#' correlation function.
#' @param getQvalues Logical: should q-values (FDRs) be calculated?
#' @param getAreaUnderROC Logical: should area under the ROC be calculated?
#' Caution, enabling the calculation will slow the function down considerably
#' for large data sets.
#' @param useRankPvalue Logical: should the \code{\link{rankPvalue}} function
#' be used to obtain alternative meta-analysis statistics?
#' @param rankPvalueOptions Additional options for function
#' \code{\link{rankPvalue}}. These include \code{na.last} (default
#' \code{"keep"}), \code{ties.method} (default \code{"average"}),
#' \code{calculateQvalue} (default copied from input \code{getQvalues}), and
#' \code{pValueMethod} (default \code{"all"}).  See the help file for
#' \code{\link{rankPvalue}} for full details.
#' @param setNames Optional specification of set names (labels). These are used
#' to label the corresponding components of the output. If not given, will be
#' taken from the \code{names} attribute of \code{multiExpr}. If
#' \code{names(multiExpr)} is \code{NULL}, generic names of the form
#' \code{Set_1, Set2, ...} will be used.
#' @param kruskalTest Logical: should the Kruskal test be performed in addition
#' to t-test? Only applies to binary traits.
#' @param var.equal Logical: should the t-test assume equal variance in both
#' groups? If \code{TRUE}, the function will warn the user that the returned
#' test statistics will be different from the results of the standard
#' \code{\link[stats]{t.test}} function.
#' @param metaKruskal Logical: should the meta-analysis be based on the results
#' of Kruskal test (\code{TRUE}) or Student t-test (\code{FALSE})?
#' @param na.action Specification of what should happen to missing values in
#' \code{\link[stats]{t.test}}.
#' @return Data frame with the following components: \item{ID}{ Identifier of
#' the input genes (or other variables) }
#'
#' \item{Z.equalWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with equal weights} \item{p.equalWeights}{ p-values corresponding to
#' \code{Z.Stouffer.equalWeights} } \item{q.equalWeights}{ q-values
#' corresponding to \code{p.Stouffer.equalWeights}, only present if
#' \code{getQvalues} is \code{TRUE}.}
#'
#' \item{Z.RootDoFWeights}{ Meta-analysis Z statistics obtained using
#' Stouffer's method with weights given by the square root of the number of
#' (non-missing) samples in each data set} \item{p.RootDoFWeights}{ p-values
#' corresponding to \code{Z.DoFWeights} } \item{q.RootDoFWeights}{ q-values
#' corresponding to \code{p.DoFWeights}, only present if \code{getQvalues} is
#' \code{TRUE}. }
#'
#' \item{Z.DoFWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with weights given by the number of (non-missing) samples in each
#' data set} \item{p.DoFWeights}{ p-values corresponding to \code{Z.DoFWeights}
#' } \item{q.DoFWeights}{ q-values corresponding to \code{p.DoFWeights}, only
#' present if \code{getQvalues} is \code{TRUE}. }
#'
#' \item{Z.userWeights}{ Meta-analysis Z statistics obtained using Stouffer's
#' method with user-defined weights. Only present if input
#' \code{metaAnalysisWeights} are present.} \item{p.userWeights}{ p-values
#' corresponding to \code{Z.userWeights} } \item{q.userWeights}{ q-values
#' corresponding to \code{p.userWeights}, only present if \code{getQvalues} is
#' \code{TRUE}. }
#'
#' The next set of columns is present only if input \code{useRankPvalue} is
#' \code{TRUE} and contain the output of the function \code{\link{rankPvalue}}
#' with the same column weights as the above meta-analysis. Depending on the
#' input options \code{calculateQvalue} and \code{pValueMethod} in
#' \code{rankPvalueOptions}, some columns may be missing. The following columns
#' are calculated using equal weights for each data set.
#'
#' \item{pValueExtremeRank.equalWeights}{This is the minimum between
#' pValueLowRank and pValueHighRank, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowRank.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the rank method.}
#'
#' \item{pValueHighRank.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the rank method.}
#'
#' \item{pValueExtremeScale.equalWeights}{This is the minimum between
#' pValueLowScale and pValueHighScale, i.e. min(pValueLow, pValueHigh)}
#'
#' \item{pValueLowScale.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the Scale
#' method.}
#'
#' \item{pValueHighScale.equalWeights}{Asymptotic p-value for observing a
#' consistently low value across the columns of datS based on the Scale
#' method.}
#'
#' \item{qValueExtremeRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueExtremeRank}
#'
#' \item{qValueLowRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueLowRank}
#'
#' \item{qValueHighRank.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueHighRank}
#'
#' \item{qValueExtremeScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueExtremeScale}
#'
#' \item{qValueLowScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueLowScale}
#'
#' \item{qValueHighScale.equalWeights}{local false discovery rate (q-value)
#' corresponding to the p-value pValueHighScale}
#'
#' \item{...}{Analogous columns calculated by weighting each input set using
#' the square root of the number of samples, number of samples, and user
#' weights (if given). The corresponding column names carry the suffixes
#' \code{RootDofWeights}, \code{DoFWeights}, \code{userWeights}.}
#'
#' The following columns contain results returned by
#' \code{\link{standardScreeningBinaryTrait}} or
#' \code{\link{standardScreeningNumericTrait}} (depending on whether the input
#' trait is binary or continuous).
#'
#' For binary traits, the following information is returned for each set:
#'
#' \item{corPearson.Set_1, corPearson.Set_2,...}{Pearson correlation with a
#' binary numeric version of the input variable. The numeric variable equals 1
#' for level 1 and 2 for level 2. The levels are given by levels(factor(y)).}
#'
#' \item{t.Student.Set_1, t.Student.Set_2, ...}{Student t-test statistic}
#'
#' \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{two-sided Student
#' t-test p-value.}
#'
#' \item{qvalueStudent.Set_1, qvalueStudent.Set_2, ...}{(if input
#' \code{qValues==TRUE}) q-value (local false discovery rate) based on the
#' Student T-test p-value (Storey et al 2004).}
#'
#' \item{foldChange.Set_1, foldChange.Set_2, ...}{a (signed) ratio of mean
#' values. If the mean in the first group (corresponding to level 1) is larger
#' than that of the second group, it equals meanFirstGroup/meanSecondGroup.
#' But if the mean of the second group is larger than that of the first group
#' it equals -meanSecondGroup/meanFirstGroup (notice the minus sign).}
#'
#' \item{meanFirstGroup.Set_1, meanSecondGroup.Set_2, ...}{means of columns in
#' input \code{datExpr} across samples in the second group.}
#'
#' \item{SE.FirstGroup.Set_1, SE.FirstGroup.Set_2, ...}{standard errors of
#' columns in input \code{datExpr} across samples in the first group.  Recall
#' that SE(x)=sqrt(var(x)/n) where n is the number of non-missing values of x.
#' }
#'
#' \item{SE.SecondGroup.Set_1, SE.SecondGroup.Set_2, ...}{standard errors of
#' columns in input \code{datExpr} across samples in the second group.}
#'
#' \item{areaUnderROC.Set_1, areaUnderROC.Set_2, ...}{the area under the ROC,
#' also known as the concordance index or C.index.  This is a measure of
#' discriminatory power. The measure lies between 0 and 1 where 0.5 indicates
#' no discriminatory power. 0 indicates that the "opposite" predictor has
#' perfect discriminatory power. To compute it we use the function
#' \link[Hmisc]{rcorr.cens} with \code{outx=TRUE} (from Frank Harrel's package
#' Hmisc).}
#'
#' \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples
#' with finite measurements for each gene.}
#'
#' If input \code{kruskalTest} is \code{TRUE}, the following columns further
#' summarize results of Kruskal-Wallis test:
#'
#' \item{stat.Kruskal.Set_1, stat.Kruskal.Set_2, ...}{Kruskal-Wallis test
#' statistic.}
#'
#' \item{stat.Kruskal.signed.Set_1, stat.Kruskal.signed.Set_2,...}{(Warning:
#' experimental) Kruskal-Wallis test statistic including a sign that indicates
#' whether the average rank is higher in second group (positive) or first group
#' (negative).  }
#'
#' \item{pvaluekruskal.Set_1, pvaluekruskal.Set_2, ...}{Kruskal-Wallis test
#' p-value.}
#'
#' \item{qkruskal.Set_1, qkruskal.Set_2, ...}{q-values corresponding to the
#' Kruskal-Wallis test p-value (if input \code{qValues==TRUE}).}
#'
#' \item{Z.Set1, Z.Set2, ...}{Z statistics obtained from
#' \code{pvalueStudent.Set1, pvalueStudent.Set2, ...} or from
#' \code{pvaluekruskal.Set1, pvaluekruskal.Set2, ...}, depending on input
#' \code{metaKruskal}.}
#'
#' For numeric traits, the following columns are returned:
#'
#' \item{cor.Set_1, cor.Set_2, ...}{correlations of all genes with the trait}
#'
#' \item{Z.Set1, Z.Set2, ...}{Fisher Z statistics corresponding to the
#' correlations}
#'
#' \item{pvalueStudent.Set_1, pvalueStudent.Set_2, ...}{Student p-values of the
#' correlations}
#'
#' \item{qvalueStudent.Set_1, qvalueStudent.Set_1, ...}{(if input
#' \code{qValues==TRUE}) q-values of the correlations calculated from the
#' p-values}
#'
#' \item{AreaUnderROC.Set_1, AreaUnderROC.Set_2, ...}{area under the ROC}
#'
#' \item{nPresentSamples.Set_1, nPresentSamples.Set_2, ...}{number of samples
#' present for the calculation of each association. }
#' @author Peter Langfelder
#' @seealso \code{\link{standardScreeningBinaryTrait}},
#' \code{\link{standardScreeningNumericTrait}} for screening functions for
#' individual data sets
#' @references For Stouffer's method, see
#'
#' Stouffer, S.A., Suchman, E.A., DeVinney, L.C., Star, S.A. & Williams, R.M.
#' Jr. 1949. The American Soldier, Vol. 1: Adjustment during Army Life.
#' Princeton University Press, Princeton.
#'
#' A discussion of weighted Stouffer's method can be found in
#'
#' Whitlock, M. C., Combining probability from independent tests: the weighted
#' Z-method is superior to Fisher's approach, Journal of Evolutionary Biology
#' 18:5 1368 (2005)
#' @keywords misc
metaAnalysis <- function(multiExpr,
                         multiTrait,
                         binary = NULL,
                         #consensusQuantile = 0,
                         metaAnalysisWeights = NULL,
                         corFnc = cor,
                         corOptions = list(use = 'p'),
                         getQvalues = FALSE,
                         getAreaUnderROC = FALSE,
                         useRankPvalue = TRUE,
                         rankPvalueOptions = list(),
                         setNames = NULL,
                         kruskalTest = FALSE,
                         var.equal = FALSE,
                         metaKruskal = kruskalTest,
                         na.action = "na.exclude")
{
    size = checkSets(multiExpr)
    nSets = size$nSets

    for (set in 1:nSets)
        multiTrait[[set]]$data = as.matrix(multiTrait[[set]]$data)

    tSize = checkSets(multiTrait)
    if (tSize$nGenes != 1)
        stop("This function only works for a single trait.")

    if (size$nSets != tSize$nSets)
        stop("The number of sets in 'multiExpr' and 'multiTrait' ",
             "must be the same.")

    if (!all.equal(size$nSamples, tSize$nSamples))
        stop("Numbers of samples in each set of 'multiExpr' and ",
             "'multiTrait' must be the same.")

    #if (!is.finite(consensusQuantile) || consensusQuantile < 0 ||
    # consensusQuantile > 1)
    #   stop("'consensusQuantile' must be between 0 and 1.")

    if (is.null(setNames))
        setNames = names(multiExpr)

    if (is.null(setNames))
        setNames = paste0("Set_", c(1:nSets))

    if (metaKruskal && !kruskalTest)
        stop("Kruskal statistic meta - analysis requires kruskal test. ",
             "Use kruskalTest = TRUE.")

    if (is.null(binary))
        binary = .isBinary(multiTrait)

    if (!is.null(metaAnalysisWeights))
    {
        if (length(metaAnalysisWeights) != nSets)
            stop(
                "Length of 'metaAnalysisWeights' must equal the number",
                " of sets in 'multiExpr'."
            )
        if (any(!is.finite(metaAnalysisWeights)) || any(metaAnalysisWeights < 0))
            stop("All weights in 'metaAnalysisWeights' must be positive.")
    }

    setResults = list()

    for (set in 1:size$nSets) {
        if (binary) {
            setResults[[set]] = standardScreeningBinaryTrait(
                multiExpr[[set]]$data,
                as.vector(multiTrait[[set]]$data),
                kruskalTest = kruskalTest,
                qValues = getQvalues,
                var.equal = var.equal,
                na.action = na.action,
                corFnc = corFnc,
                corOptions = corOptions
            )
            trafo = TRUE
            if (metaKruskal) {
                metaStat = "stat.Kruskal.signed"
                metaP = "pvaluekruskal"
            } else {
                metaStat = "t.Student"
                metaP = "pvalueStudent"
            }
        } else {
            setResults[[set]] = standardScreeningNumericTrait(
                multiExpr[[set]]$data,
                as.vector(multiTrait[[set]]$data),
                qValues = getQvalues,
                corFnc = corFnc,
                corOptions = corOptions,
                areaUnderROC = getAreaUnderROC
            )
            metaStat = "Z"
            trafo = FALSE
        }
    }

    comb = NULL
    for (set in 1:nSets) {
        if (set == 1) {
            comb = setResults[[set]] [,-1]
            ID = setResults[[set]] [, 1]
            colNames = colnames(comb)
            nColumns = ncol(comb)
            colnames(comb) = paste0("X", c(1:nColumns))
        } else {
            xx = setResults[[set]][,-1]
            colnames(xx) = paste0("X", c(1:nColumns))
            comb = rbind(comb, xx)
        }
    }

    # Re - arrange comb:

    comb = matrix(as.matrix(as.data.frame(comb)), size$nGenes, nColumns * nSets)

    colnames(comb) = paste0(rep(colNames, rep(nSets, nColumns)), ".",
                            rep(setNames, nColumns))

    # Find the columns from which to do meta - analysis
    statCols = grep(paste0("^", metaStat), colnames(comb))
    if (length(statCols) == 0)
        stop("Internal error: no columns for meta - analysis found. Sorry!")
    setStats = comb[, statCols]

    if (trafo) {
        # transform p-values to Z statistics
        # Find the pvalue columns
        pCols = grep(paste0("^", metaP), colnames(comb))
        if (length(pCols) == 0)
            stop("Internal error: no columns for meta - analysis found. Sorry!")
        setP = comb[, pCols]
        # Caution: I assume here that the returned p-values are two - sided.
        setZ = sign(setStats) * qnorm(setP / 2, lower.tail = FALSE)
    } else {
        setZ = setStats
    }

    colnames(setZ) = paste0("Z.", setNames)
    nObsCols = grep("nPresentSamples", colnames(comb))
    nObs = comb[, nObsCols]

    powers = c(0, 0.5, 1)
    nPowers = 3

    metaNames = c("equalWeights", "RootDoFWeights", "DoFWeights")
    if (is.null(metaAnalysisWeights)) {
        nMeta = nPowers
    } else {
        nMeta = nPowers + 1
        metaNames = c(metaNames, "userWeights")
    }
    metaResults = NULL
    for (m in 1:nMeta) {
        if (m <= nPowers) {
            weights = nObs ^ powers[m]
        } else
            weights = matrix(metaAnalysisWeights, size$nGenes, nSets,
                             byrow = TRUE)

        metaZ = rowSums(setZ * weights, na.rm = TRUE) / sqrt(rowSums(weights ^
                                                                         2, na.rm = TRUE))
        p.meta = 2 * pnorm(abs(metaZ), lower.tail = FALSE)
        if (getQvalues) {
            q.meta = qvalue.restricted(p.meta)
            meta1 = cbind(metaZ, p.meta, q.meta)
        } else {
            q.meta = NULL
            meta1 = cbind(metaZ, p.meta)
        }
        colnames(meta1) = paste0(c("Z.", "p.", "q.")[1:ncol(meta1)],
                                 metaNames[m])
        metaResults = cbind(metaResults, meta1)
    }

    # Use rankPvalue to produce yet another meta - analysis

    rankMetaResults = NULL
    if (useRankPvalue) {
        rankPvalueOptions$datS = as.data.frame(setZ)
        if (is.na(match("calculateQvalue", names(rankPvalueOptions))))
            rankPvalueOptions$calculateQvalue = getQvalues
        for (m in 1:nMeta) {
            if (m <= nPowers) {
                weights = nObs ^ powers[m]
            } else
                weights = matrix(metaAnalysisWeights, size$nGenes, nSets,
                                 byrow = TRUE)

            # rankPvalue requires a vector of weights... so compress the weights to a vector.
            # Output a warning if the compression loses information.
            nDifferent = apply(weights, 2, function(x) {
                length(unique(x))
            })
            if (any(nDifferent) > 1)
                printFlush(
                    paste(
                        "Warning in metaAnalysis: rankPvalue requires compressed ",
                        "weights.\nSome weights may not be entirely accurate."
                    )
                )
            rankPvalueOptions$columnweights = colMeans(weights, na.rm = TRUE)
            rankPvalueOptions$columnweights = rankPvalueOptions$columnweights /
                sum(rankPvalueOptions$columnweights)
            rp = do.call(rankPvalue, rankPvalueOptions)
            colnames(rp) = paste0(colnames(rp), ".", metaNames[m])
            rankMetaResults = cbind(rankMetaResults, as.matrix(rp))
        }
    }

    # Put together the output

    out = list(ID = ID,
               metaResults,
               rankMetaResults,
               comb,
               if (trafo)
                   setZ
               else
                   NULL,
               NULL)
    # The last NULL is necessary so the line below works even if nothing else is
    #  NULL

    out = as.data.frame(out[-(which(sapply(out, is.null), arr.ind = TRUE))])

    out
}


#===============================================================================
#
# multiUnion and multiIntersect
#
#===============================================================================



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
#' @seealso The "standard" functions \code{\link{union}} and
#' \code{\link{intersect}}.
#' @keywords misc
multiUnion <- function(setList) {
    len = length(setList)
    if (len == 0)
        return(NULL)
    if (len == 1)
        return(setList[[1]])

    out = setList[[1]]
    for (elem in 2:len)
        out = union(out, setList[[elem]])

    out
}

multiIntersect <- function(setList) {
    len = length(setList)
    if (len == 0)
        return(NULL)
    if (len == 1)
        return(setList[[1]])

    out = setList[[1]]
    for (elem in 2:len)
        out = intersect(out, setList[[elem]])

    out
}

#===============================================================================
#
# prependZeros
#
#===============================================================================
# prepend as many zeros as necessary to fill number to a certain width.
# Assumes an integer input.



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
    for (i in 1:n)
        if (lengths[i] < width)
            out[i] = paste0(paste(rep("0", width - lengths[i]), collapse = ""),
                            x[i])

    out
}

#===============================================================================
#
# Text formatting
#
#===============================================================================



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
formatLabels <-function(labels, maxCharPerLine = 14, split = " ", fixed = TRUE,
             newsplit = split, keepSplitAtEOL = TRUE) {
        n = length(labels)
        splitX = strsplit(labels, split = split, fixed = fixed)
        newLabels = rep("", n)
        for (l in 1:n) {
            nl = ""
            line = ""
            if (nchar(labels[l]) > 0)
                for (s in 1:length(splitX[[l]])) {
                    newLen = nchar(line) + nchar(splitX [[l]] [s])
                    if (nchar(line) < 5 | newLen <= maxCharPerLine) {
                        nl = paste(nl, splitX[[l]] [s], sep = newsplit)
                        line = paste(line, splitX[[l]] [s], sep = newsplit)
                    } else {
                        nl = paste(nl, splitX[[l]] [s], sep = paste0(if (keepSplitAtEOL)
                            newsplit
                            else
                                "",
                            "\n"))
                        line = splitX[[l]] [s]
                    }
                }
            newLabels[l] = nl
        }
        substring(newLabels, nchar(newsplit) + 1)
    }

#===============================================================================
#
# shortenStrings
#
#===============================================================================

.listRep <- function(data, n) {
    out = list()
    if (n >  0)
        for (i in 1:n)
            out[[i]] = data
        out
}

# Truncate labels at the last 'split' before given maximum length, add ...
# if the label is shortened.



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
#' @seealso \code{\link{gregexpr}}, the workhorse pattern matching function
#' \code{\link{formatLabels}} for splitting strings into multiple lines
#' @keywords misc
shortenStrings <-
    function(strings,
             maxLength = 25,
             minLength = 10,
             split = " ",
             fixed = TRUE,
             ellipsis = "...",
             countEllipsisInLength = FALSE) {
        dims = dim(strings)
        dnames = dimnames(strings)
        if (is.data.frame(strings))
        {
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
        if (length(split) > 0)
        {
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
                newLabels[l] = paste0(substring(strings[l], 1, splitPosition - 1), ellipsis)
            }
        }

        dim(newLabels) = dims
        dimnames(newLabels) = dnames
        if (outputDF)
            as.data.frame(newLabels)
        else
            newLabels
    }
