
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
#' @param greyLast Normally the color grey is reserved for unassigned genes
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
