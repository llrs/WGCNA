# Convenience functions for manipulating multiData structures

# Note: many of these function would be simpler to use if I used some sort of
# class/method technique to keep track of the class of each object internally.
# For example, I could then write a generic function "subset" that would work
# consistently on lists and multiData objects. Similarly, multiData2list would
# simply become a method of as.list, and as.list would be safe to use both on
# lists and on multiData objects.

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

# mtd.subset ####
#' Subset rows and columns in a multiData structure
#'
#' The function restricts each \code{data} component to the given columns and
#' rows.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function assumes a "strict" multiData structure unless
#' \code{permissive} is \code{TRUE}.
#'
#' @param multiData A multiData structure.
#' @param rowIndex A list in which each component corresponds to a set and is a
#' vector giving the rows to be retained in that set. All indexing methods
#' recognized by R can be used (numeric, logical, negative indexing, etc). If
#' \code{NULL}, all columns will be retained in each set. Note that setting
#' individual elements of \code{rowIndex} to \code{NULL} will lead to errors.
#' @param colIndex A vector giving the columns to be retained. All indexing
#' methods recognized by R can be used (numeric, logical, negative indexing,
#' etc). In addition, column names of the retained columns may be given; if a
#' given name cannot be matched to a column, an error will be thrown. If
#' \code{NULL}, all columns will be retained.
#' @param permissive logical: should the function tolerate "loose"
#' \code{multiData} input? Note that the subsetting may lead to cryptic errors
#' if the input \code{multiData} does not follow the "strict" format.
#' @param drop logical: should dimensions with extent 1 be dropped?
#' @return A multiData structure containing the selected rows and columns. Note
#' that result always retaines it dimension and other attributes.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure.
#' @keywords misc
mtd.subset = function(multiData, rowIndex = NULL, colIndex = NULL,
                      permissive = FALSE, drop = FALSE) {
  size = checkSets(multiData, checkStructure = permissive)
  if (!size$structureOK && !is.null(colIndex))
    warning(immediate. = TRUE,
            paste("mtd.subset: applying column selection on data sets that do not have\n",
                  " the same number of columns. This is treacherous territory; proceed with caution."))
  if (is.null(colIndex)) colIndex.1 = c(1:size$nGenes) else colIndex.1 = colIndex
  if (is.null(rowIndex)) rowIndex = lapply(size$nSamples, function(n) {c(1:n)})
  if (length(rowIndex)!=size$nSets)
    stop("If given, 'rowIndex' must be a list of the same length as 'multiData'.")
  out = list()
  for (set in 1:size$nSets)
  {
    if (permissive)
      if (is.null(colIndex)) colIndex.1 = c(1:ncol(multiData[[set]]$data)) else colIndex.1 = colIndex
    if (is.character(colIndex.1))
    {
      colIndex.1 = match(colIndex.1, colnames(multiData[[set]]$data))
      n1 = length(colIndex.1)
      if (any(is.na(colIndex.1)))
        stop("Cannot match the following entries in 'colIndex' to column names in set ", set, ":\n",
             paste( colIndex[is.na(colIndex.1)] [1:min(n1, 5)], collapse = ", "),
             if (n1>5) ", ... [output truncated]" else "")
    }
    out[[set]] = list(data = multiData[[set]]$data[rowIndex[[set]], colIndex.1, drop = drop])
  }
  names(out) = names(multiData)
  out
}

# list2multiData ####
#' Convert a list to a multiData structure and vice-versa.
#'
#' \code{list2multiData} converts a list to a multiData structure;
#' \code{multiData2list} does the inverse.
#'
#' A multiData structure is a vector of lists (one list for each set) where
#' each list has a component \code{data} containing some useful information.
#'
#' @aliases list2multiData multiData2list
#' @param data A list to be converted to a multiData structure.
#' @return For \code{list2multiData}, a multiData structure; for
#' \code{multiData2list}, the corresponding list.
#' @author Peter Langfelder
#' @keywords misc
list2multiData = function(data) {
  out = list()
  for (set in 1:length(data))
    out[[set]] = list(data = data[[set]])
  names(out) = names(data)
  out
}

#' @rdname list2multiData
#' @aliases multiData2list
#' @param multiData A multiData structure to be converted to a list.
#' @export
multiData2list = function(multiData) {
    lapply(multiData, getElement, 'data')
}

.calculateIndicator = function(nSets, mdaExistingResults, mdaUpdateIndex) {
  if (length(mdaUpdateIndex)==0) {
      mdaUpdateIndex = NULL
  }
  calculate = rep(TRUE, nSets)
  if (!is.null(mdaExistingResults)) {
    nSets.existing = length(mdaExistingResults)
    if (nSets.existing>nSets){
      stop("Number of sets in 'mdaExistingResults' is higher than the number ",
           "of sets in 'multiData'.\n  Please supply a valid ",
           "'mdaExistingResults' or NULL to recalculate all results.")
    }
    if (nSets.existing==0) {
      stop("Number of sets in 'mdaExistingResults' is zero.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to ",
           "recalculate all results.")
    }
    if (is.null(mdaUpdateIndex)) {
      calculate[1:length(mdaExistingResults)] = FALSE
    } else {
      if (any(! mdaUpdateIndex %in% c(1:nSets))) {
        stop("All entries in 'mdaUpdateIndex' must be between 1 and the number",
             " of sets in 'multiData'.")
      }
      calculateIndex = sort(unique(c(mdaUpdateIndex,
                                      ifelse(nSets.existing<nSets,
                                             c((nSets.existing+1):nSets),
                                             NULL))))
      calculate[ c(1:nSets)[-calculateIndex] ] = FALSE
    }
  }

  calculate
}

# mtd.apply ####
#' Apply a function to each set in a multiData structure.
#'
#' Inspired by \code{\link{lapply}}, these functions apply a given function to
#' each \code{data} component in the input \code{multiData} structure, and
#' optionally simplify the result to an array if possible.
#'
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' \code{mtd.apply} works on any "loose" multiData structure;
#' \code{mtd.applyToSubset} assumes (and checks for) a "strict" multiData
#' structure.
#'
#' @aliases mtd.apply mtd.applyToSubset
#' @param multiData A multiData structure to apply the function over
#' @param FUN Function to be applied.
#' @param \dots Other arguments to the function \code{FUN}.
#' @param mdaExistingResults Optional list that contains previously calculated
#' results. This can be useful if only a few sets in \code{multiData} have
#' changed and recalculating the unchanged ones is computationally expensive.
#' If not given, all calculations will be performed. If given, components of
#' this list are copied into the output. See \code{mdmUpdateIndex} for which
#' components are re-calculated by default.
#' @param mdaUpdateIndex Optional specification of which sets in
#' \code{multiData} the calculation should actually be carried out. This
#' argument has an effect only if \code{mdaExistingResults} is non-NULL. If the
#' length of \code{mdaExistingResults} (call the length `k') is less than the
#' number of sets in \code{multiData}, the function assumes that the existing
#' results correspond to the first `k' sets in \code{multiData} and the rest of
#' the sets are automatically calculated, irrespective of the setting of
#' \code{mdaUpdateIndex}. The argument \code{mdaUpdateIndex} can be used to
#' specify re-calculation of some (or all) of the results that already exist in
#' \code{mdaExistingResults}.
#' @param mdaCopyNonData Logical: should non-data components of
#' \code{multiData} be copied into the output? Note that the copying is
#' incompatible with simplification; enabling both will trigger an error.
#' @param mdaSimplify Logical: should the result be simplified to an array, if
#' possible? Note that this may lead to errors; if so, disable simplification.
#' @param returnList Logical: should the result be turned into a list (rather
#' than a multiData structure)? Note that this is incompatible with
#' simplification: if \code{mdaSimplify} is \code{TRUE}, this argument is
#' ignored.
#' @param mdaVerbose Integer specifying whether progress diagnistics should be
#' printed out. Zero means silent, increasing values will lead to more
#' diagnostic messages.
#' @param mdaIndent Integer specifying the indentation of the printed progress
#' messages. Each unit equals two spaces.
#' @return A multiData structure containing the results of the supplied
#' function on each \code{data} component in the input multiData structure.
#' Other components are simply copied.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure;
#' \code{\link{mtd.applyToSubset}} for applying a function to a subset of a
#' multiData structure; \code{\link{mtd.mapply}} for vectorizing over several
#' arguments.
#' @keywords misc
mtd.apply = function(
    # What to do
    multiData, FUN, ...,

    # Pre-existing results and update options
    mdaExistingResults = NULL, mdaUpdateIndex = NULL,
    mdaCopyNonData = FALSE,

    # Output formatting options
    mdaSimplify = FALSE,
    returnList = FALSE,

    # Internal behaviour options
    mdaVerbose = 0, mdaIndent = 0
) {
  printSpaces = indentSpaces(mdaIndent)

  if (!isMultiData(multiData, strict = FALSE))
    stop("Supplied 'multiData' is not a valid multiData structure.")

  if (mdaSimplify && mdaCopyNonData)
    stop("Non-data copying is not compatible with simplification.")

  nSets = length(multiData)
  if (mdaCopyNonData) out = multiData else out = vector(mode = "list", length = nSets)

  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex)

  FUN = match.fun(FUN)
  for (set in 1:nSets) {
    if (calculate[set]) {
      if (mdaVerbose > 0) {
        printFlush(paste0(printSpaces, "mtd.apply: working on set ", set))
      }
      out[[set]] = list(data = FUN(multiData[[set]]$data, ...))
    } else {
      out[set] = mdaExistingResults[set]
    }
  }

  names(out) = names(multiData)

  if (mdaSimplify) {
    if (mdaVerbose > 0) {
      printFlush(paste0(printSpaces, "mtd.apply: attempting to simplify..."))
    }
    return (mtd.simplify(out))
  } else if (returnList) {
    return(multiData2list(out))
  }

  out
}

#' @rdname mtd.apply
#' @param mdaRowIndex Index of rows to keep
#' @param mdaColIndex Index of columns to keep
#' @export
mtd.applyToSubset = function(
    # What to do
    multiData, FUN, ...,

    # Which rows and cols to keep
    mdaRowIndex = NULL, mdaColIndex = NULL,

    # Pre-existing results and update options
    mdaExistingResults = NULL, mdaUpdateIndex = NULL,
    mdaCopyNonData = FALSE,

    # Output formatting options
    mdaSimplify = FALSE,
    returnList = FALSE,

    # Internal behaviour options
    mdaVerbose = 0, mdaIndent = 0) {
  printSpaces = indentSpaces(mdaIndent)

  size = checkSets(multiData)
  if (mdaSimplify && mdaCopyNonData)
    stop("Non-data copying is not compatible with simplification.")

  if (mdaCopyNonData) res = multiData else res = vector(mode = "list", length = size$nSets)

  doSelection = FALSE
  if (!is.null(mdaColIndex))
  {
    doSelection = TRUE
    if (any(mdaColIndex < 0 | mdaColIndex > size$nGenes))
      stop("Some of the indices in 'mdaColIndex' are out of range.")
  } else {
    mdaColIndex = c(1:size$nGenes)
  }

  if (!is.null(mdaRowIndex)) {
    if (!is.list(mdaRowIndex)) {
       stop("mdaRowIndex must be a list, with one component per set.")
    }
    if (length(mdaRowIndex)!=size$nSets) {
       stop("Number of components in 'mdaRowIndex' must equal number of sets.")
    }
    doSelection = TRUE
  } else {
    mdaRowIndex = lapply(size$nSamples, function(n) { c(1:n) })
  }

  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex)

  fun = match.fun(FUN)
  for (set in 1:size$nSets) {
    if (calculate[set]) {
       if (mdaVerbose > 0) {
         printFlush(paste0(printSpaces, "mtd.applyToSubset: working on set ",
                           set))
       }
       res[[set]] = list(data = fun(
                ifelse(doSelection,
                       multiData[[set]]$data[mdaRowIndex[[set]],
                                             mdaColIndex, drop = FALSE],
                                 multiData[[set]]$data, ...)))
    } else {
       res[set] = mdaExistingResults[set]
    }
  }

  names(res) = names(multiData)

  if (mdaSimplify) {
    if (mdaVerbose > 0) {
      printFlush(paste0(printSpaces,
                        "mtd.applyToSubset: attempting to simplify..."))
        }
    return (mtd.simplify(res))
  } else if (returnList) {
    return (multiData2list(res))
  }

  return(res)
}

# mtd.simplify ####
#' If possible, simplify a multiData structure to a 3-dimensional array.
#'
#' This function attempts to put all \code{data} components into a
#' 3-dimensional array, with the last dimension corresponding to the sets. This
#' is only possible if all \code{data} components are matrices or data frames
#' with the same dimensiosn.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function assumes a "strict" multiData structure.
#'
#' @param multiData A multiData structure in the "strict" sense (see below).
#' @return A 3-dimensional array collecting all \code{data} components.
#' @note The function is relatively fragile and may fail. Use at your own risk.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure;
#'
#' \code{\link{multiData2list}} for converting multiData structures to plain
#' lists.
#' @keywords misc
mtd.simplify = function(multiData) {
  len = length(multiData[[1]]$data)
  dim = dim(multiData[[1]]$data)
  simplifiable = TRUE
  nSets = length(multiData)
  for (set in 1:nSets){
    if (len!=length(multiData[[set]]$data)) {
        simplifiable = FALSE
    }
    if (!isTRUE(all.equal( dim, dim(multiData[[set]]$data)))) {
        simplifiable = FALSE
    }
  }
  if (simplifiable) {
    if (is.null(dim)) {
       innerDim = len
       innerNames = names(multiData[[1]]$data)
       if (is.null(innerNames)) {
           innerNames = paste0("X", c(1:len))
       }
    } else {
       innerDim = dim
       innerNames = dimnames(multiData[[1]]$data)
       if (is.null(innerNames)) {
         innerNames = lapply(innerDim, function(x) {paste0("X", 1:x)})
       }
       nullIN = sapply(innerNames, is.null)
       if (any(nullIN)) {
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {
             paste0("X", 1:x)
         })
       }
    }
    setNames = names(multiData)
    if (is.null(setNames)) {
        setNames = paste0("Set_", 1:nSets)
    }
    mtd.s = matrix(NA, prod(innerDim), nSets)
    for (set in 1:nSets) {
      mtd.s[, set] = as.vector(multiData[[set]]$data)}

    dim(mtd.s) = c(innerDim, nSets)
    if (!is.null(innerNames))
      dimnames(mtd.s) = c(ifelse(is.list(innerNames), innerNames,
                                 list(innerNames)),
                          list(setNames))
    return(mtd.s)
  }
  return(multiData)
}

# isMultiData ####
#' Determine whether the supplied object is a valid multiData structure
#'
#' Attempts to determine whether the supplied object is a valid multiData
#' structure (see Details).
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function checks whether the supplied \code{x} is a multiData structure
#' in the "strict" (when \code{strict = TRUE} or "loose" \code{strict = FALSE}
#' sense.
#'
#' @param x An object.
#' @param strict Logical: should the structure of multiData be checked for
#' "strict" compliance?
#' @return Logical: \code{TRUE} if the input \code{x} is a multiData structure,
#' \code{FALSE} otherwise.
#' @author Peter Langfelder
#' @seealso Other multiData handling functions whose names start with
#' \code{mtd.}
#' @keywords misc
isMultiData = function(x, strict = TRUE) {
  if (strict) {
     !inherits(try(checkSets(x), silent = TRUE), 'try-error')
  } else {
    hasData = sapply(x, function(l) { "data" %in% names(l) })
    all(hasData)
  }
}

# mtd.mapply ####
#' Apply a function to elements of given multiData structures.
#'
#' Inspired by \code{\link{mapply}}, this function applies a given function to
#' each \code{data} component in the input multiData arguments, and optionally
#' simplify the result to an array if possible.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function applies the function \code{FUN} to each \code{data} component
#' of those arguments in \code{...} that are multiData structures in the
#' "loose" sense, and to each component of those arguments in \code{...} that
#' are not multiData structures.
#'
#' @param FUN Function to be applied.
#' @param \dots Arguments to be vectorized over. These can be multiData
#' structures or simple vectors (e.g., lists).
#' @param MoreArgs A named list that specifies the scalar arguments (if any) to
#' \code{FUN}.
#' @param mdma.argIsMultiData Optional specification whether arguments are
#' multiData structures. A logical vector where each component corresponds to
#' one entry of \code{...}. If not given, multiData status will be determined
#' using \code{\link{isMultiData}} with argument \code{strict=FALSE}.
#' @param mdmaExistingResults Optional list that contains previously calculated
#' results. This can be useful if only a few sets in \code{multiData} have
#' changed and recalculating the unchanged ones is computationally expensive.
#' If not given, all calculations will be performed. If given, components of
#' this list are copied into the output. See \code{mdmUpdateIndex} for which
#' components are re-calculated by default.
#' @param mdmaUpdateIndex Optional specification of which sets in
#' \code{multiData} the calculation should actually be carried out. This
#' argument has an effect only if \code{mdmaExistingResults} is non-NULL. If
#' the length of \code{mdmaExistingResults} (call the length `k') is less than
#' the number of sets in \code{multiData}, the function assumes that the
#' existing results correspond to the first `k' sets in \code{multiData} and
#' the rest of the sets are automatically calculated, irrespective of the
#' setting of \code{mdmaUpdateIndex}. The argument \code{mdmaUpdateIndex} can
#' be used to specify re-calculation of some (or all) of the results that
#' already exist in \code{mdmaExistingResults}.
#' @param mdmaSimplify Logical: should simplification of the result to an array
#' be attempted? The simplification is fragile and can produce unexpected
#' errors; use the default \code{FALSE} if that happens.
#' @param returnList Logical: should the result be turned into a list (rather
#' than a multiData structure)? Note that this is incompatible with
#' simplification: if \code{mdaSimplify} is \code{TRUE}, this argument is
#' ignored.
#' @param mdma.doCollectGarbage Should garbage collection be forced after each
#' application of \code{FUN}?
#' @param mdmaVerbose Integer specifying whether progress diagnistics should be
#' printed out. Zero means silent, increasing values will lead to more
#' diagnostic messages.
#' @param mdmaIndent Integer specifying the indentation of the printed progress
#' messages. Each unit equals two spaces.
#' @return A multiData structure containing (as the \code{data} components) the
#' results of \code{FUN}. If simplification is successful, an array instead.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure;
#'
#' \code{multiData.apply} for application of a function to a single multiData
#' structure.
#' @keywords misc
mtd.mapply = function(
  # What to do
  FUN, ..., MoreArgs = NULL,

  # How to interpret the input
  mdma.argIsMultiData = NULL,

  # Copy previously known results?
  mdmaExistingResults = NULL, mdmaUpdateIndex = NULL,

  # How to format output
  mdmaSimplify = FALSE,
  returnList = FALSE,

  # Options controlling internal behaviour
  mdma.doCollectGarbage = FALSE,
  mdmaVerbose = 0, mdmaIndent = 0) {
  printSpaces = indentSpaces(mdmaIndent)
  dots = list(...)
  if (length(dots)==0)
    stop("No arguments were specified. Please type ?mtd.mapply to see the help page.")
  dotLengths = sapply(dots, length)
  if (any(dotLengths!=dotLengths[1]))
    stop(paste0("All arguments to vectorize over must have the same length.\n",
                "Scalar arguments should be put into the 'MoreArgs' argument.\n",
                "Note: lengths of '...' arguments are: ", paste(dotLengths, collapse = ", ")))
  nArgs = length(dots)
  res = list()
  if (is.null(mdma.argIsMultiData)) mdma.argIsMultiData = sapply(dots, isMultiData, strict = FALSE)

  nSets = dotLengths[1]

  calculate = .calculateIndicator(nSets, mdmaExistingResults, mdmaUpdateIndex)

  FUN = match.fun(FUN)
  for (set in 1:nSets)
  {
    if (calculate[set])
    {
      if (mdmaVerbose > 0)
        printFlush(paste0(printSpaces, "mtd.mapply: working on set ", set))

      localArgs = list()
      for (arg in 1:nArgs)
        localArgs[[arg]] = if (mdma.argIsMultiData[arg]) dots[[arg]] [[set]] $ data else dots[[arg]] [[set]]
      names(localArgs) = names(dots)
      res[[set]] = list(data = do.call(FUN, c(localArgs, MoreArgs)))
      if (mdma.doCollectGarbage) collectGarbage()
    } else
      res[set] = mdmaExistingResults[set]
  }

  names(res) = names(dots[[1]])

  if (mdmaSimplify)
  {
    if (mdmaVerbose > 0)
      printFlush(paste0(printSpaces, "mtd.mapply: attempting to simplify..."))
    return (mtd.simplify(res))
  } else if (returnList) {
    return (multiData2list(res))
  }

  return(res)
}

# mtd.rbindSelf ####
#' Turn a multiData structure into a single matrix or data frame.
#'
#' This function "rbinds" the \code{data} components of all sets in the input
#' into a single matrix or data frame.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function requires a "strict" multiData structure.
#'
#' @param multiData A multiData structure.
#' @return A single matrix or data frame containing the "rbinded" result.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure;
#'
#' \code{\link{rbind}} for various subtleties of the row binding operation.
#' @keywords misc
mtd.rbindSelf = function(multiData){
  size = checkSets(multiData)
  out = NULL
  colnames = mtd.colnames(multiData)
  for (set in 1:size$nSets)
  {
    if (!is.null(colnames(multiData[[set]]$data)) &&
        !isTRUE(all.equal(colnames, colnames(multiData[[set]]$data))) )
          colnames(multiData[[set]]$data) = colnames
    out = rbind(out, multiData[[set]]$data)
  }
  out
}

# mtd.setAttr ####
#' Set attributes on each component of a multiData structure
#'
#' Set attributes on each \code{data} component of a multiData structure
#'
#'
#' @param multiData A multiData structure.
#' @param attribute Name for the attribute to be set
#' @param valueList List that gives the attribute value for each set in the
#' multiData structure.
#' @return The input multiData with the attribute set on each \code{data}
#' component.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure;
#'
#' \code{isMultiData} for a description of the multiData structure.
#' @keywords misc
mtd.setAttr = function(multiData, attribute, valueList) {
  size = checkSets(multiData)
  ind = 1
  for (set in 1:size$nSets){
    attr(multiData[[set]]$data, attribute) = valueList[[ind]]
    ind = ind + 1
    if (ind > length(valueList)) ind = 1
  }
  multiData
}

# mtd.setColnames ####
#' Get and set column names in a multiData structure.
#'
#' Get and set column names on each \code{data} component in a multiData
#' structure.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' The \code{mtd.colnames} and \code{mtd.setColnames} assume (and checks for) a
#' "strict" multiData structure.
#'
#' @aliases mtd.setColnames mtd.colnames
#' @param multiData A multiData structure
#' @param colnames A vector (coercible to character) of column names.
#' @return \code{mtd.colnames} returns the vector of column names of the
#' \code{data} component. The function assumes the column names in all sets are
#' the same.
#'
#' \code{mtd.setColnames} returns the multiData structure with the column names
#' set in all \code{data} components.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData}} to create a multiData structure.
#' @keywords misc
mtd.setColnames = function(multiData, colnames) {
  size = checkSets(multiData)
  for (set in 1:size$nSets)
    colnames(multiData[[set]]$data) = colnames
  multiData
}

#' @rdname mtd.setColnames
#' @export
mtd.colnames = function(multiData) {
    colnames(multiData[[1]]$data)
}

# multiData ####
#' Create a multiData structure.
#'
#' This function creates a multiData structure by storing its input arguments
#' as the 'data' components.
#'
#' A multiData structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiData structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiData structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' @param \dots Arguments to be stored in the multiData structure.
#' @return The resulting multiData structure.
#' @author Peter Langfelder
#' @seealso \code{\link{multiData2list}} for converting a multiData structure
#' to a list; \code{\link{list2multiData}} for an alternative way of creating a
#' multiData structure; \code{\link{mtd.apply}, \link{mtd.applyToSubset},
#' \link{mtd.mapply}} for ways of applying a function to each component of a
#' multiData structure.
#' @keywords misc
#' @examples
#'
#' data1 = matrix(rnorm(100), 20, 5);
#' data2 = matrix(rnorm(50), 10, 5);
#'
#' md = multiData(Set1 = data1, Set2 = data2);
#'
#' checkSets(md)
#'
multiData = function(...) {
  list2multiData(list(...))
}

# KeepCommonProbes ####
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
    if (nSets <= 0) {
        stop("No expression data given!")
    }

    Names = data.frame(Names = names(multiExpr[[orderBy]]$data))

    if (nSets > 1) {
        for (set in (1:nSets)) {
            SetNames = data.frame(Names = names(multiExpr[[set]]$data),
                                  index = c(1:dim(multiExpr[[set]]$data)[2]))
            Names = merge(Names, SetNames, by.x = "Names", by.y = "Names",
                          all = FALSE, sort = FALSE)
        }
    }

    for (set in 1:nSets) {
        multiExpr[[set]]$data = multiExpr[[set]]$data[, Names[, set + 1]]
    }

    multiExpr
}
