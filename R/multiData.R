# Convenience functions for manipulating multiSet structures

# Note: many of these function would be simpler to use if I used some sort of
# class/method technique to keep track of the class of each object internally.
# For example, I could then write a generic function "subset" that would work
# consistently on lists and multiSet objects. Similarly, multiSet2list would
# simply become a method of as.list, and as.list would be safe to use both on
# lists and on multiSet objects.
# Help on this: http://www.programiz.com/r-programming/S3-class
# https://abhishek-tiwari.com/hacking/class-and-objects-in-r-s3-style

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
#' @examples
#' data1 <- matrix(rnorm(100L), 20L, 5L)
#' data2 <- matrix(rnorm(50L), 10L, 5L)
#' colnames(data1) <- LETTERS[1:5]
#' colnames(data2) <- LETTERS[2:6]
#' multiUnion(list(data1, data2))
#' multiIntersect(list(data1, data2))
multiUnion <- function(setList) {
    len = length(setList)
    if (len == 0) {
        return(NULL)
    }
    if (len == 1) {
        return(setList[[1]])
    }

    out <- setList[[1]]
    for (elem in 2:len) {
        out <- union(out, setList[[elem]])
    }
    out
}

#' @export
#' @rdname multiUnion
multiIntersect <- function(setList) {
    len <- length(setList)
    if (len == 0) {
        return(NULL)
    }
    if (len == 1) {
        return(setList[[1]])
    }

    out <- setList[[1]]
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
#' singleSetData <- matrix(rnorm(100), 10,10)
#' encapsData <- fixDataStructure(singleSetData)
#' length(encapsData)
#' names(encapsData[[1]])
#' dim(encapsData[[1]]$data)
#' all.equal(encapsData[[1]]$data, singleSetData)
#'
#' @export
fixDataStructure <- function(data, verbose = 0, indent = 0) {
    spaces = indentSpaces(indent)
    if ((class(data) != "list") || (class(data[[1]]) != "list")) {
        if (verbose > 0) {
            printFlush(paste(spaces,
                             "fixDataStructure: data is not a vector of lists:",
                             "converting it into one."))
        }
        x <- data
        data <- vector(mode = "list", length = 1)
        data[[1]] <- list(data = x)
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
#' @return A list with components
#' \item{nSets}{Number of sets (length of the vector \code{data}).}
#' \item{nGenes}{Number of columns in the \code{data} components in the lists.
#' This number must be the same for all sets.}
#' \item{nSamples}{A vector of length \code{nSets} giving the number of rows in
#' the \code{data} components.}
#' \item{structureOK}{Only set if the argument \code{checkStructure} equals
#' \code{TRUE}.  The value is \code{TRUE} if the paramter \code{data} passes a
#' few tests of its structure, and \code{FALSE} otherwise. The tests are not
#' exhaustive and are meant to catch obvious user errors rather than be
#' bulletproof.}
#' @examples
#' data1 <- matrix(rnorm(100L), 20L, 5L)
#' data2 <- matrix(rnorm(50L), 10L, 5L)
#' md <- multiSet(Set1 = data1, Set2 = data2)
#' checkSets(md)
#' @author Peter Langfelder, \email{Peter.Langfelder@@gmail.com}
#' @keywords misc
checkSets <- function(data, checkStructure = FALSE, useSets = NULL) {
    nSets = length(data)
    if (is.null(useSets)) {
        useSets <- c(1:nSets)
    }
    if (nSets <= 0) {
        stop("No data given.")
    }
    structureOK <- TRUE
    if ("multiSet" %in% class(data)) {
        if (checkStructure) {
            structureOK <- FALSE
            nGenes <- 0
            nSamples <- 0
        } else {

            nSamples <- vector(length = nSets)
            nGenes <- dim(data[[useSets[1]]]$data)[2]
            for (set in useSets) {
                if (nGenes != dim(data[[set]]$data)[2]) {
                    if (checkStructure) {
                        structureOK <- FALSE
                    } else {
                        stop("Incompatible number of genes in set 1 and ", set)
                    }
                }
                nSamples[set] <- dim(data[[set]]$data)[1]
            }
        }
        list(nSets = nSets, nGenes = nGenes, nSamples = nSamples,
             structureOK = structureOK)
    } else {
        stop("Not in multiSet format")
    }
}

# subset.multiSet ####


#' Subset rows and columns in a multiSet structure
#'
#' The function restricts each \code{data} component to the given columns and
#' rows.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function assumes a "strict" multiSet structure unless \code{permissive}
#' is \code{TRUE}.
#'
#' @aliases subset.multiSet subset subset.default [.multiSet subset.multiSet
#' @param multiSet A multiSet structure.
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
#' \code{multiSet} input? Note that the subsetting may lead to cryptic errors
#' if the input \code{multiSet} does not follow the "strict" format.
#' @param drop logical: should dimensions with extent 1 be dropped?
#' @param \dots Other arguments
##' @param x MultiSet Object to subset
##' @param i Columns to subset of the object
##' @param j Rows to subset of the object
#' @return A multiSet structure containing the selected rows and columns. Note
#' that result always retaines it dimension and other attributes.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet}} to create a multiSet structure.
#' @keywords misc
#' @examples
#'
#'
#' data1 <- matrix(rnorm(100), 20, 5)
#' data2 <- matrix(rnorm(50), 10, 5)
#'
#' md <- multiSet(Set1 = data1, Set2 = data2)
#' md[, c(1,3)]
#' md[c(1,3), ]
#'
subset.multiSet <- function(multiSet, rowIndex = NULL, colIndex = NULL,
                      permissive = FALSE, drop = FALSE, ...) {
  size <- checkSets(multiSet, checkStructure = permissive)
  if (!size$structureOK && !is.null(colIndex)) {
    warning(immediate. = TRUE,
            "subset.multiSet: applying column selection on data sets that do ",
            "not have\n the same number of columns. This is treacherous ",
            "territory; proceed with caution.")
  }
  if (is.null(colIndex)) {
      colIndex.1 <- c(1:size$nGenes)
  } else {
      colIndex.1 <- colIndex
  }
  if (is.null(rowIndex)) {
      rowIndex <- lapply(size$nSamples, function(n) {c(1:n)})
  }
  if (length(rowIndex)!=size$nSets) {
    stop("If given, 'rowIndex' must be a list of the same length as 'multiSet'")
  }
  out <- list()
  for (set in 1:size$nSets) {
    if (permissive)
      if (is.null(colIndex)) {
          colIndex.1 <- c(1:ncol(multiSet[[set]]$data))
      } else {
          colIndex.1 <- colIndex
      }
    if (is.character(colIndex.1)) {
      colIndex.1 <- match(colIndex.1, colnames(multiSet[[set]]$data))
      n1 <- length(colIndex.1)
      if (any(is.na(colIndex.1))) {
        stop("Cannot match the following entries in 'colIndex' to column ",
             "names in set ", set, ":\n",
             paste(colIndex[is.na(colIndex.1)][1:min(n1, 5)],
                    collapse = ", "),
             ifelse(n1 > 5, ", ... [output truncated]", ""))
          }
    }
    out[[set]] <- list(data = multiSet[[set]]$data[rowIndex[[set]],
                                                  colIndex.1, drop = drop])
  }
  names(out) <- names(multiSet)
  out
}

#' @export
#' @rdname subset.multiSet
subset <- function(multiSet, ...) {
    UseMethod("subset", multiSet)
}

#' @export
#' @rdname subset.multiSet
subset.default <- function(multiSet, ...) {
    subset(multiSet, ...)
}

# Method to allow subseting the data of a multiSet
#' @export
#' @rdname subset.multiSet
#' @aliases  subset.multiSet
#' @param i Columns to subset of the object
#' @param j Rows to subset of the object
#' @param x MultiSet Object to subset
`[.multiSet` <- function(x, i, j, drop = TRUE) {
    y <- list(names(x))
    for (set in 1:length(x)){
        y[[set]] <- list(data = x[[set]]$data[i, j, drop = drop])
    }
    names(y) <- names(x)
    class(y) <- c("multiSet", "list")
    y
}

# list2multiSet ####


#' Convert a list to a multiSet structure and vice-versa.
#'
#' \code{list2multiSet} converts a list to a multiSet structure
#' \code{multiSet2list} does the inverse.
#'
#' A multiSet structure is a vector of lists (one list for each set) where each
#' list has a component \code{data} containing some useful information.
#'
#' @aliases list2multiSet multiSet2list as.multiSet as.multiSet.default.list
#' as.multiSet.default multiSet2list
#' @param data A list to be converted to a multiSet structure.
#' @param x Object to convert to multiSet or to list
#' @param \dots Other arguments
#' @param multiSet A multiSet structure to be converted to a list.
#' @return For \code{list2multiSet}, a multiSet structure; for
#' \code{multiSet2list}, the corresponding list.
#' @author Peter Langfelder
#' @keywords misc
#' @export list2multiSet
list2multiSet <- function(data) {
  out = list()
  for (set in 1:length(data)) {
    out[[set]] <- list(data = data[[set]])
  }
  names(out) <- names(data)
  class(out) <- c("multiSet", "list")
  out
}



#' Convert a list to a multiSet
#'
#' Convert a list to a multiSet
#'
#'
#' @param multiSet list to be converted to a multiSet object, or a multiSet to
#' be converted to a list
#' @param \dots other arguments to be passed to \code{list2multiSet} or
#' \code{list2multiset}
#' @return An object of class multiSet or a list.
#' @seealso \code{\link{list2multiSet}}, \code{\link{multiSet2list}}
as.multiSet.list <- function(multiSet) {
    list2multiSet(multiSet)
}

#' @export
#' @rdname list2multiSet
as.multiSet <- list2multiSet

#' @export
#' @method as.multiSet.default list
#' @rdname list2multiSet
as.multiSet.default.list <- function(x) {
    list2multiSet(x)
}

#' @export
#' @rdname list2multiSet
as.multiSet.default <- function(x, ...) {
    list2multiSet(x, ...)
}


#' Convert a multiSet object to a list
#'
#' Each dataset is a list with a data element with the expression.
#'
#'
#' @param x muiltiSet object to be converted to list
#' @param \dots other arguments passed to multiSet2list
as.list.multiSet <- function(x, ...) {
    multiSet2list(x, ...)
}

#' @rdname list2multiSet
#' @aliases multiSet2list
#' @param multiSet A multiSet structure to be converted to a list.
#' @export
multiSet2list <- function(multiSet) {
    out <- vector("list", length(multiSet))
    names(out) <- names(multiSet)
    for (i in 1:length(multiSet)) {
        out[[i]] <- multiSet[[i]]$data
    }
    return(out)
}

# apply.multiSet ####
.calculateIndicator <- function(nSets, mdaExistingResults, mdaUpdateIndex) {
  if (length(mdaUpdateIndex) == 0) {
      mdaUpdateIndex <- NULL
  }
  calculate = rep(TRUE, nSets)
  if (!is.null(mdaExistingResults)) {
    nSets.existing = length(mdaExistingResults)
    if (nSets.existing>nSets){
      stop("Number of sets in 'mdaExistingResults' is higher than the number ",
           "of sets in 'multiSet'.\n  Please supply a valid ",
           "'mdaExistingResults' or NULL to recalculate all results.")
    }
    if (nSets.existing == 0) {
      stop("Number of sets in 'mdaExistingResults' is zero.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to ",
           "recalculate all results.")
    }
    if (is.null(mdaUpdateIndex)) {
      calculate[1:length(mdaExistingResults)] = FALSE
    } else {
      if (any(! mdaUpdateIndex %in% c(1:nSets))) {
        stop("All entries in 'mdaUpdateIndex' must be between 1 and the number",
             " of sets in 'multiSet'.")
      }
      calculateIndex <- sort(unique(c(mdaUpdateIndex,
                                      ifelse(nSets.existing<nSets,
                                             c((nSets.existing+1):nSets),
                                             NULL))))
      calculate[ c(1:nSets)[-calculateIndex] ] = FALSE
    }
  }

  calculate
}

#' Apply a function to each set in a multiSet structure.
#'
#' Inspired by \code{\link{lapply}}, these functions apply a given function to
#' each \code{data} component in the input \code{multiSet} structure, and
#' optionally simplify the result to an array if possible.
#'
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' \code{apply.multiSet} works on any "loose" multiSet structure
# #' \code{apply.multiSetToSubset} assumes (and checks for) a "strict" multiSet
# #' structure.
#'
# #' @aliases apply
#' @param multiSet A multiSet structure to apply the function over
#' @param FUN Function to be applied.
#' @param \dots Other arguments to the function \code{FUN}.
# #' @param mdaExistingResults Optional list that contains previously calculated
# #' results. This can be useful if only a few sets in \code{multiSet} have
# #' changed and recalculating the unchanged ones is computationally expensive.
# #' If not given, all calculations will be performed. If given, components of
# #' this list are copied into the output. See \code{mdmUpdateIndex} for which
# #' components are re-calculated by default.
# #' @param mdaUpdateIndex Optional specification of which sets in
# #' \code{multiSet} the calculation should actually be carried out. This
# #' argument has an effect only if \code{mdaExistingResults} is non-NULL. If the
# #' length of \code{mdaExistingResults} (call the length `k') is less than the
# #' number of sets in \code{multiSet}, the function assumes that the existing
# #' results correspond to the first `k' sets in \code{multiSet} and the rest of
# #' the sets are automatically calculated, irrespective of the setting of
# #' \code{mdaUpdateIndex}. The argument \code{mdaUpdateIndex} can be used to
# #' specify re-calculation of some (or all) of the results that already exist in
# #' \code{mdaExistingResults}.
# #' @param mdaCopyNonData Logical: should non-data components of
# #' \code{multiSet} be copied into the output? Note that the copying is
# #' incompatible with simplification; enabling both will trigger an error.
# #' @param mdaSimplify Logical: should the result be simplified to an array, if
# #' possible? Note that this may lead to errors; if so, disable simplification.
# #' @param returnList Logical: should the result be turned into a list (rather
# #' than a multiSet structure)? Note that this is incompatible with
# #' simplification: if \code{mdaSimplify} is \code{TRUE}, this argument is
# #' ignored.
# #' @param mdaVerbose Integer specifying whether progress diagnistics should be
# #' printed out. Zero means silent, increasing values will lead to more
# #' diagnostic messages.
# #' @param mdaIndent Integer specifying the indentation of the printed progress
# #' messages. Each unit equals two spaces.


#' Apply a function to each set in a multiSet structure.
#'
#' Inspired by \code{\link{lapply}}, these functions apply a given function to
#' each \code{data} component in the input \code{multiSet} structure, and
#' optionally simplify the result to an array if possible.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' \code{apply.multiSet} works on any "loose" multiSet structure
#'
#' @aliases apply.multiSet apply apply.default
#' @param multiSet A multiSet structure to apply the function over
#' @param MARGIN Column or row you want to use
#' @param FUN Function to be applied.
#' @param \dots Other arguments to the function \code{FUN}.
#' @return A multiSet structure containing the results of the supplied function
#' on each \code{data} component in the input multiSet structure. Other
#' components are simply copied.
#' @author Lluís Revilla
#' @seealso \code{\link{multiSet}} to create a multiSet structure
#' \code{\link{multiSet.mapply}} for vectorizing over several arguments.
#' @keywords misc
apply.multiSet <- function(multiSet, MARGIN, FUN, ...) {

  if (!is.multiSet(multiSet, strict = FALSE)) {
    stop("Supplied 'multiSet' is not a valid multiSet structure.")
  }

  nSets <- length(multiSet)
  out <- vector(mode = "list", length = nSets)

  FUN <- match.fun(FUN)
  for (set in 1:nSets) {
      data <- base::apply(multiSet[[set]]$data, MARGIN = MARGIN,
                    FUN = FUN, ...)
      out[[set]] <- list(data = data)
  }

  names(out) <- names(multiSet)

  out
}

#' @export
#' @rdname apply.multiSet
apply <- function(multiSet, ...) {
    UseMethod("apply", multiSet)
}

#' @export
#' @rdname apply.multiSet
apply.default <- function(multiSet, MARGIN, FUN, ...) {
    base::apply(X = multiSet, MARGIN = MARGIN, FUN = FUN, ...)
}

multiSet.simplify = function(multiSet) {
  len = length(multiSet[[1]]$data)
  dim = dim(multiSet[[1]]$data)
  simplifiable = TRUE
  nSets = length(multiSet)
  for (set in 1:nSets) {
    if (len!=length(multiSet[[set]]$data)) simplifiable = FALSE
    if (!isTRUE(all.equal( dim, dim(multiSet[[set]]$data)))) simplifiable = FALSE
  }
  if (simplifiable) {
    if (is.null(dim)) {
       innerDim = len
       innerNames = names(multiSet[[1]]$data)
       if (is.null(innerNames)) innerNames = paste0("X", c(1:len))
    } else {
       innerDim = dim
       innerNames = dimnames(multiSet[[1]]$data)
       if (is.null(innerNames))
         innerNames = lapply(innerDim, function(x) {paste0("X", 1:x)})
       nullIN = sapply(innerNames, is.null)
       if (any(nullIN))
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {paste0("X", 1:x)})
    }
    setNames = names(multiSet)
    if (is.null(setNames)) setNames = paste0("Set_", 1:nSets)
    mtd.s = matrix(NA, prod(innerDim), nSets)
    for (set in 1:nSets)
      mtd.s[, set] = as.vector(multiSet[[set]]$data)

    dim(mtd.s) = c(innerDim, nSets)
    if (!is.null(innerNames))
      dimnames(mtd.s) = c (if (is.list(innerNames)) innerNames else list(innerNames), list(setNames))
    return(mtd.s)
  }
  return(multiSet)
}

#' @export
#' @rdname apply.multiSet
apply.default <- function(multiSet, MARGIN, FUN, ...) {
    base::apply(X = multiSet, MARGIN = MARGIN, FUN = FUN, ...)
}


# is.multiSet ####


#' Determine whether the supplied object is a valid multiSet structure
#'
#' Attempts to determine whether the supplied object is a valid multiSet
#' structure (see Details).
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function checks whether the supplied \code{x} is a multiSet structure
#' in the "strict" (when \code{strict = TRUE} or "loose" \code{strict = FALSE}
#' sense.
#'
#' @param x An object.
#' @param strict Logical: should the structure of multiSet be checked for
#' "strict" compliance?
#' @return Logical: \code{TRUE} if the input \code{x} is a multiSet structure,
#' \code{FALSE} otherwise.
#' @author Peter Langfelder
#' @seealso Other multiSet handling functions whose names start with
#' \code{multiSet.}
#' @keywords misc
#' @export is.multiSet
is.multiSet <- function(x, strict = TRUE) {
    if ("multiSet" %in% class(x)) {
        return(TRUE)
    }
    if (strict) {
        !inherits(try(checkSets(x), silent = TRUE), 'try-error')
    } else {
        hasData <- sapply(x, function(l) { "data" %in% names(l) })
        all(hasData)
    }
}

# multiSet.mapply ####


#' Apply a function to elements of given multiSet structures.
#'
#' Inspired by \code{\link{mapply}}, this function applies a given function to
#' each \code{data} component in the input multiSet arguments, and optionally
#' simplify the result to an array if possible.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function applies the function \code{FUN} to each \code{data} component
#' of those arguments in \code{...} that are multiSet structures in the "loose"
#' sense, and to each component of those arguments in \code{...} that are not
#' multiSet structures.
#'
#' @param FUN Function to be applied.
#' @param \dots Arguments to be vectorized over. These can be multiSet
#' structures or simple vectors (e.g., lists).
#' @param MoreArgs A named list that specifies the scalar arguments (if any) to
#' \code{FUN}.
#' @param mdma.argIsMultiData Optional specification whether arguments are
#' multiSet structures. A logical vector where each component corresponds to
#' one entry of \code{...}. If not given, multiSet status will be determined
#' using \code{\link{is.multiSet}} with argument \code{strict=FALSE}.
#' @param mdmaExistingResults Optional list that contains previously calculated
#' results. This can be useful if only a few sets in \code{multiSet} have
#' changed and recalculating the unchanged ones is computationally expensive.
#' If not given, all calculations will be performed. If given, components of
#' this list are copied into the output. See \code{mdmUpdateIndex} for which
#' components are re-calculated by default.
#' @param mdmaUpdateIndex Optional specification of which sets in
#' \code{multiSet} the calculation should actually be carried out. This
#' argument has an effect only if \code{mdmaExistingResults} is non-NULL. If
#' the length of \code{mdmaExistingResults} (call the length `k') is less than
#' the number of sets in \code{multiSet}, the function assumes that the
#' existing results correspond to the first `k' sets in \code{multiSet} and the
#' rest of the sets are automatically calculated, irrespective of the setting
#' of \code{mdmaUpdateIndex}. The argument \code{mdmaUpdateIndex} can be used
#' to specify re-calculation of some (or all) of the results that already exist
#' in \code{mdmaExistingResults}.
#' @param mdmaSimplify Logical: should simplification of the result to an array
#' be attempted? The simplification is fragile and can produce unexpected
#' errors; use the default \code{FALSE} if that happens.
#' @param returnList Logical: should the result be turned into a list (rather
#' than a multiSet structure)? Note that this is incompatible with
#' simplification: if \code{mdaSimplify} is \code{TRUE}, this argument is
#' ignored.
#' @param mdma.doCollectGarbage Should garbage collection be forced after each
#' application of \code{FUN}?
#' @param mdmaVerbose Integer specifying whether progress diagnistics should be
#' printed out. Zero means silent, increasing values will lead to more
#' diagnostic messages.
#' @param mdmaIndent Integer specifying the indentation of the printed progress
#' messages. Each unit equals two spaces.
#' @return A multiSet structure containing (as the \code{data} components) the
#' results of \code{FUN}. If simplification is successful, an array instead.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet}} to create a multiSet structure
#'
#' \code{apply.multiSet} for application of a function to a single multiSet
#' structure.
#' @keywords misc
#' @examples
#'
#' data1 <- matrix(rnorm(100L), 20L, 5L)
#' data2 <- matrix(rnorm(50L), 10L, 5L)
#' colnames(data1) <- LETTERS[1:5]
#' colnames(data2) <- LETTERS[2:6]
#' md <- list2multiSet(list(Set1 = data1, Set2 = data2))
#' multiSet.mapply(md, FUN = sum)
#'
#' @export multiSet.mapply
multiSet.mapply <- function(FUN, ..., MoreArgs = NULL,
                            # How to interpret the input
                            mdma.argIsMultiData = NULL,
                            # Copy previously known results?
                            mdmaExistingResults = NULL, mdmaUpdateIndex = NULL,
                            # How to format output
                            mdmaSimplify = FALSE, returnList = FALSE,
                            # Options controlling internal behaviour
                            mdma.doCollectGarbage = FALSE, mdmaVerbose = 0,
                            mdmaIndent = 0) {
    printSpaces <- indentSpaces(mdmaIndent)
    dots <- list(...)
    if (length(dots) == 0) {
        stop("No arguments were specified. Please type ?multiSet.mapply to see",
             " the help page.")
    }
    dotLengths <- sapply(dots, length)
    if (any(dotLengths != dotLengths[1])) {
        stop("All arguments to vectorize over must have the same length.\n",
             "Scalar arguments should be put into the 'MoreArgs' argument.\n",
             "Note: lengths of '...' arguments are: ",
             paste(dotLengths, collapse = ", "))
    }
    nArgs <- length(dots)
    res <- list()
    if (is.null(mdma.argIsMultiData)) {
        mdma.argIsMultiData <- sapply(dots, is.multiSet, strict = FALSE)
    }
    nSets <- dotLengths[1]

    calculate <- .calculateIndicator(nSets, mdmaExistingResults, mdmaUpdateIndex)

    FUN <- match.fun(FUN)
    for (set in 1:nSets) {
        if (calculate[set]) {
            if (mdmaVerbose > 0) {
                printFlush(paste0(printSpaces, "working on set ", set))
            }
            localArgs <- list()
            for (arg in 1:nArgs) {
                localArgs[[arg]] <- ifelse(mdma.argIsMultiData[arg],
                                          dots[[arg]][[set]]$data,
                                          dots[[arg]][[set]])
            }
            names(localArgs) <- names(dots)
            res[[set]] <- list(data = do.call(FUN, c(localArgs, MoreArgs)))
            if (mdma.doCollectGarbage) {
                collectGarbage()
            }
        } else
            res[set] <- mdmaExistingResults[set]
    }

    names(res) <- names(dots[[1]])

    if (mdmaSimplify) {
        if (mdmaVerbose > 0)
            printFlush(paste0(printSpaces,
                              "multiSet.mapply: attempting to simplify..."))
    return(simplify(res))
  } else if (returnList) {
    return (multiSet2list(res))
  }
  return(res)
}



#' If possible, simplify a multiSet structure to a 3-dimensional array
#'
#' This function attempts to put all \code{data} components into a
#' 3-dimensional array, with the last dimension corresponding to the sets. This
#' is only possible if all \code{data} components are matrices or data frames
#' with the same dimension.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' @param multiSet Object of class multiSet
#' @return A 3-dimensional array collecting all \code{data} components.
#' @note The function is relatively fragile and may fail. Use at your own risk.
#' @seealso \code{\link{multiSet}} to create a multiSet structure
#' \code{\link{multiSet2list}} for converting multiSet structures to plain
#' lists.
simplify <- function(multiSet) {
        len <- length(multiSet[[1]]$data)
        dim <- dim(multiSet[[1]]$data)
        simplifiable <- TRUE
        nSets <- length(multiSet)
        for (set in 1:nSets) {
            if (len != length(multiSet[[set]]$data)) {
                simplifiable <- FALSE
                }
            if (!isTRUE(all.equal( dim, dim(multiSet[[set]]$data)))) {
                simplifiable <- FALSE
            }
        }
        if (simplifiable) {
            if (is.null(dim)) {
                innerDim <- len
                innerNames <- names(multiSet[[1]]$data)
                if (is.null(innerNames)) {
                    innerNames <- paste0("X", c(1:len))
                }
            } else {
                innerDim <- dim
                innerNames <- dimnames(multiSet[[1]]$data)
                if (is.null(innerNames)) {
                    innerNames <- lapply(innerDim, function(x) {
                        paste0("X", 1:x)
                    })
                }
                nullIN <- sapply(innerNames, is.null)
                if (any(nullIN)) {
                    innerNames[nullIN] <- lapply(innerDim[nullIN], function(x) {
                        paste0("X", 1:x)
                    })
                }
            }
            setNames <- names(multiSet)
            if (is.null(setNames)) {
                setNames <- paste0("Set_", 1:nSets)
            }
            mtd.s <- matrix(NA, prod(innerDim), nSets)
            for (set in 1:nSets) {
                mtd.s[, set] <- as.vector(multiSet[[set]]$data)
            }

            dim(mtd.s) <- c(innerDim, nSets)
            if (!is.null(innerNames))
                dimnames(mtd.s) <- c(ifelse(is.list(innerNames),
                                           innerNames,
                                           list(innerNames)),
                                    list(setNames))
            return(mtd.s)
        }
        return(multiSet)

}
# multiSet.rbindSelf ####


#' Turn a multiSet structure into a single matrix or data frame.
#'
#' This function "rbinds" the \code{data} components of all sets in the input
#' into a single matrix or data frame.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' This function requires a "strict" multiSet structure.
#'
#' @param multiSet A multiSet structure.
#' @return A single matrix or data frame containing the "rbinded" result.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet}} to create a multiSet structure
#'
#' \code{\link{rbind}} for various subtleties of the row binding operation.
#' @keywords misc
#' @aliases multiSet.rbindSelf
#' @export multiSet.rbindSelf
multiSet.rbindSelf <- function(multiSet){
  size <- checkSets(multiSet)
  out <- NULL
  colnames <- colnames(multiSet)
  for (set in 1:size$nSets) {
    if (!is.null(colnames(multiSet[[set]]$data)) &&
        !isTRUE(all.equal(colnames, colnames(multiSet[[set]]$data))) ) {
          colnames(multiSet[[set]]$data) <- colnames
    }
    out <- rbind(out, multiSet[[set]]$data)
  }
  out
}

# multiSet.setAttr ####


#' Set attributes on each component of a multiSet structure
#'
#' Set attributes on each \code{data} component of a multiSet structure
#'
#'
#' @param multiSet A multiSet structure.
#' @param attribute Name for the attribute to be set
#' @param valueList List that gives the attribute value for each set in the
#' multiSet structure.
#' @return The input multiSet with the attribute set on each \code{data}
#' component.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet}} to create a multiSet structure
#'
#' \code{is.multiSet} for a description of the multiSet structure.
#' @keywords misc
#' @export multiSet.setAttr
multiSet.setAttr <- function(multiSet, attribute, valueList) {
  size = checkSets(multiSet)
  ind <- 1
  for (set in 1:size$nSets){
    attr(multiSet[[set]]$data, attribute) <- valueList[[ind]]
    ind <- ind + 1
    if (ind > length(valueList)) {
        ind <- 1
    }
  }
  multiSet
}

# multiSet.colnames ####


#' Get and set column names in a multiSet structure.
#'
#' Get and set column names on each \code{data} component in a multiSet
#' structure.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' The \code{colnames.multiSet} assume (and checks for) a"strict" multiSet
#' structure.
#'
#' @param multiSet A multiSet structure
#' @param value A vector (coercible to character) of column names.
#' @param \dots Other arguments passed
#' @return The multiSet structure with the column names set in all \code{data}
#' components.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet}} to create a multiSet structure.
#' @keywords misc
`colnames<-.multiSet` <- function(multiSet, value) {
  size <- checkSets(multiSet)
  for (set in 1:size$nSets) {
      # warning(dim(multiSet[[set]]$data), " dimensions on data")
      base::colnames(multiSet[[set]]$data) <- value
  }
  multiSet
}

#' @export
#' @rdname colnames.multiSet
`colnames<-` <- function(multiSet, value) {
    UseMethod("colnames<-", multiSet)
}

#' @export
#' @rdname colnames.multiSet
`colnames<-.default` <- function(multiSet, value) {
    base::`colnames<-`(multiSet, value)
}

#' @export
#' @rdname colnames.multiSet
colnames <- function(multiSet) {
    UseMethod("colnames", multiSet)
}

#' @export
#' @rdname colnames.multiSet
#' @param \dots other arguments passed
colnames.default <- function(multiSet) {
    base::colnames(multiSet)
}

#' Extract the colnames of the multiSet
#' @export
#' @param multiSet multiSet object to retrieve the columns of the first data
#' @return \code{colnames.multiSet} returns the vector of column names of the
#' \code{data} component. The function assumes the column names in all sets are
#' the same.
colnames.multiSet <- function(multiSet) {
    base::colnames(multiSet[[1]]$data)
}


# multiSet ####
#' Create a multiSet structure.
#'
#' This function creates a multiSet structure by storing its input arguments
#' as the 'data' components.
#'
#' A multiSet structure is intended to store (the same type of) data for
#' multiple, possibly independent, realizations (for example, expression data
#' for several independent experiments). It is a list where each component
#' corresponds to an (independent) data set. Each component is in turn a list
#' that can hold various types of information but must have a \code{data}
#' component. In a "strict" multiSet structure, the \code{data} components are
#' required to each be a matrix or a data frame and have the same number of
#' columns. In a "loose" multiSet structure, the \code{data} components can be
#' anything (but for most purposes should be of comparable type and content).
#'
#' @param \dots Arguments to be stored in the multiSet structure.
#' @return The resulting multiSet structure.
#' @author Peter Langfelder
#' @seealso \code{\link{multiSet2list}} for converting a multiSet structure
#' to a list; \code{\link{list2multiSet}} for an alternative way of creating a
#' multiSet structure; \code{\link{apply.multiSet}, \link{multiSet.mapply}} for
#' ways of applying a function to each component of a multiSet structure.
#' @keywords misc
#' @examples
#'
#' data1 <- matrix(rnorm(100), 20, 5)
#' data2 <- matrix(rnorm(50), 10, 5)
#'
#' md <- multiSet(Set1 = data1, Set2 = data2)
#'
#' checkSets(md)
#' @export
multiSet <- function(...) {
  list2multiSet(list(...))
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
#' @examples
#' data1 <- matrix(rnorm(100L), 20L, 5L)
#' data2 <- matrix(rnorm(50L), 10L, 5L)
#' colnames(data1) <- LETTERS[1:5]
#' colnames(data2) <- LETTERS[2:6]
#' md <- multiSet(Set1 = data1, Set2 = data2)
#' keepCommonProbes(md)
#' @author Peter Langfelder
#' @seealso \code{\link{checkSets}}
#' @keywords misc
#' @export
keepCommonProbes <- function(multiExpr, orderBy = 1) {
    size <- checkSets(multiExpr)
    nSets <- size$nSets
    if (nSets <= 0) {
        stop("No expression data given!")
    } else if (nSets == 1) {
        warning("A single data set was provided")
    } else if (nSets > 1) {
        Names <- data.frame(Names = colnames(multiExpr[[orderBy]]$data))
        for (set in (1:nSets)) {
            SetNames <- data.frame(Names = colnames(multiExpr[[set]]$data),
                                  index = c(1:dim(multiExpr[[set]]$data)[2]))
            Names <- merge(Names, SetNames, by.x = "Names", by.y = "Names",
                          all = FALSE, sort = FALSE)
        }
        for (set in 1:nSets) {
            multiExpr[[set]]$data <- multiExpr[[set]]$data[, Names[, set + 1]]
        }
    }
    multiExpr
}

