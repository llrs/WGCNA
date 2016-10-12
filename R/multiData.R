#========================================================================================================
#
# Convenience functions for manipulating multiData structures
#
#========================================================================================================

# Note: many of these function would be simpler to use if I used some sort of class/method technique to keep
# track of the class of each object internally. For example, I could then write a generic function "subset" that
# would work consistently on lists and multiData objects. Similarly, multiData2list would simply become a
# method of as.list, and as.list would be safe to use both on lists and on multiData objects. 



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
mtd.subset = function(multiData, rowIndex = NULL, colIndex = NULL, permissive = FALSE, drop = FALSE)
{
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

multiData2list = function(multiData)
{
  lapply(multiData, getElement, 'data')
}



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
#' @param multiData A multiData structure to be converted to a list.
#' @return For \code{list2multiData}, a multiData structure; for
#' \code{multiData2list}, the corresponding list.
#' @author Peter Langfelder
#' @keywords misc
list2multiData = function(data)
{
  out = list()
  for (set in 1:length(data))
    out[[set]] = list(data = data[[set]])
  names(out) = names(data)
  out
}

mtd.colnames = function(multiData)
{
  colnames(multiData[[1]]$data)
}

.calculateIndicator = function(nSets, mdaExistingResults, mdaUpdateIndex)
{
  if (length(mdaUpdateIndex)==0) mdaUpdateIndex = NULL
  calculate = rep(TRUE, nSets)
  if (!is.null(mdaExistingResults))
  {
    nSets.existing = length(mdaExistingResults)
    if (nSets.existing>nSets)
      stop("Number of sets in 'mdaExistingResults' is higher than the number of sets in 'multiData'.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to recalculate all results.")
    if (nSets.existing==0)
      stop("Number of sets in 'mdaExistingResults' is zero.\n",
           "  Please supply a valid 'mdaExistingResults' or NULL to recalculate all results.")
    if (is.null(mdaUpdateIndex))
    {
      calculate[1:length(mdaExistingResults)] = FALSE
    } else {
      if (any(! mdaUpdateIndex %in% c(1:nSets)))
        stop("All entries in 'mdaUpdateIndex' must be between 1 and the number of sets in 'multiData'.")
      calculateIndex = sort(unique(c(mdaUpdateIndex,
                                      if (nSets.existing<nSets) c((nSets.existing+1):nSets) else NULL)))
      calculate[ c(1:nSets)[-calculateIndex] ] = FALSE
    }
  }

  calculate
}

  



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
#' @param mdaRowIndex If given, must be a list of the same length as
#' \code{multiData}. Each element must be a logical or numeric vector that
#' specifies rows in each \code{data} component to select before applying the
#' function.
#' @param mdaColIndex A logical or numeric vector that specifies columns in
#' each \code{data} component to select before applying the function.
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
)
{
  printSpaces = indentSpaces(mdaIndent)

  if (!isMultiData(multiData, strict = FALSE))
    stop("Supplied 'multiData' is not a valid multiData structure.")

  if (mdaSimplify && mdaCopyNonData) 
    stop("Non-data copying is not compatible with simplification.")

  nSets = length(multiData)
  if (mdaCopyNonData) out = multiData else out = vector(mode = "list", length = nSets)

  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex)

  FUN = match.fun(FUN)
  for (set in 1:nSets) 
  {
    if (calculate[set])
    {
      if (mdaVerbose > 0)
        printFlush(paste0(printSpaces, "mtd.apply: working on set ", set)); 
      out[[set]] = list(data = FUN(multiData[[set]]$data, ...))
    } else
      out[set] = mdaExistingResults[set]
  }

  names(out) = names(multiData)

  if (mdaSimplify) 
  {
    if (mdaVerbose > 0)
      printFlush(paste0(printSpaces, "mtd.apply: attempting to simplify...")); 
    return (mtd.simplify(out))
  } else if (returnList) {
    return (multiData2list(out))
  }

  out
}

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
    mdaVerbose = 0, mdaIndent = 0
)
{
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

  if (!is.null(mdaRowIndex))
  {
    if (!is.list(mdaRowIndex))
       stop("mdaRowIndex must be a list, with one component per set.")
    if (length(mdaRowIndex)!=size$nSets)
       stop("Number of components in 'mdaRowIndex' must equal number of sets.")
    doSelection = TRUE
  } else {
    mdaRowIndex = lapply(size$nSamples, function(n) { c(1:n) })
  }
    
  calculate = .calculateIndicator(nSets, mdaExistingResults, mdaUpdateIndex)

  fun = match.fun(FUN) 
  for (set in 1:size$nSets)
  {
    if (calculate[set])
    {
       if (mdaVerbose > 0)
         printFlush(paste0(printSpaces, "mtd.applyToSubset: working on set ", set))
       res[[set]] = list(data = fun( 
                if (doSelection) multiData[[set]] $ data[mdaRowIndex[[set]], mdaColIndex, drop = FALSE] else
                                 multiData[[set]] $ data, ...))
    } else
       res[set] = mdaExistingResults[set]
  }

  names(res) = names(multiData)

  if (mdaSimplify) 
  {
    if (mdaVerbose > 0)
      printFlush(paste0(printSpaces, "mtd.applyToSubset: attempting to simplify..."))
    return (mtd.simplify(res))
  } else if (returnList) {
    return (multiData2list(res))
  }

  return(res)
}



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
mtd.simplify = function(multiData)
{
  len = length(multiData[[1]]$data)
  dim = dim(multiData[[1]]$data)
  simplifiable = TRUE
  nSets = length(multiData)
  for (set in 1:nSets)
  {
    if (len!=length(multiData[[set]]$data)) simplifiable = FALSE
    if (!isTRUE(all.equal( dim, dim(multiData[[set]]$data)))) simplifiable = FALSE
  }
  if (simplifiable)
  {
    if (is.null(dim)) {
       innerDim = len
       innerNames = names(multiData[[1]]$data)
       if (is.null(innerNames)) innerNames = paste0("X", c(1:len))
    } else {
       innerDim = dim
       innerNames = dimnames(multiData[[1]]$data)
       if (is.null(innerNames)) 
         innerNames = lapply(innerDim, function(x) {paste0("X", 1:x)})
       nullIN = sapply(innerNames, is.null)
       if (any(nullIN))
         innerNames[nullIN] = lapply(innerDim[nullIN], function(x) {paste0("X", 1:x)})
    }
    setNames = names(multiData)
    if (is.null(setNames)) setNames = paste0("Set_", 1:nSets)
    mtd.s = matrix(NA, prod(innerDim), nSets)
    for (set in 1:nSets)
      mtd.s[, set] = as.vector(multiData[[set]]$data)

    dim(mtd.s) = c(innerDim, nSets)
    if (!is.null(innerNames))
      dimnames(mtd.s) = c (if (is.list(innerNames)) innerNames else list(innerNames), list(setNames))
    return(mtd.s)
  }
  return(multiData)
}



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
isMultiData = function(x, strict = TRUE)
{
  if (strict) {
     !inherits(try(checkSets(x), silent = TRUE), 'try-error')
  } else {
    hasData = sapply(x, function(l) { "data" %in% names(l) })
    all(hasData)
  }
}



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
  mdmaVerbose = 0, mdmaIndent = 0)

{
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
mtd.rbindSelf = function(multiData)
{
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
mtd.setAttr = function(multiData, attribute, valueList)
{
  size = checkSets(multiData)
  ind = 1
  for (set in 1:size$nSets)
  {
    attr(multiData[[set]]$data, attribute) = valueList[[ind]]
    ind = ind + 1
    if (ind > length(valueList)) ind = 1
  }
  multiData
}



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
mtd.setColnames = function(multiData, colnames)
{
  size = checkSets(multiData)
  for (set in 1:size$nSets)
    colnames(multiData[[set]]$data) = colnames
  multiData
}




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
multiData = function(...)
{
  list2multiData(list(...))
}  


