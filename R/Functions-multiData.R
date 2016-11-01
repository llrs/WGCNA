#======================================================================================================
#
# multiData.eigengeneSignificance
#
#======================================================================================================



#' Eigengene significance across multiple sets
#'
#' This function calculates eigengene significance and the associated
#' significance statistics (p-values, q-values etc) across several data sets.
#'
#' This is a convenience function that calculates module eigengene
#' significances (i.e., correlations of module eigengenes with a given trait)
#' across all sets in a multi-set analysis. Also returned are p-values, Z
#' scores, numbers of present (i.e., non-missing) observations for each
#' significance, and optionally the q-values (false discovery rates)
#' corresponding to the p-values.
#'
#' The function \code{corAndPvalueFnc} is currently is expected to accept
#' arguments \code{x} (gene expression profiles) and \code{y} (eigengene
#' expression profiles).  Any additional arguments can be passed via
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
#' @param multiData Expression data (or other data) in multi-set format (see
#' \code{\link{checkSets}}). A vector of lists; in each list there must be a
#' component named \code{data} whose content is a matrix or dataframe or array
#' of dimension 2.
#' @param multiTrait Trait or ourcome data in multi-set format. Only one trait
#' is allowed; consequesntly, the \code{data} component of each component list
#' can be either a vector or a data frame (matrix, array of dimension 2).
#' @param moduleLabels Module labels: one label for each gene in
#' \code{multiExpr}.
#' @param multiEigengenes Optional eigengenes of modules specified in
#' \code{moduleLabels}. If not given, will be calculated from \code{multiExpr}.
#' @param useModules Optional specification of module labels to which the
#' analysis should be restricted. This could be useful if there are many
#' modules, most of which are not interesting. Note that the "grey" module
#' cannot be used with \code{useModules}.
#' @param corAndPvalueFnc Function that calculates associations between
#' expression profiles and eigengenes. See details.
#' @param corOptions List giving additional arguments to function
#' \code{corAndPvalueFnc}. See details.
#' @param corComponent Name of the component of output of
#' \code{corAndPvalueFnc} that contains the actual correlation.
#' @param getQvalues logical: should q-values (estimates of FDR) be calculated?
#' @param setNames names for the input sets. If not given, will be taken from
#' \code{names(multiExpr)}. If those are \code{NULL} as well, the names will be
#' \code{"Set_1", "Set_2", ...}.
#' @param excludeGrey logical: should the grey module be excluded from the kME
#' tables? Since the grey module is typically not a real module, it makes
#' little sense to report kME values for it.
#' @param greyLabel label that labels the grey module.
#' @return A list containing the following components. Each component is a
#' matrix in which the rows correspond to module eigengenes and columns to data
#' sets. Row and column names are set appropriately.
#' \item{eigengeneSignificance}{Module eigengene significance.}
#' \item{p.value}{p-values (returned by \code{corAndPvalueFnc}). }
#' \item{q.value}{q-values corresponding to the p-values above. Only returned
#' in input \code{getWvalues} is \code{TRUE}. } \item{Z}{Z statistics (if
#' returned by \code{corAndPvalueFnc}). } \item{nObservations}{Number of
#' non-missing observations in each correlation/p-value.}
#' @author Peter Langfelder
#' @keywords misc
multiData.eigengeneSignificance = function(multiData, multiTrait, moduleLabels,
                        multiEigengenes = NULL,
                        useModules = NULL,
                        corAndPvalueFnc = corAndPvalue, corOptions = list(),
                        corComponent = "cor", getQvalues = FALSE,
                        setNames = NULL, excludeGrey = TRUE,
                        greyLabel = ifelse(is.numeric(moduleLabels), 0, "grey"))
{
  corAndPvalueFnc = match.fun(corAndPvalueFnc)

  size = checkSets(multiData)
  nSets = size$nSets
  nGenes = size$nGenes
  nSamples = size$nSamples

  if (is.null(multiEigengenes)) {
    multiEigengenes = multiSetMEs(multiData, universalColors = moduleLabels, verbose = 0,
                                  excludeGrey = excludeGrey, grey = greyLabel)
  } else {
    eSize = checkSets(multiEigengenes)
    if (!isTRUE(all.equal(eSize$nSamples, nSamples)))
      stop("Numbers of samples in multiData and multiEigengenes must agree.")
  }

  if (!is.null(useModules))
  {
    keep = substring(colnames(multiEigengenes[[1]]$data), 3) %in% useModules
    if (sum(keep)==0)
      stop("Incorrectly specified 'useModules': no such module(s).")
    if (any( ! (useModules %in% substring(colnames(multiEigengenes[[1]]$data), 3))))
      stop("Some entries in 'useModules' do not exist in the module labels or eigengenes.")
    for (set in 1:nSets)
      multiEigengenes[[set]]$data = multiEigengenes[[set]]$data[, keep, drop = FALSE]
  }

  modLevels = substring(colnames(multiEigengenes[[1]]$data), 3)
  nModules = length(modLevels)

  MES = p = Z = nObs = array(NA, dim = c(nModules, nSets))

  haveZs = FALSE
  for (set in 1:nSets) {
    corOptions$x = multiEigengenes[[set]]$data
    corOptions$y = multiTrait[[set]]$data
    cp = do.call(corAndPvalueFnc, args = corOptions)
    corComp = grep(corComponent, names(cp))
    pComp = match("p", names(cp))
    if (is.na(pComp)) pComp = match("p.value", names(cp))
    if (is.na(pComp)) stop("Function `corAndPvalueFnc' did not return a p-value.")
    MES[, set] = cp[[corComp]]
    p[, set] = cp[[pComp]]
    if (!is.null(cp$Z)) { Z[, set] = cp$Z; haveZs = TRUE}
    if (!is.null(cp$nObs))
    {
       nObs[, set] = cp$nObs
    } else
       nObs[, set] = t(is.na(multiEigengenes[[set]]$data)) %*% (!is.na(multiTrait[[set]]$data))
  }

  if (is.null(setNames))
     setNames = names(multiData)

  if (is.null(setNames))
     setNames = paste0("Set_", c(1:nSets))

  colnames(MES) = colnames(p) = colnames(Z) = colnames(nObs) = setNames
  rownames(MES) = rownames(p) = rownames(Z) = rownames(nObs) = colnames(multiEigengenes[[1]]$data)

  if (getQvalues)
  {
    q = apply(p, 2, qvalue.restricted)
    dim(q) = dim(p)
    dimnames(q) = dimnames(p)
  } else q = NULL

  if (!haveZs) Z = NULL

  list(eigengeneSignificance = MES,
       p.value = p,
       q.value = q,
       Z = Z,
       nObservations = nObs)
}

#' Number of sets in a multi-set variable
#'
#' A convenience function that returns the number of sets in a multi-set
#' variable.
#'
#'
#' @param multiData vector of lists; in each list there must be a component
#' named \code{data} whose content is a matrix or dataframe or array of
#' dimension 2.
#' @param \dots Other arguments to function \code{\link{checkSets}}.
#' @return A single integer that equals the number of sets given in the input
#' \code{multiData}.
#' @author Peter Langfelder
#' @seealso \code{\link{checkSets}}
#' @keywords misc
nSets = function(multiData, ...) {
  size = checkSets(multiData, ...)
  size$nSets
}
