#============================================================================================================
#
# mdx version of collapseRows
#
#============================================================================================================


.cr.absMaxMean <- function(x, robust)
{
   if (robust)
   {
     colMedians(abs(x), na.rm = TRUE)
   } else
     colMeans(abs(x), na.rm = TRUE)
}

.cr.absMinMean <- function(x, robust)
{
   if (robust)
   {
     -colMedians(abs(x), na.rm = TRUE)
   } else
     -colMeans(abs(x), na.rm = TRUE)
}


.cr.MaxMean <- function(x, robust)
{
   if (robust)
   {
     colMedians(x, na.rm = TRUE)
   } else
     colMeans(x, na.rm = TRUE)
}

.cr.MinMean <- function(x, robust)
{
   if (robust)
   {
     -colMedians(x, na.rm = TRUE)
   } else
     -colMeans(x, na.rm = TRUE)
}

.cr.maxVariance <- function(x, robust)
{
   if (robust) 
   {
     colMads(x, na.rm = TRUE)
   } else 
     colSds(x, na.rm = TRUE)

}

.checkConsistencyOfGroupAndColID = function(mdx, colID, group)
{
  colID  = as.character(colID)
  group = as.character(group)
  if (length(colID)!=length(group))
     stop("'group' and 'colID' must have the same length.")

  if (any(duplicated(colID)))
     stop("'colID' contains duplicate entries.")

  rnDat = mtd.colnames(mdx)

  if ( sum(is.na(colID))>0 )
      warning(paste0("The argument colID contains missing data. It is recommend you choose non-missing,\n",
              "unique values for colID, e.g. character strings."))

   if ( sum(group=="",na.rm=TRUE)>0 ){
      warning(paste("group contains blanks. It is strongly recommended that you remove",
                    "these rows before calling the function.\n",
                    "   But for your convenience, the collapseRow function will remove these rows"))
      group[group==""]=NA
   }

   if ( sum(is.na(group))>0 ){
       warning(paste("The argument group contains missing data. It is strongly recommended\n",
              "   that you remove these rows before calling the function. Or redefine group\n",
           "   so that it has no missing data. But for convenience, we remove these data."))
   }   

  if ((is.null(rnDat))&(checkSets(mdx)$nGenes==length(colID)))
  {
    write("Warning: mdx does not have column names. Using 'colID' as column names.","")
    rnDat = colID
    mdx = mtd.setColnames(mdx, colID)
  }
  if (is.null(rnDat))
     stop("'mdx' does not have row names and \n",
          "length of 'colID' is not the same as # variables in mdx.")

  keepProbes = rep(TRUE, checkSets(mdx)$nGenes)

  if (sum(is.element(rnDat,colID))!=length(colID)){
      write("Warning: row names of input data and probes not identical...","")
      write("... Attempting to proceed anyway. Check results carefully.","")
      keepProbes = is.element(colID, rnDat)
      colID = colID[keepProbes]
      mdx = mtd.subset(mdx, , colID)
      group = group[colID]
  }

  restCols = (group!="" & !is.na(group))
  if (any(!restCols))
  {
    keepProbes[keepProbes] = restCols
    mdx = mtd.subset(mdx, , restCols)
    group = group[restCols]
    colID = colID[restCols]
    rnDat = rnDat[restCols]
  }

  list(mdx = mdx, group = group, colID = colID, keepProbes = keepProbes)
}





#' Select columns with the lowest consensus number of missing data
#' 
#' Given a \code{\link{multiData}} structure, this function calculates the
#' consensus number of present (non-missing) data for each variable (column)
#' across the data sets, forms the consensus and for each group selects
#' variables whose consensus proportion of present data is at least
#' \code{selectFewestMissing} (see usage below).
#' 
#' A 'consensus' of a vector (say 'x') is simply defined as the quantile with
#' probability \code{consensusQuantile} of the vector x. This function
#' calculates, for each variable in \code{mdx}, its proportion of present
#' (i.e., non-NA and non-NaN) values in each of the data sets in \code{mdx},
#' and forms the consensus. Only variables whose consensus proportion of
#' present data is at least \code{selectFewestMissing} are retained.
#' 
#' @param mdx A \code{\link{multiData}} structure. All sets must have the same
#' columns.
#' @param colID Character vector of column identifiers.  This must include all
#' the column names from \code{mdx}, but can include other values as well. Its
#' entries must be unique (no duplicates) and no missing values are permitted.
#' @param group Character vector whose components contain the group label (e.g.
#' a character string) for each entry of \code{colID}. This vector must be of
#' the same length as the vector \code{colID}. In gene expression applications,
#' this vector could contain the gene symbol (or a co-expression module label).
#' @param minProportionPresent A numeric value between 0 and 1 (logical values
#' will be coerced to numeric).  Denotes the minimum consensus fraction of
#' present data in each column that will result in the column being retained.
#' @param consensusQuantile A number between 0 and 1 giving the quantile
#' probability for consensus calculation. 0 means the minimum value (true
#' consensus) will be used.
#' @param verbose Level of verbosity; 0 means silent, larger values will cause
#' progress messages to be printed.
#' @param ... Other arguments that should be considered undocumented and
#' subject to change.
#' @return A logical vector with one element per variable in \code{mdx}, giving
#' \code{TRUE} for the retained variables.
#' @author Jeremy Miller and Peter Langfelder
#' @seealso \code{\link{multiData}}
#' @keywords misc
selectFewestConsensusMissing <- function(mdx, colID, group, 
                                 minProportionPresent = 1,
                                 consensusQuantile = 0,
                                 verbose = 0, ...)

{
## For each gene, select the gene with the fewest missing probes, and return the results.
#   If there is a tie, keep all probes involved in the tie.
#   The main part of this function is run only if omitGroups=TRUE
   
   otherArgs = list(...)
   nVars = checkSets(mdx)$nGenes
   nSamples = checkSets(mdx)$nSamples
   nSets = length(mdx)

   if ((!"checkConsistency" %in% names(otherArgs)) || otherArgs$checkConsistency)
   {
      cd = .checkConsistencyOfGroupAndColID(mdx, colID, group)
      mdx = cd$mdx
      group = cd$group
      colID = cd$colID
      keep = cd$keepProbes
   } else
      keep = rep(TRUE, nVars)

   # First, return datET if there is no missing data, otherwise run the function
   if (sum(mtd.apply(mdx, function(x) sum(is.na(x)), mdaSimplify = TRUE))==0) 
     return(rep(TRUE, nVars))
   
   # Set up the variables.
   names(group)     = colID
   probes              = mtd.colnames(mdx)
   genes               = group[probes]
   keepGenes           = rep(TRUE, nVars)
   tGenes              = table(genes)
   checkGenes          = sort(names(tGenes)[tGenes>1])
   presentData         = as.matrix(mtd.apply(mdx, function(x) colSums(is.finite(x)), mdaSimplify = TRUE))

   presentFrac = presentData/matrix(nSamples, nVars, nSets, byrow = TRUE)
   consensusPresentFrac = .consensusCalculation(setTomMat = presentFrac, useMean = FALSE,
                            setWeightMat = NULL,
                            consensusQuantile = consensusQuantile)$consensus
   
   # Omit all probes with at least omitFrac genes missing
   #keep = consensusPresentFrac > omitFraction
   minProportionPresent = as.numeric(minProportionPresent)

   # Omit relevant genes and return results
   if (minProportionPresent > 0)
   {
      if (verbose) pind = initProgInd()
      for (gi in 1:length(checkGenes))
      {
         g = checkGenes[gi]
         gn            = which(genes==g)
         keepGenes[gn] = (consensusPresentFrac[gn] >= minProportionPresent * max(consensusPresentFrac[gn]))
         if (verbose) pind = updateProgInd(gi/length(checkGenes), pind)
      }
      if (verbose) printFlush("")
   }

   keep[keep] = keepGenes
   return (keep)
}

# ----------------- Main Function ------------------- #



#' Consensus selection of group representatives
#' 
#' Given multiple data sets corresponding to the same variables and a grouping
#' of variables into groups, the function selects a representative variable for
#' each group using a variety of possible selection approaches. Typical uses
#' include selecting a representative probe for each gene in microarray data.
#' 
#' This function was inspired by \code{\link{collapseRows}}, but there are also
#' important differences. This function focuses on selecting representatives;
#' when summarization is more important, \code{collapseRows} provides more
#' flexibility since it does not require that a single representative be
#' selected.
#' 
#' This function and \code{collapseRows} use different input and ouput
#' conventions; user-specified functions need to be tailored differently for
#' \code{collapseRows} than for \code{consensusRepresentatives}.
#' 
#' Missing data are allowed and are treated as missing at random. If
#' \code{rowID} is \code{NULL}, it is replaced by the variable names in
#' \code{mdx}.
#' 
#' All groups with a single variable are represented by that variable, unless
#' the consensus proportion of present data in the variable is lower than
#' \code{minProportionPresent}, in which case the variable and the group are
#' excluded from the output.
#' 
#' For all variables belonging to groups with 2 variables (when
#' \code{useGroupHubs=TRUE}) or with at least 2 variables (when
#' \code{useGroupHubs=FALSE}), selection statistics are calculated in each set
#' (e.g., the selection statistic may be the mean, variance, etc). This results
#' in a matrix of selection statistics (one entry per variable per data set).
#' The selection statistics are next optionally calibrated (normalized) between
#' sets to make them comparable; currently the only implemented calibration
#' method is quantile normalization.
#' 
#' For each variable, the consensus selection statistic is defined as the
#' consensus of the (calibrated) selection statistics across the data sets is
#' calculated. The 'consensus' of a vector (say 'x') is simply defined as the
#' quantile with probability \code{consensusQuantile} of the vector x.
#' Important exception: for the \code{"MinMean"} and \code{"absMinMean"}
#' methods, the consensus is the quantile with probability
#' \code{1-consensusQuantile}, since the idea of the consensus is to select the
#' worst (or close to worst) value across the data sets.
#' 
#' For each group, the representative is selected as the variable with the best
#' (typically highest, but for \code{"MinMean"} and \code{"absMinMean"} methods
#' the lowest) consensus selection statistic.
#' 
#' If \code{useGroupHubs=TRUE}, the intra-group connectivity is calculated for
#' all variables in each set. The intra-group connectivities are optionally
#' calibrated (normalized) between sets, and consensus intra-group connectivity
#' is calculated similarly to the consensus selection statistic above. In each
#' group, the variable with the highest consensus intra-group connectivity is
#' chosen as the representative.
#' 
#' @param mdx A \code{\link{multiData}} structure. All sets must have the same
#' columns.
#' @param group Character vector whose components contain the group label (e.g.
#' a character string) for each entry of \code{colID}. This vector must be of
#' the same length as the vector \code{colID}. In gene expression applications,
#' this vector could contain the gene symbol (or a co-expression module label).
#' @param colID Character vector of column identifiers.  This must include all
#' the column names from \code{mdx}, but can include other values as well. Its
#' entries must be unique (no duplicates) and no missing values are permitted.
#' @param consensusQuantile A number between 0 and 1 giving the quantile
#' probability for consensus calculation. 0 means the minimum value (true
#' consensus) will be used.
#' @param method character string for determining which method is used to
#' choose the representative (when \code{useGroupHubs} is \code{TRUE}, this
#' method is only used for groups with 2 variables). The following values can
#' be used: "MaxMean" (default) or "MinMean" return the variable with the
#' highest or lowest mean value, respectively; "maxRowVariance" return the
#' variable with the highest variance; "absMaxMean" or "absMinMean" return the
#' variable with the highest or lowest mean absolute value; and "function" will
#' call a user-input function (see the description of the argument
#' \code{selectionStatisticFnc}). The built-in functions can be instructed to
#' use robust analogs (median and median absolute deviation) by also specifying
#' \code{statisticFncArguments=list(robust = TRUE)}.
#' @param useGroupHubs Logical: if \code{TRUE}, groups with 3 or more variables
#' will be represented by the variable with the highest connectivity according
#' to a signed weighted correlation network adjacency matrix among the
#' corresponding rows. The connectivity is defined as the row sum of the
#' adjacency matrix. The signed weighted adjacency matrix is defined as
#' A=(0.5+0.5*COR)^power where power is determined by the argument
#' \code{connectivityPower} and COR denotes the matrix of pairwise correlation
#' coefficients among the corresponding rows. Additional arguments to the
#' underlying function \code{\link{adjacency}} can be specified using the
#' argument \code{adjacencyArguments} below.
#' @param calibration Character string describing the method of calibration of
#' the selection statistic among the data sets. Recognized values are
#' \code{"none"} (no calibration) and \code{"full quantile"} (quantile
#' normalization).
#' @param selectionStatisticFnc User-supplied function used to calculate the
#' selection statistic when \code{method} above equals \code{"function"}.  The
#' function must take argumens \code{x} (a matrix) and possibly other arguments
#' that can be specified using \code{statisticFncArguments} below. The return
#' value must be a vector with one component per column of \code{x} giving the
#' selection statistic for each column.
#' @param connectivityPower Positive number (typically integer) for specifying
#' the soft-thresholding power used to construct the signed weighted adjacency
#' matrix, see the description of \code{useGroupHubs}. This option is only used
#' if \code{useGroupHubs} is \code{TRUE}.
#' @param minProportionPresent A number between 0 and 1 specifying a filter of
#' candidate probes. Specifically, for each group, the variable with the
#' maximum consensus proportion of present data is found. Only variables whose
#' consensus proportion of present data is at least \code{minProportionPresent}
#' times the maximum consensus proportion are retained as candidates for being
#' a representative.
#' @param getRepresentativeData Logical: should the representative data, i.e.,
#' \code{mdx} restricted to the representative variables, be returned?
#' @param statisticFncArguments A list giving further arguments to the
#' selection statistic function. Can be used to supply additional arguments to
#' the user-specified \code{selectionStatisticFnc}; the value \code{list(robust
#' = TRUE)} can be used with the built-in functions to use their robust
#' variants.
#' @param adjacencyArguments Further arguments to the function
#' \code{adjacency}, e.g. \code{adjacencyArguments=list(corFnc = "bicor",
#' corOptions = "use = 'p', maxPOutliers = 0.05")} will select the robust
#' correlation \code{\link{bicor}} with a good set of options. Note that the
#' \code{\link{adjacency}} arguments \code{type} and \code{power} cannot be
#' changed.
#' @param verbose Level of verbosity; 0 means silent, larger values will cause
#' progress messages to be printed.
#' @param indent Indent for the diagnostic messages; each unit equals two
#' spaces.
#' @return \item{representatives}{A named vector giving, for each group, the
#' selected representative (input \code{rowID} or the variable (column) name in
#' \code{mdx}). Names correspond to groups.} \item{varSelected}{A logical
#' vector with one entry per variable (column) in input \code{mdx} (possibly
#' after restriction to variables occurring in \code{colID}), \code{TRUE} if
#' the column was selected as a representative.} \item{representativeData}{Only
#' present if \code{getRepresentativeData} is \code{TRUE}; the input \code{mdx}
#' restricted to the representative variables, with column names changed to the
#' corresponding groups.}
#' @author Peter Langfelder, based on code by Jeremy Miller
#' @seealso \code{\link{multiData}} for a description of the \code{multiData}
#' structures; \code{\link{collapseRows}} that solves a related but different
#' problem. Please note the differences in input and output!
#' @keywords misc
consensusRepresentatives = function(mdx, 
                            group, colID, 
                            consensusQuantile = 0,
                            method = "MaxMean", 
                            useGroupHubs = TRUE,
                            calibration = c("none", "full quantile"),
                            selectionStatisticFnc = NULL, 
                            connectivityPower=1, 
                            minProportionPresent=1,
                            getRepresentativeData = TRUE,
                            statisticFncArguments = list(),
                            adjacencyArguments = list(),
                            verbose = 2, indent = 0)

# Change in methodFunction: if the methodFunction picks a single representative, it should return it in
# attribute "selectedRepresentative".
# minProportionPresent now gives the fraction of the maximum of present values that will still be included.
# minProportionPresent=1 corresponds to minProportionPresent=TRUE in original collapseRows.

# In connectivity-based collapsing, use simple connectivity, do not normalize. This way the connectivities
# retain a larger spread which should prevent quantile normalization from making big changes and potentially
# suprious changes. 

{

   if (!is.null(dim(mdx)))
   {
     warning("consensusRepresentatives: wrapping matrix-like input into a mdx structure.")
     mdx = multiData(mdx)
   }
   spaces = indentSpaces(indent)
   nSamples= checkSets(mdx)$nSamples
   nSets = length(mdx)

   colnames.in = mtd.colnames(mdx)

   calibration = match.arg(calibration)
   
    ## Test to make sure the variables are the right length.
    #     if not, fix it if possible, or stop.

   cd = .checkConsistencyOfGroupAndColID(mdx, colID, group)
   colID = cd$colID
   group = cd$group
   mdx = cd$mdx
   keepVars = cd$keepProbes

   rnDat = mtd.colnames(mdx)

      
## For each gene, select the gene with the fewest missing probes (if minProportionPresent==TRUE)
##  Also, remove all probes with more than 90% missing data
   
   if (verbose > 0)
     printFlush(paste0(spaces, "..selecting variables with lowest numbers of missing data.."))
   keep = selectFewestConsensusMissing(mdx, colID, group, minProportionPresent, 
                                 consensusQuantile = consensusQuantile, verbose = verbose -1)
   mdx = mtd.subset(mdx, , keep)
   keepVars[keepVars] = keep

   group = group[keep]
   colID = colID[keep]

   rnDat = mtd.colnames(mdx)
   
##   If method="function", use the function "methodFunction" as a way of combining genes
#    Alternatively, use one of the built-in functions 
#    Note: methodFunction must be a function that takes a vector of numbers as input and
#     outputs a single number. This function will return(0) or crash otherwise.

   recMethods = c("function","MaxMean","maxVariance","MinMean","absMinMean","absMaxMean")
   imethod = pmatch(method, recMethods)
        
   if (is.na(imethod)) 
      stop("Error: entered method is not a legal option. Recognized options are\n",
           "       *maxVariance*, *MaxMean*, *MinMean*, *absMaxMean*, *absMinMean*\n",
           "       or *function* for a user-defined function.")

   if (imethod > 1) 
   {
     selectionStatisticFnc = paste0(".cr.", method)
     selStatFnc = get(selectionStatisticFnc, mode = "function")
   } else {
     selStatFnc = match.fun(selectionStatisticFnc)
     if((!is.function(selStatFnc))&(!is.null(selStatFnc)))
            stop("Error: 'selectionStatisticFnc must be a function... please read the help file.")
   }
      
## Format the variables for use by this function
   colID[is.na(colID)] = group[is.na(colID)]    # Use group if row is missing
   rnDat[is.na(rnDat)]   = group[is.na(rnDat)]
   mdx = mtd.setColnames(mdx, rnDat)

   remove       = (is.na(colID))|(is.na(group)) # Omit if both gene and probe are missing
   colID  = colID[!remove]
   group = group[!remove]
   names(group) = colID
   colID = sort(intersect(rnDat,colID))
   if (length(colID)<=1)
      stop("None of the variable names in 'mdx' are in 'colID'.")

   group = group[colID]
   mdx  = mtd.apply(mdx, as.matrix)
   keepVars[keepVars] =  mtd.colnames(mdx) %in% colID
   mdx = mtd.subset(mdx, , colID)

   probes = mtd.colnames(mdx)
   genes  = group[probes]
   tGenes = table(genes)
   colnames.out = sort(names(tGenes))
    
   if (getRepresentativeData)
   {
     mdxOut = mtd.apply(mdx, function(x) 
     {
       out = matrix(0, nrow(x), length(tGenes))
       rownames(out) = rownames(x)
       colnames(out) = colnames.out
       out
     })
     names(mdxOut) = names(mdx)
   }

   representatives = rep("", length(colnames.out))
   names(representatives) = colnames.out; 
   
##  If !is.null(connectivityPower), default to the connectivity method with power=method
#      Collapse genes with multiple probe sets together using the following algorthim:
#      1) If there is one ps/g = keep
#      2) If there are 2 ps/g = (use "method" or "methodFunction")
#      3) If there are 3+ ps/g = take the max connectivity
#   Otherwise, use "method" if there are 3+ ps/g as well. 
   if(!is.null(connectivityPower)){
     if(!is.numeric(connectivityPower))
        stop("Error: if entered, connectivityPower must be numeric.")
     if(connectivityPower<=0)
       stop("Warning: connectivityPower must be >= 0.")

     if(any(nSamples<=5)){
       write("Warning: 5 or fewer samples, this method of probe collapse is unreliable...","")
       write("...Running anyway, but we suggest trying another method (for example, *mean*).","")
     }
   }
   
   # Run selectionStatisticFnc on all data; if quantile normalization is requested, normalize the selection
   # statistics across data sets.

   selectionStatistics = mtd.apply(mdx, function(x) 
             do.call(selStatFnc, c(list(x), statisticFncArguments)),
              mdaSimplify = TRUE)

   #if (FALSE) xxx = selectionStatistics

   if (is.null(dim(selectionStatistics)))
      stop("Calculation of selection statistics produced results of zero or unqual lengths.")

   if (calibration=="full quantile")
      selectionStatistics = normalize.quantiles(selectionStatistics)

   #if (FALSE)
   #{
   #   sizeGrWindow(14, 5)
   #   par(mfrow = c(1,3))
   #   for (set in 1:nSets)
   #     hist(xxx[, set], breaks = 200)
#
##      for (set in 1:nSets)
 #       verboseScatterplot(xxx[, set], selectionStatistics[, set], samples = 10000)
 #  }

   consensusSelStat = .consensusCalculation(selectionStatistics, useMean = FALSE, setWeightMat = NULL,
                             consensusQuantile = consensusQuantile)$consensus

   # Actually run the summarization.

   ones = sort(names(tGenes)[tGenes==1])
   if(useGroupHubs){
      twos = sort(names(tGenes)[tGenes==2]) # use "method" and connectivity
      more = sort(names(tGenes)[tGenes>2])
   } else { 
      twos = sort(names(tGenes)[tGenes>1]) # only use "method"
      more = character(0)
   }
   ones2genes =  match(ones, genes)
   if (getRepresentativeData) for (set in 1:nSets)
       mdxOut[[set]]$data[,ones] = mdx[[set]]$data[, ones2genes]
   representatives[ones] = probes[ones2genes]
   count = 0

   if (length(twos) > 0)
   {
     if (verbose > 0)
       printFlush(paste0(spaces, "..selecting representatives for 2-variable groups.."))
     if (verbose > 1) pind = initProgInd(paste(spaces, ".."))
     repres = rep(NA, length(twos))
     for (ig in 1:length(twos))
     {
        g = twos[ig]
        probeIndex = which(genes==g)
        repres[ig] = probeIndex[which.max(consensusSelStat[probeIndex])]
        if (verbose > 1) pind = updateProgInd(ig/length(twos), pind)
     }
     if (verbose > 1) printFlush("")
     if (getRepresentativeData) for (set in 1:nSets)
        mdxOut[[set]]$data[, twos] = mdx[[set]]$data[, repres]
     representatives[twos] = probes[repres]
   }
   if (length(more) > 0)
   {
     if (verbose > 0)
       printFlush(paste0(spaces, "..selecting representatives for 3-variable groups.."))
     if (verbose > 1) pind = initProgInd(paste(spaces, ".."))
     genes.more = genes[genes %in% more]
     nAll = length(genes.more)
     connectivities = matrix(NA, nAll, nSets)
     for (ig in 1:length(more))
     {
        g = more[ig]
        keepProbes1 = which(genes==g)
        keep.inMore = which(genes.more==g)
        mdxTmp = mtd.subset(mdx, , keepProbes1)
        adj = mtd.apply(mdxTmp, function(x) do.call(adjacency, 
                c(list(x, type = "signed", power = connectivityPower), adjacencyArguments)))
        connectivities[keep.inMore, ] = mtd.apply(adj, colSums, mdaSimplify = TRUE)
        count = count + 1
        if (count %% 50000 == 0) collectGarbage()
        if (verbose > 1) pind = updateProgInd(ig/(2*length(more)), pind)
     }

     if (calibration=="full quantile")
       connectivities = normalize.quantiles(connectivities)

     consConn = .consensusCalculation(connectivities, useMean = FALSE, setWeightMat = NULL,
                                               consensusQuantile = consensusQuantile)$consensus
     repres.inMore = rep(0, length(more))
     for (ig in 1:length(more))
     {
        probeIndex = which(genes.more==more[ig])
        repres.inMore[ig] = probeIndex[which.max(consConn[probeIndex])]
        if (verbose > 1) pind = updateProgInd(ig/(2*length(more)) + 0.5, pind)
     }
     repres = which(genes %in% more)[repres.inMore]
     if (verbose > 1) printFlush("")
     if (getRepresentativeData) for (set in 1:nSets)
       mdxOut[[set]]$data[, more] = as.numeric(mdx[[set]]$data[, repres]); 
     representatives[more] = probes[repres]
   }
      
   # Retreive the information about which probes were saved, and include that information
   #   as part of the output.  

   out2 = cbind(colnames.out, representatives)
   colnames(out2) = c("group","selectedColID")

   reprIndicator = keepVars
   reprIndicator[keepVars] [match(representatives, mtd.colnames(mdx))] = TRUE
   reprIndicator = colnames.in
   out = list(representatives = out2, 
              varSelected = reprIndicator,
              representativeData = if (getRepresentativeData) mdxOut else NULL)
   return(out)
      
} 

