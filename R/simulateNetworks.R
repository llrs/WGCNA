# Functions to simulate networks

# simulateEigengeneNetwork ####
# Note: The returned vector may contain multiple occurences of the same child.
.causalChildren <- function(parents, causeMat) {
  nNodes <- dim(causeMat)[[1]]

  # print(paste("Length of parents: ", length(parents)))
  if (length(parents) == 0)
    return(NULL)

  Child_ind <- apply(as.matrix(abs(causeMat[, parents])), 1, sum) > 0
  if (sum(Child_ind) > 0) {
    children <- c(1:nNodes)[Child_ind]
  } else {
    children <- NULL
  }
  children
}

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
simulateEigengeneNetwork <- function(causeMat, anchorIndex, anchorVectors,
                                     noise = 1, verbose = 0, indent = 0) {
  spaces <- indentSpaces(indent)

  if (verbose > 0){
    printFlush(paste(spaces, "Creating seed vectors..."))
  }
  nNodes <- dim(causeMat)[[1]]
  nSamples <- dim(anchorVectors)[[1]]

  if (length(anchorIndex) != dim(anchorVectors)[[2]]) {
    stop("Length of anchorIndex must equal the number of vectors in",
         "anchorVectors.")
  }

  if (length(noise) == 1) {
    noise <- rep(noise, nNodes)
  }
  if (length(noise) != nNodes) {
    stop("Length of noise must equal the number of nodes as given by the ",
         "dimension of the causeMat matrix."
    )
  }

  # Initialize all node vectors to noise with given standard deviation
  NodeVectors <- matrix(0, nrow = nSamples, ncol = nNodes)
  for (i in 1:nNodes) {
    NodeVectors[, i] <- rnorm(n = nSamples, mean = 0,
                              sd = noise[i])}

  Levels <- rep(0, times = nNodes)

  # Calculate levels for all nodes: start from anchors and go through each
  # successive level of children
  level <- 0
  parents <- anchorIndex
  Children <- .causalChildren(parents = parents, causeMat = causeMat)
  if (verbose > 1) {
    printFlush(spaces, "..Determining level structure...")
  }
  while (!is.null(Children)) {
    # print(paste("level:", level))
    # print(paste("   parents:", parents))
    # print(paste("   Children:", Children))
    level <- level + 1
    if ((verbose > 1) & (level / 10 == as.integer(level / 10))) {
      printFlush(spaces, "  ..Detected level", level)
    }
    #printFlush(paste("Detected level", level))
    Levels[Children] <- level
    parents <- Children
    Children <- .causalChildren(parents = parents, causeMat = causeMat)
  }

  HighestLevel <- level

  # Generate the whole network
  if (verbose > 1) {
    printFlush(spaces, "..Calculating network...")
  }
  NodeVectors[, anchorIndex] = NodeVectors[, anchorIndex] + anchorVectors
  for (level in (1:HighestLevel)) {
    if ((verbose > 1) & (level / 10 == as.integer(level / 10))) {
      printFlush(spaces, " .Working on level", level)
    }
    #printFlush(paste("Working on level", level))
    LevelChildren = c(1:nNodes)[Levels == level]
    for (child in LevelChildren) {
      LevelParents <- c(1:nNodes)[causeMat[child,] != 0]
      for (parent in LevelParents) {
        NodeVectors[, child] <- scale(NodeVectors[, child] +
                                        causeMat[child, parent]  *
                                        NodeVectors[, parent])
      }
    }
  }

  Nodes <- list(eigengenes = NodeVectors, causeMat = causeMat,
                levels = Levels, anchorIndex = anchorIndex)
  Nodes
}

# simulateModule ####
# The resulting data is normalized.
# Attributes contain the component trueKME giving simulated correlation with
# module eigengene for both module genes and near - module genes.
# corPower controls how fast the correlation drops with index i in the module
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
#' given causal structure
#'
#' \code{\link{simulateDatExpr}} for simulations of whole datasets consisting
#' of multiple modules
#'
#' \code{\link{simulateDatExpr5Modules}} for a simplified interface to
#' expression simulations
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
#' \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/EigengeneNetwork/SupplementSimulations.pdf}.
#' @keywords misc
#' @examples
#' nSamples <- 50
#' nGenes <- 30
#' y <- sample(c(1,2), nSamples, replace = TRUE)
#' module1 <- simulateModule(scale(y), nGenes)
#' module2 <- simulateModule(scale(y), nGenes, signed = TRUE)
#' geneMeans <- sample(c(4, 5, 10), nGenes, replace = TRUE)
#' module3 <- simulateModule(scale(y), nGenes, signed = TRUE,
#'                           geneMeans = geneMeans)
simulateModule <- function(ME, nGenes, nNearGenes = 0, minCor = 0.3, maxCor = 1,
                           corPower = 1, signed = FALSE, propNegativeCor = 0.3,
                           geneMeans = NULL, verbose = 0, indent = 0) {
  ME <- as.vector(ME)
  nSamples <- length(ME)

  datExpr <- matrix(rnorm((nGenes + nNearGenes) * nSamples),
                    nrow = nSamples,
                    ncol = nGenes + nNearGenes)

  VarME <- var(ME)

  # generate the in - module genes
  CorME <- maxCor - (c(1:nGenes) / nGenes) ^ (1 / corPower) * (maxCor - minCor)
  noise <- sqrt(VarME * (1 - CorME ^ 2) / CorME ^ 2)
  sign <- rep(1, nGenes)
  if (!signed) {
    negGenes <- as.integer(seq(
      from = 1 / propNegativeCor,
      by = 1 / propNegativeCor,
      length.out = nGenes * propNegativeCor)
    )
    negGenes <- negGenes[negGenes <= nGenes]
    sign[negGenes] <- - 1
  }
  for (gene in 1:nGenes) {
    datExpr[, gene] <- sign[gene] * (ME + rnorm(nSamples, sd = noise[gene]))
  }

  trueKME <- CorME
  # generate the near - module genes
  if (nNearGenes > 0) {
    CorME <- c(1:nNearGenes) / nNearGenes * minCor
    noise <- sqrt(VarME * (1 - CorME ^ 2) / CorME ^ 2)
    sign <- rep(1, nNearGenes)
    if (!signed) {
      negGenes <- as.integer(
        seq(from = 1 / propNegativeCor,
            by = 1 / propNegativeCor,
            length.out = nNearGenes * propNegativeCor)
      )
      negGenes <- negGenes[negGenes <= nNearGenes]
      sign[negGenes] <- - 1
    }
    for (gene in 1:nNearGenes) {
      datExpr[, nGenes + gene] <- ME + sign[gene] * rnorm(
        nSamples, sd = noise[gene])
    }
    trueKME <- c(trueKME, CorME)
  }

  datExpr <- scale(datExpr)
  if (!is.null(geneMeans)) {
    if (any(is.na(geneMeans))) {
      stop("All entries of 'geneMeans' must be finite.")
    }
    if (length(geneMeans) != nGenes + nNearGenes) {
      stop("The lenght of 'geneMeans' must equal nGenes + nNearGenes.")
    }
    datExpr <- datExpr + matrix(geneMeans, nSamples, nGenes + nNearGenes,
                                byrow = TRUE)
  }

  attributes(datExpr)$trueKME <- trueKME

  datExpr
}

# simulateSmallLayer ####
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
#' \code{\link{simulateModule}} for simulation of individual modules
#'
#' \code{\link{simulateDatExpr}} for the main gene expression simulation
#' function.
#' @keywords misc
simulateSmallLayer <- function(order, nSamples, minCor = 0.3, maxCor = 0.5,
                               corPower = 1, averageModuleSize, averageExpr,
                               moduleSpacing, verbose = 4, indent = 0) {
  spaces <- indentSpaces(indent)
  nGenes <- length(order)
  datExpr <- matrix(0, nrow = nSamples, ncol = nGenes)

  maxCorN0 <- averageModuleSize

  if (verbose > 0) {
    printFlush(spaces, "simulateSmallLayer: simulating modules with min ",
               "corr ", minCor, ", average expression ", averageExpr,
               ", average module size ", averageModuleSize,
               ", inverse density ",moduleSpacing)
  }
  index <- 0
  while (index < nGenes) {
    ModSize <- as.integer(rexp(1, 1 / averageModuleSize))
    if (ModSize < 3) {
      ModSize <- 3
    }
    if (index + ModSize > nGenes) {
      ModSize <- nGenes - index
    }
    if (ModSize > 2) {
      ModuleExpr = rexp(1, 1 / averageExpr)
      if (verbose > 4) {
        printFlush(spaces,"  Module of size ", ModSize, ", expression ",
                   ModuleExpr, ", min corr ", minCor," inserted at index ",
                   index + 1)
      }
      ME <- rnorm(nSamples, sd = ModuleExpr)
      NInModule <- as.integer(ModSize * 2 / 3)
      nNearModule <- ModSize - NInModule
      EffMinCor <- minCor * maxCor
      datExpr[, order[(index + 1):(index + ModSize)]] <- ModuleExpr *
        simulateModule(ME, NInModule, nNearModule, EffMinCor, maxCor,
                       corPower)
    }
    index <- index + ModSize * moduleSpacing
  }
  datExpr
}

# simulateDatExpr ####
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
#' given causal structure
#'
#' \code{\link{simulateModule}} for simulations of individual modules
#'
#' \code{\link{simulateDatExpr5Modules}} for a simplified interface to
#' expression simulations
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
simulateDatExpr <- function(eigengenes, nGenes, modProportions, minCor = 0.3,
                            maxCor = 1, corPower = 1, signed = FALSE,
                            propNegativeCor = 0.3, geneMeans = NULL,
                            backgroundNoise = 0.1, leaveOut = NULL,
                            nSubmoduleLayers = 0, nScatteredModuleLayers = 0,
                            averageNGenesInSubmodule = 10,
                            averageExprInSubmodule = 0.2, submoduleSpacing = 2,
                            verbose = 1, indent = 0) {
  spaces <- indentSpaces(indent)

  nMods <- length(modProportions) - 1

  nSamples <- dim(eigengenes)[[1]]

  if (length(minCor) == 1) {
    minCor <- rep(minCor, nMods)
  }
  if (length(maxCor) == 1) {
    maxCor <- rep(maxCor, nMods)
  }

  if (length(minCor) != nMods) {
    stop("Input error: minCor is an array of different lentgh than",
         " the length - 1 of modProportions array.")
  }

  if (length(maxCor) != nMods) {
    stop("Input error: maxCor is an array of different lentgh than",
         "the length - 1 of modProportions array.")
  }

  if (dim(eigengenes)[[2]] != nMods) {
    stop("Input error: Number of seed vectors must equal the",
         "length of modProportions.")
  }

  if (is.null(geneMeans)) {
    geneMeans <- rep(0, nGenes)
  }
  if (length(geneMeans) != nGenes) {
    stop("Length of 'geneMeans' must equal 'nGenes'.")
  }
  if (any(is.na(geneMeans))) {
    stop("All entries of 'geneMeans' must be finite.")
  }
  grey <- 0
  moduleLabels <- c(1:nMods)

  if (sum(modProportions) > 1) {
    stop("Input error: the sum of Mod.Props must be less than 1")
  }

  no.in.modules <- as.integer(nGenes * modProportions)
  no.in.proper.modules <- no.in.modules[c(1:(length(modProportions) - 1))]
  no.near.modules <- as.integer((nGenes - sum(no.in.modules))  *
                                  no.in.proper.modules / sum(
                                    no.in.proper.modules))

  simulate.module <- rep(TRUE, times = nMods)
  if (!is.null(leaveOut)){
    simulate.module[leaveOut] <- FALSE
  }
  no.in.modules[nMods + 1] <- nGenes - sum(
    no.in.proper.modules[simulate.module]) -
    sum(no.near.modules[simulate.module])

  labelOrder <- moduleLabels[rank(-modProportions[-length(modProportions)],
                                  ties.method = "first")]
  labelOrder <- c(labelOrder, grey)

  if (verbose > 0) {
    printFlush(spaces, "simulateDatExpr: simulating ", nGenes, " genes in ",
               nMods, " modules. ")
  }

  if (verbose > 1) {
    printFlush(spaces, "    Simulated labels:", paste(labelOrder[1:nMods],
                                                      collapse = ", "),
               " and ", grey)
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
      if (verbose > 1) {
        printFlush(spaces, "Simulating unordereded extra layer ",
                   layer)
      }
      OrderVector <- sample(nGenes)
      datExpr <- datExpr + simulateSmallLayer(OrderVector, nSamples,
                                              minCor[1], maxCor[1],
                                              corPower,
                                              averageNGenesInSubmodule,
                                              averageExprInSubmodule,
                                              submoduleSpacing,
                                              verbose = verbose - 1,
                                              indent = indent + 1)
    }
  collectGarbage()
  if (verbose > 1) {
    printFlush(spaces, "  Adding background noise with amplitude ",
               backgroundNoise)
  }
  datExpr <- datExpr + rnorm(n = nGenes * nSamples, sd = backgroundNoise)
  means <- colMeans(datExpr)

  datExpr <- datExpr + matrix(geneMeans - means, nSamples, nGenes,
                              byrow = TRUE)

  colnames(datExpr) <- paste0("Gene.", c(1:nGenes))
  rownames(datExpr) <- paste0("Sample.", c(1:nSamples))

  list(datExpr = datExpr,
       setLabels = truemodule,
       allLabels = allLabels,
       labelOrder = labelOrder,
       trueKME = trueKME,
       trueKME.whichMod = trueKME.whichMod
  )
} # end of function

# simulateMultiExpr ####
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
#' given causal structure
#'
#' \code{\link{simulateDatExpr}} for simulation of individual data sets
#'
#' \code{\link{simulateDatExpr5Modules}} for a simple simulation of a data set
#' consisting of 5 modules
#'
#' \code{\link{simulateModule}} for simulations of individual modules
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

# simulateDatExpr5Modules ####
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
#' \code{\link{simulateModule}} for simulation of individual modules
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
