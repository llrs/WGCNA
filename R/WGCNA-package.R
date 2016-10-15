#' Blood Cell Types with Corresponding Gene Markers
#'
#' This matrix gives a predefined set of marker genes for many blood cell
#' types, as reported in several previously-published studies.  It is used with
#' userListEnrichment to search user-defined gene lists for enrichment.
#'
#'
#' @name BloodLists
#' @docType data
#' @format A 2048 x 2 matrix of characters containing Gene / Category pairs.
#' The first column (Gene) lists genes corresponding to a given category
#' (second column).  Each Category entry is of the form <Blood cell
#' type>__<reference>, where the references can be found at
#' \code{\link{userListEnrichment}}.  Note that the matrix is sorted first by
#' Category and then by Gene, such that all genes related to the same category
#' are listed sequentially.
#' @source For references used in this variable, please see
#' \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(BloodLists)
#' head(BloodLists)
#'
NULL

#' Brain-Related Categories with Corresponding Gene Markers
#'
#' This matrix gives a predefined set of marker genes for many brain-related
#' categories (ie., cell type, organelle, changes with disease, etc.), as
#' reported in several previously-published studies.  It is used with
#' userListEnrichment to search user-defined gene lists for enrichment.
#'
#'
#' @name BrainLists
#' @docType data
#' @format A 48319 x 2 matrix of characters containing Gene / Category pairs.
#' The first column (Gene) lists genes corresponding to a given category
#' (second column).  Each Category entry is of the form <Brain
#' descriptor>__<reference>, where the references can be found at
#' \code{\link{userListEnrichment}}.  Note that the matrix is sorted first by
#' Category and then by Gene, such that all genes related to the same category
#' are listed sequentially.
#' @source For references used in this variable, please see
#' \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(BrainLists)
#' head(BrainLists)
#'
NULL

#' Gene Markers for Regions of the Human Brain
#'
#' This matrix gives a predefined set of marker genes for many regions of the
#' human brain, using data from the Allen Human Brain Atlas
#' (http://human.brain-map.org/) as reported in: Hawrylycz MJ, Lein ES,
#' Guillozet-Bongaarts AL, Shen EH, Ng L, Miller JA, et al. (2012) An
#' Anatomically Comprehensive Atlas of the Adult Human Brain Transcriptome.
#' Nature (in press).  It is used with userListEnrichment to search
#' user-defined gene lists for enrichment.
#'
#'
#' @name BrainRegionMarkers
#' @docType data
#' @format A 28477 x 2 matrix of characters containing Gene / Category pairs.
#' The first column (Gene) lists genes corresponding to a given category
#' (second column).  Each Category entry is of the form <Brain Region>_<Marker
#' Type>__HBA.  Note that the matrix is sorted first by Category and then by
#' Gene, such that all genes related to the same category are listed
#' sequentially.
#' @source For references used in this variable, or other information, please
#' see \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(BrainRegionMarkers)
#' head(BrainRegionMarkers)
#'
NULL

#' Immune Pathways with Corresponding Gene Markers
#'
#' This matrix gives a predefined set of marker genes for many immune response
#' pathways, as assembled by Brian Modena (a member of Daniel R Salomon's lab
#' at Scripps Research Institute), and colleagues.  It is used with
#' userListEnrichment to search user-defined gene lists for enrichment.
#'
#'
#' @name ImmunePathwayLists
#' @docType data
#' @format A 3597 x 2 matrix of characters containing Gene / Category pairs.
#' The first column (Gene) lists genes corresponding to a given category
#' (second column).  Each Category entry is of the form <Immune
#' Pathway>__ImmunePathway.  Note that the matrix is sorted first by Category
#' and then by Gene, such that all genes related to the same category are
#' listed sequentially.
#' @source For more information about this list, please see
#' \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(ImmunePathwayLists)
#' head(ImmunePathwayLists)
#'
NULL

#' Inline display of progress
#'
#' These functions provide an inline display of pregress.
#'
#' A progress indicator is a simple inline display of progress intended to
#' satisfy impatient users during lengthy operations. The function
#' \code{initProgInd} initializes a progress indicator (at zero);
#' \code{updateProgInd} updates it to a specified fraction.
#'
#' Note that excessive use of \code{updateProgInd} may lead to a performance
#' penalty (see examples).
#'
#' @aliases initProgInd updateProgInd
#' @param leadStr character string that will be printed before the actual
#' progress number.
#' @param trailStr character string that will be printed after the actual
#' progress number.
#' @param quiet can be used to silence the indicator for non-interactive
#' sessions whose output is typically redirected to a file.
#' @param newFrac new fraction of progress to be displayed.
#' @param progInd an object of class \code{progressIndicator} that encodes
#' previously printed message.
#' @return Both functions return an object of class \code{progressIndicator}
#' that holds information on the last printed value and should be used for
#' subsequent updates of the indicator.
#' @author Peter Langfelder
#' @keywords misc
#' @examples
#'
#' max = 10;
#' prog = initProgInd("Counting: ", "done");
#' for (c in 1:max)
#' {
#'   Sys.sleep(0.10);
#'   prog = updateProgInd(c/max, prog);
#' }
#' printFlush("");
#'
#' printFlush("Example 2:");
#' prog = initProgInd();
#' for (c in 1:max)
#' {
#'   Sys.sleep(0.10);
#'   prog = updateProgInd(c/max, prog);
#' }
#' printFlush("");
#'
#' ## Example of a significant slowdown:
#'
#' ## Without progress indicator:
#'
#' system.time( {a = 0; for (i in 1:10000) a = a+i; } )
#'
#' ## With progress indicator, some 50 times slower:
#'
#' system.time(
#'   {
#'     prog = initProgInd("Counting: ", "done");
#'     a = 0;
#'     for (i in 1:10000)
#'     {
#'       a = a+i;
#'       prog = updateProgInd(i/10000, prog);
#'     }
#'   }
#' )
#'
NULL

#' Pathways with Corresponding Gene Markers - Compiled by Mike Palazzolo and
#' Jim Wang from CHDI
#'
#' This matrix gives a predefined set of marker genes for many immune response
#' pathways, as assembled by Mike Palazzolo and Jim Wang from CHDI, and
#' colleagues.  It is used with userListEnrichment to search user-defined gene
#' lists for enrichment.
#'
#'
#' @name PWLists
#' @docType data
#' @format A 124350 x 2 matrix of characters containing 2724 Gene / Category
#' pairs.  The first column (Gene) lists genes corresponding to a given
#' category (second column).  Each Category entry is of the form <gene
#' set>__<reference>.
#' @source For more information about this list, please see
#' \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(PWLists)
#' head(PWLists)
#'
NULL

#' Stem Cell-Related Genes with Corresponding Gene Markers
#'
#' This matrix gives a predefined set of genes related to several stem cell
#' (SC) types, as reported in two previously-published studies.  It is used
#' with userListEnrichment to search user-defined gene lists for enrichment.
#'
#'
#' @name SCsLists
#' @docType data
#' @format A 14003 x 2 matrix of characters containing Gene / Category pairs.
#' The first column (Gene) lists genes corresponding to a given category
#' (second column).  Each Category entry is of the form <Stem cell-related
#' category>__<reference>, where the references can be found at
#' \code{\link{userListEnrichment}}.  Note that the matrix is sorted first by
#' Category and then by Gene, such that all genes related to the same category
#' are listed sequentially.
#' @source For references used in this variable, please see
#' \code{\link{userListEnrichment}}
#' @keywords datasets
#' @examples
#'
#' data(SCsLists)
#' head(SCsLists)
#'
NULL

#' Weighted Gene Co-Expression Network Analysis
#'
#' Functions necessary to perform Weighted Correlation Network Analysis. WGCNA
#' is also known as weighted gene co-expression network analysis when dealing
#' with gene expression data. Many functions of WGCNA can also be used for
#' general association networks specified by a symmetric adjacency matrix.
#'
#' \tabular{ll}{ Package: \tab WGCNA\cr Version: \tab 1.51\cr Date: \tab
#' 2016-03-08\cr Depends: \tab R (>= 3.0), dynamicTreeCut (>= 1.62),
#' fastcluster, Hmisc \cr Imports: \tab stats, impute, grDevices, utils,
#' splines, reshape, foreach, doParallel, matrixStats (>= 0.8.1), GO.db,
#' AnnotationDbi \cr Suggests: \tab org.Hs.eg.db, org.Mm.eg.db, infotheo,
#' entropy, minet, survival \cr ZipData: \tab no\cr License: \tab GPL (>= 2)\cr
#' URL: \tab
#' \url{http://www.genetics.ucla.edu/labs/horvath/CoexpressionNetwork/Rpackages/WGCNA/}\cr
#' }
#'
#' Index: \preformatted{ GTOMdist Generalized Topological Overlap Measure
#' TOMdist Topological overlap matrix dissimilarity TOMplot Graphical
#' representation of the Topological Overlap Matrix TOMsimilarity Topological
#' overlap matrix similarity TOMsimilarityFromExpr Topological overlap matrix
#' similarity WGCNA-package Weighted Gene Co-Expression Network Analysis
#' accuracyMeasures Accuracy measures for a 2x2 confusion matrix addErrorBars
#' Add error bars to a barplot. addGrid Add grid lines to an existing plot.
#' addGuideLines Add vertical "guide lines" to a dendrogram plot addTraitToMEs
#' Add trait information to multi-set module eigengene structure adjacency
#' Calculate network adjacency adjacency.fromSimilarity Calculate network
#' adjacency from a similarity matrix adjacency.polyReg Adjacency based on
#' polynomial regression adjacency.splineReg Adjacency based on natural cubic
#' spline regression alignExpr Align expression data with given vector
#' automaticNetworkScreening One-step automatic network gene screening
#' automaticNetworkScreeningGS One-step automatic network gene screening with
#' external gene significance AFcorMI Prediction of weighted mutual information
#' adjacency matrix by correlation bicor Biweight Midcorrelation bicorAndPvalue
#' Biweight Midcorrelation and the associated p-value blockwiseConsensusModules
#' Find consensus modules across several datasets. blockwiseIndividualTOMs
#' Calculate individual topological overlaps across multi-set data
#' blockwiseModules Automatic network construction and module detection
#' BloodLists (data) Gene sets for user enrichment analysis blueWhiteRed
#' Blue-white-red color sequence BrainLists (data) Gene sets for user
#' enrichment analysis BrainRegionMarkers (data) Gene Markers for Regions of
#' the Human Brain checkAdjMat Check adjacency matrix checkSets Check structure
#' and retrieve sizes of a group of datasets checkSimilarity Check a similarity
#' matrix chooseOneHubInEachModule Choose a hub gene in each module
#' chooseTopHubInEachModule Choose the top hub gene in each module clusterCoef
#' Clustering coefficient calculation coClustering Cluster preservation based
#' on co-clustering coClustering.permutationTest Permutation test for
#' co-clustering collapseRows Collapse Rows in Numeric Matrix
#' collapseRowsUsingKME Selects one representative row per group based on kM
#' collectGarbage Iterative garbage collection colQuantileC Fast colunm-wise
#' quantile of a matrix conformityBasedNetworkConcepts Calculation of
#' conformity-based network concepts conformityDecomposition Conformity
#' vector(s) and factorizability measure(s) of a network
#' consensusDissTOMandTree Consensus TOM-based dissimilarity and clustering
#' tree consensusKME Consensus eigengene-based connectivity
#' consensusMEDissimilarity Consensus dissimilarity of module eigengenes.
#' consensusOrderMEs Put close eigenvectors next to each other in several sets.
#' consensusProjectiveKMeans Consensus projective K-means (pre-)clustering of
#' expression data cor Faster calculation of Pearson correlations corAndPvalue
#' Correlation and the associated p-value cor1 Faster calculation of column
#' correlations of a matrix corFast Faster calculation of Pearson correlations
#' corPredictionSuccess ~~function to do ... ~~ corPvalueFisher Fisher's
#' asymptotic p-value for correlation corPvalueStudent Student asymptotic
#' p-value for correlation correlationPreservation Preservation of eigengene
#' correlations coxRegressionResiduals Deviance- and martingale residuals from
#' a Cox regression model cutreeStatic Constant height tree cut
#' cutreeStaticColor Constant height tree cut using color labels displayColors
#' Show colors used to label modules dynamicMergeCut Threshold for module
#' merging exportNetworkToVisANT Export network data in format readable by
#' VisANT exportNetworkToCytoscape Export network data in format readable by
#' Cytoscape fixDataStructure Put single-set data into a form useful for
#' multiset calculations fundamentalNetworkConcepts Calculation of fundamental
#' network concepts GOenrichmentAnalysis Calculate enrichment p-values of
#' clusters in GO terms goodGenes Filter genes with too many missing entries
#' goodGenesMS Filter genes with too many missing entries across multiple sets
#' goodSamples Filter samples with too many missing entries goodSamplesGenes
#' Iterative filtering of samples and genes with too many missing entries
#' goodSamplesGenesMS Iterative filtering of samples and genes with too many
#' missing entries across multiple data sets goodSamplesMS Filter samples with
#' too many missing entries across multiple data sets greenBlackRed
#' Green-black-red color sequence greenWhiteRed Green-white-red color sequence
#' hubGeneSignificance Hubgene significance ImmunePathwayLists (data) Immune
#' Pathways with Corresponding Gene Markers initProgInd Inline display of
#' progress intramodularConnectivity intramodularConnectivity.fromExpr
#' Calculation of intramodular connectivity keepCommonProbes Keep probes that
#' are shared among given data sets kMEcomparisonScatterplot Scatterplots of
#' eigengene-based connectivity labeledBarplot Barplot with text or color
#' labels labeledHeatmap Produce a labeled heatmap plot labelPoints Attempt to
#' intelligently label points in a scatterplot labels2colors Convert numerical
#' labels to colors lowerTri2matrix Reconstruct a symmetric matrix from a
#' distance (lower-triangular) representation matchLabels Relabel modules to
#' best approximate a reference labeling mergeCloseModules Merge close modules
#' in gene expression data metaAnalysis Meta-analysis significance statistics
#' metaZfunction Meta-analysis Z statistic moduleColor.getMEprefix Get the
#' prefix used to label module eigengenes moduleEigengenes Calculate module
#' eigengenes moduleMergeUsingKME Merge modules and reassign genes using kME
#' moduleNumber Fixed-height cut of a dendrogram modulePreservation Calculation
#' of module preservation statistics multiSetMEs Calculate module eigengenes
#' multiData.eigengeneSignificance Calculate eigengene significance for
#' multiple data sets mutualInfoAdjacency Calculate weighted adjacency matrices
#' based on mutual information nPresent Number of present data entries
#' nearestNeighborConnectivity Connectivity to a constant number of nearest
#' neighbors nearestNeighborConnectivityMS Connectivity to a constant number of
#' nearest neighbors across multiple data sets nearestCentroidPredictor Nearest
#' centroid predictor for two-class problems networkConcepts Calculations of
#' network concepts networkScreening Network screening networkScreeningGS
#' Network screening with external gene significance normalizeLabels Transform
#' numerical labels into normal order numbers2colors Color representation for a
#' numeric variable orderBranchesUsingHubGenes Optimize dendrogram using branch
#' swaps and reflections orderMEs Put close eigenvectors next to each other
#' overlapTable Overlap counts and Fisher exact tests for two sets of module
#' labels overlapTableUsingKME Determines significant overlap between modules
#' in two networks based on kME tables pickHardThreshold Analysis of scale free
#' topology for hard-thresholding pickHardThreshold.fromSimilarity Analysis of
#' scale free topology for hard-thresholding pickSoftThreshold Analysis of
#' scale free topology for soft-thresholding pickSoftThreshold.fromSimilarity
#' Analysis of scale free topology for soft-thresholding plotClusterTreeSamples
#' Annotated clustering dendrogram of microarray samples plotColorUnderTree
#' Plot color rows under a dendrogram plotDendroAndColors Dendrogram plot with
#' color annotation of objects plotEigengeneNetworks Eigengene network plot
#' plotMEpairs Pairwise scatterplots of eigengenes plotModuleSignificance
#' Barplot of module significance plotNetworkHeatmap Network heatmap plot pmean
#' Parallel mean pmedian Parallel median populationMeansInAdmixture Estimation
#' of population-specific mean values in an admixed population pquantile
#' Parallel quantile preservationNetworkConnectivity Network preservation
#' calculations projectiveKMeans Projective K-means (pre-)clustering of
#' expression data propVarExplained Proportion of variance explained by
#' eigengenes proportionsInAdmixture Estimation of proportion of pure
#' populations in an admixture qvalue q-value calculation from package qvalue
#' randIndex Calculation of (adjusted) Rand index randomGLMpredictor Ensemble
#' predictor based on bagging of generalized linear models rankPvalue Estimate
#' the p-value for ranking consistently high (or low) on multiple lists
#' recutBlockwiseTrees Repeat blockwise module detection from pre-calculated
#' data recutConsensusTrees Repeat blockwise consensus module detection from
#' pre-calculated data redWhiteGreen Red-white-green color sequence
#' reflectTwoBranches Reflect branches in a dendrogram
#' relativeCorPredictionSuccess Compare prediction success removeGreyME Removes
#' the grey eigengene from a given collection of eigengenes.
#' removePrincipalComponents Remove leading principal components from data
#' returnGeneSetsAsLists Return pre-defined gene lists in several biomedical
#' categories. scaleFreeFitIndex Calculation of fitting statistics for
#' evaluating scale free topology fit. scaleFreePlot Visual check of scale-free
#' topology SCsLists (data) Stem Cell-Related Genes with Corresponding Gene
#' Markers selectBranch Find a branch in a dendrogram
#' setCorrelationPreservation Summary correlation preservation measure
#' sigmoidAdjacencyFunction Sigmoid-type adacency function signedKME Signed
#' eigengene-based connectivity signumAdjacencyFunction Hard-thresholding
#' adjacency function simulateDatExpr Simulation of expression data
#' simulateDatExpr5Modules
#'
#' simulateEigengeneNetwork Simulate eigengene network from a causal model
#' simulateModule Simulate a gene co-expression module simulateMultiExpr
#' Simulate multi-set expression data simulateSmallLayer Simulate small modules
#' sizeGrWindow Open a graphics window of given width and height
#' softConnectivity Calculation of soft (weighted) connectevity
#' softConnectivity.fromSimilarity Calculation of soft (weighted) connectevity
#' modules standardScreeningBinaryTrait Standard screening for a binary trait
#' standardScreeningCensoredTime Standard screening with regard to a Censored
#' Time Variable stdErr Standard error stratifiedBarplot Bar plots of data
#' across two splitting parameters swapTwoBranches Swap branches in a
#' dendrogram TrueTrait Estimate the true trait underlying a list of surrogate
#' markers transposeBigData Block-by-block transpose of large matrices
#' unsignedAdjacency Calculation of unsigned adjacency userListEnrichment
#' Measure enrichment between inputted and user-defined lists vectorTOM
#' Topological overlap for a subset of the whole set of genes vectorizeMatrix
#' Turn a matrix into a vector of non-redundant components verboseBarplot
#' Barplot with error bars, annotated by Kruskal-Wallis p-value verboseBoxplot
#' Boxplot annotated by a Kruskal-Wallis p-value verboseScatterplot Scatterplot
#' annotated by regression line and p-value verboseIplot Scatterplot annotated
#' by regression line, p-value, and color for density votingLinearPredictor
#' Voting linear predictor }
#'
#' @name WGCNA-package
#' @aliases WGCNA-package WGCNA
#' @docType package
#' @author Peter Langfelder <Peter.Langfelder@@gmail.com> and Steve Horvath
#' <SHorvath@@mednet.ucla.edu>, with contributions by Jun Dong, Jeremy Miller,
#' Lin Song, Andy Yip, and Bin Zhang
#'
#' Maintainer: Peter Langfelder <Peter.Langfelder@@gmail.com>
#' @references Peter Langfelder and Steve Horvath (2008) WGCNA: an R package
#' for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559
#'
#' Peter Langfelder, Steve Horvath (2012) Fast R Functions for Robust
#' Correlations and Hierarchical Clustering. Journal of Statistical Software,
#' 46(11), 1-17. \url{http://www.jstatsoft.org/v46/i11/}
#'
#' Bin Zhang and Steve Horvath (2005) "A General Framework for Weighted Gene
#' Co-Expression Network Analysis", Statistical Applications in Genetics and
#' Molecular Biology: Vol. 4: No. 1, Article 17
#'
#' Dong J, Horvath S (2007) Understanding Network Concepts in Modules, BMC
#' Systems Biology 2007, 1:24
#'
#' Horvath S, Dong J (2008) Geometric Interpretation of Gene Coexpression
#' Network Analysis. PLoS Comput Biol 4(8): e1000117
#'
#' Yip A, Horvath S (2007) Gene network interconnectedness and the generalized
#' topological overlap measure. BMC Bioinformatics 2007, 8:22
#'
#' Langfelder P, Zhang B, Horvath S (2007) Defining clusters from a
#' hierarchical cluster tree: the Dynamic Tree Cut library for R.
#' Bioinformatics. November/btm563
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co-expression modules. BMC Systems Biology 2007, 1:54
#' @keywords package
NULL

#' Fast calculations of Pearson correlation.
#'
#' These functions implements a faster calculation of Pearson correlation.
#'
#' The speedup against the R's standard \code{\link[stats]{cor}} function will
#' be substantial particularly if the input matrix only contains a small number
#' of missing data. If there are no missing data, or the missing data are
#' numerous, the speedup will be smaller but still present.
#'
#' The fast calculations are currently implemented only for
#' \code{method="pearson"} and \code{use} either \code{"all.obs"} or
#' \code{"pairwise.complete.obs"}.  The \code{corFast} function is a wrapper
#' that calls the function \code{cor}. If the combination of \code{method} and
#' \code{use} is implemented by the fast calculations, the fast code is
#' executed; otherwise, R's own correlation \code{\link[stats]{cor}} is
#' executed.
#'
#' The argument \code{quick} specifies the precision of handling of missing
#' data. Zero will cause all calculations to be executed precisely, which may
#' be significantly slower than calculations without missing data.
#' Progressively higher values will speed up the calculations but introduce
#' progressively larger errors. Without missing data, all column means and
#' variances can be pre-calculated before the covariances are calculated. When
#' missing data are present, exact calculations require the column means and
#' variances to be calculated for each covariance. The approximate calculation
#' uses the pre-calculated mean and variance and simply ignores missing data in
#' the covariance calculation. If the number of missing data is high, the
#' pre-calculated means and variances may be very different from the actual
#' ones, thus potentially introducing large errors.  The \code{quick} value
#' times the number of rows specifies the maximum difference in the number of
#' missing entries for mean and variance calculations on the one hand and
#' covariance on the other hand that will be tolerated before a recalculation
#' is triggered. The hope is that if only a few missing data are treated
#' approximately, the error introduced will be small but the potential speedup
#' can be significant.
#'
#' @aliases cor1 corFast cor
#' @param x a numeric vector or a matrix. If \code{y} is null, \code{x} must be
#' a matrix.
#' @param y a numeric vector or a matrix. If not given, correlations of columns
#' of \code{x} will be calculated.
#' @param use a character string specifying the handling of missing data. The
#' fast calculations currently support \code{"all.obs"} and
#' \code{"pairwise.complete.obs"}; for other options, see R's standard
#' correlation function \code{\link[stats]{cor}}.  Abbreviations are allowed.
#' @param method a character string specifying the method to be used. Fast
#' calculations are currently available only for \code{"pearson"}.
#' @param quick real number between 0 and 1 that controls the precision of
#' handling of missing data in the calculation of correlations. See details.
#' @param cosine logical: calculate cosine correlation? Only valid for
#' \code{method="pearson"}. Cosine correlation is similar to Pearson
#' correlation but the mean subtraction is not performed. The result is the
#' cosine of the angle(s) between (the columns of) \code{x} and \code{y}.
#' @param cosineX logical: use the cosine calculation for \code{x}? This
#' setting does not affect \code{y} and can be used to give a hybrid
#' cosine-standard correlation.
#' @param cosineY logical: use the cosine calculation for \code{y}? This
#' setting does not affect \code{x} and can be used to give a hybrid
#' cosine-standard correlation.
#' @param drop logical: should the result be turned into a vector if it is
#' effectively one-dimensional?
#' @param nThreads non-negative integer specifying the number of parallel
#' threads to be used by certain parts of correlation calculations. This option
#' only has an effect on systems on which a POSIX thread library is available
#' (which currently includes Linux and Mac OSX, but excludes Windows). If zero,
#' the number of online processors will be used if it can be determined
#' dynamically, otherwise correlation calculations will use 2 threads.
#' @param verbose Controls the level of verbosity. Values above zero will cause
#' a small amount of diagnostic messages to be printed.
#' @param indent Indentation of printed diagnostic messages. Each unit above
#' zero adds two spaces.
#' @return The matrix of the Pearson correlations of the columns of \code{x}
#' with columns of \code{y} if \code{y} is given, and the correlations of the
#' columns of \code{x} if \code{y} is not given.
#' @note The implementation uses the BLAS library matrix multiplication
#' function for the most expensive step of the calculation. Using a tuned,
#' architecture-specific BLAS may significantly improve the performance of this
#' function.
#'
#' The values returned by the corFast function may differ from the values
#' returned by R's function \code{\link[stats]{cor}} by rounding errors on the
#' order of 1e-15.
#' @author Peter Langfelder
#' @seealso R's standard Pearson correlation function \code{\link{cor}}.
#' @references Peter Langfelder, Steve Horvath (2012) Fast R Functions for
#' Robust Correlations and Hierarchical Clustering.  Journal of Statistical
#' Software, 46(11), 1-17.  \url{http://www.jstatsoft.org/v46/i11/}
#' @keywords misc
#' @examples
#'
#'
#' ## Test the speedup compared to standard function cor
#'
#' # Generate a random matrix with 200 rows and 1000 columns
#'
#' set.seed(10)
#' nrow = 100;
#' ncol = 500;
#' data = matrix(rnorm(nrow*ncol), nrow, ncol);
#'
#' ## First test: no missing data
#'
#' system.time( {corStd = stats::cor(data)} );
#'
#' system.time( {corFast = cor(data)} );
#'
#' all.equal(corStd, corFast)
#'
#' # Here R's standard correlation performs very well.
#'
#' # We now add a few missing entries.
#'
#' data[sample(nrow, 10), 1] = NA;
#'
#' # And test the correlations again...
#'
#' system.time( {corStd = stats::cor(data, use ='p')} );
#'
#' system.time( {corFast = cor(data, use = 'p')} );
#'
#' all.equal(corStd, corFast)
#'
#' # Here the R's standard correlation slows down considerably
#' # while corFast still retains it speed. Choosing
#' # higher ncol above will make the difference more pronounced.
#'
#'
#'
NULL
