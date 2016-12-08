
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
    dim <- dim(adjMat)
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
        stop("some entries are not between ", min, " and ", max)
}

# unsignedAdjacency ####
#' Calculation of unsigned adjacency
#'
#' Calculation of the unsigned network adjacency from expression data. The
#' restricted set of parameters for this function should allow a faster and
#' less memory-hungry calculation.
#'
#' The correlation function will be called with arguments \code{datExpr,
#' datExpr2} plus any extra arguments given in \code{corOptions}. If
#' \code{datExpr2} is \code{NULL}, the standard correlation functions will
#' calculate the corelation of columns in \code{datExpr}.
#'
#' @param datExpr expression data. A data frame in which columns are genes and
#' rows ar samples. Missing values are ignored.
#' @param datExpr2 optional specification of a second set of expression data.
#' See details.
#' @param power soft-thresholding power for network construction.
#' @param corFnc character string giving the correlation function to be used
#' for the adjacency calculation. Recommended choices are \code{"cor"} and
#' \code{"bicor"}, but other functions can be used as well.
#' @param corOptions character string giving further options to be passed to
#' the correlation function
#' @return Adjacency matrix of dimensions \code{n*n}, where \code{n} is the
#' number of genes in \code{datExpr}.
#' @author Steve Horvath and Peter Langfelder
#' @seealso \code{\link{adjacency}}
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
#' @examples
#' datExpr <- matrix(rnorm(150), ncol = 5)
#' unsignedAdjacency(datExpr)
unsignedAdjacency <- function(datExpr, datExpr2 = NULL, power = 6,
                              corFnc = "cor", corOptions = "use = 'p'") {
    corExpr <- parse(text = paste(corFnc, "(datExpr, datExpr2 ",
                                 prepComma(corOptions), ")"))
    # abs(cor(datExpr, datExpr2, use = "p"))^power
    adj <- abs(eval(corExpr))^power
    checkAdjMat(adj)
    adj
}

# adjacnecy ####
#' Calculate network adjacency
#'
#' Calculates (correlation or distance) network adjacency from given expression
#' data or from a similarity. Computes the adjacency from the expression data:
#' takes cor, transforms it as appropriate and possibly adds a sign if
#' requested.
#'
#' The argument \code{type} determines whether a correlation (\code{type} one
#' of \code{"unsigned"}, \code{"signed"}, \code{"signed hybrid"}), or a
#' distance network (\code{type} equal \code{"distance"}) will be calculated.
#' In correlation networks the adajcency is constructed from correlations
#' (values between -1 and 1, with high numbers meaning high similarity). In
#' distance networks, the adjacency is constructed from distances (non-negative
#' values, high values mean low similarity).
#'
#' The function calculates the similarity of columns (genes) in \code{datExpr}
#' by calling the function given in \code{corFnc} (for correlation networks) or
#' \code{distFnc} (for distance networks), transforms the similarity according
#' to \code{type} and raises it to \code{power}, resulting in a weighted
#' network adjacency matrix. If \code{selectCols} is given, the \code{corFnc}
#' function will be given arguments \code{(datExpr, datExpr[selectCols], ...)}
#' hence the returned adjacency will have rows corresponding to all genes and
#' columns corresponding to genes selected by \code{selectCols}.
#'
#' Correlation and distance are transformed as follows: for \code{type =
#' "unsigned"}, adjacency = |cor|^power; for \code{type = "signed"} , adjacency
#' = (0.5 * (1 + cor))^power; for \code{type = "signed hybrid"}, adjacency =
#' cor^power if cor > 0 and 0 otherwise; and for \code{type = "distance"},
#' adjacency = (1 - (dist/max(dist))^2)^power.
#'
#' The function \code{adjacency.fromSimilarity} inputs a similarity matrix,
#' that is it skips the correlation calculation step but is otherwise
#' identical.
#'
#' @aliases adjacency adjacency.fromSimilarity
#' @param similarity A (signed) similarity matrix: square, symmetric matrix
#' with entries between - 1 and 1.
#' @param type Network type, allowed values are \code{"unsigned"},
#' \code{"signed"}, \code{"signed hybrid"}, \code{"distance"}.
#' @param power Soft thresholding power, integer used in the function of the
#' network type.
#' @param datExpr data.frame containing expression data. Columns correspond to
#' genes and rows to samples.
#' @param selectCols For correlation networks only (see below) can be used to
#' select genes whose adjacencies will be calculated. Should be either a
#' numeric vector giving the indices of the genes to be used, or a boolean
#' vector indicating which genes are to be used.
#' @param corFnc Character string specifying the function to be used to
#' calculate co - expression similarity for correlation networks. Defaults to
#' Pearson correlation. Any function returning values between - 1 and 1 can be
#' used.
#' @param corOptions Character string specifying additional arguments to be
#' passed to the function given by \code{corFnc}. Use \code{"use = 'p', method
#' = 'spearman'"} to obtain Spearman correlation.
#' @param distFnc Character string specifying the function to be used to
#' calculate co - expression similarity for distance networks. Defaults to the
#' function \code{\link{dist}}. Any function returning non - negative values
#' can be used.
#' @param distOptions Character string specifying additional arguments to be
#' passed to the function given by \code{distFnc}. For example, when the
#' function \code{\link{dist}} is used, the argument \code{method} can be used
#' to specify various ways of computing the distance..
#' @return Adjacency matrix of dimensions \code{ncol(datExpr)} times
#' \code{ncol(datExpr)} (or the same dimensions as \code{similarity}). If
#' \code{selectCols} was given, the number of columns will be the length (if
#' numeric) or sum (if boolean) of \code{selectCols}.
#' @note When calculated from the \code{datExpr}, the network is always
#' calculated among the columns of \code{datExpr} irrespective of whether a
#' correlation or a distance network is requested.
#' @author Peter Langfelder and Steve Horvath
#' @references Bin Zhang and Steve Horvath (2005) A General Framework for
#' Weighted Gene Co-Expression Network Analysis, Statistical Applications in
#' Genetics and Molecular Biology, Vol. 4 No. 1, Article 17
#'
#' Langfelder P, Horvath S (2007) Eigengene networks for studying the
#' relationships between co - expression modules. BMC Systems Biology 2007,
#' 1:54
#' @keywords array
#' @examples
#' similarity <- matrix(seq(-1, 1, length.out = 25), 5)
#' similarity[lower.tri(similarity)] <- t(similarity)[lower.tri(similarity)]
#' adj <-  adjacency.fromSimilarity(similarity)
#'
#' datExpr <- matrix(seq(-1, 1, length.out = 25), 5)
#' datExpr[lower.tri(datExpr)] <- t(datExpr)[lower.tri(datExpr)]
#' adj <- adjacency(datExpr)
#'
#' @export adjacency
adjacency <- function(datExpr, selectCols = NULL, type = "unsigned",
                      power = if (type == "distance") 1 else 6,
                      corFnc = "cor", corOptions = "use = 'p'",
                      distFnc = "dist", distOptions = "method = 'euclidean'") {
    intType = charmatch(type, .adjacencyTypes)
    if (is.na(intType)) {
        stop(paste("Unrecognized 'type'. Recognized values are",
                   paste(.adjacencyTypes, collapse = ", ")))
    }

    if (intType < 4) {
        if (is.null(selectCols)) {
            corExpr <- parse(text = paste(corFnc, "(datExpr ",
                                         prepComma(corOptions), ")"))
            # cor_mat = cor(datExpr, use = "p")
            cor_mat <- eval(corExpr)
        } else {
            corExpr <- parse(text = paste(corFnc,
                                         "(datExpr, datExpr[, selectCols] ",
                                         prepComma(corOptions), ")"))
            #cor_mat = cor(datExpr, datExpr[, selectCols], use = "p")
            cor_mat <- eval(corExpr)
        }
    } else {
        if (!is.null(selectCols)) {
            stop("The argument 'selectCols' cannot be used for distance ",
                 "adjacency.")
        }
        corExpr <- parse(text = paste(distFnc, "(t(datExpr) ",
                                     prepComma(distOptions), ")"))
        # cor_mat = cor(datExpr, use = "p")
        d <- eval(corExpr)
        if (any(d < 0)) {
            warning("Distance function returned (some) negative values.")
        }
        cor_mat <- 1 - as.matrix((d/max(d, na.rm = TRUE))^2)
    }

    if (intType == 1) {
        cor_mat <- abs(cor_mat)
    } else if (intType == 2) {
        cor_mat <- (1 + cor_mat)/2
    } else if (intType == 3) {
        cor_mat[cor_mat < 0] <- 0
    }
    cor_mat^power
}

#' @rdname adjacency
#' @inheritParams adjacency
#' @description
#' similarity calculates the pairwise correlation of the expression data.
#' @export
similarity <- function(datExpr, selectCols = NULL, corFnc = "cor",
                       corOptions = "use = 'p'", type = "cor",
                       distFnc = "dist", distOptions = "method = 'euclidean'") {
    intType <- ifelse(type == "cor", 1, 3)
    if (intType < 2) {
        if (is.null(selectCols)) {
            corExpr = parse(text = paste(corFnc,
                                         "(datExpr ",
                                         prepComma(corOptions), ")"))
            # cor_mat = cor(datExpr, use = "p")
            cor_mat <- eval(corExpr)
        } else {
            corExpr <- parse(text = paste(
                corFnc, "(datExpr, datExpr[, selectCols] ",
                prepComma(corOptions), ")"))
            #cor_mat = cor(datExpr, datExpr[, selectCols], use = "p")
            cor_mat <- eval(corExpr)
        }
    } else {
        if (!is.null(selectCols)) {
            stop("The argument 'selectCols' cannot be used for distance ",
                 "adjacency.")
        }
        corExpr <- parse(text = paste(distFnc, "(t(datExpr) ",
                                     prepComma(distOptions), ")"))
        # cor_mat = cor(datExpr, use = "p")
        d <- eval(corExpr)
        if (any(d < 0)) {
            warning("Distance function returned (some) negative values.")
        }
        cor_mat <- 1 - as.matrix((d/max(d, na.rm = TRUE))^2)
    }

    checkAdjMat(cor_mat, min = -1)
    cor_mat
}

# adjacency.polyReg ####
#' Adjacency matrix based on polynomial regression
#'
#' adjacency.polyReg calculates a network adjacency matrix by fitting
#' polynomial regression models to pairs of variables (i.e. pairs of columns
#' from \code{datExpr}). Each polynomial fit results in a model fitting index
#' R.squared.  Thus, the n columns of \code{datExpr} result in an n x n
#' dimensional matrix whose entries contain R.squared measures. This matrix is
#' typically non-symmetric. To arrive at a (symmetric) adjacency matrix, one
#' can specify different symmetrization methods with
#' \code{symmetrizationMethod}.
#'
#' A network adjacency matrix is a symmetric matrix whose entries lie between 0
#' and 1. It is a special case of a similarity matrix. Each variable (column of
#' \code{datExpr}) is regressed on every other variable, with each model
#' fitting index recorded in a square matrix. Note that the model fitting index
#' of regressing variable x and variable y is usually different from that of
#' regressing y on x.  From the polynomial regression model glm(y ~
#' poly(x,degree)) one can calculate the model fitting index R.squared(y,x).
#' R.squared(y,x) is a number between 0 and 1. The closer it is to 1, the
#' better the polynomial describes the relationship between x and y and the
#' more significant is the pairwise relationship between the 2 variables. One
#' can also reverse the roles of x and y to arrive at a model fitting index
#' R.squared(x,y). If \code{degree}>1 then R.squared(x,y) is typically
#' different from R.squared(y,x). Assume a set of n variables x1,...,xn
#' (corresponding to the columns of \code{datExpr} then one can define
#' R.squared(xi,xj). The model fitting indices for the elements of an n x n
#' dimensional matrix (R.squared(ij)).  \code{symmetrizationMethod} implements
#' the following symmetrization methods:
#' A.min(ij)=min(R.squared(ij),R.squared(ji)),
#' A.ave(ij)=(R.squared(ij)+R.squared(ji))/2,
#' A.max(ij)=max(R.squared(ij),R.squared(ji)).
#'
#' @inheritParams adjacency
#' @param degree the degree of the polynomial. Must be less than the number of
#' unique points.
#' @param symmetrizationMethod character string (eg "none", "min","max","mean")
#' that specifies the method used to symmetrize the pairwise model fitting
#' index matrix (see details).
#' @return An adjacency matrix of dimensions ncol(datExpr) times ncol(datExpr).
#' @author Lin Song, Steve Horvath
#' @seealso For more information about polynomial regression, please refer to
#' functions \code{\link{poly}} and \code{\link{glm}}
#' @references Song L, Langfelder P, Horvath S Avoiding mutual information
#' based co-expression measures (to appear).
#'
#' Horvath S (2011) Weighted Network Analysis. Applications in Genomics and
#' Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8
#' @keywords misc
#' @examples
#' #Simulate a data frame datE which contains 5 columns and 50 observations
#' m <- 50
#' x1 <- rnorm(m)
#' r <- 0.5
#' x2 <- r * x1 + sqrt(1 - r^2) * rnorm(m)
#' r <- 0.3
#' x3 <- r * (x1 - 0.5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
#' x4 <- rnorm(m)
#' r <- 0.3
#' x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
#' datE <- data.frame(x1, x2, x3, x4, x5)
#' #calculate adjacency by symmetrizing using max
#' A.max <- adjacency.polyReg(datE, symmetrizationMethod = "max")
#' A.max
#' #calculate adjacency by symmetrizing using mean
#' A.mean <- adjacency.polyReg(datE, symmetrizationMethod = "mean")
#' A.mean
#' #calculate adjacency by symmetrizing using min
#' A.min <- adjacency.polyReg(datE, symmetrizationMethod = "min")
#' A.min
#' # output the unsymmetrized pairwise model fitting indices R.squared
#' R.squared <- adjacency.polyReg(datE, symmetrizationMethod = "none")
#' R.squared
adjacency.polyReg <- function(datExpr, degree = 3,
                              symmetrizationMethod = "mean") {

    if (!is.element(symmetrizationMethod, c("none", "min", "max", "mean"))) {
        stop("Unrecognized symmetrization method.")
    }

    datExpr <- matrix(as.numeric(as.matrix(datExpr)), nrow(datExpr),
                      ncol(datExpr))
    n <- ncol(datExpr)
    polyRsquare <- matrix(NA, n,n)

    for (i in 2:n) {
        for (j in 1:(i-1)) {
            del <- is.na(datExpr[, i] + datExpr[,j])
            if (sum(del) >= (n-1) | var(datExpr[, i], na.rm = T) == 0 |
                var(datExpr[, j], na.rm = T) == 0) {
                polyRsquare[i, j] = polyRsquare[j, i] = NA
            } else {
                dati <- datExpr[!del, i]
                datj <- datExpr[!del, j]
                lmPij <- glm(dati ~ poly( datj, degree))
                polyRsquare[i, j] <- cor(dati, predict(lmPij))^2
                lmPji <- glm( datj ~ poly( dati, degree))
                polyRsquare[j, i] <- cor( datj, predict(lmPji))^2
                rm(dati, datj, lmPij, lmPji)
            }
        }
    }

    diag(polyRsquare) = rep(1,n)

    if (symmetrizationMethod =="none") {
        adj <- polyRsquare
    } else {
        adj <- switch(symmetrizationMethod,
                     min = pmin(polyRsquare, t(polyRsquare)),
                     max = pmax(polyRsquare, t(polyRsquare)),
                     mean = (polyRsquare + t(polyRsquare))/2)
        checkAdjMat(adj)
    }
    adj

}

# adjacency.splineReg ####
#' Calculate network adjacency based on natural cubic spline regression
#'
#' adjacency.splineReg calculates a network adjacency matrix by fitting spline
#' regression models to pairs of variables (i.e. pairs of columns from
#' \code{datExpr}). Each spline regression model results in a fitting index
#' R.squared.  Thus, the n columns of \code{datExpr} result in an n x n
#' dimensional matrix whose entries contain R.squared measures. This matrix is
#' typically non-symmetric. To arrive at a (symmetric) adjacency matrix, one
#' can specify different symmetrization methods with
#' \code{symmetrizationMethod}.
#'
#' A network adjacency matrix is a symmetric matrix whose entries lie between 0
#' and 1. It is a special case of a similarity matrix. Each variable (column of
#' \code{datExpr}) is regressed on every other variable, with each model
#' fitting index recorded in a square matrix. Note that the model fitting index
#' of regressing variable x and variable y is usually different from that of
#' regressing y on x.  From the spline regression model glm( y ~ ns( x, df))
#' one can calculate the model fitting index R.squared(y,x).  R.squared(y,x) is
#' a number between 0 and 1. The closer it is to 1, the better the spline
#' regression model describes the relationship between x and y and the more
#' significant is the pairwise relationship between the 2 variables. One can
#' also reverse the roles of x and y to arrive at a model fitting index
#' R.squared(x,y). R.squared(x,y) is typically different from R.squared(y,x).
#' Assume a set of n variables x1,...,xn (corresponding to the columns of
#' \code{datExpr}) then one can define R.squared(xi,xj). The model fitting
#' indices for the elements of an n x n dimensional matrix (R.squared(ij)).
#' \code{symmetrizationMethod} implements the following symmetrization methods:
#' A.min(ij)=min(R.squared(ij),R.squared(ji)),
#' A.ave(ij)=(R.squared(ij)+R.squared(ji))/2,
#' A.max(ij)=max(R.squared(ij),R.squared(ji)). For more information about
#' natural cubic spline regression, please refer to functions "ns" and "glm".
#'
#' @inheritParams adjacency
#' @param df degrees of freedom in generating natural cubic spline. The default
#' is as follows: if nrow(datExpr)>100 use 6, if nrow(datExpr)>30 use 4,
#' otherwise use 5.
#' @param symmetrizationMethod character string (eg "none", "min","max","mean")
#' that specifies the method used to symmetrize the pairwise model fitting
#' index matrix (see details).
#' @param ...  other arguments from function \code{\link[splines]{ns}}
#' @return An adjacency matrix of dimensions ncol(datExpr) times ncol(datExpr).
#' @author Lin Song, Steve Horvath
#' @seealso \code{\link[splines]{ns}}, \code{\link{glm}}
#' @references Song L, Langfelder P, Horvath S Avoiding mutual information
#' based co-expression measures (to appear).
#'
#' Horvath S (2011) Weighted Network Analysis. Applications in Genomics and
#' Systems Biology. Springer Book. ISBN: 978-1-4419-8818-8
#' @keywords misc
#' @examples
#'
#' #Simulate a data frame datE which contains 5 columns and 50 observations
#' m <- 50
#' x1 <- rnorm(m)
#' r <- 0.5
#' x2 <- r * x1 + sqrt(1 - r ^ 2) * rnorm(m)
#' r <- 0.3
#' x3 <- r *(x1 - .5) ^ 2 + sqrt(1 - r ^ 2) * rnorm(m)
#' x4 <- rnorm(m)
#' r <- 0.3
#' x5 <- r * x4 + sqrt(1 - r ^ 2) * rnorm(m)
#' datE <- data.frame(x1, x2, x3, x4, x5)
#' #calculate adjacency by symmetrizing using max
#' A.max <- adjacency.splineReg(datE, symmetrizationMethod = "max")
#' A.max
#' #calculate adjacency by symmetrizing using mean
#' A.mean <- adjacency.splineReg(datE, symmetrizationMethod = "mean")
#' A.mean
#' #calculate adjacency by symmetrizing using min
#' A.min <- adjacency.splineReg(datE, symmetrizationMethod = "min")
#' A.min
#' # output the unsymmetrized pairwise model fitting indices R.squared
#' R.squared <- adjacency.splineReg(datE, symmetrizationMethod = "none")
#' R.squared
adjacency.splineReg <- function(datExpr,
                               df = 6-(nrow(datExpr)<100)-(nrow(datExpr)<30),
                               symmetrizationMethod = "mean", ...) {

    if (!is.element(symmetrizationMethod, c("none", "min" ,"max", "mean"))) {
        stop("Unrecognized symmetrization method.")
    }

    datExpr <- matrix(as.numeric(as.matrix(datExpr)),
                     nrow(datExpr), ncol(datExpr))
    n <- ncol(datExpr)
    splineRsquare <- matrix(NA, n, n)

    for (i in 2:n) {
        for (j in 1:(i-1)) {
            del = is.na(datExpr[, i]+datExpr[,j])
            if (sum(del) >= (n-1) | var(datExpr[, i], na.rm=T) == 0 |
                var(datExpr[, j], na.rm = T) == 0) {
                splineRsquare[i, j] = splineRsquare[j, i] = NA
            } else {
                dati <- datExpr[!del, i]
                datj <- datExpr[!del, j]
                lmSij <- glm( dati ~ ns( datj, df = df, ...))
                splineRsquare[i, j] <- cor(dati, predict(lmSij))^2
                lmSji <- glm( datj ~ ns(dati, df = df, ...))
                splineRsquare[j, i] <- cor(datj, predict(lmSji))^2
                rm(dati, datj, lmSij, lmSji)
            }
        }
    }

    diag(splineRsquare) <- rep(1, n)
    if (symmetrizationMethod == "none") {
        adj <- splineRsquare
    } else {
        adj <- switch(symmetrizationMethod,
                     min = pmin(splineRsquare, t(splineRsquare)),
                     max = pmax(splineRsquare, t(splineRsquare)),
                     mean = (splineRsquare + t(splineRsquare))/2)
        checkAdjMat(adj)
    }
    adj
}

# sigmoidAdjacency ####
#' Sigmoid-type adacency function.
#'
#' Sigmoid-type function that converts a similarity to a weighted network
#' adjacency.
#'
#' The sigmoid adjacency function is defined as \eqn{1/(1+\exp[-\alpha(ss -
#' \mu)])}{1/(1 + exp(-alpha * (ss - mu)))}.
#'
#' @param ss similarity, a number between 0 and 1. Can be given as a scalar,
#' vector or a matrix.
#' @param mu shift parameter.
#' @param alpha slope parameter.
#' @return
#' Adjacencies returned in the same form as the input \code{ss}
#' @author Steve Horvath
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
sigmoidAdjacency <- function(ss, mu = 0.8, alpha = 20) {
    if (any(ss > 1 | ss < 0)) {
        stop("ss should be a between 0 and 1")
    }
    adj <- 1/(1 + exp(- alpha * (ss - mu)))
    checkAdjMat(adj)
    adj
}

# SignumAdjacency ####
#' Hard-thresholding adjacency function
#'
#' This function transforms correlations or other measures of similarity into
#' an unweighted network adjacency.
#'
#'
#' @param corMat a matrix of correlations or other measures of similarity.
#' @param threshold threshold for connecting nodes: all nodes whose
#' \code{corMat} is above the threshold will be connected in the resulting
#' network.
#' @return An unweighted adjacency matrix of the same dimensions as the input
#' \code{corMat}.
#' @author Steve Horvath
#' @seealso \code{\link{adjacency}} for soft-thresholding and creating weighted
#' networks.
#' @references Bin Zhang and Steve Horvath (2005) "A General Framework for
#' Weighted Gene Co-Expression Network Analysis", Statistical Applications in
#' Genetics and Molecular Biology: Vol. 4: No. 1, Article 17
#' @keywords misc
#' @examples
#' datExpr <- matrix(rnorm(150), ncol = 5)
#' corMat <- cor(datExpr)
#' signumAdjacency(corMat, 0.1)
signumAdjacency <- function(corMat, threshold) {
    adjmat <- as.matrix(abs(corMat) >= threshold)
    dimnames(adjmat) <- dimnames(corMat)
    diag(adjmat) <- 0
    checkAdjMat(adjmat)
    adjmat
}
