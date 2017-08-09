#  networkScreeningGS ####


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

# networkScreening ####


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
#' @param weightESy weight to use for the trait as an additional eigengene
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
                             getQValues = TRUE) {
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

# automaticNetworkScreening ####


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

# automaticNetworkScreeningGS ####


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
        hubGeneSignificance = data.frame(HGS1)
    )
    output
} # end of function automaticNetworkScreeningGS
