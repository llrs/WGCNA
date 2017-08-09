# GOenrichmentAnalysis
#' Calculation of GO enrichment (experimental)
#'
#' WARNING: This function should be considered experimental. The arguments and
#' resulting values (in particular, the enrichment p-values) are not yet
#' finalized and may change in the future. The function should only be used to
#' get a quick and rough overview of GO enrichment in the modules in a data
#' set; for a publication-quality analysis, please use an established tool.
#'
#' Using Bioconductor's annotation packages, this function calculates
#' enrichments and returns terms with best enrichment values.
#'
#' This function is basically a wrapper for the annotation packages available
#' from Bioconductor. It requires the packages GO.db, AnnotationDbi, and
#' org.xx.eg.db, where xx is the code corresponding to the organism that the
#' user wishes to analyze (e.g., Hs for human Homo Sapiens, Mm for mouse Mus
#' Musculus etc). For each cluster specified in the input, the function
#' calculates all enrichments in the specified ontologies, and collects
#' information about the terms with highest enrichment. The enrichment p-value
#' is calculated using Fisher exact test. As background we use all of the
#' supplied genes that are present in at least one term in GO (in any of the
#' ontologies).
#'
#' For best results, the newest annotation libraries should be used. Because of
#' the way Bioconductor is set up, to get the newest annotation libraries you
#' may have to use the current version of R.
#'
#' According to http://www.geneontology.org/GO.evidence.shtml, the following
#' codes are used by GO: \preformatted{ Experimental Evidence Codes EXP:
#' Inferred from Experiment IDA: Inferred from Direct Assay IPI: Inferred from
#' Physical Interaction IMP: Inferred from Mutant Phenotype IGI: Inferred from
#' Genetic Interaction IEP: Inferred from Expression Pattern
#'
#' Computational Analysis Evidence Codes ISS: Inferred from Sequence or
#' Structural Similarity ISO: Inferred from Sequence Orthology ISA: Inferred
#' from Sequence Alignment ISM: Inferred from Sequence Model IGC: Inferred from
#' Genomic Context IBA: Inferred from Biological aspect of Ancestor IBD:
#' Inferred from Biological aspect of Descendant IKR: Inferred from Key
#' Residues IRD: Inferred from Rapid Divergence RCA: inferred from Reviewed
#' Computational Analysis
#'
#' Author Statement Evidence Codes TAS: Traceable Author Statement NAS:
#' Non-traceable Author Statement
#'
#' Curator Statement Evidence Codes IC: Inferred by Curator ND: No biological
#' Data available
#'
#' Automatically-assigned Evidence Codes IEA: Inferred from Electronic
#' Annotation
#'
#' Obsolete Evidence Codes NR: Not Recorded }
#'
#' @param labels cluster (module, group) labels of genes to be analyzed. Either
#' a single vector, or a matrix. In the matrix case, each column will be
#' analyzed separately; analyzing a collection of module assignments in one
#' function call will be faster than calling the function several tinmes. For
#' each row, the labels in all columns must correspond to the same gene
#' specified in \code{entrezCodes}.
#' @param entrezCodes Entrez (a.k.a. LocusLink) codes of the genes whose labels
#' are given in \code{labels}. A single vector; the i-th entry corresponds to
#' row i of the matrix \code{labels} (or to the i-the entry if \code{labels} is
#' a vector).
#' @param yeastORFs if \code{organism == "yeast"} (below), this argument can be
#' used to input yeast open reading frame (ORF) identifiers instead of Entrez
#' codes. Since the GO mappings for yeast are provided in terms of ORF
#' identifiers, this may lead to a more accurate GO enrichment analysis. If
#' given, the argument \code{entrezCodes} is ignored.
#' @param organism character string specifying the organism for which to
#' perform the analysis. Recognized values are (unique abbreviations of)
#' \code{"human", "mouse", "rat", "malaria", "yeast", "fly", "bovine", "worm",
#' "canine", "zebrafish", "chicken"}.
#' @param ontologies vector of character strings specifying GO ontologies to be
#' included in the analysis.  Can be any subset of \code{"BP", "CC", "MF"}. The
#' result will contain the terms with highest enrichment in each specified
#' category, plus a separate list of terms with best enrichment in all
#' ontologies combined.
#' @param evidence vector of character strings specifying admissible evidence
#' for each gene in its specific term, or "all" for all evidence codes. See
#' Details or http://www.geneontology.org/GO.evidence.shtml for available
#' evidence codes and their meaning.
#' @param includeOffspring logical: should genes belonging to the offspring of
#' each term be included in the term? As a default, only genes belonging
#' directly to each term are associated with the term. Note that the
#' calculation of enrichments with offspring included can be quite slow for
#' large data sets.
#' @param backgroundType specification of the background to use. Recognized
#' values are (unique abbreviations of) \code{"allGiven", "allInGO",
#' "givenInGO"}, meaning that the functions will take all genes given in
#' \code{labels} as backround (\code{"allGiven"}), all genes present in any of
#' the GO categories (\code{"allInGO"}), or the intersection of given genes and
#' genes present in GO (\code{"givenInGO"}). The default is recommended for
#' genome-wide enrichment studies.
#' @param removeDuplicates logical: should duplicate entries in
#' \code{entrezCodes} be removed? If \code{TRUE}, only the first occurence of
#' each unique Entrez code will be kept. The cluster labels \code{labels} will
#' be adjusted accordingly.
#' @param leaveOutLabel optional specifications of module labels for which
#' enrichment calculation is not desired. Can be a single label or a vector of
#' labels to be ignored. However, if in any of the sets no labels are left to
#' calculate enrichment of, the function will stop with an error.
#' @param nBestP specifies the number of terms with highest enrichment whose
#' detailed information will be returned.
#' @param pCut alternative specification of terms to be returned: all terms
#' whose enrichment p-value is more significant than \code{pCut} will be
#' returned. If \code{pCut} is given, \code{nBestP} is ignored.
#' @param nBiggest in addition to returning terms with highest enrichment,
#' terms that contain most of the genes in each cluster can be returned by
#' specifying the number of biggest terms per cluster to be returned. This may
#' be useful for development and testing purposes.
#' @param getTermDetails logical indicating whether detailed information on the
#' most enriched terms should be returned.
#' @param verbose integer specifying the verbosity of the function. Zero means
#' silent, positive values will cause the function to print progress reports.
#' @param indent integer specifying indentation of the diagnostic messages.
#' Zero means no indentation, each unit adds two spaces.
#' @return A list with the following components: \item{keptForAnalysis }{
#' logical vector with one entry per given gene. \code{TRUE} if the entry was
#' used for enrichment analysis. Depending on the setting of
#' \code{removeDuplicates} above, only a single entry per gene may be used. }
#'
#' \item{inGO }{ logical vector with one entry per given gene. \code{TRUE} if
#' the gene belongs to any GO term, \code{FALSE} otherwise. Also \code{FALSE}
#' for genes not used for the analysis because of duplication. }
#'
#' If input \code{labels} contained only one vector of labels, the following
#' components:
#'
#' \item{countsInTerms }{ a matrix whose rows correspond to given cluster, and
#' whose columns correspond to GO terms, contaning number of genes in the
#' intersection of the corresponding module and GO term. Row and column names
#' are set appropriately.}
#'
#' \item{enrichmentP}{a matrix whose rows correspond to given cluster, and
#' whose columns correspond to GO terms, contaning enrichment p-values of each
#' term in each cluster. Row and column names are set appropriately.}
#'
#' \item{bestPTerms}{a list of lists with each inner list corresponding to an
#' ontology given in \code{ontologies} in input, plus one component
#' corresponding to all given ontologies combined.  The name of each component
#' is set appropriately. Each inner list contains two components:
#' \code{enrichment} is a data frame containing the highest enriched terms for
#' each module; and \code{forModule} is a list of lists with one inner list per
#' module, appropriately named. Each inner list contains one component per
#' term. If input \code{getTermDeyails} is \code{TRUE}, this component is yet
#' another list and contains components \code{termName} (term name),
#' \code{enrichmentP} (enrichment P value), \code{termDefinition} (GO term
#' definition), \code{termOntology} (GO term ontology), \code{geneCodes}
#' (Entrez codes of module genes in this term), \code{genePositions} (indices
#' of the genes listed in \code{geneCodes} within the given \code{labels}).
#' Thus, to obtain information on say the second term of the 5th module in
#' ontology BP, one can look at the appropriate row of
#' \code{bestPTerms$BP$enrichment}, or one can reference
#' \code{bestPTerms$BP$forModule[[5]][[2]]}. The author of the function
#' apologizes for any confusion this structure of the output may cause. }
#'
#' \item{biggestTerms}{a list of the same format as \code{bestPTerms},
#' containing information about the terms with most genes in the module for
#' each supplied ontology. }
#'
#' If input \code{labels} contained more than one vector, instead of the above
#' components the return value contains a list named \code{setResults} that has
#' one component per given set; each component is a list containing the above
#' components for the corresponding set.
#' @author Peter Langfelder
#' @seealso Bioconductor's annotation packages such as GO.db and
#' organism-specific annotation packages such as org.Hs.eg.db.
#' @keywords misc
GOenrichmentAnalysis <- function(labels, entrezCodes, yeastORFs = NULL,
               organism = "human", ontologies = c("BP", "CC", "MF"),
               evidence = "all", includeOffspring = TRUE,
               backgroundType = "givenInGO", removeDuplicates = TRUE,
               leaveOutLabel = NULL, nBestP = 10, pCut = NULL, nBiggest = 0,
               getTermDetails = TRUE, verbose = 2, indent = 0 ) {

   sAF <- options("stringsAsFactors")
   options(stringsAsFactors = FALSE)
   on.exit(options(stringsAsFactors = sAF[[1]]), TRUE)

   organisms <- c("human", "mouse", "rat", "malaria", "yeast", "fly", "bovine",
                 "worm", "canine", "zebrafish", "chicken")
   allEvidence <-  c("EXP", "IDA", "IPI", "IMP", "IGI", "IEP", "ISS", "ISO",
                    "ISA", "ISM", "IGC", "IBA", "IBD", "IKR", "IRD", "RCA",
                    "TAS", "NAS", "IC", "ND", "IEA", "NR")
   allOntologies <- c("BP", "CC", "MF")

   backgroundTypes <- c("allGiven", "allInGO", "givenInGO")

   spaces <- indentSpaces(indent)
   orgInd <- pmatch(organism, organisms)
   if (is.na(orgInd)) {
     stop("Unrecognized 'organism' given. Recognized values are ",
                paste(organisms, collapse = ", "))
   }
   if (length(evidence) == 0) {
     stop("At least one valid evidence code must be given in 'evidence'.")
   }
   if (length(ontologies) == 0) {
     stop("At least one valid ontology code must be given in 'ontology'.")
   }
   if (evidence == "all") {
      evidence <- allEvidence
   }
   evidInd <- pmatch(evidence, allEvidence)
   if (sum(is.na(evidInd)) != 0) {
     stop("Unrecognized 'evidence' given. Recognized values are ",
                paste(allEvidence, collapse = ", "))
   }
   ontoInd <- pmatch(ontologies, allOntologies)
   if (sum(is.na(ontoInd)) != 0) {
     stop("Unrecognized 'ontologies' given. Recognized values are ",
                paste(allEvidence, collapse = ", "))
   }

   backT <- pmatch(backgroundType, backgroundTypes)
   if (is.na(backT)) {
     stop("Unrecognized 'backgroundType' given. Recognized values are ",
           paste(backgroundTypes,  collapse = ", "))
   }
   orgCodes <- c("Hs", "Mm", "Rn", "Pf", "Sc", "Dm", "Bt", "Ce", "Cf", "Dr",
                "Gg")
   orgExtensions <- c(rep(".eg", 4), ".sgd", rep(".eg", 6))
   reverseMap <- c(rep(".egGO2EG", 4), ".sgdGO2ORF", rep(".egGO2EG", 6))

   missingPacks <- NULL
   packageName <- paste0("org.", orgCodes[orgInd], orgExtensions[orgInd],
                         ".db")
   if (!eval(parse(text = "require(packageName, character.only = TRUE)"))){
     missingPacks <- c(missingPacks, packageName)
   }

   if (!eval(parse(text="require(GO.db)"))) {
     missingPacks <- c(missingPacks, "GO.db")
   }
   if (!is.null(missingPacks)) {
     stop("Could not load the requisite package(s) ",
           paste(missingPacks, collapse = ", "),
          ". Please install the package(s).")
   }
   if (verbose > 0) {
     printFlush(paste(spaces,
                      "GOenrichmentAnalysis: loading annotation data..."))
   }

   if (orgInd == 5) {
      # Yeast needs special care.
      if (!is.null(yeastORFs)) {
        entrezCodes <- yeastORFs
      } else {
        # Map the entrez IDs to yeast ORFs
        x <- eval(parse(text = "org.Sc.sgd:::org.Sc.sgdENTREZID"))
        # x <- org.Sc.sgd:::org.Sc.sgdENTREZID
        xx <- as.list(x[mapped_genes])
        allORFs <- names(xx)
        mappedECs <- as.character(sapply(xx, as.character))
        entrez2orf <- match(entrezCodes, mappedECs)
        fin <- is.finite(entrez2orf)
        newCodes <- paste("InvalidCode", c(1:length(entrezCodes)), sep = ".")
        newCodes[fin] <- allORFs[entrez2orf[fin]]
        entrezCodes <- newCodes
     }
   }

   labels <- as.matrix(labels)

   nSets <- ncol(labels)
   nGivenRaw <- nrow(labels)

   if (removeDuplicates) {
     # Restrict given entrezCodes such that each code is unique
     ECtab <- table(entrezCodes)
     uniqueEC <- names(ECtab)
     keepEC <- match(uniqueEC, entrezCodes)
     entrezCodes <- entrezCodes[keepEC]
     labels <- labels[keepEC, , drop = FALSE]
   } else {
     keepEC <- c(1:nGivenRaw)
   }
   egGO <- eval(parse(text = paste(packageName, ":::org.", orgCodes[orgInd],
                                   orgExtensions[orgInd],
                                  "GO", sep = "")))

   if (orgInd == 5) {
      mapped_genes <- as.character(do.call(match.fun("mappedkeys"),
                                           list(egGO)))
      encodes2mapped <- match(as.character(entrezCodes), mapped_genes)
   } else {
      mapped_genes <- as.numeric(as.character(do.call(match.fun("mappedkeys"),
                                                      list(egGO))))
      encodes2mapped <- match(as.numeric(entrezCodes), mapped_genes)
   }
   encMapped <- is.finite(encodes2mapped)
   nAllIDsInGO <- sum(encMapped)

   mapECodes <- entrezCodes[encMapped]
   mapLabels <- labels[encMapped, , drop = FALSE]
   nMappedGenes <- nrow(mapLabels)

   if (nMappedGenes == 0){
     stop("None of the supplied gene identifiers map to the GO database.\n",
          "Please make sure you have specified the correct organism (default",
          " is human).")
   }
   Go2eg <- eval(parse(text = paste("AnnotationDbi::as.list(", packageName,
                                    ":::org.", orgCodes[orgInd],
                                           reverseMap[orgInd],")", sep = "")))
   nTerms <- length(Go2eg)

   goInfo <- as.list(GO.db::GOTERM)
   if (length(goInfo) > 0) {
      orgGoNames <- names(Go2eg)
      dbGoNames <- as.character(sapply(goInfo, GOID))
      dbGoOntologies <- as.character(sapply(goInfo, Ontology))
   } else {
      dbGoNames <- ""
   }
   goOffSpr <- list()
   if (includeOffspring) {
     goOffSpr[[1]] <- as.list(GOBPOFFSPRING)
     goOffSpr[[2]] <- as.list(GOCCOFFSPRING)
     goOffSpr[[3]] <- as.list(GOMFOFFSPRING)
   }
   term2info <- match(names(Go2eg), names(goInfo))
   termOntologies <- dbGoOntologies[term2info]

   if (backT == 1) {
     nBackgroundGenes <- sum(!is.na(entrezCodes))
   } else {
     if (backT == 2) {
       nBackgroundGenes <- length(mapped_genes)
     } else {
       nBackgroundGenes <- nMappedGenes
     }
   }

   termCodes <- vector(mode="list", length = nTerms)
   collectGarbage()
   nExpandLength <- 0
   blockSize <- 3000; # For a more efficient concatenating of offspring genes
   nAllInTerm <- rep(0, nTerms)
   if (verbose > 0) {
      printFlush(paste(spaces, " ..of the", length(entrezCodes),
                      " Entrez identifiers submitted,", sum(encMapped),
                      "are mapped in current GO categories."))
      printFlush(paste(spaces, " ..will use", nBackgroundGenes,
                       "background genes for enrichment calculations."))
      cat(paste(spaces, " ..preparing term lists (this may take a while).."))
      pind <- initProgInd()
   }

   for (c in 1:nTerms)
       if (!is.na(Go2eg[[c]][[1]])) {
           {
               te <- as.character(names(Go2eg[[c]])); # Term evidence codes
               tc <- Go2eg[[c]]
               if (includeOffspring) {
                   termOffspring <- NULL
                   for (ont in 1:length(goOffSpr)) {
                       term2off <- match(names(Go2eg)[c],
                                         names(goOffSpr[[ont]]))
                       if (!is.na(term2off)){
                           termOffspring <- c(termOffspring,
                                              goOffSpr[[ont]][[term2off]])
                       }
                   }
        if (length(termOffspring) > 0){
           maxLen <- blockSize
           tex <- rep("", maxLen)
           tcx <- rep("", maxLen)
           ind <- 1
           len <- length(te)
           tex[ ind:len ] <- te
           tcx[ ind:len ] <- tc
           ind <- len + 1
           o2go <- match(termOffspring, as.character(names(Go2eg)))
           o2go <- o2go[is.finite(o2go)]
           if (length(o2go) > 0) for (o in 1:length(o2go)) if (!is.na(Go2eg[[o2go[o]]][[1]]))
           {
             #printFlush(paste("Have offspring for term", c, ": ", names(Go2eg)[c],
             #           Term(goInfo[[term2info[c]]])))
             newc <- Go2eg[[o2go[o]]]
             newe <- names(newc)
             newl <- length(newe)
             if ((len + newl) > maxLen)
             {
               nExpand <- ceiling( (len + newl - maxLen)/blockSize)
               maxLen <- maxLen + blockSize * nExpand
               tex <- c(tex, rep("", maxLen - length(tex)))
               tcx <- c(tcx, rep("", maxLen - length(tex)))
               nExpandLength <- nExpandLength + 1
             }
             tex[ind:(len + newl)] <- newe
             tcx[ind:(len + newl)] <- newc
             ind <- ind + newl
             len <- len + newl
           }
           te <- tex[1:len]
           tc <- tcx[1:len]
        }
      }
      use <- is.finite(match(te, evidence))
      if (orgInd == 5) {
          if (backT == 2) {
              termCodes[[c]] <- unique(as.character(tc[use]))
          } else
              termCodes[[c]] <- as.character(intersect(tc[use], mapECodes))
      } else {
          if (backT == 2)
          {
              termCodes[[c]] <- unique(as.character(tc[use]))
          } else
              termCodes[[c]] <- as.numeric(as.character(intersect(tc[use], mapECodes)))
      }
      nAllInTerm[c] <- length(termCodes[[c]])
      if ( (c %%50  == 0) & (verbose > 0)) pind <- updateProgInd(c/nTerms, pind)
   }
       }
   if (verbose > 0) {
      pind <- updateProgInd(1, pind)
      printFlush("")
   }
   if ((verbose > 5) & (includeOffspring))
      printFlush(paste(spaces, " ..diagnostic for the developer: offspring buffer was expanded",
                       nExpandLength, "times."))

   ftp <- function(...) { fisher.test(...)$p.value }

   setResults <- list()

   for (set in 1:nSets) {
      if (verbose > 0)
        printFlush(paste(spaces, " ..working on label set", set, ".."))
      labelLevels <- levels(factor(labels[, set]))
      if (!is.null(leaveOutLabel)) {
        keep <- !(labelLevels %in% as.character(leaveOutLabel))
        if (sum(keep) == 0)
          stop("No labels were kept after removing labels that are supposed to be ignored.")
        labelLevels <- labelLevels[keep]
      }
      nLabelLevels <- length(labelLevels)

      modCodes <- list()
      nModCodes <- rep(0, nLabelLevels)
      if (backT == 1) {
        for (ll in 1:nLabelLevels) {
           modCodes[[ll]] <- entrezCodes[labels[, set] == labelLevels[ll]]
           nModCodes[ll] <- length(modCodes[[ll]])
        }
      } else {
        for (ll in 1:nLabelLevels) {
           modCodes[[ll]] <- mapECodes[mapLabels[, set] == labelLevels[ll]]
           nModCodes[ll] <- length(modCodes[[ll]])
        }
      }

      countsInTerm <- matrix(0, nLabelLevels, nTerms)
      enrichment <- matrix(1, nLabelLevels, nTerms)

      for (ll in 1:nLabelLevels)
        countsInTerm[ll, ] <- sapply(lapply(termCodes, intersect, modCodes[[ll]]), length)

      nAllInTermMat <- matrix(nAllInTerm, nLabelLevels, nTerms, byrow = TRUE)
      nModCodesMat <- matrix(nModCodes, nLabelLevels, nTerms)
      tabArr <- array(c(countsInTerm, nAllInTermMat - countsInTerm,
                       nModCodesMat - countsInTerm,
                       nBackgroundGenes - nModCodesMat - nAllInTermMat + countsInTerm),
                     dim = c(nLabelLevels * nTerms, 2, 2))

      if (verbose > 0)
         printFlush(paste(spaces, "   ..calculating enrichments (this may also take a while).."))
      calculate <- c(countsInTerm) > 0
      enrichment[calculate] <- apply(tabArr[calculate, , ], 1, ftp, alternative = "g")

      dimnames(enrichment) <- list (labelLevels, names(Go2eg))
      dimnames(countsInTerm) <- list (labelLevels, names(Go2eg))

      bestPTerms <- list()
      modSizes <- table(labels[ !(labels[, set] %in% leaveOutLabel), set])

      if (!is.null(pCut) || nBestP > 0) {
         printFlush(paste(spaces, "   ..putting together terms with highest enrichment significance.."))
         nConsideredOntologies <- length(ontologies) + 1
         for (ont in 1:nConsideredOntologies) {
            if (ont == nConsideredOntologies) {
               ontTerms <- is.finite(match(termOntologies, ontologies))
               bestPTerms[[ont]] <- list(ontology = ontologies)
               names(bestPTerms)[ont] <- paste(ontologies, collapse = ", ")
            } else {
               ontTerms <- termOntologies == ontologies[ont]
               bestPTerms[[ont]] <- list(ontology = ontologies[ont])
               names(bestPTerms)[ont] <- ontologies[ont]
            }
            bestPTerms[[ont]]$enrichment <- NULL
            bestPTerms[[ont]]$termInfo <- list()
            nOntTerms <- sum(ontTerms)
            ontEnr <- enrichment[, ontTerms, drop = FALSE]
            order <- apply(ontEnr, 1, order)
            for (ll in 1:nLabelLevels) {
              if (!is.null(pCut)) {
                 reportTerms <- c(1:nTerms)[ontTerms][ontEnr[ll, ] < pCut]
                 reportTerms <- reportTerms[order(ontEnr[ll, ][reportTerms])]
              } else
                 reportTerms <-  c(1:nTerms)[ontTerms][order[1:nBestP, ll]]
              nRepTerms <- length(reportTerms)
              enrTab <- data.frame(module = rep(labelLevels[ll], nRepTerms),
                                  modSize = rep(modSizes[ll], nRepTerms),
                                  bkgrModSize = rep(nModCodes[ll], nRepTerms),
                                  rank = seq(length.out = nRepTerms),
                                  enrichmentP = enrichment[ll, reportTerms],
                                  BonferoniP = pmin(rep(1, nRepTerms),
                                                    enrichment[ll, reportTerms] * nOntTerms),
                                  nModGenesInTerm = countsInTerm[ll, reportTerms],
                                  fracOfBkgrModSize = countsInTerm[ll, reportTerms]/nModCodes[ll],
                                  fracOfBkgrTermSize =
                                       countsInTerm[ll, reportTerms]/nAllInTerm[reportTerms],
                                  bkgrTermSize = nAllInTerm[reportTerms],
                                  termID = names(Go2eg)[reportTerms],
                                  termOntology =  rep("NA", nRepTerms),
                                  termName = rep("NA", nRepTerms),
                                  termDefinition =  rep("NA", nRepTerms))
              bestPTerms[[ont]]$forModule[[ll]] <- list()
              for (rci in seq(length.out = nRepTerms)) {
                 term <- reportTerms[rci]
                 termID <- names(Go2eg)[term]
                 dbind <- match(termID, dbGoNames)
                 if (is.finite(dbind)) {
                   enrTab$termName[rci] <- eval(parse(text = "AnnotationDbi:::Term(goInfo[[dbind]])"))
                   enrTab$termDefinition[rci] =
                         eval(parse(text = "AnnotationDbi:::Definition(goInfo[[dbind]])"))
                   enrTab$termOntology[rci] =
                         eval(parse(text = "AnnotationDbi:::Ontology(goInfo[[dbind]])"))
                 }
                 if (getTermDetails) {
                   geneCodes <- intersect(modCodes[[ll]], termCodes[[term]])
                   bestPTerms[[ont]]$forModule[[ll]][[rci]] <- list(termID = termID,
                           termName = enrTab$termName[rci],
                           enrichmentP = enrTab$enrichmentP[rci],
                           termDefinition = enrTab$termDefinition[rci],
                           termOntology = enrTab$termOntology[rci],
                           geneCodes = geneCodes,
                           genePositions = keepEC[match(geneCodes, entrezCodes)])
                }
              }
              if (ll == 1) {
                 bestPTerms[[ont]]$enrichment <- enrTab
              } else {
                 bestPTerms[[ont]]$enrichment <- rbind(bestPTerms[[ont]]$enrichment, enrTab)
              }
            }
         }
      }

      biggestTerms <- list()

      if (nBiggest > 0) {
         printFlush(paste(spaces, "   ..putting together terms with largest number of genes in modules.."))
         nConsideredOntologies <- length(ontologies) + 1
         for (ont in 1:nConsideredOntologies) {
            if (ont == nConsideredOntologies) {
               ontTerms <- is.finite(match(termOntologies, ontologies))
               biggestTerms[[ont]] <- list(ontology = ontologies)
               names(biggestTerms)[ont] <- paste(ontologies, collapse = ", ")
            } else {
               ontTerms <- termOntologies == ontologies[ont]
               biggestTerms[[ont]] <- list(ontology = ontologies[ont])
               names(biggestTerms)[ont] <- ontologies[ont]
            }
            biggestTerms[[ont]]$enrichment <- NULL
            biggestTerms[[ont]]$termInfo <- list()
            nOntTerms <- sum(ontTerms)
            ontNGenes <- countsInTerm[, ontTerms, drop = FALSE]
            order <- apply(-ontNGenes, 1, order)
            for (ll in 1:nLabelLevels) {
              reportTerms <- c(1:nTerms)[ontTerms][order[1:nBiggest, ll]]
              nRepTerms <- length(reportTerms)
              enrTab <- data.frame(module = rep(labelLevels[ll], nRepTerms),
                                  modSize = rep(modSizes[ll], nRepTerms),
                                  bkgrModSize = rep(nModCodes[ll], nRepTerms),
                                  rank = c(1:nRepTerms),
                                  enrichmentP = enrichment[ll, reportTerms],
                                  BonferoniP = pmin(rep(1, nRepTerms),
                                                    enrichment[ll, reportTerms] * nOntTerms),
                                  nModGenesInTerm = countsInTerm[ll, reportTerms],
                                  fracOfModSize = countsInTerm[ll, reportTerms]/nModCodes[ll],
                                  fracOfBkgrTermSize =
                                       countsInTerm[ll, reportTerms]/nAllInTerm[reportTerms],
                                  bkgrTermSize = nAllInTerm[reportTerms],
                                  termID = names(Go2eg)[reportTerms],
                                  termOntology =  rep("NA", nRepTerms),
                                  termName = rep("NA", nRepTerms),
                                  termDefinition =  rep("NA", nRepTerms))
              biggestTerms[[ont]]$forModule[[ll]] <- list()
              for (rci in seq(length.out = nRepTerms)) {
                 term <- reportTerms[rci]
                 termID <- names(Go2eg)[term]
                 dbind <- match(termID, dbGoNames)
                 if (is.finite(dbind)) {
                   enrTab$termName[rci] <- eval(parse(text="AnnotationDbi:::Term(goInfo[[dbind]])"))
                   enrTab$termDefinition[rci] =
                       eval(parse(text="AnnotationDbi:::Definition(goInfo[[dbind]])"))
                   enrTab$termOntology[rci] <- eval(parse(text="AnnotationDbi:::Ontology(goInfo[[dbind]])"))
                 }
                 if (getTermDetails) {
                   geneCodes <- intersect(modCodes[[ll]], termCodes[[term]])
                   biggestTerms[[ont]]$forModule[[ll]][[rci]] <- list(termID = termID,
                           termName = enrTab$termName[rci],
                           enrichmentP = enrTab$enrichmentP[rci],
                           termDefinition = enrTab$termDefinition[rci],
                           termOntology = enrTab$termOntology[rci],
                           geneCodes = geneCodes,
                           genePositions = keepEC[match(geneCodes, entrezCodes)])
                }
              }
              if (ll == 1) {
                 biggestTerms[[ont]]$enrichment <- enrTab
              } else {
                 biggestTerms[[ont]]$enrichment <- rbind(
                     biggestTerms[[ont]]$enrichment, enrTab)
              }
           }
         }
      }
      setResults[[set]] <- list(countsInTerm = countsInTerm,
                               enrichmentP = enrichment,
                               bestPTerms = bestPTerms,
                               biggestTerms = biggestTerms)
   }

   inGO <- rep(FALSE, nGivenRaw)
   inGO[keepEC] <- encMapped
   kept <- rep(FALSE, nGivenRaw)
   kept[keepEC] <- TRUE
   if (nSets == 1) {
      list(keptForAnalysis = kept,
           inGO = inGO,
           countsInTerms = setResults[[1]]$countsInTerm,
           enrichmentP = setResults[[1]]$enrichmentP,
           bestPTerms = setResults[[1]]$bestPTerms,
           biggestTerms = setResults[[1]]$biggestTerms)
   } else {
      list(keptForAnalysis = kept,
           inGO = inGO,
           setResults = setResults)
   }
}
