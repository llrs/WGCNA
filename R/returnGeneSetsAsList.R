#' Return pre-defined gene lists in several biomedical categories.
#' 
#' This function returns gene sets for use with other R functions.  These gene
#' sets can include inputted lists of genes and files containing user-defined
#' lists of genes, as well as a pre-made collection of brain, blood, and other
#' biological lists.  The function returns gene lists associated with each
#' category for use with other enrichment strategies (i.e., GSVA).
#' 
#' User-inputted files for fnIn can be in one of three formats:
#' 
#' 1) Text files (must end in ".txt") with one list per file, where the first
#' line is the list descriptor and the remaining lines are gene names
#' corresponding to that list, with one gene per line.  For example Ribosome
#' RPS4 RPS8 ...
#' 
#' 2) Gene / category files (must be csv files), where the first line is the
#' column headers corresponding to Genes and Lists, and the remaining lines
#' correspond to the genes in each list, for any number of genes and lists.
#' For example: Gene, Category RPS4, Ribosome RPS8, Ribosome ... NDUF1,
#' Mitohcondria NDUF3, Mitochondria ... MAPT, AlzheimersDisease PSEN1,
#' AlzheimersDisease PSEN2, AlzheimersDisease ...
#' 
#' 3) Module membership (kME) table in csv format.  Currently, the module
#' assignment is the only thing that is used, so as long as the Gene column is
#' 2nd and the Module column is 3rd, it doesn't matter what is in the other
#' columns.  For example, PSID, Gene, Module, <other columns> <psid>, RPS4,
#' blue, <other columns> <psid>, NDUF1, red, <other columns> <psid>, RPS8,
#' blue, <other columns> <psid>, NDUF3, red, <other columns> <psid>, MAPT,
#' green, <other columns> ...
#' 
#' @param fnIn A vector of file names containing user-defined lists.  These
#' files must be in one of three specific formats (see details section).  The
#' default (NULL) may only be used if one of the "use_____" parameters is TRUE.
#' @param catNmIn A vector of category names corresponding to each fnIn.  This
#' name will be appended to each overlap corresponding to that filename.  The
#' default sets the category names as the corresponding file names.
#' @param useBrainLists If TRUE, a pre-made set of brain-derived enrichment
#' lists will be added to any user-defined lists for enrichment comparison.
#' The default is FALSE.  See references section for related references.
#' @param useBloodAtlases If TRUE, a pre-made set of blood-derived enrichment
#' lists will be added to any user-defined lists for enrichment comparison.
#' The default is FALSE.  See references section for related references.
#' @param useStemCellLists If TRUE, a pre-made set of stem cell (SC)-derived
#' enrichment lists will be added to any user-defined lists for enrichment
#' comparison.  The default is FALSE.  See references section for related
#' references.
#' @param useBrainRegionMarkers If TRUE, a pre-made set of enrichment lists for
#' human brain regions will be added to any user-defined lists for enrichment
#' comparison.  The default is FALSE.  These lists are derived from data from
#' the Allen Human Brain Atlas (http://human.brain-map.org/).  See references
#' section for more details.
#' @param useImmunePathwayLists If TRUE, a pre-made set of enrichment lists for
#' immune system pathways will be added to any user-defined lists for
#' enrichment comparison.  The default is FALSE.  These lists are derived from
#' the lab of Daniel R Saloman.  See references section for more details.
#' @param geneSubset A vector of gene (or other) identifiers.  If entered, only
#' genes in this list will be returned in the output, otherwise all genes in
#' each category will be returned (default, geneSubset=NULL).
#' @return \item{geneSets}{ A list of categories in alphabetical order, where
#' each compnent of the list is a character vector of all genes corresponding
#' to the named category.  For example: geneSets =
#' list(category1=c("gene1","gene2"),category2=c("gene3","gene4","gene5")) }
#' @author Jeremy Miller
#' @references Please see the help file for userListEnrichment in the WGCNA
#' library for references for the pre-defined lists.
#' @keywords misc
#' @examples
#' 
#' # Example: Return a list of genes for various immune pathways
#' geneSets   = returnGeneSetsAsList(useImmunePathwayLists=TRUE)
#' geneSets[7:8]
#' 
returnGeneSetsAsList <- function (fnIn = NULL, catNmIn = fnIn, useBrainLists = FALSE, useBloodAtlases = FALSE, 
    useStemCellLists = FALSE, useBrainRegionMarkers = FALSE, useImmunePathwayLists = FALSE, geneSubset=NULL) 
{
    if (length(catNmIn) < length(fnIn)) {
        catNmIn = c(catNmIn, fnIn[(length(catNmIn) + 1):length(fnIn)])
        write("WARNING: not enough category names.  \n\t\t\t   Naming remaining categories with file names.", 
            "")
    }
    if (is.null(fnIn) & (!(useBrainLists | useBloodAtlases | 
        useStemCellLists | useBrainRegionMarkers | useImmunePathwayLists))) 
        stop("Either enter user-defined lists or set one of the use_____ parameters to TRUE.")
    glIn = NULL
    if (length(fnIn) > 0) {
        for (i in 1:length(fnIn)) {
            ext = substr(fnIn[i], nchar(fnIn[i]) - 2, nchar(fnIn[i]))
            if (ext == "csv") {
                datIn = read.csv(fnIn[i])
                if (colnames(datIn)[2] == "Gene") {
                  datIn = datIn[, 2:3]
                }
                else {
                  datIn = datIn[, 1:2]
                }
            }
            else {
                datIn = scan(fnIn[i], what = "character", sep = "\n")
                datIn = cbind(datIn[2:length(datIn)], datIn[1])
            }
            colnames(datIn) = c("Gene", "Category")
            datIn[, 2] = paste(datIn[, 2], catNmIn[i], sep = "__")
            glIn = rbind(glIn, datIn)
        }
        glIn = cbind(glIn, Type = rep("User", nrow(glIn)))
    }
    if (useBrainLists) {
        if (!(exists("BrainLists"))) 
            BrainLists = NULL
        data("BrainLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding brain list references.", 
            "")
        glIn = rbind(glIn, cbind(BrainLists, Type = rep("Brain", 
            nrow(BrainLists))))
    }
    if (useBloodAtlases) {
        if (!(exists("BloodLists"))) 
            BloodLists = NULL
        data("BloodLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding blood atlas references.", 
            "")
        glIn = rbind(glIn, cbind(BloodLists, Type = rep("Blood", 
            nrow(BloodLists))))
    }
    if (useStemCellLists) {
        if (!(exists("SCsLists"))) 
            SCsLists = NULL
        data("SCsLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding stem cell list references.", 
            "")
        glIn = rbind(glIn, cbind(SCsLists, Type = rep("StemCells", 
            nrow(SCsLists))))
    }
    if (useBrainRegionMarkers) {
        if (!(exists("BrainRegionMarkers"))) 
            BrainRegionMarkers = NULL
        data("BrainRegionMarkers", envir = sys.frame(sys.nframe()))
        write("Brain region markers from http://human.brain-map.org/ -- See userListEnrichment help file for details.", 
            "")
        glIn = rbind(glIn, cbind(BrainRegionMarkers, Type = rep("HumanBrainRegions", 
            nrow(BrainRegionMarkers))))
    }
    if (useImmunePathwayLists) {
        if (!(exists("ImmunePathwayLists"))) 
            ImmunePathwayLists = NULL
        data("ImmunePathwayLists", envir = sys.frame(sys.nframe()))
        write("See userListEnrichment help file for details regarding immune pathways.", 
            "")
        glIn = rbind(glIn, cbind(ImmunePathwayLists, Type = rep("Immune", 
            nrow(ImmunePathwayLists))))
    }
    removeDups = unique(paste(as.character(glIn[, 1]), as.character(glIn[, 
        2]), as.character(glIn[, 3]), sep = "@#$%"))
    if (length(removeDups) < length(glIn[, 1])) 
        glIn = t(as.matrix(as.data.frame(strsplit(removeDups, 
            "@#$%", fixed = TRUE))))
    geneIn = as.character(glIn[, 1])
    labelIn = paste(as.character(glIn[, 2]),as.character(glIn[, 3]),sep="__")
	if(!is.null(geneSubset)){
      keep = is.element(geneIn, geneSubset)
      geneIn = geneIn[keep]
      labelIn = labelIn[keep]
	}
	if(length(geneIn)<2)
	   stop("Please include a larger geneSubset, or set geneSubset=NULL.")
    allLabels <- sort(unique(labelIn))
	geneSet <- list()
	for (i in 1:length(allLabels))  geneSet[[i]] = geneIn[labelIn==allLabels[i]]
	names(geneSet) = allLabels
    return(geneSet)
}
