.getAssayName <- function(pat) {
    dplyr::case_when(
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "LFQ.intensity." ~ "LFQ.intensity",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Intensity." ~ "Intensity",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "iBAQ." ~ "iBAQ",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Abundance.F.+.Sample." ~ "Abundance",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Abundances.Count.F.+.Sample." ~ "Abundances.count",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "MS.MS.Count." ~ "MS.MS.Count",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Sequence.coverage." ~ "Sequence.coverage",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Unique.peptides." ~ "Unique.peptides",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Razor.+unique.peptides." ~ "Razor.unique.peptides",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Identification.type." ~ "Identification.type",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Peptides." ~ "Peptides",
    )
}

#' Import a protein_groups.txt from MaxQuant or Proteins.txt file from PD
#'
#' @param inFile The path to an input file (e.g. MaxQuant
#'     peptideGroups.txt or PD Proteins.txt).
#' @param iColPattern Character scalar defining a regular expression to
#'     identify sample columns.
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude.
#' @param ... Additional arguments that will be passed on to
#'     \code{QFeatures::readSummarizedExperiment}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with two elements: a \code{SummarizedExperiment} and
#' a character scalar with the main assay name.
#'
#' @importFrom QFeatures readSummarizedExperiment
#' @importFrom SummarizedExperiment rowData assay assayNames
#' @importFrom methods as
#'
importExperiment <- function(inFile, iColPattern, includeOnlySamples = "",
                             excludeSamples = "", ...) {
    ## List of assays to create/allowed column patterns
    pats <- c("^MS\\.MS\\.Count\\.", "^LFQ\\.intensity\\.",
              "^Intensity\\.", "^Sequence\\.coverage\\.",
              "^Unique\\.peptides\\.", "^Razor\\.+unique\\.peptides\\.",
              "^Peptides\\.", "^iBAQ\\.", "^Identification\\.type\\.",
              "^Abundances\\.Count\\.F.+\\.Sample\\.",
              "^Abundance\\.F.+\\.Sample\\.")

    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = iColPattern, type = "character",
                  validValues = pats)
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")

    ## Put the iColPattern as the first assay
    pats <- unique(c(iColPattern, pats))
    names(pats) <- vapply(pats, .getAssayName, "ERROR")
    if (any(names(pats) == "ERROR")) {
        stop("Unsupported assay")
    }

    coln <- NULL
    ## Create a list of SummarizedExperiment objects
    assayList <- lapply(pats, function(pat) {
        icols <- getIntensityColumns(
            inFile = mqFile, iColPattern = pat,
            includeOnlySamples = includeOnlySamples,
            excludeSamples = excludeSamples, stopIfEmpty = FALSE)$iCols
        ## Don't consider summary columns (just the column pattern
        ## and one or more final periods) - keep these in rowData
        icols <- icols[!grepl(paste0(pat, "+$"), icols)]

        if (length(icols) > 0) {
            se <- QFeatures::readSummarizedExperiment(
                inFile, ecol = icols, sep = "\t", ...
            )
            ## Remove trailing periods from colnames
            colnames(se) <- gsub("\\.+$", "", gsub(pat, "", colnames(se)))

            ## Remove columns from rowData
            findCol <- grepl(paste(pats, collapse = "|"),
                             colnames(SummarizedExperiment::rowData(se))) &
                !grepl(paste(paste0(pats, "+$"), collapse = "|"),
                       colnames(SummarizedExperiment::rowData(se)))
            remCol <- colnames(SummarizedExperiment::rowData(se))[findCol]
            SummarizedExperiment::rowData(se) <-
                SummarizedExperiment::rowData(se)[, !findCol]

            ## Return SE if it is compatible with previous ones
            if (is.null(coln)) {
                coln <- colnames(se)
                rown <- rownames(se)
                return(se)
            } else {
                ## Only return SE of colnames agree
                if (all(colnames(se) == coln) && all(rownames(se) == rown)) {
                    return(se)
                } else {
                    return(NULL)
                }
            }
        } else {
            return(NULL)
        }
    })

    assayList <- assayList[!vapply(assayList, is.null, TRUE)]
    aName <- .getAssayName(iColPattern)
    sce <- assayList[[aName]]
    SummarizedExperiment::assayNames(sce) <- aName
    for (a in names(assayList)) {
        if (a != aName) {
            SummarizedExperiment::assay(sce, a) <-
                SummarizedExperiment::assay(assayList[[a]])
        }
    }

    ## Clean up column names of row data
    colnames(SummarizedExperiment::rowData(sce)) <-
        gsub("\\.$", "",
             gsub("\\.+", ".", colnames(SummarizedExperiment::rowData(sce))))

    sce <- methods::as(sce, "SingleCellExperiment")

    return(list(sce = sce, aName = aName))
}
