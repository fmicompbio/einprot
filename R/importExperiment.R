#' @keywords internal
#' @noRd
#' @importFrom dplyr case_when
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
            "Abundances.Normalized.F.+.Sample." ~ "Abundances.normalized",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Abundances.Grouped.Count." ~ "Abundances.grouped.count",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Abundances.Grouped.CV.in.Percent." ~ "Abundances.grouped.CV",
        gsub("\\^", "", gsub("\\\\", "", pat)) ==
            "Abundances.Grouped." ~ "Abundances.grouped",
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
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Unique.Spectral.Count" ~ "Unique.spectral.count",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Total.Spectral.Count" ~ "Total.spectral.count",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Spectral.Count" ~ "Spectral.count",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Unique.Intensity" ~ "Unique.intensity",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Total.Intensity" ~ "Total.intensity",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".MaxLFQ.Unique.Intensity" ~ "MaxLFQ.unique.intensity",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".MaxLFQ.Total.Intensity" ~ "MaxLFQ.total.intensity",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".MaxLFQ.Intensity" ~ "MaxLFQ.intensity",
        gsub("\\$", "", gsub("\\\\", "", pat)) ==
            ".Intensity" ~ "Intensity"
    )
}

#' Import an abundance file
#'
#' @param inFile The path to an input text file (e.g. MaxQuant
#'     peptideGroups.txt, PD Proteins.txt or FragPipe combined_protein.tsv).
#' @param iColPattern Character scalar defining a regular expression to
#'     identify sample columns. For MaxQuant output, this is typically
#'     one of "^iBAQ\\.", "^LFQ\\.intensity\\." or "^Intensity\\.". For PD,
#'     it is typically "^Abundance\\.F.+\\.Sample\\.". For FragPipe,
#'     it is typically "\\.MaxLFQ\\.Intensity$". Columns matching the
#'     given pattern will form the first assay in the output object.
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude samples.
#' @param ... Additional arguments that will be passed on to
#'     \code{QFeatures::readSummarizedExperiment} (e.g., the number of rows
#'     to import).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with two elements: a \code{SingleCellExperiment} object and
#' a character scalar with the main assay name (containing the columns
#' matching the provided \code{iColPattern}).
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")
#'
#' @importFrom QFeatures readSummarizedExperiment
#' @importFrom SummarizedExperiment rowData assay assayNames
#' @importFrom methods as
#'
importExperiment <- function(inFile, iColPattern, includeOnlySamples = "",
                             excludeSamples = "", ...) {
    ## List of assays to create/allowed column patterns
    pats <- c(## MaxQuant
        "^MS\\.MS\\.Count\\.", "^LFQ\\.intensity\\.",
        "^Intensity\\.", "^Sequence\\.coverage\\.",
        "^Unique\\.peptides\\.", "^Razor\\.+unique\\.peptides\\.",
        "^Peptides\\.", "^iBAQ\\.", "^Identification\\.type\\.",
        ## ProteomeDiscoverer
        "^Abundances\\.Count\\.F.+\\.Sample\\.",
        "^Abundance\\.F.+\\.Sample\\.",
        "^Abundances\\.Normalized\\.F.+\\.Sample\\.",
        "^Abundances\\.Grouped\\.Count\\.",
        "^Abundances\\.Grouped\\.CV\\.in\\.Percent\\.",
        "^Abundances\\.Grouped\\.",
        ## FragPipe
        "\\.Unique\\.Spectral\\.Count$",
        "\\.Total\\.Spectral\\.Count$",
        "\\.Spectral\\.Count$",
        "\\.Unique\\.Intensity$",
        "\\.Total\\.Intensity$",
        "\\.MaxLFQ\\.Unique\\.Intensity$",
        "\\.MaxLFQ\\.Total\\.Intensity$",
        "\\.MaxLFQ\\.Intensity$",
        "\\.Intensity$")

    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = iColPattern, type = "character",
                  validValues = pats)
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")

    ## Currently, don't allow Abundances.Grouped and similar columns -
    ## regular expression matches also other columns
    if (iColPattern %in% c("^Abundances\\.Grouped\\.", "\\.Spectral\\.Count$",
                           "\\.Total\\.Intensity$", "\\.Intensity$")) {
        stop("Specifying ", iColPattern, " as the main assay is currently ",
             "not supported.")
    }

    ## Put the iColPattern as the first assay
    pats <- unique(c(iColPattern, pats))
    names(pats) <- vapply(pats, .getAssayName, "ERROR")
    if (any(names(pats) == "ERROR")) {
        stop("Unsupported assay")
    }

    ## Create a list of SummarizedExperiment objects
    assayList <- lapply(pats, function(pat) {
        icols <- getIntensityColumns(
            inFile = inFile, iColPattern = pat,
            includeOnlySamples = includeOnlySamples,
            excludeSamples = excludeSamples, stopIfEmpty = FALSE)$iCols
        ## Don't consider summary columns (just the column pattern
        ## and one or more final periods) - keep these in rowData
        icols <- icols[!grepl(paste0("^", pat, "+$"), icols)]
        ## For FragPipe, don't consider the Combined column
        icols <- icols[!grepl(paste0("Combined", pat, "$"), icols)]

        if (length(icols) > 0) {
            se <- QFeatures::readSummarizedExperiment(
                inFile, ecol = icols, sep = "\t", ...
            )
            ## Remove column pattern and trailing periods from colnames
            colnames(se) <- gsub("\\.+$", "", gsub(pat, "", colnames(se)))

            ## Remove columns from rowData
            findCol <- grepl(paste(pats, collapse = "|"),
                             colnames(SummarizedExperiment::rowData(se))) &
                !grepl(paste(c(paste0("^", pats, "+$"),
                               paste0("Combined", pats, "$")), collapse = "|"),
                       colnames(SummarizedExperiment::rowData(se)))
            remCol <- colnames(SummarizedExperiment::rowData(se))[findCol]
            SummarizedExperiment::rowData(se) <-
                SummarizedExperiment::rowData(se)[, !findCol]

            return(se)
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
            if (all(colnames(sce) %in% colnames(assayList[[a]])) &&
                all(rownames(sce) == rownames(assayList[[a]]))) {
                SummarizedExperiment::assay(sce, a) <-
                    SummarizedExperiment::assay(assayList[[a]][, colnames(sce)])
            }
        }
    }

    ## Clean up column names of row data
    colnames(SummarizedExperiment::rowData(sce)) <-
        gsub("\\.$", "",
             gsub("\\.+", ".", colnames(SummarizedExperiment::rowData(sce))))

    sce <- methods::as(sce, "SingleCellExperiment")

    return(list(sce = sce, aName = aName))
}
