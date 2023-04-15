# Determine the name to use for the assays
#' @keywords internal
#' @noRd
.getAssayName <- function(pat) {
    patmatch1 <- gsub("^", "", gsub("\\", "", pat, fixed = TRUE), fixed = TRUE)
    patmatch2 <- gsub("$", "", gsub("\\", "", pat, fixed = TRUE), fixed = TRUE)
    if (patmatch1 == "LFQ.intensity.") {
        "LFQ.intensity"
    } else if (patmatch1 == "Intensity.") {
        "Intensity"
    } else if (patmatch1 == "iBAQ.") {
        "iBAQ"
    } else if (patmatch1 == "Abundance.") {
        "Abundance"
    } else if (patmatch1 == "Abundance.F[0-9]+.") {
        "Abundance"
    } else if (patmatch1 == "Abundance.F.+.Sample.") {
        "Abundance"
    } else if (patmatch1 == "Abundances.Count.") {
        "Abundances.count"
    } else if (patmatch1 == "Abundances.Count.F[0-9]+.") {
        "Abundances.count"
    } else if (patmatch1 == "Abundances.Count.F.+.Sample.") {
        "Abundances.count"
    } else if (patmatch1 == "Abundances.Normalized.") {
        "Abundances.normalized"
    } else if (patmatch1 == "Abundances.Normalized.F[0-9]+.") {
        "Abundances.normalized"
    } else if (patmatch1 == "Abundances.Normalized.F.+.Sample.") {
        "Abundances.normalized"
    } else if (patmatch1 == "Abundances.Grouped.Count.") {
        "Abundances.grouped.count"
    } else if (patmatch1 == "Abundances.Grouped.CV.in.Percent.") {
        "Abundances.grouped.CV"
    } else if (patmatch1 == "Abundances.Grouped.") {
        "Abundances.grouped"
    } else if (patmatch1 == "MS.MS.Count.") {
        "MS.MS.Count"
    } else if (patmatch1 == "Sequence.coverage.") {
        "Sequence.coverage"
    } else if (patmatch1 == "Unique.peptides.") {
        "Unique.peptides"
    } else if (patmatch1 == "Razor.+unique.peptides.") {
        "Razor.unique.peptides"
    } else if (patmatch1 == "Identification.type.") {
        "Identification.type"
    } else if (patmatch1 == "Peptides.") {
        "Peptides"
    } else if (patmatch2 == ".Unique.Spectral.Count") {
        "Unique.spectral.count"
    } else if (patmatch2 == ".Total.Spectral.Count") {
        "Total.spectral.count"
    } else if (patmatch2 == ".Spectral.Count") {
        "Spectral.count"
    } else if (patmatch2 == ".Unique.Intensity") {
        "Unique.intensity"
    } else if (patmatch2 == ".Total.Intensity") {
        "Total.intensity"
    } else if (patmatch2 == ".MaxLFQ.Unique.Intensity") {
        "MaxLFQ.unique.intensity"
    } else if (patmatch2 == ".MaxLFQ.Total.Intensity") {
        "MaxLFQ.total.intensity"
    } else if (patmatch2 == ".MaxLFQ.Intensity") {
        "MaxLFQ.intensity"
    } else if (patmatch2 == ".Intensity") {
        "Intensity"
    } else {
        "ERROR"
    }
}

#' Import an abundance file
#'
#' @param inFile The path to an input text file (e.g. MaxQuant
#'     peptideGroups.txt, PD Proteins.txt or FragPipe combined_protein.tsv).
#' @param iColPattern Character scalar defining a regular expression to
#'     identify sample columns. For MaxQuant output, this is typically
#'     one of "^iBAQ\\.", "^LFQ\\.intensity\\." or "^Intensity\\.". For PD,
#'     it is typically "^Abundance\\.", "^Abundance\\.F[0-9]+\\." or
#'     "^Abundance\\.F.+\\.Sample\\.". For FragPipe,
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
#' a character scalar with the main assay name (this assay contains the values
#' from the columns matching the provided \code{iColPattern}).
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")
#'
#' @importFrom QFeatures readSummarizedExperiment
#' @importFrom SummarizedExperiment rowData assay assayNames
#' @importFrom S4Vectors metadata
#' @importFrom methods as
#' @importFrom stringdist amatch
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
        "^Abundance\\.",
        "^Abundance\\.F[0-9]+\\.",
        "^Abundance\\.F.+\\.Sample\\.",
        "^Abundances\\.Count\\.",
        "^Abundances\\.Count\\.F[0-9]+\\.",
        "^Abundances\\.Count\\.F.+\\.Sample\\.",
        "^Abundances\\.Normalized\\.",
        "^Abundances\\.Normalized\\.F[0-9]+\\.",
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

    ## Without the escaped characters
    patsexp <- gsub("\\", "", pats, fixed = TRUE)

    ## Check input arguments
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = iColPattern, type = "character")
    if (!(iColPattern %in% c(pats, patsexp))) {
        closest_match <- stringdist::amatch(iColPattern, c(pats, patsexp),
                                            method = "lv", maxDist = 5)
        if (!is.na(closest_match)) {
            stop("Invalid iColPattern. Did you mean ",
                 c(pats, patsexp)[closest_match], "?\nValid values: ",
                 paste(c(pats, patsexp), collapse = ", "))
        } else {
            stop("Invalid iColPattern. Valid values: ",
                 paste(c(pats, patsexp), collapse = ", "))
        }
    }
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")

    ## Currently, don't allow columns where the regular expression can match
    ## also other columns
    if (iColPattern %in% c("\\.Spectral\\.Count$",
                           "\\.Total\\.Intensity$", "\\.Intensity$",
                           ".Spectral.Count$",
                           ".Total.Intensity$", ".Intensity$")) {
        stop("Specifying ", iColPattern, " as the main assay is currently ",
             "not supported.")
    }
    ## The exception is Abundances.Grouped - in this case,
    ## allow it but give a warning
    if (iColPattern %in% c("^Abundances\\.Grouped\\.",
                           "^Abundances.Grouped.")) {
        warning("Note that the specified iColPattern may match different ",
                "types of columns - please check the selected samples ",
                "carefully.")
    }

    ## Put the iColPattern as the first assay
    ## If provided without escaped periods, find the corresponding pattern with
    ## escaped periods
    if (iColPattern %in% patsexp) {
        pos <- which(patsexp == iColPattern)[1]
        if (gsub("\\", "", pats[pos], fixed = TRUE) != iColPattern) {
            ## Should not end up here
            #nocov start
            stop("An error occurred - unmatched pats/patsexp")
            #nocov end
        }
        iColPattern <- pats[pos]
    }
    pats <- unique(c(iColPattern, pats))

    ## For PD, if the pattern without 'Sample' is chosen, remove the
    ## ones with 'Sample' from the list of patterns (otherwise there would be
    ## multiple patterns corresponding to the same assay)
    if (iColPattern %in% c("^Abundance\\.F[0-9]+\\.", "^Abundances\\.Count\\.F[0-9]+\\.",
                           "^Abundances\\.Normalized\\.F[0-9]+\\.")) {
        pats <- pats[!(pats %in% c("^Abundance\\.Sample\\.",
                                   "^Abundance\\.F.+\\.Sample\\.",
                                   "^Abundances\\.Count\\.F.+\\.Sample\\.",
                                   "^Abundances\\.Normalized\\.F.+\\.Sample\\.",
                                   "^Abundance\\.",
                                   "^Abundances\\.Count\\.",
                                   "^Abundances\\.Normalized\\."))]
    } else if (iColPattern %in% c("^Abundance\\.", "^Abundances\\.Count\\.",
                                  "^Abundances\\.Normalized\\.")) {
        pats <- pats[!(pats %in% c("^Abundance\\.Sample\\.",
                                   "^Abundance\\.F.+\\.Sample\\.",
                                   "^Abundances\\.Count\\.F.+\\.Sample\\.",
                                   "^Abundances\\.Normalized\\.F.+\\.Sample\\.",
                                   "^Abundance\\.F[0-9]+\\.",
                                   "^Abundances\\.Count\\.F[0-9]+\\.",
                                   "^Abundances\\.Normalized\\.F[0-9]+\\."))]
    } else {
        pats <- pats[!(pats %in% c("^Abundance\\.",
                                   "^Abundance\\.F[0-9]+\\.",
                                   "^Abundances\\.Count\\.",
                                   "^Abundances\\.Count\\.F[0-9]+\\.",
                                   "^Abundances\\.Normalized\\.",
                                   "^Abundances\\.Normalized\\.F[0-9]+\\."))]
    }

    names(pats) <- vapply(pats, .getAssayName, "ERROR")

    if (any(names(pats) == "ERROR")) {
        ## Should never end up in here, as we check the validity of the
        ## specified assay above
        #nocov start
        stop("Unsupported assay")
        #nocov end
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
        ## If iColPattern is "^Abundances\\.Grouped\\.", remove grouped,
        ## grouped cv columns (also matched by the regex)
        if (pat %in% c("^Abundances\\.Grouped\\.",
                       "^Abundances.Grouped.")) {
            icols <- icols[!grepl(
                paste0("Abundances\\.Grouped\\.CV\\.in\\.Percent", "|",
                       "Abundances\\.Grouped\\.Count"), icols)]
        }

        if (length(icols) > 0) {
            se <- QFeatures::readSummarizedExperiment(
                inFile, ecol = icols, sep = "\t", ...
            )
            ## Add list of columns to metadata
            S4Vectors::metadata(se)$cols <- icols

            ## Remove column pattern and trailing periods from colnames
            colnames(se) <- sub("\\.+$", "", sub(pat, "", colnames(se)))

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

    ## Put together into a single SingleCellExperiment object
    assayList <- assayList[!vapply(assayList, is.null, TRUE)]

    if (any(duplicated(names(assayList)))) {
        ## Should never end up in here
        #nocov start
        warning("Multiple column patterns corresponding to the same assay name: ",
                paste(names(assayList)[duplicated(names(assayList))],
                      collapse = "; "))
        #nocov end
    }

    ## Get the SCE corresponding to the main assay
    aName <- .getAssayName(iColPattern)
    sce <- assayList[[aName]]
    SummarizedExperiment::assayNames(sce) <- aName
    S4Vectors::metadata(sce)$colList <- list()
    S4Vectors::metadata(sce)$colList[[aName]] <-
        S4Vectors::metadata(assayList[[aName]])$cols
    S4Vectors::metadata(sce)$cols <- NULL
    ## Add the other assays
    for (a in names(assayList)) {
        if (a != aName) {
            if (all(colnames(sce) %in% colnames(assayList[[a]])) &&
                all(rownames(sce) == rownames(assayList[[a]]))) {
                SummarizedExperiment::assay(sce, a) <-
                    SummarizedExperiment::assay(assayList[[a]][, colnames(sce)])
                S4Vectors::metadata(sce)$colList[[a]] <-
                    S4Vectors::metadata(assayList[[a]])$cols
            }
        }
    }

    ## Clean up column names of row data
    colnames(SummarizedExperiment::rowData(sce)) <-
        sub("\\.$", "",
             gsub("\\.+", ".", colnames(SummarizedExperiment::rowData(sce))))

    sce <- methods::as(sce, "SingleCellExperiment")

    return(list(sce = sce, aName = aName))
}
