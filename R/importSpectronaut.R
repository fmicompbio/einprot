#' Import data from Spectronaut into a SingleCellExperiment object
#'
#' Import data from a Spectronaut quantification file into a
#' \code{SingleCellExperiment} object. Note that Spectronaut support in einprot
#' is currently experimental - please be aware that the interface may change,
#' and interpret results with caution.
#'
#' @param inFile Path to a tab-delimited input text file from Spectronaut.
#' @param outLevel Character scalar indicating the desired output level;
#'     currently only \code{"pg"} (protein group) is supported.
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude samples.
#' @param stopIfEmpty Logical scalar, whether to raise an error if no
#'     columns matching the patterns are found.
#' @param aName Character scalar giving the
#'     column from which to get the values for the main assay.
#' @param ... Additional arguments that will be passed on to
#'     \code{QFeatures::readSummarizedExperiment} (e.g., the number of rows
#'     to import).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A list with two elements: a \code{SingleCellExperiment} object and
#' a character scalar with the main assay name.
#'
importSpectronaut <- function(inFile, outLevel = "pg",
                              includeOnlySamples = "",
                              excludeSamples = "", stopIfEmpty = FALSE,
                              aName = "PG.Quantity", ...) {
    ## Check input arguments
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = outLevel, type = "character",
                  validValues = c("pg"))
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    .assertScalar(x = stopIfEmpty, type = "logical")
    .assertScalar(x = aName, type = "character")

    tmp <- read.delim(inFile, header = TRUE, sep = "\t")
    stopifnot(aName %in% colnames(tmp))
    iColsAll <- unique(tmp$R.FileName)
    iCols <- .getiCols(iColsAll = iColsAll,
                       includeOnlySamples = includeOnlySamples,
                       excludeSamples = excludeSamples,
                       stopIfEmpty = stopIfEmpty)
    tmp <- tmp[tmp$R.FileName %in% iCols, ]

    if (outLevel == "pg") {
        tmp <- tmp[, unique(c("R.FileName", "PG.ProteinGroups", "PG.ProteinNames",
                              aName,
                              grep("PG\\.", colnames(tmp), value = TRUE)))]
        # tmp <- tmp %>%
        #     dplyr::group_by(dplyr::across(c(-PG.ProteinNames))) %>%
        #     dplyr::summarize(PG.ProteinNames =
        #                          paste(unique(unlist(strsplit(.data$PG.ProteinNames, ";"))), collapse = ";"),
        #                      .groups = "drop")
        tmp <- unique(tmp)
        rd <- DataFrame(PG.ProteinGroups = unique(tmp$PG.ProteinGroups))
        aL <- list()
        for (nm in setdiff(colnames(tmp), c("R.FileName", "PG.ProteinGroups"))) {
            ## Check if it's a sample-specific value or not
            tmpsub <- tmp %>%
                dplyr::select(dplyr::all_of(c("PG.ProteinGroups", nm)))
            tmpcount <- tmpsub %>%
                dplyr::count(.data$PG.ProteinGroups)
            if (all(tmpcount$n == 1)) {
                ## One value per protein group -> annotation
                tmpsub <- tmpsub %>%
                    dplyr::distinct()
                rd[[nm]] <- tmpsub[[nm]][match(rd$PG.ProteinGroups, tmpsub$PG.ProteinGroups)]
            } else {
                ## One value per protein group/sample -> assay
                fv <- NA
                # fv <- ifelse(is.numeric(tmp[[nm]]), 0,
                #              ifelse(is.logical(tmp[[nm]]), NA, NA_character_)
                tmpsub <- tmp %>%
                    dplyr::select(dplyr::all_of(c("R.FileName", "PG.ProteinGroups", nm))) %>%
                    tidyr::pivot_wider(names_from = "R.FileName",
                                       values_from = dplyr::all_of(nm),
                                       values_fill = fv) %>%
                    as.data.frame()
                rownames(tmpsub) <- tmpsub$PG.ProteinGroups
                tmpsub$PG.ProteinGroups <- NULL
                tmpsub <- as.matrix(tmpsub)
                aL[[nm]] <- tmpsub[match(rd$PG.ProteinGroups, rownames(tmpsub)),
                                   match(iCols, colnames(tmpsub)), drop = FALSE]
            }
        }
        ## Add annotation columns corresponding to specific assays
        for (aNm in c("PG.Genes", "PG.ProteinNames", "PG.ProteinAccessions",
                        "PG.ProteinDescriptions", "PG.UniProtIds")) {
            if (aNm %in% names(aL)) {
                ## Aggregate across columns for each row
                rd[[aNm]] <- apply(aL[[aNm]], 1, function(x) {
                    paste(setdiff(unique(unlist(strsplit(x, ";"))), NA), collapse = ";")
                })
            }
        }

        sce <- SingleCellExperiment::SingleCellExperiment(
            assays = aL,
            rowData = rd,
            colData = DataFrame(sampleName = iCols)
        )

        ## Add list of columns to metadata
        S4Vectors::metadata(sce)$colList <- list()
        return(list(sce = sce, aName = aName))
    } else {
        stop("Not yet implemented")
    }
}
