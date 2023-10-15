#' Import data from DIA-NN into a SingleCellExperiment object
#'
#' Import data from a DIA-NN quantification file into a
#' \code{SingleCellExperiment} object.
#'
#' @param inFile Path to a tab-delimited input text file from DIA-NN; either
#'     \code{pg_matrix.tsv}, \code{pr_matrix.tsv} or the main \code{report.tsv}.
#' @param fileType Character scalar indicating the type of input file; either
#'     \code{"pg_matrix"}, \code{"pr_matrix"} or \code{"main_report"}.
#' @param outLevel Character scalar indicating the desired output level;
#'     either \code{"pg"} (protein group) or \code{"pr"} (precursor).
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude samples.
#' @param stopIfEmpty Logical scalar, whether to raise an error if no
#'     columns matching the patterns are found.
#' @param aName Character scalar giving the name of the main assay (if
#'     \code{fileType} is \code{"pg_matrix"} or \code{"pr_matrix"}), or the
#'     column from which to get the values for the main assay (if
#'     \code{fileType} is \code{"main_report"}).
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
#' @examples
#' sceL <- importDIANN(system.file("extdata", "diann_example",
#'                                 "PXD028735.pg_matrix.tsv",
#'                                 package = "einprot"),
#'                     fileType = "pg_matrix", outLevel = "pg",
#'                     aName = "MaxLFQ")
#' sceL
#'
importDIANN <- function(inFile, fileType = "pg_matrix", outLevel = "pg",
                        includeOnlySamples = "",
                        excludeSamples = "", stopIfEmpty = FALSE,
                        aName = "MaxLFQ", ...) {
    ## Check input arguments
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = fileType, type = "character",
                  validValues = c("pg_matrix", "pr_matrix", "main_report"))
    .assertScalar(x = outLevel, type = "character",
                  validValues = c("pg", "pr"))
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    .assertScalar(x = stopIfEmpty, type = "logical")
    .assertScalar(x = aName, type = "character")

    if (fileType == "pg_matrix" && outLevel == "pr") {
        stop("To obtain precursor output, please provide either the ",
             "pr_matrix or the main report.")
    }
    if (fileType == "pr_matrix" && outLevel == "pg") {
        stop("To obtain protein group output, please provide either the ",
             "pg_matrix or the main report.")
    }

    if (fileType %in% c("pg_matrix", "pr_matrix")) {
        cn <- getColumnNames(inFile, check.names = FALSE)
        iColsAll <- grep("\\.mzML$|\\.wiff$|\\.dia$|\\.d$|\\.raw$|\\|/", cn, value = TRUE)
        iCols <- .getiCols(iColsAll = iColsAll,
                           includeOnlySamples = includeOnlySamples,
                           excludeSamples = excludeSamples,
                           stopIfEmpty = stopIfEmpty)

        sce <- QFeatures::readSummarizedExperiment(
            inFile, ecol = iCols, sep = "\t", check.names = FALSE, ...
        )

        SummarizedExperiment::assayNames(sce) <- aName

        ## Add list of columns to metadata
        S4Vectors::metadata(sce)$colList <- list()
        S4Vectors::metadata(sce)$colList[[aName]] <- iCols

        ## Remove directory name and extension from colnames
        colnames(sce) <- tools::file_path_sans_ext(basename(
            gsub("\\\\", "/", colnames(sce))))

        sce <- methods::as(sce, "SingleCellExperiment")
        return(list(sce = sce, aName = aName))

    } else if (fileType == "main_report") {
        tmp <- read.delim(inFile, header = TRUE, sep = "\t")
        stopifnot(aName %in% colnames(tmp))
        iColsAll <- unique(tmp$Run)
        iCols <- .getiCols(iColsAll = iColsAll,
                           includeOnlySamples = includeOnlySamples,
                           excludeSamples = excludeSamples,
                           stopIfEmpty = stopIfEmpty)
        tmp <- tmp[tmp$Run %in% iCols, ]

        if (outLevel == "pg") {
            tmp <- tmp[, unique(c("Run", "Protein.Group", "Protein.Ids",
                                  aName,
                                  grep("PG\\.", colnames(tmp), value = TRUE)))]
            tmp <- tmp %>%
                dplyr::group_by(dplyr::across(c(-Protein.Ids))) %>%
                dplyr::summarize(Protein.Ids =
                    paste(unique(unlist(strsplit(Protein.Ids, ";"))), collapse = ";"),
                    .groups = "drop")
            tmp <- unique(tmp)
            rd <- DataFrame(Protein.Group = unique(tmp$Protein.Group))
            aL <- list()
            for (nm in setdiff(colnames(tmp), c("Run", "Protein.Group"))) {
                ## Check if it's a sample-specific value or not
                tmpsub <- tmp %>%
                    dplyr::select(dplyr::all_of(c("Protein.Group", nm)))
                tmpcount <- tmpsub %>%
                    dplyr::count(Protein.Group)
                if (all(tmpcount$n == 1)) {
                    ## One value per protein group -> annotation
                    tmpsub <- tmpsub %>%
                        dplyr::distinct()
                    rd[[nm]] <- tmpsub[[nm]][match(rd$Protein.Group, tmpsub$Protein.Group)]
                } else {
                    ## One value per protein group/sample -> assay
                    fv <- ifelse(is.numeric(tmp[[nm]]), 0, "0")
                    tmpsub <- tmp %>%
                        dplyr::select(dplyr::all_of(c("Run", "Protein.Group", nm))) %>%
                        tidyr::pivot_wider(names_from = Run,
                                           values_from = .data[[nm]],
                                           values_fill = fv) %>%
                        as.data.frame()
                    rownames(tmpsub) <- tmpsub$Protein.Group
                    tmpsub$Protein.Group <- NULL
                    tmpsub <- as.matrix(tmpsub)
                    aL[[nm]] <- tmpsub[match(rd$Protein.Group, rownames(tmpsub)),
                                       match(iCols, colnames(tmpsub)), drop = FALSE]
                }
            }
            sce <- SingleCellExperiment::SingleCellExperiment(
                assays = aL,
                rowData = rd,
                colData = DataFrame(sample = iCols)
            )

            ## Add list of columns to metadata
            S4Vectors::metadata(sce)$colList <- list()
            return(list(sce = sce, aName = aName))
        } else {
            stop("Not yet implemented")
        }
    }
}
