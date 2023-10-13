#' @export
#'
importDIANN <- function(inFile, fileType = "pg_matrix", outLevel = "pg",
                        includeOnlySamples = "", excludeSamples = "",
                        stopIfEmpty = FALSE, aName = "MaxLFQ", ...) {
    if (fileType == "pg_matrix") {
        if (outLevel == "pr") {
            stop("To obtain precursor output, please provide either the ",
                 "pr_matrix or the main report.")
        }
        cn <- getColumnNames(inFile, check.names = FALSE)
        iColsAll <- grep("\\.mzML$|\\.wiff$|\\.dia$|\\.d$|\\.raw$|\\|/", cn, value = TRUE)
        if (length(includeOnlySamples) > 1 || includeOnlySamples != "") {
            ## Specify samples to include
            iCols <- iColsAll[!is.na(stringr::str_extract(
                iColsAll, paste(includeOnlySamples,
                                collapse = "|")))]
        } else if (length(excludeSamples) > 1 || excludeSamples != "") {
            ## Specify samples to exclude
            iCols <- iColsAll[is.na(stringr::str_extract(
                iColsAll, paste(excludeSamples,
                                collapse = "|")))]
        } else {
            iCols <- iColsAll
        }
        if (stopIfEmpty && length(iCols) == 0) {
            stop("No samples were retained - please check the specification ",
                 "of includeOnlySamples/excludeSamples.")
        }

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

    } else if (fileType == "pr_matrix") {
        stop("Not yet implemented")
    } else if (fileType == "main_report") {
        stop("Not yet implemented")
    }
}
