#' Get column names
#'
#' Utility function to retrieve column names from quantification text file.
#'
#' @param inFile Path to a tab-delimited input text file (e.g. \code{MaxQuant}
#'     peptideGroups.txt or \code{ProteomeDiscoverer} Proteins.txt).
#' @param ... Additional arguments passed on to \code{utils::read.delim()}.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A character vector with column names (converted to
#'     valid names) for the input file.
#'
#' @examples
#' colnm <- getColumnNames(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"))
#' colnm
#'
#' @importFrom utils read.delim
#'
getColumnNames <- function(inFile, ...) {
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))

    names(utils::read.delim(inFile, nrows = 2, ...))
}

# Helper function to get iCols from iColsAll
#' @importFrom stringr str_extract
.getiCols <- function(iColsAll, includeOnlySamples, excludeSamples,
                      stopIfEmpty) {
    .assertScalar(x = stopIfEmpty, type = "logical")
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    ## Can only specify one of includeOnlySamples and excludeSamples
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }

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
    iCols
}



#' Extract intensity columns from intensity file
#'
#' Extract all column names in the \code{inFile} that match the provided
#' \code{iColPattern}. Also allows specifying regular expressions to
#' define specific columns to retain or exclude. All column names in the
#' file can be listed using the \code{getColumnNames} function.
#'
#' @param inFile Path to a tab-delimited input text file (e.g. \code{MaxQuant}
#'     peptideGroups.txt or \code{ProteomeDiscoverer} Proteins.txt).
#' @param iColPattern Character scalar defining a regular expression to
#'     identify intensity columns.
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude samples.
#' @param stopIfEmpty Logical scalar, whether to raise an error if no
#'     columns matching the patterns are found.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A list with two character vectors: \code{iColsAll}, which contains
#'     all column names matching the \code{iColPattern}, and \code{iCols},
#'     which contains the remaining column names after applying the
#'     selection imposed by \code{includeOnlySamples} or
#'     \code{excludeSamples}.
#'
#' @examples
#' icols <- getIntensityColumns(system.file("extdata", "mq_example",
#'                                          "1356_proteinGroups.txt",
#'                                          package = "einprot"),
#'                              iColPattern = "^LFQ\\.intensity\\.",
#'                              excludeSamples = "Adnp")
#' icols
#'
getIntensityColumns <- function(inFile, iColPattern,
                                includeOnlySamples = "",
                                excludeSamples = "", stopIfEmpty = FALSE) {
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = iColPattern, type = "character")
    .assertScalar(x = stopIfEmpty, type = "logical")

    ## -------------------------------------------------------------------------
    ## Columns representing intensities
    ## -------------------------------------------------------------------------
    iColsAll <- grep(iColPattern, getColumnNames(inFile), value = TRUE)
    if (stopIfEmpty && length(iColsAll) == 0) {
        stop("No samples were found matching the specified iColPattern.")
    }

    ## -------------------------------------------------------------------------
    ## Samples to include/exclude
    ## -------------------------------------------------------------------------
    iCols <- .getiCols(iColsAll = iColsAll,
                       includeOnlySamples = includeOnlySamples,
                       excludeSamples = excludeSamples,
                       stopIfEmpty = stopIfEmpty)

    list(iColsAll = iColsAll, iCols = iCols)
}
