#' Get column names from MaxQuant peptideGroups file
#'
#' @param mqFile The path to a MaxQuant peptideGroups.txt file.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A character vector with column names (converted to
#'     valid names) for the input file.
#'
#' @examples
#' colnm <- getColumnNames(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"))
#'
#' @importFrom utils read.delim
#'
getColumnNames <- function(mqFile) {
    .assertScalar(x = mqFile, type = "character")
    stopifnot(file.exists(mqFile))

    names(utils::read.delim(mqFile, nrows = 2))
}

#' Extract intensity columns from MaxQuant peptideGroups file
#'
#' Extract all column names in the \code{mqFile} that match the provided
#' \code{iColPattern}. Also allows specifying regular expressions to
#' define specific columns to retain or exclude. All column names in the
#' file can be listed using the \code{getColumnNames} function.
#'
#' @param mqFile The path to a MaxQuant peptideGroups.txt file.
#' @param iColPattern Character scalar defining a regular expression to
#'     identify sample columns.
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude.
#' @param stopIfEmpty Logical scalar, whether to raise an error if no
#'     columns matching the patterns are found.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A list with two character vectors: \code{iColsAll}, which contains
#'     all column names matching the \code{iColPattern}, and \code{iCols},
#'     which contains the remaining column names after applying the
#'     selection imposed by \code{includeOnlySamples} and
#'     \code{excludeSamples}.
#'
#' @examples
#' icols <- getIntensityColumns(system.file("extdata", "mq_example",
#'                                          "1356_proteinGroups.txt",
#'                                          package = "einprot"),
#'                              iColPattern = "^LFQ\\.intensity\\.")
#'
#' @importFrom utils read.delim
#' @importFrom stringr str_extract
#'
getIntensityColumns <- function(mqFile, iColPattern,
                                includeOnlySamples = "",
                                excludeSamples = "", stopIfEmpty = FALSE) {
    .assertScalar(x = mqFile, type = "character")
    stopifnot(file.exists(mqFile))
    .assertScalar(x = iColPattern, type = "character")
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    .assertScalar(x = stopIfEmpty, type = "logical")
    ## Can only specify one of includeOnlySamples and excludeSamples
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }

    ## --------------------------------------------------------------------- ##
    ## Columns representing intensities
    ## --------------------------------------------------------------------- ##
    iColsAll <- grep(iColPattern, getColumnNames(mqFile),
                     value = TRUE)
    if (stopIfEmpty && length(iColsAll) == 0) {
        stop("No samples were found matching the specified iColPattern.")
    }

    ## --------------------------------------------------------------------- ##
    ## Samples to include/exclude
    ## --------------------------------------------------------------------- ##
    if (length(includeOnlySamples) > 1 || includeOnlySamples != "") {
        ## Specify samples to include
        iCols <- iColsAll[!is.na(stringr::str_extract(
            iColsAll, paste(includeOnlySamples, collapse = "|")))]
    } else if (length(excludeSamples) > 1 || excludeSamples != "") {
        ## Specify samples to exclude
        iCols <- iColsAll[is.na(stringr::str_extract(
            iColsAll, paste(excludeSamples, collapse = "|")))]
    } else {
        iCols <- iColsAll
    }

    if (stopIfEmpty && length(iCols) == 0) {
        stop("No samples were retained - please check the specification ",
             "of includeOnlySamples/excludeSamples.")
    }

    list(iColsAll = iColsAll, iCols = iCols)
}
