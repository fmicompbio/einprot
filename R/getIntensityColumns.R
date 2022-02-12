#' @export
#' @author Charlotte Soneson
#'
#' @importFrom utils read.delim
#' @importFrom stringr str_extract
#'
getIntensityColumns <- function(mqFile, iColPattern,
                                includeOnlySamples,
                                excludeSamples, stopIfEmpty = FALSE) {
    ## Columns representing intensities
    iColsAll <- grep(iColPattern, names(utils::read.delim(mqFile, nrows = 2)),
                     value = TRUE)
    if (stopIfEmpty && length(iColsAll) == 0) {
        stop("No samples were found matching the specified iColPattern.")
    }

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
             "of includeOnlySamples/excludeSamples (note that these should ",
             "specify individual samples, not group names).")
    }

    list(iColsAll = iColsAll, iCols = iCols)
}
