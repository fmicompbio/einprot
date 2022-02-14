#' Make feature IDs (row names) unique
#'
#' Make the feature IDs (row names) of \code{qft} unique by first extracting
#' the first entry in the \code{geneIdCol} and \code{proteinIdCol} columns
#' (multiple entries for each row are separated by semicolons), and then
#' using the gene ID as the feature ID if it exists and is unique, and
#' otherwise appending the protein ID. If it's still not unique, append an
#' integer to the name.
#'
#' @param qft A \code{QFeatures} object.
#' @param geneIdCol,proteinIdCol Character scalars indicating which columns of
#'     \code{rowData(qft[[1]])} should be used as gene/protein identifiers,
#'     respectively.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{QFeatures} object with modified, unique row names (see
#'     description for how these are generated).
#'
#' @importFrom SummarizedExperiment rowData
#'
fixFeatureIds <- function(qft, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs") {
    ## Extract the first annotated gene name
    gName <- vapply(strsplit(
        SummarizedExperiment::rowData(qft[[1]])[[geneIdCol]], ";"),
        .subset, 1, FUN.VALUE = "NA")
    rowData(qft[[1]])$geneIdSingle <- gName

    ## Extract the first annotated majority protein ID
    majProtID <- vapply(strsplit(
        SummarizedExperiment::rowData(qft[[1]])[[proteinIdCol]], ";"),
        .subset, 1, FUN.VALUE = "NA")
    rowData(qft[[1]])$proteinIdSingle <- majProtID

    ## Find features with missing gene name and replace with the corresponding
    ## majority protein ID
    idxna <- which(is.na(gName))
    gName[idxna] <- majProtID[idxna]

    ## Add the derived gNames to the rowData - these will be used for STRING
    rowData(qft[[1]])$IDsForSTRING <- gName

    ## If there are duplicated gene IDs, make them unique by appending the
    ## respective majority protein ID
    idxdup <- which(duplicated(gName))
    idxdup <- which(gName %in% gName[idxdup])
    gName[idxdup] <- paste0(gName[idxdup], ".", majProtID[idxdup])

    ## Check that there are no duplicated IDs and set as row names
    stopifnot(all(!duplicated(gName)))
    rownames(qft[[1]]) <- gName

    qft
}
