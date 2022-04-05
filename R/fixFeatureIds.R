#' Make feature IDs (row names) unique
#'
#' Make the feature IDs (row names) of \code{sce} unique by first extracting
#' the first entry in the \code{primaryIdCol} and \code{secondaryIdCol} columns
#' (multiple entries for each row are separated by semicolons), and then
#' using the primary ID as the feature ID if it exists and is unique, and
#' otherwise appending the secondary ID. If it's still not unique, append an
#' integer to the name.
#'
#' @param sce A \code{SummarizedExperiment} object (or derivative).
#' @param primaryIdCol,secondaryIdCol Character scalars indicating which columns of
#'     \code{rowData(sce)} should be used as primary/secondary identifiers,
#'     respectively.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return An object of the same type as \code{sce} with modified, unique row
#'     names (see description for how these are generated).
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#' sce <- fixFeatureIds(sce, primaryIdCol = "Gene.names",
#'                      secondaryIdCol = "Majority.protein.IDs")
#' head(rownames(sce))
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
fixFeatureIds <- function(sce, primaryIdCol = "Gene.names",
                          secondaryIdCol = "Majority.protein.IDs") {
    .assertVector(x = sce, type = "SummarizedExperiment")
    vvs <- colnames(SummarizedExperiment::rowData(sce))
    .assertScalar(x = primaryIdCol, type = "character", validValues = vvs)
    .assertScalar(x = secondaryIdCol, type = "character", validValues = vvs)

    ## Extract the first annotated gene name
    gName <- vapply(strsplit(
        SummarizedExperiment::rowData(sce)[[primaryIdCol]], ";"),
        .subset, 1, FUN.VALUE = "NA")
    rowData(sce)$primaryIdSingle <- gName

    ## Extract the first annotated majority protein ID
    majProtID <- vapply(strsplit(
        SummarizedExperiment::rowData(sce)[[secondaryIdCol]], ";"),
        .subset, 1, FUN.VALUE = "NA")
    rowData(sce)$secondaryIdSingle <- majProtID

    ## Generate IDs for STRING
    stringIDs <- gName
    idxna <- which(is.na(stringIDs))
    stringIDs[idxna] <- majProtID[idxna]
    SummarizedExperiment::rowData(sce)$IDsForSTRING <- stringIDs

    ## If there are duplicated gene IDs, make them unique by appending the
    ## respective majority protein ID
    idxdup <- which(duplicated(gName) | is.na(gName))
    idxdup <- which(gName %in% gName[idxdup])
    gName[idxdup] <- paste0(gName[idxdup], ".", majProtID[idxdup])
    gName[idxdup] <- sub("^NA\\.", "", gName[idxdup])

    ## Check that there are no duplicated IDs and set as row names
    gName <- make.unique(gName, sep = ".")
    stopifnot(all(!duplicated(gName)))
    rownames(sce) <- gName

    sce
}
