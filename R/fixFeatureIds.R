#' @export
#' @author Charlotte Soneson
#'
#' @importFrom SummarizedExperiment rowData
#'
fixFeatureIds <- function(qft, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs") {
    ## Extract the first annotated gene name
    gName <- vapply(strsplit(rowData(qft[[1]])[[geneIdCol]], ";"),
                    .subset, 1, FUN.VALUE = "NA")
    rowData(qft[[1]])$geneIdSingle <- gName

    ## Extract the first annotated majority protein ID
    majProtID <- vapply(strsplit(rowData(qft[[1]])[[proteinIdCol]], ";"),
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
