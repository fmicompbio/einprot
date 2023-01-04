#' Extract feature identifiers
#'
#' Extract feature names by splitting a given column by a separator and
#' keeping the first entry.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @param df A \code{data.frame}.
#' @param colName A character scalar indicating which of the columns in
#'     \code{df} to consider.
#' @param separator A character scalar giving the separator to split the
#'     column entries by.
#'
#' @return A vector with extracted feature identifiers
#'
getFirstId <- function(df, colName, separator = ";") {
    df <- as.data.frame(df)
    vvs <- colnames(df)
    .assertScalar(x = colName, type = "character", validValues = vvs)
    .assertScalar(x = separator, type = "character")

    vapply(strsplit(as.character(df[[colName]]), separator),
           .subset, 1, FUN.VALUE = "NA")
}

#' Combine multiple columns into an ID
#'
#' @export
#' @author Charlotte Soneson
#'
#' @param df A \code{data.frame}.
#' @param combineCols A character vector giving the columns of \code{df} to
#'     combine.
#' @param combineWhen A character scalar indicating when to combine columns.
#'     Must be either 'always' (which always combines the columns),
#'     'nonunique' (which only combines the columns if it's necessary to
#'     obtain unique names), or 'missing' (which uses subsequent columns if
#'     all previous columns have missing values in a given position).
#' @param splitSeparator A character scalar, a character vector of length
#'     equal to the length of \code{combineCols}, or \code{NULL}. If not
#'     \code{NULL}, indicates the separator by which to split the entries in
#'     the corresponding column before combining columns.
#' @param joinSeparator A character scalar giving the separator to use when
#'     combining columns.
#' @param makeUnique Logical scalar, indicating whether or not the feature IDs
#'     should be guaranteed to be unique.
#'
#' @return A vector with combined feature identifiers
#'
#' @importFrom dplyr bind_cols
#' @importFrom rlang .data
#'
combineIds <- function(df, combineCols, combineWhen = "nonunique",
                       splitSeparator = ";", joinSeparator = ".",
                       makeUnique = TRUE) {
    df <- as.data.frame(df)
    vvs <- colnames(df)
    .assertVector(x = combineCols, type = "character", validValues = vvs,
                  allowNULL = TRUE)
    .assertScalar(x = combineWhen, type = "character",
                  validValues = c("always", "nonunique", "missing"))
    .assertVector(x = splitSeparator, type = "character", allowNULL = TRUE)
    .assertScalar(x = joinSeparator, type = "character")
    if (!is.null(splitSeparator) && length(splitSeparator) == 1) {
        splitSeparator <- rep(splitSeparator, length(combineCols))
    }
    if (!is.null(splitSeparator) && length(splitSeparator) > 1) {
        stopifnot(length(splitSeparator) == length(combineCols))
    }
    .assertScalar(x = makeUnique, type = "logical")

    if (is.null(combineCols)) {
        return(rep(NA_character_, nrow(df)))
    }

    colExtr <- do.call(
        dplyr::bind_cols,
        lapply(structure(seq_along(combineCols),
                         names = combineCols), function(i) {
            if (!is.null(splitSeparator)) {
                ## Empty cells become NA here
                getFirstId(df, colName = combineCols[i],
                           separator = splitSeparator[i])
            } else {
                df[[combineCols[i]]]
            }
        }))

    if (combineWhen == "always") {
        finalIds <- colExtr %>% tidyr::unite(col = "finalId", everything(),
                                             sep = joinSeparator) %>%
            dplyr::pull(.data$finalId)
    } else if (combineWhen == "nonunique") {
        finalIds <- colExtr[[1]]
        j <- 1
        while (any(duplicated(finalIds) | is.na(finalIds) | finalIds == "") &&
               j < ncol(colExtr)) {
            idxdup <- which(duplicated(finalIds) | is.na(finalIds) | finalIds == "")
            idxdup <- which(finalIds %in% finalIds[idxdup])
            finalIds[idxdup] <- paste0(finalIds[idxdup],
                                       joinSeparator, colExtr[[j + 1]][idxdup])
            finalIds[idxdup] <- sub(paste0("^NA", ifelse(joinSeparator == ".",
                                                         "\\.", joinSeparator)),
                                    "", finalIds[idxdup])
            finalIds[idxdup] <- sub(paste0("^", ifelse(joinSeparator == ".",
                                                       "\\.", joinSeparator)),
                                    "", finalIds[idxdup])
            j <- j + 1
        }
    } else if (combineWhen == "missing") {
        finalIds <- colExtr[[1]]
        j <- 1
        while (any(is.na(finalIds) | finalIds == "") && j < ncol(colExtr)) {
            idxdup <- which(is.na(finalIds) | finalIds == "")
            finalIds[idxdup] <- colExtr[[j + 1]][idxdup]
            j <- j + 1
        }
    } else {
        stop("Unknown value of 'combineWhen'")
    }

    if (makeUnique) {
        finalIds <- make.unique(finalIds, sep = joinSeparator)
    }
    finalIds
}


#' Define various ID columns and set row names
#'
#' Define various types of feature IDs based on the information in the
#' \code{rowData} of \code{sce}.
#'
#' @param sce A \code{SummarizedExperiment} object (or derivative).
#' @param idCol,labelCol,geneIdCol,proteinIdCol,stringIdCol Arguments defining
#'     the feature identifiers (row names, should be unique),
#'     feature labels (for plots, can be anything),
#'     gene IDs (single gene symbols, will be matched against complexes and
#'     GO terms, can be \code{NULL}),
#'     protein IDs (UniProt IDs, will be used to create automatic URLs and
#'     match to species-specific identifiers, each entry can consist of
#'     multiple UniProt IDs separated by semicolons), and
#'     stringIdCol (single gene or protein ID, will be used to retrieve
#'     STRING networks, can be \code{NULL}).
#'     Each of these arguments can be either a character vector of column
#'     names in \code{rowData(sce)}, in which case the corresponding feature ID
#'     is generated by simply concatenating
#'     the values in these columns, or a function with one input argument
#'     (a data.frame, corresponding to \code{rowData(sce)}), returning a
#'     character vector corresponding to the desired feature IDs.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return An object of the same type as \code{sce} with modified row names
#' and additional columns \code{"einprotId"}, \code{"einprotLabel"},
#' \code{"einprotGene"} and \code{"einprotProtein"} in \code{rowData(sce)}.
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#' sce <- fixFeatureIds(
#'     sce,
#'     idCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs")),
#'     labelCol = c("Gene.names", "Majority.protein.IDs"),
#'     geneIdCol = function(df) getFirstId(df, "Gene.names"),
#'     proteinIdCol = "Majority.protein.IDs",
#'     stringIdCol = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
#'                                           combineWhen = "missing", makeUnique = FALSE)
#' )
#' head(rownames(sce))
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
fixFeatureIds <- function(
        sce,
        idCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                        combineWhen = "nonunique",
                                        splitSeparator = ";", joinSeparator = ".",
                                        makeUnique = TRUE),
        labelCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                           combineWhen = "nonunique",
                                           splitSeparator = ";", joinSeparator = ".",
                                           makeUnique = FALSE),
        geneIdCol = function(df) getFirstId(df, colName = "Gene.names",
                                            separator = ";"),
        proteinIdCol = "Majority.protein.IDs",
        stringIdCol = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                              combineWhen = "missing",
                                              splitSeparator = ";", joinSeparator = ".",
                                              makeUnique = FALSE)) {

    .assertVector(x = sce, type = "SummarizedExperiment")
    vvs <- colnames(SummarizedExperiment::rowData(sce))

    if (is.function(idCol)) {
        ## Should be a function with a single argument
        stopifnot(length(formals(idCol)) == 1)
        idColFun <- idCol
    } else {
        .assertVector(x = idCol, type = "character", validValues = vvs,
                      allowNULL = FALSE)
        ## If it's a vector of column names, create a function that
        ## pastes them together
        idColFun <- function(df) combineIds(df, combineCols = idCol,
                                            combineWhen = "always",
                                            splitSeparator = NULL,
                                            joinSeparator = ".",
                                            makeUnique = TRUE)
    }

    if (is.function(labelCol)) {
        ## Should be a function with a single argument
        stopifnot(length(formals(labelCol)) == 1)
        labelColFun <- labelCol
    } else {
        .assertVector(x = labelCol, type = "character", validValues = vvs,
                      allowNULL = FALSE)
        ## If it's a vector of column names, create a function that
        ## pastes them together
        labelColFun <- function(df) combineIds(df, combineCols = labelCol,
                                               combineWhen = "always",
                                               splitSeparator = NULL,
                                               joinSeparator = ".",
                                               makeUnique = FALSE)
    }

    if (is.function(geneIdCol)) {
        ## Should be a function with a single argument
        stopifnot(length(formals(geneIdCol)) == 1)
        geneIdColFun <- geneIdCol
    } else {
        .assertVector(x = geneIdCol, type = "character", validValues = vvs,
                      allowNULL = TRUE)
        ## If it's a vector of column names, create a function that
        ## pastes them together
        geneIdColFun <- function(df) combineIds(df, combineCols = geneIdCol,
                                                combineWhen = "always",
                                                splitSeparator = NULL,
                                                joinSeparator = ".",
                                                makeUnique = FALSE)
    }

    if (is.function(proteinIdCol)) {
        ## Should be a function with a single argument
        stopifnot(length(formals(proteinIdCol)) == 1)
        proteinIdColFun <- proteinIdCol
    } else {
        .assertVector(x = proteinIdCol, type = "character", validValues = vvs,
                      allowNULL = FALSE)
        ## If it's a vector of column names, create a function that
        ## pastes them together
        proteinIdColFun <- function(df) combineIds(df, combineCols = proteinIdCol,
                                                   combineWhen = "always",
                                                   splitSeparator = NULL,
                                                   joinSeparator = ".",
                                                   makeUnique = FALSE)
    }

    if (is.function(stringIdCol)) {
        ## Should be a function with a single argument
        stopifnot(length(formals(stringIdCol)) == 1)
        stringIdColFun <- stringIdCol
    } else {
        .assertVector(x = stringIdCol, type = "character", validValues = vvs,
                      allowNULL = TRUE)
        ## If it's a vector of column names, create a function that
        ## pastes them together
        stringIdColFun <- function(df) combineIds(df, combineCols = stringIdCol,
                                                  combineWhen = "always",
                                                  splitSeparator = NULL,
                                                  joinSeparator = ".",
                                                  makeUnique = FALSE)
    }

    ## Get new columns
    SummarizedExperiment::rowData(sce)[["einprotLabel"]] <-
        labelColFun(SummarizedExperiment::rowData(sce))
    SummarizedExperiment::rowData(sce)[["einprotGene"]] <-
        geneIdColFun(SummarizedExperiment::rowData(sce))
    SummarizedExperiment::rowData(sce)[["einprotProtein"]] <-
        proteinIdColFun(SummarizedExperiment::rowData(sce))
    SummarizedExperiment::rowData(sce)[["IDsForSTRING"]] <-
        stringIdColFun(SummarizedExperiment::rowData(sce))
    SummarizedExperiment::rowData(sce)[["einprotId"]] <-
        idColFun(SummarizedExperiment::rowData(sce))
    rownames(sce) <- SummarizedExperiment::rowData(sce)[["einprotId"]]

    sce
}
