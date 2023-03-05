#' @export
#' @rdname getNthId
getFirstId <- function(df, colName, separator = ";") {
    getNthId(df = df, colName = colName, N = 1,
             separator = separator)
}

#' Extract feature identifiers
#'
#' Extract feature names by splitting a given column by a separator and
#' keeping the Nth entry.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @param df A \code{data.frame}.
#' @param colName A character scalar indicating which of the columns in
#'     \code{df} to consider.
#' @param N Numeric scalar indicating which part of the column elements to
#'     extract, after separating by \code{separator}.
#' @param separator A character scalar giving the separator to split the
#'     column entries by.
#'
#' @name getNthId
#'
#' @return A vector with extracted feature identifiers
#'
#' @examples
#' df <- data.frame(x = c("g1;p1;h2", "g2;p1;h3"))
#' getFirstId(df, colName = "x", separator = ";")
#' getNthId(df, colName = "x", N = 2, separator = ";")
#'
#' ## Return value will be NA if the field doesn't exist
#' df <- data.frame(x = c("g1;p1;h2", "g2"))
#' getFirstId(df, colName = "x", separator = ";")
#' getNthId(df, colName = "x", N = 2, separator = ";")
#'
getNthId <- function(df, colName, N, separator = ";") {
    df <- as.data.frame(df)
    vvs <- colnames(df)
    .assertScalar(x = colName, type = "character", validValues = vvs)
    .assertScalar(x = N, type = "numeric", rngExcl = c(0, Inf))
    .assertScalar(x = separator, type = "character")

    vapply(strsplit(as.character(df[[colName]]), separator),
           .subset, N, FUN.VALUE = "NA")
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
        #nocov start
        stop("Unknown value of 'combineWhen'")
        #nocov end
    }

    if (makeUnique) {
        finalIds <- make.unique(finalIds, sep = joinSeparator)
    }
    finalIds
}


#' Define various ID columns
#'
#' Define various types of feature IDs based on the information in the
#' \code{rowData} of \code{sce}.
#'
#' @param sce A \code{SummarizedExperiment} object (or derivative).
#' @param colDefs A named list defining how each new column should be defined.
#'     The names will be used as the column names. Each entry can be either
#'     a character vector of column names in \code{rowData(sce)}, in which case
#'     the corresponding feature ID is generated by simply concatenating
#'     the values in these columns, or a function with one input argument
#'     (a data.frame, corresponding to \code{rowData(sce)}), returning a
#'     character vector corresponding to the desired feature IDs.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return An object of the same type as \code{sce} with additional columns
#' in \code{rowData(sce)}.
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#' sce <- fixFeatureIds(
#'     sce,
#'     colDefs = list(
#'         einprotId = function(df) combineIds(df, combineCols = c("Gene.names",
#'                                            "Majority.protein.IDs")),
#'         einprotLabel = c("Gene.names", "Majority.protein.IDs"),
#'         einprotGene = function(df) getFirstId(df, "Gene.names"),
#'         einprotProtein = "Majority.protein.IDs",
#'         IDsForSTRING = function(df) combineIds(df, c("Gene.names",
#'                                                      "Majority.protein.IDs"),
#'                                                combineWhen = "missing",
#'                                                makeUnique = FALSE))
#' )
#' head(SummarizedExperiment::rowData(sce)$einprotId)
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
fixFeatureIds <- function(
        sce,
        colDefs = list(
            einprotId = function(df) combineIds(df, combineCols = c("Gene.names",
                                                                    "Majority.protein.IDs"),
                                                combineWhen = "nonunique",
                                                splitSeparator = ";", joinSeparator = ".",
                                                makeUnique = TRUE),
            einprotLabel = function(df) combineIds(df, combineCols = c("Gene.names",
                                                                       "Majority.protein.IDs"),
                                                   combineWhen = "nonunique",
                                                   splitSeparator = ";", joinSeparator = ".",
                                                   makeUnique = FALSE),
            einprotGene = function(df) getFirstId(df, colName = "Gene.names",
                                                  separator = ";"),
            einprotProtein = "Majority.protein.IDs",
            IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                   combineWhen = "missing",
                                                   splitSeparator = ";", joinSeparator = ".",
                                                   makeUnique = FALSE))
) {

    .assertVector(x = sce, type = "SummarizedExperiment")
    vvs <- colnames(SummarizedExperiment::rowData(sce))

    for (nm in names(colDefs)) {
        if (is.function(colDefs[[nm]])) {
            ## Should be a function with a single argument
            stopifnot(length(formals(colDefs[[nm]])) == 1)
            colFun <- colDefs[[nm]]
        } else {
            .assertVector(x = colDefs[[nm]], type = "character",
                          validValues = vvs, allowNULL = TRUE)
            ## If it's a vector of column names, create a function that
            ## pastes them together
            colFun <- function(df) combineIds(df, combineCols = colDefs[[nm]],
                                              combineWhen = "always",
                                              splitSeparator = NULL,
                                              joinSeparator = ".",
                                              makeUnique = FALSE)
        }
        SummarizedExperiment::rowData(sce)[[nm]] <-
            colFun(SummarizedExperiment::rowData(sce))
    }

    sce
}

