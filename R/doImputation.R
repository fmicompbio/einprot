#' @importFrom stats quantile median rnorm
#' @keywords internal
#' @noRd
#'
.minProbGlobalFun <- function(mat, lowQuantile = 0.01, multSigma = 1) {
    .assertVector(x = mat, type = "matrix")
    .assertScalar(x = lowQuantile, type = "numeric", rngExcl = c(0, 1))
    .assertScalar(x = multSigma, type = "numeric")

    nMissing <- sum(is.na(c(mat)))
    if (nMissing == 0) {
        return(mat)
    } else {
        ## Get a global low quantile
        globalMin <- stats::quantile(c(mat), prob = lowQuantile, na.rm = TRUE)

        ## Only use features observed in at least half the samples to calculate
        ## the standard deviation
        keepRows <- which(apply(!is.na(mat), 1, sum) / ncol(mat) > 0.5)
        globalSd <- stats::median(apply(mat[keepRows, ], 1, sd, na.rm = TRUE),
                                  na.rm = TRUE) * multSigma
        print(globalSd)

        ## Generate random values and replace NAs
        randValues <- stats::rnorm(nMissing, mean = globalMin, sd = globalSd)
        mat[is.na(mat)] <- randValues
        return(mat)
    }
}

#' @importFrom stats quantile median rnorm
#' @keywords internal
#' @noRd
#'
.minProbGlobalOffsetFun <- function(mat, lowQuantile = 0.01, multSigma = 1,
                                    constOffset = NULL, relOffset = NULL) {
    .assertVector(x = mat, type = "matrix")
    .assertScalar(x = lowQuantile, type = "numeric", rngExcl = c(0, 1))
    .assertScalar(x = multSigma, type = "numeric")
    .assertScalar(x = constOffset, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = relOffset, type = "numeric", allowNULL = TRUE)
    if (!is.null(constOffset) && !is.null(relOffset)) {
        stop("At most one of constOffset and relOffset can be specified")
    }

    nMissing <- sum(is.na(c(mat)))
    if (nMissing == 0) {
        return(mat)
    } else {
        ## Get a global low quantile
        globalMin <- stats::quantile(c(mat), prob = lowQuantile, na.rm = TRUE)

        ## Only use features observed in at least half the samples to calculate
        ## the standard deviation
        keepRows <- which(apply(!is.na(mat), 1, sum) / ncol(mat) > 0.5)
        globalSd <- stats::median(apply(mat[keepRows, ], 1, sd, na.rm = TRUE),
                                  na.rm = TRUE) * multSigma
        print(globalSd)

        if (!is.null(constOffset)) {
            globalMin <- globalMin + constOffset
        } else if (!is.null(relOffset)) {
            globalMin <- globalMin + relOffset * globalSd
        }

        ## Generate random values and replace NAs
        randValues <- stats::rnorm(nMissing, mean = globalMin, sd = globalSd)
        mat[is.na(mat)] <- randValues
        return(mat)
    }
}

#' @importFrom stats quantile median rnorm
#' @keywords internal
#' @noRd
#'
.minProbOffsetFun <- function(mat, lowQuantile = 0.01, multSigma = 1,
                              margin = 2L, constOffset = NULL, relOffset = NULL) {
    .assertVector(x = mat, type = "matrix")
    .assertScalar(x = lowQuantile, type = "numeric", rngExcl = c(0, 1))
    .assertScalar(x = multSigma, type = "numeric")
    .assertScalar(x = margin, type = "numeric", validValues = c(1, 2))
    .assertScalar(x = constOffset, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = relOffset, type = "numeric", allowNULL = TRUE)
    if (!is.null(constOffset) && !is.null(relOffset)) {
        stop("At most one of constOffset and relOffset can be specified")
    }

    nMissing <- sum(is.na(c(mat)))
    if (nMissing == 0) {
        return(mat)
    } else {
        if (margin == 1) {
            mat <- t(mat)
        }

        ## Get a low quantile per sample
        minPerSample <- apply(mat, MARGIN = 2, FUN = stats::quantile,
                              prob = lowQuantile, na.rm = TRUE)

        ## Only use features observed in at least half the samples to calculate
        ## the standard deviation
        keepRows <- which(apply(!is.na(mat), 1, sum) / ncol(mat) > 0.5)
        globalSd <- stats::median(apply(mat[keepRows, ], 1, sd, na.rm = TRUE),
                                  na.rm = TRUE) * multSigma
        print(globalSd)

        if (!is.null(constOffset)) {
            minPerSample <- minPerSample + constOffset
        } else if (!is.null(relOffset)) {
            minPerSample <- minPerSample + relOffset * globalSd
        }

        ## Generate random values and replace NAs
        for (i in seq_len(ncol(mat))) {
            randValues <- stats::rnorm(nrow(mat), mean = minPerSample[i],
                                       sd = globalSd)
            naidx <- which(is.na(mat[, i]))
            if (length(naidx) > 0) {
                mat[naidx, i] <- randValues[naidx]
            }
        }

        if (margin == 1) {
            mat <- t(mat)
        }

        return(mat)
    }
}

#' Perform imputation of NA values
#'
#' Perform imputation of missing values (represented by \code{NA}) in one assay
#' in a \code{SummarizedExperiment}, and generate a new assay containing the
#' complete data (including imputed values).
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param imputeMethod Character scalar giving the imputation method. Currently,
#'     \code{"MinProb"} (provided in the \code{MsCoreUtils} package),
#'     \code{"impSeqRob"} (provided in the \code{rrcovNA} package),
#'     \code{"MinProbGlobal"} (a reimplementation of the MinProb algorithm
#'     using a global mean value rather than sample-specific ones),
#'     \code{"MinProbGlobalOffset"} (like \code{"MinProbGlobal"}, but a
#'     constant or relative offset can be specified to shift the mean of the
#'     imputed values), and \code{"MinProbOffset"} (like
#'     \code{"MinProbGlobalOffset"} but imputing separately for each row or
#'     column) are supported. In addition, a custom function can be used by
#'     setting \code{imputeMethod = "custom"} and providing the function in the
#'     \code{FUN} argument (this will be passed on to
#'     \code{MsCoreUtils::impute_matrix}).
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be imputed. The matrix should have missing values represented as
#'     \code{NA}.
#' @param imputedAssayName Character scalar providing the name that will be
#'     given to the assay containing the imputed values.
#' @param ... Additional arguments that will be passed on to the imputation
#'     function.
#'
#' @details For the pre-defined imputation functions defined by einprot,
#' these are the available arguments:
#' \describe{
#'     \item{"MinProbGlobal":}{lowQuantile = 0.01, multSigma = 1 (similar to
#'         \code{q} and \code{sigma} in
#'         \code{MsCoreUtils::MsCoreUtils::impute_MinProb}).}
#'     \item{"MinProbOffset":}{lowQuantile = 0.01, multSigma = 1,
#'         margin = 2L, constOffset = NULL, relOffset = NULL. \code{margin}
#'         indicates whether to impute along rows or columns.
#'         \code{constOffset} defines a constant shift of the mean of the
#'         imputed values relative to the low quantile. \code{relOffset}
#'         instead defines the shift as a multiple of the standard deviation
#'         of the distribution of imputed values.}
#'     \item{"MinProbGlobalOffset":}{ lowQuantile = 0.01, multSigma = 1,
#'         constOffset = NULL, relOffset = NULL. Similar interpretation as
#'         for \code{"MinProbOffset"}.}
#' }
#' If \code{imputeMethod = "custom"}, any suitable function that accepts a
#' matrix as input and returns a matrix of the same shape can be provided to
#' the \code{FUN} argument. In addition, any argument to the function can
#' also be specified. The function and all arguments will be passed on to
#' \code{\link[MsCoreUtils]{impute_matrix}}.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns An object of the same type as \code{sce} with an additional assay
#'     named \code{imputedAssayName}.
#'
#' @examples
#' ## Import data
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#'
#' ## Log-transform iBAQ values
#' SummarizedExperiment::assay(sce, "log2_iBAQ") <-
#'     log2(SummarizedExperiment::assay(sce, "iBAQ"))
#'
#' ## Replace non-finite values by NA
#' SummarizedExperiment::assay(sce, "log2_iBAQ")[!is.finite(
#'     SummarizedExperiment::assay(sce, "log2_iBAQ"))] <- NA
#'
#' ## Impute missing values
#' sce <- doImputation(sce, imputeMethod = "MinProb", assayName = "log2_iBAQ",
#'                     imputedAssayName = "imputed_iBAQ")
#' SummarizedExperiment::assayNames(sce)
#'
#' ## Use another method from MsCoreUtils (specified as a custom function)
#' sce2 <- doImputation(sce, imputeMethod = "custom", assayName = "log2_iBAQ",
#'                      imputedAssayName = "imputed_iBAQ_knn",
#'                      FUN = MsCoreUtils::impute_knn)
#' ## Alternatively, just specify one of the methods supported by MsCoreUtils::impute_matrix
#' sce3 <- doImputation(sce, imputeMethod = "custom", assayName = "log2_iBAQ",
#'                      imputedAssayName = "imputed_iBAQ_knn",
#'                      method = "knn")
#' identical(SummarizedExperiment::assay(sce2, "imputed_iBAQ_knn"),
#'           SummarizedExperiment::assay(sce3, "imputed_iBAQ_knn"))
#'
#' @importFrom MsCoreUtils impute_matrix
#' @importFrom rrcovNA impSeqRob
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom imputeLCMD impute.MinProb
#'
doImputation <- function(sce, imputeMethod, assayName, imputedAssayName, ...) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = imputeMethod, type = "character",
                  validValues = c("MinProb", "impSeqRob", "MinProbGlobal",
                                  "MinProbOffset", "MinProbGlobalOffset",
                                  "custom"))
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = imputedAssayName, type = "character")

    if (imputeMethod == "MinProb") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                method = "MinProb", ...
            )
    } else if (imputeMethod == "impSeqRob") {
        tmp <- rrcovNA::impSeqRob(
            x = SummarizedExperiment::assay(sce, assayName), ...
        )
        SummarizedExperiment::assay(sce, imputedAssayName) <- tmp$x
    } else if (imputeMethod == "MinProbGlobal") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                FUN = .minProbGlobalFun, ...
            )
    } else if (imputeMethod == "MinProbGlobalOffset") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                FUN = .minProbGlobalOffsetFun, ...
            )
    } else if (imputeMethod == "MinProbOffset") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                FUN = .minProbOffsetFun, ...
            )
    } else if (imputeMethod == "custom") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                ...
            )
    }
    sce
}
