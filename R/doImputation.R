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

#' Perform imputation of NA values
#'
#' Perform imputation of missing values (represented by \code{NA}) in one assay
#' in a \code{SummarizedExperiment}, and generate a new assay containing the
#' complete data (including imputed values).
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param method Character scalar giving the imputation method. Currently,
#'     \code{"MinProb"} (provided in the \code{MsCoreUtils} package),
#'     \code{"impSeqRob"} (provided in the \code{rrcovNA} package), and
#'     \code{"MinProbGlobal"} (a reimplementation of the MinProb algorithm
#'     using a global mean value rather than sample-specific ones) are
#'     supported.
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be imputed. The matrix should have missing values represented as
#'     \code{NA}.
#' @param imputedAssayName Character scalar providing the name that will be
#'     given to the assay containing the imputed values.
#' @param ... Additional arguments that will be passed on to the imputation
#'     function.
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
#' sce <- doImputation(sce, method = "MinProb", assayName = "log2_iBAQ",
#'                     imputedAssayName = "imputed_iBAQ")
#' SummarizedExperiment::assayNames(sce)
#'
#' @importFrom MsCoreUtils impute_matrix
#' @importFrom rrcovNA impSeqRob
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom imputeLCMD impute.MinProb
#'
doImputation <- function(sce, method, assayName, imputedAssayName, ...) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = method, type = "character",
                  validValues = c("MinProb", "impSeqRob", "MinProbGlobal"))
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = imputedAssayName, type = "character")

    if (method == "MinProb") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                method = "MinProb", ...
            )
    } else if (method == "impSeqRob") {
        tmp <- rrcovNA::impSeqRob(SummarizedExperiment::assay(sce, assayName))
        SummarizedExperiment::assay(sce, imputedAssayName) <- tmp$x
    } else if (method == "MinProbGlobal") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName),
                FUN = .minProbGlobalFun, ...
            )
    }
    sce
}
