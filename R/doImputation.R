#' Perform imputation of NA values
#'
#' Perform imputation of missing values (represented by \code{NA}) in one assay
#' in a \code{SummarizedExperiment}, and generate a new assay containing the
#' complete data (including imputed values).
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param method Character scalar giving the imputation method. Currently,
#'     \code{"MinProb"} (provided in the \code{MsCoreUtils} package) and
#'     \code{"impSeqRob"} (provided in the \code{rrcovNA} package) are
#'     supported.
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be imputed. The matrix should have missing values represented as
#'     \code{NA}.
#' @param imputedAssayName Character scalar providing the name that will be
#'     given to the assay containing the imputed values.
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
doImputation <- function(sce, method, assayName, imputedAssayName) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = method, type = "character",
                  validValues = c("MinProb", "impSeqRob"))
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = imputedAssayName, type = "character")

    if (method == "MinProb") {
        SummarizedExperiment::assay(sce, imputedAssayName) <-
            MsCoreUtils::impute_matrix(
                SummarizedExperiment::assay(sce, assayName), method = "MinProb"
            )
    } else if (method == "impSeqRob") {
        tmp <- rrcovNA::impSeqRob(SummarizedExperiment::assay(sce, assayName))
        SummarizedExperiment::assay(sce, imputedAssayName) <- tmp$x
    }
    sce
}
