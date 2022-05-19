#' Perform imputation of NA values
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param method Character scalar giving the imputation method.
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be imputed.
#' @param imputedAssayName Character scalar providing the name that will be
#'     given to the assay containing imputed values.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{SummarizedExperiment} object with an additional assay
#'     named \code{imputedAssayName}.
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#' SummarizedExperiment::assay(sce, "log2_iBAQ") <-
#'     log2(SummarizedExperiment::assay(sce, "iBAQ"))
#' SummarizedExperiment::assay(sce, "log2_iBAQ")[!is.finite(
#'     SummarizedExperiment::assay(sce, "log2_iBAQ"))] <- NA
#' sce <- doImputation(sce, method = "MinProb", assayName = "log2_iBAQ",
#'                     imputedAssayName = "imputed_iBAQ")
#'
#' @importFrom MsCoreUtils impute_matrix
#' @importFrom rrcovNA impSeqRob
#' @importFrom SummarizedExperiment assay assay<-
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
