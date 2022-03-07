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
