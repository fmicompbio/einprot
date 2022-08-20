#' Define assay names to use in workflow
#'
#' @param aName Base assay name, typically obtained from
#'     \code{importExperiment}.
#' @param normMethod Character scalar, indicating the normalization method.
#' @param doBatchCorr Logical scalar, whether or not batch correction will
#'     be performed.
#'
#' @return A list with assay names that will be used for assays created at
#' different steps in the einprot workflows.
#'
#' @author Charlotte Soneson
#' @export
defineAssayNames <- function(aName, normMethod, doBatchCorr) {
    .assertScalar(x = aName, type = "character")
    .assertScalar(x = normMethod, type = "character")
    .assertScalar(x = doBatchCorr, type = "logical")

    if (normMethod == "none") {
        nsuf <- ""
    } else {
        nsuf <- "_norm"
    }

    if (doBatchCorr) {
        bsuf <- "_batchCorr"
    } else {
        bsuf <- ""
    }

    list(
        assayInput = aName,
        assayLog2WithNA = paste0("log2_", aName, "_withNA"),
        assayImputIndic = paste0("imputed_", aName),
        assayLog2NormWithNA = paste0("log2_", aName, "_withNA", nsuf),
        assayImputed = paste0("log2_", aName, nsuf),
        assayBatchCorr = paste0("log2_", aName, nsuf, bsuf)
    )
}
