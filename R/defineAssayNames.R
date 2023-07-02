#' Define assay names
#'
#' Starting from a base assay name and a few decisions about the workflow,
#' define the names of assays that will be generated in the \code{einprot}
#' workflow and included in the final \code{SingleCellExperiment} object.
#'
#' @param aName Base assay name, typically obtained from
#'     \code{importExperiment}.
#' @param normMethod Character scalar, indicating the normalization method.
#'     See \code{doNormalization} for available options. Set to \code{"none"} if
#'     no between-sample normalization will be performed.
#' @param doBatchCorr Logical scalar indicating whether or not batch correction
#'     will be performed.
#'
#' @returns A list with assay names that will be used for assays created at
#' different steps in the \code{einprot} workflows.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")
#' defineAssayNames(sce$aName, normMethod = "none", doBatchCorr = FALSE)
#'
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
