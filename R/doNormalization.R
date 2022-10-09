#' Apply normalization to an assay in a SummarizedExperiment object
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param method Character scalar giving the normalization method. Currently,
#'     the methods from \code{MsCoreUtils::normalizeMethods()} are supported.
#'     If \code{spikeFeatures} is not \code{NULL}, only
#'     \code{center.mean}, \code{center.median}, \code{div.mean} and
#'     \code{div.median} are supported.
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be normalized
#' @param normalizedAssayName Character scalar providing the name that will be
#'     given to the assay containing normalized values.
#' @param spikeFeatures Character vector of feature IDs (rownames of sce)
#'     that will be used to calculate normalization factors. If \code{NULL}
#'     (default), all features are used.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{SummarizedExperiment} object with an additional assay
#'     named \code{normalizedAssayName}.
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
#' sce <- doNormalization(sce, method = "center.median", assayName = "log2_iBAQ",
#'                        normalizedAssayName = "normalized_iBAQ")
#'
#' @importFrom MsCoreUtils normalizeMethods impute_matrix
#' @importFrom SummarizedExperiment assay assay<-
#'
doNormalization <- function(sce, method, assayName, normalizedAssayName,
                            spikeFeatures = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = method, type = "character",
                  validValues = MsCoreUtils::normalizeMethods())
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = normalizedAssayName, type = "character")
    .assertVector(x = spikeFeatures, type = "character", allowNULL = TRUE,
                  validValues = rownames(sce))

    assayIn <- SummarizedExperiment::assay(sce, assayName)
    if (is.null(spikeFeatures)) {
        ## No spike features - apply one of the methods in
        ## MsCoreUtils::normalizeMethods()
        if (method %in% MsCoreUtils::normalizeMethods()) {
            assayOut <-
                MsCoreUtils::normalize_matrix(assayIn,
                                              method = method)
        } else {
            stop("Unrecognized normalization method: ", method)
        }
    } else {
        ## Spike features - only a subset of the methods available
        if (method == "center.mean") {
            cvec <- colMeans(assayIn[spikeFeatures, , drop = FALSE],
                             na.rm = TRUE)
            cvec <- cvec - mean(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "-", check.margin = FALSE)
        } else if (method == "center.median") {
            cvec <- apply(assayIn[spikeFeatures, , drop = FALSE], 2L, stats::median,
                          na.rm = TRUE)
            cvec <- cvec - median(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "-", check.margin = FALSE)
        } else if (method == "div.mean") {
            cvec <- colMeans(assayIn[spikeFeatures, , drop = FALSE],
                             na.rm = TRUE)
            cvec <- cvec/mean(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "/", check.margin = FALSE)
        } else if (method == "div.median") {
            cvec <- apply(assayIn[spikeFeatures, , drop = FALSE], 2L, stats::median,
                          na.rm = TRUE)
            cvec <- cvec/median(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "/", check.margin = FALSE)
        } else {
            stop("Unsupported normalization method with spike features: ", method)
        }
        rownames(assayOut) <- rownames(assayIn)
        colnames(assayOut) <- colnames(assayIn)
    }

    SummarizedExperiment::assay(sce, normalizedAssayName) <- assayOut
    sce
}
