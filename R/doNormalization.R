#' Apply normalization
#'
#' Apply normalization to an assay in a \code{SummarizedExperiment} object and
#' add a new assay containing the normalized values.
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param method Character scalar giving the normalization method. Currently,
#'     the methods from \code{MsCoreUtils::normalizeMethods()} are supported.
#'     If \code{spikeFeatures} is not \code{NULL}, only
#'     \code{"center.mean"}, \code{"center.median"}, \code{"div.mean"} and
#'     \code{"div.median"} are supported.
#' @param assayName Character scalar giving the name of the assay in \code{sce}
#'     to be normalized.
#' @param normalizedAssayName Character scalar providing the name that will be
#'     given to the assay containing normalized values.
#' @param spikeFeatures Character vector of feature IDs (rownames of sce)
#'     that will be used to calculate normalization factors. If \code{NULL}
#'     (default), all features are used.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return An object of the same type as \code{sce} with an additional assay
#'     named \code{normalizedAssayName}.
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
#' ## Normalize between samples using median centering
#' sce <- doNormalization(sce, method = "center.median",
#'                        assayName = "log2_iBAQ",
#'                        normalizedAssayName = "normalized_iBAQ")
#' SummarizedExperiment::assayNames(sce)
#'
#' ## Check that the median is zero for all samples in the normalized data
#' apply(SummarizedExperiment::assay(sce, "normalized_iBAQ"), 2, median,
#'       na.rm = TRUE)
#'
#' @importFrom MsCoreUtils normalizeMethods impute_matrix
#' @importFrom SummarizedExperiment assay assay<- assayNames
#'
doNormalization <- function(sce, method, assayName, normalizedAssayName,
                            spikeFeatures = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = method, type = "character",
                  validValues = MsCoreUtils::normalizeMethods())
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = normalizedAssayName, type = "character")
    .assertVector(x = spikeFeatures, type = "character", allowNULL = TRUE)

    ## Find the spike features that are also present in the data
    if (!is.null(spikeFeatures)) {
        nsp <- length(spikeFeatures)
        spikeFeatures <- intersect(spikeFeatures, rownames(sce))
        if (length(spikeFeatures) > 0) {
            ## Only keep spike features that have non-missing values in
            ## all samples - otherwise normalization becomes difficult to
            ## interpret
            tmpmat <-
                SummarizedExperiment::assay(sce, assayName)[spikeFeatures, ,
                                                            drop = FALSE]
            tmpmat <- tmpmat[rowSums(is.na(tmpmat)) == 0, , drop = FALSE]
            spikeFeatures <- rownames(tmpmat)
        }
        message(length(spikeFeatures), "/", nsp, " spike feature",
                ifelse(length(spikeFeatures) == 1, "", "s"), " found ",
                "with no missing values in the data set.")
        if (length(spikeFeatures) == 0) {
            spikeFeatures <- NULL
            warning("No valid spike features found - falling back to ",
                    "normalizing based on all features.")
        }
    }

    assayIn <- SummarizedExperiment::assay(sce, assayName)
    if (is.null(spikeFeatures)) {
        ## No spike features - apply one of the methods in
        ## MsCoreUtils::normalizeMethods()
        if (method %in% MsCoreUtils::normalizeMethods()) {
            assayOut <-
                MsCoreUtils::normalize_matrix(assayIn,
                                              method = method)
        } else {
            ## Should never end up here as we check the validity of method above
            #nocov start
            stop("Unrecognized normalization method: ", method)
            #nocov end
        }
    } else {
        ## Spike features - only a subset of the methods available
        if (method == "center.mean") {
            cvec <- colMeans(assayIn[spikeFeatures, , drop = FALSE],
                             na.rm = TRUE)
            cvec <- cvec - mean(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "-",
                              check.margin = FALSE)
        } else if (method == "center.median") {
            cvec <- apply(assayIn[spikeFeatures, , drop = FALSE], 2L,
                          stats::median, na.rm = TRUE)
            cvec <- cvec - median(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "-",
                              check.margin = FALSE)
        } else if (method == "div.mean") {
            cvec <- colMeans(assayIn[spikeFeatures, , drop = FALSE],
                             na.rm = TRUE)
            cvec <- cvec/mean(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "/",
                              check.margin = FALSE)
        } else if (method == "div.median") {
            cvec <- apply(assayIn[spikeFeatures, , drop = FALSE], 2L,
                          stats::median, na.rm = TRUE)
            cvec <- cvec/median(cvec, na.rm = TRUE)
            assayOut <- sweep(assayIn, 2L, cvec, FUN = "/",
                              check.margin = FALSE)
        } else {
            stop("Unsupported normalization method with spike features: ",
                 method)
        }
        rownames(assayOut) <- rownames(assayIn)
        colnames(assayOut) <- colnames(assayIn)
    }

    SummarizedExperiment::assay(sce, normalizedAssayName) <- assayOut
    sce
}
