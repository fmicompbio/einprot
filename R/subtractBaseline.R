#' Extract a matrix and subtract the background/baseline for each feature
#'
#' Extract an assay from a SummarizedExperiment, and for each feature subtract
#' the average value across a set of reference/background/baseline samples,
#' and add the average across all the baseline samples to retain information
#' about overall abundance. Typically used to adjust for varying baselines
#' between batches, when a reference sample is included in each batch.
#'
#' @param sce A \code{SummarizedExperiment} object. The colData of the object
#'     must have at least two columns named 'group' and 'batch'.
#' @param assayName The name of the assay to extract from \code{sce}.
#' @param baselineGroup The value of the 'group' column of
#'     \code{colData(sceFull)} that corresponds to the baseline/reference
#'     samples.
#' @param sceFull A \code{SummarizedExperiment} object containing at least
#'     the reference samples for each batch in \code{sce}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @importFrom SummarizedExperiment assay assayNames
getMatSubtractedBaseline <- function(sce, assayName, baselineGroup, sceFull) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertVector(x = sceFull, type = "SummarizedExperiment")
    stopifnot(all(!is.null(sce$group), !is.null(sce$batch),
                  !is.null(sceFull$group), !is.null(sceFull$batch)))
    .assertScalar(x = baselineGroup, type = "character",
                  validValues = unique(sceFull$group))
    stopifnot(all(sce$batch %in% sceFull$batch[sceFull$group == baselineGroup]))
    .assertScalar(x = assayName, type = "character",
                  validValues = intersect(SummarizedExperiment::assayNames(sce),
                                          SummarizedExperiment::assayNames(sceFull)))

    ## Extract matrix and batch info from sce
    mat <- SummarizedExperiment::assay(sce, assayName,
                                       withDimnames = TRUE)
    batch <- sce$batch

    ## Extract matrix and batch/group info from sceFull
    matFull <- SummarizedExperiment::assay(sceFull, assayName,
                                           withDimnames = TRUE)
    batchFull <- sceFull$batch
    groupFull <- sceFull$group

    ## Get median abundance across reference columns. This will be added
    ## to all columns after the adjustment, to retain overall abundance information
    meanRef <- rowMeans(matFull[, which(groupFull == baselineGroup &
                                            batchFull %in% batch), drop = FALSE],
                        na.rm = TRUE)

    ## Adjust each column of mat
    for (i in seq_len(ncol(mat))) {
        idx <- which(groupFull == baselineGroup & batchFull == batch[i])
        mat[, i] <- mat[, i] -
            rowMeans(matFull[, idx, drop = FALSE], na.rm = TRUE) +
            meanRef
    }
    mat
}
