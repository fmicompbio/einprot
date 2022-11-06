#' Add sample annotations to SummarizedExperiment object
#'
#' @param sce A \code{SummarizedExperiment} object.
#' @param sampleAnnot A \code{data.frame} with sample annotations. Must
#'     have at least columns named \code{sample} (which must contain
#'     all the column names of \code{sce}) and \code{group} (which contains
#'     the group assignment for each sample).
#'
#' @return A \code{SummarizedExperiment} object with additional sample
#'     annotations.
#'
#' @examples
#' mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
#'                       package = "einprot")
#' samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
#'              "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
#'              "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
#' out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
#'                         includeOnlySamples = samples)
#' sampleAnnot <- data.frame(sample = samples,
#'                           group = gsub("_IP.*", "", samples))
#' sce <- addSampleAnnots(out$sce, sampleAnnot = sampleAnnot)
#' SummarizedExperiment::colData(sce)  ## group information added to sce
#'
#' @export
#' @author Charlotte Soneson
#'
#' @importFrom SummarizedExperiment colData
addSampleAnnots <- function(sce, sampleAnnot) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertVector(x = sampleAnnot, type = "data.frame")
    stopifnot(all(c("sample", "group") %in% colnames(sampleAnnot)))
    cex <- c("sample", "group") %in%
        colnames(SummarizedExperiment::colData(sce))
    if (any(cex)) {
        stop("'sce' already have column(s) named ",
             paste(c("sample", "group")[cex], collapse = ", "))
    }
    stopifnot(all(!duplicated(sampleAnnot$sample)))

    sce$sample <- colnames(sce)

    if (!all(sce$sample %in% sampleAnnot$sample)) {
        stop("Some samples are missing from the sample annotation: ",
             paste(setdiff(sce$sample, sampleAnnot$sample), collapse = ", "))
    }

    sce$group <- sampleAnnot$group[match(sce$sample, sampleAnnot$sample)]
    for (cn in setdiff(colnames(sampleAnnot), c("sample", "group"))) {
        if (cn %in% colnames(SummarizedExperiment::colData(sce))) {
            stop("Column already exists in SummarizedExperiment: ", cn)
        }
        SummarizedExperiment::colData(sce)[[cn]] <-
            sampleAnnot[[cn]][match(sce$sample, sampleAnnot$sample)]
    }

    sce
}
