#' Add sample annotations to SummarizedExperiment object
#'
#' @param sce A \code{SummarizedExperiment} object.
#' @param sampleAnnot \code{data.frame} with sample annotations. Must
#'     have at least columns named \code{sample} (which must contain
#'     all the column names of \code{sce}) and \code{group} (which contains
#'     the group assignment for each sample). Can also have a column named
#'     \code{batch}; in this case this will be used as a covariate in the
#'     limma model.
#' @param mergeGroups Named list defining groups to merge. Each entry of
#'     the list corresponds to a new group, and consists of a vector
#'     with the original group names to merge.
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
#' sce <- addSampleAnnots(out$sce, sampleAnnot = sampleAnnot, mergeGroups = list())
#' SummarizedExperiment::colData(sce)  ## group information added to sce
#'
#' @export
#' @author Charlotte Soneson
#'
addSampleAnnots <- function(sce, sampleAnnot,
                            mergeGroups = list()) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertVector(x = sampleAnnot, type = "data.frame")
    .assertVector(x = mergeGroups, type = "list")
    stopifnot(all(c("sample", "group") %in% colnames(sampleAnnot)))
    if (length(mergeGroups) > 0) {
        .assertVector(x = names(mergeGroups), type = "character")
        if (is.null(names(mergeGroups)) || any(names(mergeGroups) == "") ||
            any(duplicated(names(mergeGroups)))) {
            stop("'mergeGroups' must be a named list, without duplicated names")
        }
        if (any(duplicated(unlist(mergeGroups)))) {
            stop("A given name can just be part of one merged group")
        }
    }
    stopifnot(all(!duplicated(sampleAnnot$sample)))

    sce$sample <- colnames(sce)

    if (!all(sce$sample %in% sampleAnnot$sample)) {
        stop("Some samples are missing from the sample annotation: ",
             paste(setdiff(sce$sample, sampleAnnot$sample), collapse = ", "))
    }

    sce$group_orig <- sampleAnnot$group[match(sce$sample, sampleAnnot$sample)]
    if ("batch" %in% colnames(sampleAnnot)) {
        sce$batch <- sampleAnnot$batch[match(sce$sample, sampleAnnot$sample)]
    }

    ## Define a new grouping based on the defined mergeGroups
    sce$group <- sce$group_orig
    for (nm in names(mergeGroups)) {
        sce$group[sce$group_orig %in% mergeGroups[[nm]]] <- nm
    }

    sce
}
