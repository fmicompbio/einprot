#' @export
#' @author Charlotte Soneson
#'
addSampleAnnots <- function(qft, iColPattern, sampleAnnot, mergeGroups) {
    ## Get sample ID by removing the iColPattern from the colnames
    cdn <- sub(iColPattern, "", colnames(qft)[[1]])
    qft$sample <- cdn

    qft$group_orig <- sampleAnnot$group[match(qft$sample, sampleAnnot$sample)]
    if ("batch" %in% colnames(sampleAnnot)) {
        qft$batch <- sampleAnnot$batch[match(qft$sample, sampleAnnot$sample)]
    }

    ## Define a new grouping based on the defined mergeGroups
    qft$group <- qft$group_orig
    for (nm in names(mergeGroups)) {
        qft$group[qft$group_orig %in% mergeGroups[[nm]]] <- nm
    }

    qft
}
