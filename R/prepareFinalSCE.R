#' Prepare "final" SingleCellExperiment for use with iSEE
#'
#' This function is intended to use within the einprot workflows and assumes
#' that the data has been processed as outlined in there.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param baseFileName Character scalar or \code{NULL}, the base file name of
#'     the output files. If \code{NULL}, no result files are generated.
#' @param featureCollections List of \code{CharacterList}s with results
#'     from gene set testing.
#' @param expType Character scalar, either "MaxQuant" or "ProteomeDiscoverer"
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{SingleCellExperiment} object. If \code{baseFileName} is
#'     not \code{NULL}, also save a text file with feature annotation columns
#'     that will not be included in the final SCE object.
#'
#' @importFrom SummarizedExperiment assayNames assay rowData colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom iSEEu registerLogFCFields registerAveAbFields
#'     registerPValueFields registerFeatureSetCollections
#' @importFrom utils write.table
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
prepareFinalSCE <- function(sce, baseFileName, featureCollections, expType) {
    .assertVector(x = sce, type = "SingleCellExperiment")
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = expType, type = "character",
                  validValues = c("MaxQuant", "ProteomeDiscoverer"))

    if (expType == "MaxQuant") {
        ## If not already there, also include log-transformed iBAQ values
        ## for use in heatmaps
        if (("iBAQ" %in% SummarizedExperiment::assayNames(sce)) &&
            !("log2_iBAQ_withNA" %in% SummarizedExperiment::assayNames(sce))) {
            tmplogibaq <- log2(as.matrix(SummarizedExperiment::assay(sce, "iBAQ")))
            tmplogibaq[!is.finite(tmplogibaq)] <- NA
            SummarizedExperiment::assay(sce, "log2_iBAQ_withNA") <- tmplogibaq
        }

        ## Remove some columns from rowData (save to text file)
        colsToRemove <- c("Peptide.counts.all", "Peptide.counts.razor.unique",
                          "Peptide.counts.unique", "Fasta.headers",
                          "Peptide.IDs", "Peptide.is.razor", "Mod.peptide.IDs",
                          "Evidence.IDs", "MS.MS.IDs", "Best.MS.MS",
                          "Sequence.lengths", "Oxidation.M.site.IDs",
                          "Oxidation.M.site.positions")
        if (!is.null(baseFileName)) {
            utils::write.table(as.data.frame(
                SummarizedExperiment::rowData(sce)[, colsToRemove]) %>%
                    tibble::rownames_to_column("ID"),
                file = paste0(baseFileName, "_sce_extra_annots.tsv"),
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        }
        SummarizedExperiment::rowData(sce) <-
            SummarizedExperiment::rowData(sce)[, !colnames(
                SummarizedExperiment::rowData(sce)) %in% colsToRemove]
    } else if (expType == "ProteomeDiscoverer") {
        ## Remove some columns from rowData (save to text file)
        colsToRemove <- c(grep(paste(c(
            "Found.in.Fraction", "Abundances.Grouped.",
            "Abundances.Grouped.CV.in.Percent.", "Abundances.Grouped.Count.",
            "Abundances.Normalized.", "Found.in.Sample.Group.",
            "Found.in.Sample", "Proteins.Unique.Sequence.ID"), collapse = "|"),
            colnames(rowData(sce)), value = TRUE), "Sequence", "GO.Accessions")
        write.table(as.data.frame(rowData(sce)[, colsToRemove]) %>%
                        tibble::rownames_to_column("ID"),
                    file = sub("\\.Rmd$", paste0("_sce_extra_annots.tsv"),
                               knitr::current_input(dir = TRUE)),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        rowData(sce) <- SummarizedExperiment::rowData(sce)[, !colnames(
            SummarizedExperiment::rowData(sce)) %in% colsToRemove]
    }

    ## Register logFC/AveAb/pvalue fields for use in iSEE
    sce <- iSEEu::registerLogFCFields(
        sce, grep("\\.logFC$", colnames(SummarizedExperiment::rowData(sce)),
                  value = TRUE)
    )
    sce <- iSEEu::registerAveAbFields(
        sce, grep("\\.AveExpr$", colnames(SummarizedExperiment::rowData(sce)),
                  value = TRUE)
    )
    sce <- iSEEu::registerPValueFields(
        sce, grep("\\.P.Value$", colnames(SummarizedExperiment::rowData(sce)),
                  value = TRUE)
    )

    if (length(featureCollections) > 0) {
        sce <- iSEEu::registerFeatureSetCollections(sce, featureCollections)
    }

    sce
}