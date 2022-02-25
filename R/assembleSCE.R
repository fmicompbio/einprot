#' Assemble SingleCellExperiment object from QFeatures object
#'
#' This function is intended to use within the einprot workflows and assumes
#' that the data has been processed as outlined in there.
#'
#' @param qft A \code{QFeatures} object.
#' @param aName Character scalar giving the name of the base assay.
#' @param testResults A \code{data.frame} with test results for all
#'     comparisons.
#' @param iColPattern Character scalar providing the pattern used to
#'     identify intensity columns from the input files.
#' @param iColsAll Character vector with all intensity column names.
#' @param baseFileName Character scalar or \code{NULL}, the base file name of
#'     the output files. If \code{NULL}, no result files are generated.
#' @param nbrNA A \code{list} with information about the number of NA values
#'     for rows and columns, respectively.
#' @param featureCollections List of \code{CharacterList}s with results
#'     from gene set testing.
#' @param expType Character scalar, either "MaxQuant" or "ProteomeDiscoverer"
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @importFrom SummarizedExperiment assayNames assay rowData colData
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr filter
#' @importFrom iSEEu registerLogFCFields registerAveAbFields
#'     registerPValueFields registerFeatureSetCollections
#' @importFrom utils write.table
#' @importFrom SingleCellExperiment SingleCellExperiment
#'
assembleSCE <- function(qft, aName, testResults, iColPattern,
                        iColsAll, baseFileName, nbrNA,
                        featureCollections, expType) {
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = aName, type =  "character", validValues = names(qft))
    .assertVector(x = testResults, type = "data.frame")
    .assertScalar(x = iColPattern, type = "character")
    .assertVector(x = iColsAll, type = "character")
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertVector(x = nbrNA, type = "list")
    .assertVector(x = names(nbrNA), type = "character",
                  validValues = c("nNA", "nNArows", "nNAcols"))
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = expType, type = "character",
                  validValues = c("MaxQuant", "ProteomeDiscoverer"))

    ## Make 'base' SCE and add assays
    sce <- qft[[aName]]
    SummarizedExperiment::assayNames(sce) <- aName
    for (nm in names(qft)) {
        SummarizedExperiment::assay(sce, nm) <-
            SummarizedExperiment::assay(qft[[nm]])
    }

    ## Add sample and feature annotations
    SummarizedExperiment::colData(sce) <-
        SummarizedExperiment::colData(qft)
    stopifnot(rownames(sce) == rownames(testResults))
    SummarizedExperiment::rowData(sce) <-
        cbind(SummarizedExperiment::rowData(sce), testResults)

    if (expType == "MaxQuant") {
        ## Move some values to assays rather than rowData
        colnames(SummarizedExperiment::rowData(sce)) <-
            gsub("\\.+$", "", colnames(SummarizedExperiment::rowData(sce)))
        colnames(SummarizedExperiment::rowData(sce)) <-
            gsub("\\.+", ".", colnames(SummarizedExperiment::rowData(sce)))
        moveToAssay <- c("MS.MS.Count.", "LFQ.intensity.", "Intensity.",
                         "Sequence.coverage.", "Unique.peptides.",
                         "Razor.unique.peptides.", "Peptides.", "iBAQ.")
        for (mta in moveToAssay) {
            cols <- sub(iColPattern, mta, colnames(sce))
            if (all(cols %in% colnames(SummarizedExperiment::rowData(sce)))) {
                ## Use withDimnames = FALSE since colnames are not identical
                SummarizedExperiment::assay(sce, sub("\\.$", "", mta),
                                            withDimnames = FALSE) <-
                    as.matrix(SummarizedExperiment::rowData(sce)[, cols])
                ## For iBAQ, also include log-transformed values (for use in heatmaps)
                if (mta == "iBAQ." && !("log2_iBAQ_withNA" %in%
                                        SummarizedExperiment::assayNames(sce))) {
                    tmplogibaq <- log2(as.matrix(SummarizedExperiment::rowData(sce)[, cols]))
                    tmplogibaq[!is.finite(tmplogibaq)] <- NA
                    SummarizedExperiment::assay(sce, paste0("log2_", sub("\\.$", "", mta), "_withNA"),
                                                withDimnames = FALSE) <- tmplogibaq
                }
            }
            ## Remove all columns corresponding to this assay from the rowData,
            ## even if only some samples are retained
            colsall <- sub(iColPattern, mta, iColsAll)
            SummarizedExperiment::rowData(sce) <-
                SummarizedExperiment::rowData(sce)[, !colnames(rowData(sce)) %in% colsall]
        }

        ## Remove some columns from rowData (save to text file)
        colsToRemove <- c("Peptide.counts.all", "Peptide.counts.razor.unique",
                          "Peptide.counts.unique", "Fasta.headers",
                          "Peptide.IDs", "Peptide.is.razor", "Mod.peptide.IDs",
                          "Evidence.IDs", "MS.MS.IDs", "Best.MS.MS", "Sequence.lengths",
                          "Oxidation.M.site.IDs", "Oxidation.M.site.positions")
        if (!is.null(baseFileName)) {
            utils::write.table(as.data.frame(
                SummarizedExperiment::rowData(sce)[, colsToRemove]) %>%
                    tibble::rownames_to_column("ID"),
                file = paste0(baseFileName, "_sce_extra_annots.tsv"),
                row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        }
        SummarizedExperiment::rowData(sce) <-
            SummarizedExperiment::rowData(sce)[, !colnames(SummarizedExperiment::rowData(sce)) %in%
                                                   colsToRemove]
    } else if (expType == "ProteomeDiscoverer") {
        ## Move some values to assays rather than rowData
        if (iColPattern == "^Abundance\\.F.+\\.Sample\\.") {
            cols <- sub("Abundance", "Abundances.Count", colnames(sce))
            colsall <- sub("Abundance", "Abundances.Count", iColsAll)
            anm <- "Abundances.Count"
        } else if (iColPattern == "^Abundances.Count\\.F.+\\.Sample\\.") {
            cols <- sub("Abundances.Count", "Abundance", colnames(sce))
            colsall <- sub("Abundances.Count", "Abundance", iColsAll)
            anm <- ifelse(aName == "Abundance", "AbundanceRaw", "Abundance")
        }
        if (all(cols %in% colnames(rowData(sce)))) {
            ## Use withDimnames = FALSE since colnames are not identical
            assay(sce, anm, withDimnames = FALSE) <-
                as.matrix(rowData(sce)[, cols])
        }
        ## Remove all columns corresponding to this assay from the rowData,
        ## even if only some samples are retained
        rowData(sce) <- rowData(sce)[, !colnames(rowData(sce)) %in% colsall]

        ## Remove some columns from rowData (save to text file)
        colsToRemove <- c(grep("Found.in.Fraction|Abundances.Grouped.|Abundances.Grouped.CV.in.Percent.|Abundances.Grouped.Count.|Abundances.Normalized.|Found.in.Sample.Group.|Found.in.Sample|Proteins.Unique.Sequence.ID",
                               colnames(rowData(sce)), value = TRUE), "Sequence", "GO.Accessions")
        write.table(as.data.frame(rowData(sce)[, colsToRemove]) %>%
                        tibble::rownames_to_column("ID"),
                    file = sub("\\.Rmd$", paste0("_sce_extra_annots.tsv"),
                               knitr::current_input(dir = TRUE)),
                    row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
        rowData(sce) <- rowData(sce)[, !colnames(rowData(sce)) %in% colsToRemove]
    }

    ## Add information about missing values
    nacols <- as.data.frame(nbrNA$nNAcols) %>%
        dplyr::filter(assay == aName)
    stopifnot(all(sub(iColPattern, "", nacols$name) == sce$sample))
    SummarizedExperiment::colData(sce) <-
        cbind(SummarizedExperiment::colData(sce), nacols[, c("nNA", "pNA")])
    narows <- as.data.frame(nbrNA$nNArows) %>%
        dplyr::filter(assay == aName)
    stopifnot(all(narows$name == rownames(sce)))
    SummarizedExperiment::rowData(sce) <-
        cbind(SummarizedExperiment::rowData(sce), narows[, c("nNA", "pNA")])

    ## Register logFC/AveAb/pvalue fields for use in iSEE
    sce <- iSEEu::registerLogFCFields(
        sce, grep("\\.logFC$", colnames(SummarizedExperiment::rowData(sce)), value = TRUE)
    )
    sce <- iSEEu::registerAveAbFields(
        sce, grep("\\.AveExpr$", colnames(SummarizedExperiment::rowData(sce)), value = TRUE)
    )
    sce <- iSEEu::registerPValueFields(
        sce, grep("\\.P.Value$", colnames(SummarizedExperiment::rowData(sce)), value = TRUE)
    )

    if (length(featureCollections) > 0) {
        sce <- iSEEu::registerFeatureSetCollections(sce, featureCollections)
    }

    sce <- as(sce, "SingleCellExperiment")

    sce
}
