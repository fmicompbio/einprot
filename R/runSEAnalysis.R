#' Run analysis on pre-imported data
#'
#' Launch an analysis workflow on data stored in a
#' \code{\link[SummarizedExperiment]{SummarizedExperiment}} object, generated
#' from arbitrary source data.
#'
#' @inheritParams runMaxQuantAnalysis
#'
#' @param inFile Character string pointing to an \code{.rds} file containing
#'     a serialized \code{\link[SummarizedExperiment]{SummarizedExperiment}}
#'     object (or a derivative, such as
#'     \code{\link[SingleCellExperiment]{SingleCellExperiment}}, with features
#'     (e.g. proteins) in rows and samples in columns. The object
#'     must have at least one assay (\code{aName}), containing raw intensity
#'     values or equivalent (these will be automatically log-transformed by
#'     einprot).
#' @param summaryInfo Named list with information about the identification
#'     and quantification workflow. This information will be displayed in
#'     table format in the report, under the 'Analysis summary' header. einprot
#'     provides functions to generate such lists for common processing tools,
#'     such as \code{\link{readMaxQuantXML}}, \code{\link{readFragPipeInfo}},
#'     \code{\link{readProteomeDiscovererInfo}} and \code{\link{readDIANNInfo}}.
#' @param aName Character scalar providing the name of the assay of the
#'     \code{SummarizedExperiment} object in \code{inputFile} that should be
#'     used as the main assay for the analysis.
#'     The assay should contain values on the "raw intensity scale" (i.e.,
#'     they will be log-transformed as part of the analysis workflow). Missing
#'     values should be represented by either zeros or NA values.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns Invisibly, the path to the compiled html report.
#'
#' @examples
#' if (interactive()) {
#'     sampleAnnot <- read.delim(
#'         system.file("extdata/mq_example/1356_sampleAnnot.txt",
#'                     package = "einprot"))
#'     se <- importExperiment(
#'         system.file("extdata/mq_example/1356_proteinGroups.txt",
#'                     package = "einprot"),
#'         iColPattern = "^LFQ.intensity.")$sce
#'     se <- addSampleAnnots(se, sampleAnnot)
#'     tmpf <- tempfile(fileext = ".rds")
#'     saveRDS(se, file = tmpf)
#'     settings <- readMaxQuantXML(
#'         system.file("extdata/mq_example/1356_mqpar.xml",
#'                     package = "einprot"))
#'
#'     ## Basic analysis
#'     out <- runSEAnalysis(
#'         outputDir = tempdir(),
#'         outputBaseName = "MQ_LFQ_basic",
#'         species = "mouse",
#'         inFile = tmpf,
#'         summaryInfo = settings,
#'         aName = "LFQ.intensity",
#'         idCol = function(df) getFirstId(df, "Majority.protein.IDs", ";"),
#'         labelCol = function(df) getFirstId(df, "Majority.protein.IDs", ";"),
#'         filtersSE = list()
#'     )
#'     ## Output file
#'     out
#' }
#'
#' @importFrom xfun Rscript_call
#' @importFrom rmarkdown render
#' @importFrom readr read_file write_file
#' @import STRINGdb
#' @importFrom SummarizedExperiment rowData colData assay assayNames
#' @importFrom DT datatable
#' @importFrom limma removeBatchEffect
#' @importFrom ExploreModelMatrix VisualizeDesign
#' @importFrom cowplot plot_grid theme_cowplot
#' @importFrom htmltools tagList
#' @importFrom dplyr %>% select starts_with full_join filter matches everything
#'     mutate
#' @importFrom knitr current_input
#' @importFrom ComplexUpset upset
#' @importFrom ggplot2 ggplot aes geom_bar coord_flip theme_bw labs theme
#'     element_text geom_point ggtitle
#' @importFrom tibble rownames_to_column
#' @importFrom S4Vectors metadata
#' @importFrom scater runPCA
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom BiocSingular ExactParam
#' @importFrom plotly ggplotly
#' @importFrom ComplexHeatmap Heatmap columnAnnotation draw
#'
runSEAnalysis <- function(
        templateRmd = system.file("extdata/process_basic_template.Rmd",
                                  package = "einprot"),
        outputDir = ".", outputBaseName = "SEAnalysis",
        reportTitle = "Pregenerated data processing", reportAuthor = "",
        forceOverwrite = FALSE,
        experimentInfo = list(), species,
        inFile, summaryInfo = NULL, aName,
        idCol, labelCol, geneIdCol = NULL, proteinIdCol = NULL, stringIdCol = NULL,
        extraFeatureCols = NULL, filtersSE,
        imputeMethod = "MinProb", assaysForExport = c("iBAQ", "Top3"),
        addAbundanceValues = TRUE, addHeatmaps = TRUE,
        mergeGroups = list(), comparisons = list(),
        ctrlGroup = "", allPairwiseComparisons = TRUE, singleFit = TRUE,
        subtractBaseline = FALSE, baselineGroup = "", normMethod = "none",
        spikeFeatures = NULL, stattest = "limma", minNbrValidValues = 2,
        minlFC = 0, samSignificance = TRUE, nperm = 250, volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1, volcanoMaxFeatures = 25, volcanoLabelSign = "both",
        volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
        addInteractiveVolcanos = FALSE, interactiveDisplayColumns = NULL,
        interactiveGroupColumn = NULL, complexFDRThr = 0.1,
        maxNbrComplexesToPlot = Inf, maxComplexSimilarity = 1, seed = 42,
        includeFeatureCollections = c(), minSizeToKeepSet = 2,
        customComplexes = list(),
        complexSpecies = "all", complexDbPath = NULL, stringVersion = "11.5",
        stringDir = NULL, linkTableColumns = c(), customYml = NULL, doRender = TRUE
) {
    ## -------------------------------------------------------------------------
    ## Fix ctrlGroup/mergeGroups
    ## -------------------------------------------------------------------------
    ## For backward compatibility: If mergeGroups is list(), and ctrlGroup
    ## is a vector (the way things were specified before v0.3.2), add the
    ## merged ctrl group to mergeGroups. Raise an error if mergeGroups is
    ## not empty and ctrlGroup is a vector.
    if (length(mergeGroups) > 0 && length(ctrlGroup) > 1) {
        stop("If 'mergeGroups' is specified, 'ctrlGroup' should not ",
             "be a vector.")
    }
    if (length(mergeGroups) == 0 && length(ctrlGroup) > 1) {
        mergeGroups <- list()
        newCtrlName <- paste(ctrlGroup, collapse = ".")
        mergeGroups[[newCtrlName]] <- ctrlGroup
        ctrlGroup <- newCtrlName
    }

    ## -------------------------------------------------------------------------
    ## Check arguments
    ## -------------------------------------------------------------------------
    .checkArgumentsSE(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, reportTitle = reportTitle,
        reportAuthor = reportAuthor, forceOverwrite = forceOverwrite,
        experimentInfo = experimentInfo, species = species,
        inFile = inFile, summaryInfo = summaryInfo, aName = aName,
        idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
        proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
        extraFeatureCols = extraFeatureCols, filtersSE = filtersSE,
        imputeMethod = imputeMethod, assaysForExport = assaysForExport,
        addAbundanceValues = addAbundanceValues, addHeatmaps = addHeatmaps,
        mergeGroups = mergeGroups,
        comparisons = comparisons, ctrlGroup = ctrlGroup,
        allPairwiseComparisons = allPairwiseComparisons, singleFit = singleFit,
        subtractBaseline = subtractBaseline, baselineGroup = baselineGroup,
        normMethod = normMethod, spikeFeatures = spikeFeatures,
        stattest = stattest, minNbrValidValues = minNbrValidValues,
        minlFC = minlFC, samSignificance = samSignificance,
        nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
        volcanoLog2FCThr = volcanoLog2FCThr,
        volcanoMaxFeatures = volcanoMaxFeatures,
        volcanoLabelSign = volcanoLabelSign,
        volcanoS0 = volcanoS0, volcanoFeaturesToLabel = volcanoFeaturesToLabel,
        addInteractiveVolcanos = addInteractiveVolcanos,
        interactiveDisplayColumns = interactiveDisplayColumns,
        interactiveGroupColumn = interactiveGroupColumn,
        complexFDRThr = complexFDRThr,
        maxNbrComplexesToPlot = maxNbrComplexesToPlot,
        maxComplexSimilarity = maxComplexSimilarity, seed = seed,
        includeFeatureCollections = includeFeatureCollections,
        minSizeToKeepSet = minSizeToKeepSet,
        customComplexes = customComplexes, complexSpecies = complexSpecies,
        complexDbPath = complexDbPath, stringVersion = stringVersion,
        stringDir = stringDir, linkTableColumns = linkTableColumns,
        customYml = customYml, doRender = doRender)

    ## If pandoc is not available, don't run it (just generate .md file)
    ## Gives a warning if pandoc and/or pandoc-citeproc is not available
    pandocOK <- .checkPandoc(ignorePandoc = TRUE)

    ## -------------------------------------------------------------------------
    ## Normalize paths
    ## -------------------------------------------------------------------------
    inFile <- normalizePath(inFile)

    ## -------------------------------------------------------------------------
    ## Copy Rmd template and insert arguments
    ## -------------------------------------------------------------------------
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentInfo = experimentInfo, species = species,
             inFile = inFile, summaryInfo = summaryInfo, aName = aName,
             idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
             proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
             extraFeatureCols = extraFeatureCols, filtersSE = filtersSE,
             reportTitle = reportTitle, reportAuthor = reportAuthor,
             imputeMethod = imputeMethod, assaysForExport = assaysForExport,
             addAbundanceValues = addAbundanceValues,
             addHeatmaps = addHeatmaps, mergeGroups = mergeGroups,
             comparisons = comparisons, ctrlGroup = ctrlGroup,
             allPairwiseComparisons = allPairwiseComparisons,
             singleFit = singleFit,
             subtractBaseline = subtractBaseline, baselineGroup = baselineGroup,
             normMethod = normMethod, spikeFeatures = spikeFeatures,
             stattest = stattest, minNbrValidValues = minNbrValidValues,
             minlFC = minlFC, samSignificance = samSignificance,
             nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
             volcanoLog2FCThr = volcanoLog2FCThr,
             volcanoMaxFeatures = volcanoMaxFeatures,
             volcanoLabelSign = volcanoLabelSign,
             volcanoS0 = volcanoS0,
             volcanoFeaturesToLabel = volcanoFeaturesToLabel,
             addInteractiveVolcanos = addInteractiveVolcanos,
             interactiveDisplayColumns = interactiveDisplayColumns,
             interactiveGroupColumn = interactiveGroupColumn,
             complexFDRThr = complexFDRThr,
             maxNbrComplexesToPlot = maxNbrComplexesToPlot,
             maxComplexSimilarity = maxComplexSimilarity, seed = seed,
             includeFeatureCollections = includeFeatureCollections,
             minSizeToKeepSet = minSizeToKeepSet,
             customComplexes = customComplexes, complexSpecies = complexSpecies,
             complexDbPath = complexDbPath, stringVersion = stringVersion,
             stringDir = stringDir, linkTableColumns = linkTableColumns,
             expType = "SE")
    )

    ## Read Rmd
    rmd <- readr::read_file(templateRmd)

    ## Determine where to insert the config chunk
    ## From https://community.rstudio.com/t/how-to-write-r-script-into-rmd-as-functioning-code-chunk/37453/2
    # header_regex <- sprintf("\\{\\{%sStart\\}\\}(.*?)\\{\\{%sEnd\\}\\}",
    #                         confighook,
    #                         confighook)
    header_regex <- sprintf("{{%sStart}}\n\n{{%sEnd}}",
                            confighook,
                            confighook)

    ## Replace hooks with config chunk
    output <- sub(header_regex, configchunk, rmd, fixed = TRUE)

    ## Similarly, add any custom yaml
    ymlhook <- "YmlParameters"
    # header_regex_yml <- sprintf("\\{\\{%sStart\\}\\}(.*?)\\{\\{%sEnd\\}\\}",
    #                             ymlhook,
    #                             ymlhook)
    header_regex_yml <- sprintf("{{%sStart}}\n\n{{%sEnd}}",
                                ymlhook,
                                ymlhook)
    if (!is.null(customYml)) {
        customYml <- paste(readLines(customYml), collapse = "\n")
    } else {
        customYml <- ""
    }
    output <- sub(header_regex_yml, customYml, output, fixed = TRUE)

    ## Write output to file
    if (!dir.exists(outputDir)) {
        dir.create(outputDir, recursive = TRUE)
    }
    outputFile <- file.path(outputDir, paste0(outputBaseName, ".Rmd"))
    if (file.exists(outputFile) && !forceOverwrite) {
        stop(outputFile,
             " already exists and forceOverwrite = FALSE, stopping.")
    } else if (file.exists(outputFile) && forceOverwrite) {
        message(outputFile,
                " already exists but forceOverwrite = TRUE, overwriting.")
    }
    readr::write_file(output, file = outputFile)

    ## -------------------------------------------------------------------------
    ## Render the Rmd file
    ## -------------------------------------------------------------------------
    args <- list()
    args$input <- outputFile
    args$output_format <- "html_document"
    args$output_dir <- outputDir
    args$intermediates_dir <- outputDir
    args$quiet <- FALSE
    args$run_pandoc <- pandocOK

    if (doRender) {
        #nocov start
        outputReport <- xfun::Rscript_call(
            rmarkdown::render,
            args
        )
        #nocov end
    } else {
        outputReport <- outputFile
    }

    ## -------------------------------------------------------------------------
    ## Return (invisibly) the path to the rendered html file
    ## -------------------------------------------------------------------------
    invisible(outputReport)
}
