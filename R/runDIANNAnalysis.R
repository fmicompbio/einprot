#' Run analysis on DIA-NN data
#'
#' Launch an analysis workflow on quantifications obtained with \code{DIA-NN}.
#' Note that DIA-NN support in einprot is currently experimental - please be
#' aware that the interface may change, and interpret results with caution.
#'
#' @inheritParams runMaxQuantAnalysis
#'
#' @param diannFile Character string pointing to the DIA-NN
#'     \code{pg_matrix.tsv}, \code{pr_matrix.tsv} or main \code{report.tsv}
#'     file. File paths will be expressed in canonical form (using
#'     \code{normalizePath()}) before they are processed.
#' @param diannFileType Character string indicating what type of input file
#'     \code{diannFile} represents. Either \code{"pg_matrix"},
#'     \code{"pr_matrix"} or \code{"main_report"}.
#' @param outLevel Character string indicating the desired output level.
#'     Either \code{"pg"} or \code{"pr"}.
#' @param diannLogFile Character string pointing to the DIA-NN log file.
#'     File paths will be expressed in canonical form (using
#'     \code{normalizePath()}) before they are processed.
#' @param aName Character scalar indicating the desired name of the main
#'     assay (if \code{diannFileType} is \code{"pg_matrix"} or
#'     \code{"pr_matrix"}), or the column to use for the main assay (if
#'     \code{diannFileType} is \code{"main_report"}).
#' @param minScore Numeric, minimum score for a protein to be retained in the
#'     analysis. Set to \code{NULL} if no score filtering is desired.
#' @param minPeptides Numeric, minimum number of peptides for a protein to be
#'     retained in the analysis. Set to \code{NULL} if no filtering on the
#'     number of peptides is desired.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns Invisibly, the path to the compiled html report.
#'
#' @examples
#' if (interactive()) {
#'     sampleAnnot <- read.delim(
#'         system.file("extdata/diann_example/PXD028735_sampleAnnot.txt",
#'                     package = "einprot"))
#'
#'     ## Basic analysis, pg_matrix.tsv
#'     out <- runDIANNAnalysis(
#'         outputDir = tempdir(),
#'         outputBaseName = "DIANN_LFQ_basic",
#'         species = "human",
#'         diannFile = system.file("extdata/diann_example/PXD028735.pg_matrix.tsv",
#'                                 package = "einprot"),
#'         diannFileType = "pg_matrix",
#'         outLevel = "pg",
#'         diannLogFile = system.file("extdata/diann_example/diann-output.log.txt",
#'                                    package = "einprot"),
#'         sampleAnnot = sampleAnnot,
#'         includeFeatureCollections = "complexes",
#'         stringIdCol = NULL,
#'         aName = "MaxLFQ",
#'         idCol = function(df) combineIds(df, combineCols = c("Genes", "Protein.Ids")),
#'         labelCol = function(df) getFirstId(df, colName = "Protein.Names"),
#'         geneIdCol = function(df) getFirstId(df, colName = "Genes"),
#'         proteinIdCol = "Protein.Ids"
#'     )
#'     ## Output file
#'     out
#'
#'     ## Basic analysis, main report
#'     outM <- runDIANNAnalysis(
#'         outputDir = tempdir(),
#'         outputBaseName = "DIANN_LFQ_basic",
#'         species = "human",
#'         diannFile = system.file("extdata/diann_example/PXD028735.report.tsv",
#'                                 package = "einprot"),
#'         diannFileType = "main_report",
#'         outLevel = "pg",
#'         diannLogFile = system.file("extdata/diann_example/diann-output.log.txt",
#'                                    package = "einprot"),
#'         sampleAnnot = sampleAnnot,
#'         includeFeatureCollections = "complexes",
#'         stringIdCol = NULL,
#'         aName = "PG.MaxLFQ",
#'         idCol = "Protein.Group",
#'         labelCol = "Protein.Group",
#'         geneIdCol = NULL,
#'         proteinIdCol = "Protein.Group"
#'     )
#'     ## Output file
#'     outM
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
#' @importFrom ggalt geom_encircle
#' @importFrom plotly ggplotly
#' @importFrom ComplexHeatmap Heatmap columnAnnotation draw
#'
runDIANNAnalysis <- function(
    templateRmd = system.file("extdata/process_basic_template.Rmd",
                              package = "einprot"),
    outputDir = ".", outputBaseName = "DIANNAnalysis",
    reportTitle = "DIA-NN LFQ data processing", reportAuthor = "",
    forceOverwrite = FALSE,
    experimentInfo = list(), species, diannFile, diannFileType,
    outLevel, diannLogFile, aName,
    idCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                        "Majority.protein.IDs")),
    labelCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                           "Majority.protein.IDs")),
    geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
    proteinIdCol = "Majority.protein.IDs",
    stringIdCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                              "Majority.protein.IDs"),
                                          combineWhen = "missing",
                                          makeUnique = FALSE),
    extraFeatureCols = NULL,
    sampleAnnot,
    includeOnlySamples = "", excludeSamples = "",
    minScore = 10, minPeptides = 2, imputeMethod = "MinProb",
    assaysForExport = NULL, addHeatmaps = TRUE, mergeGroups = list(),
    comparisons = list(),
    ctrlGroup = "", allPairwiseComparisons = TRUE, singleFit = TRUE,
    subtractBaseline = FALSE, baselineGroup = "", normMethod = "none",
    spikeFeatures = NULL, stattest = "limma", minNbrValidValues = 2,
    minlFC = 0, samSignificance = TRUE, nperm = 250, volcanoAdjPvalThr = 0.05,
    volcanoLog2FCThr = 1, volcanoMaxFeatures = 25, volcanoLabelSign = "both",
    volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
    addInteractiveVolcanos = FALSE, interactiveDisplayColumns = NULL,
    interactiveGroupColumn = NULL, complexFDRThr = 0.1,
    maxNbrComplexesToPlot = Inf, seed = 42,
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
    .checkArgumentsDIANN(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, reportTitle = reportTitle,
        reportAuthor = reportAuthor, forceOverwrite = forceOverwrite,
        experimentInfo = experimentInfo, species = species,
        diannFile = diannFile, diannFileType = diannFileType,
        outLevel = outLevel, diannLogFile = diannLogFile, aName = aName,
        idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
        proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
        extraFeatureCols = extraFeatureCols, sampleAnnot = sampleAnnot,
        includeOnlySamples = includeOnlySamples,
        excludeSamples = excludeSamples, minScore = minScore,
        minPeptides = minPeptides, imputeMethod = imputeMethod,
        assaysForExport = assaysForExport, addHeatmaps = addHeatmaps,
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
        maxNbrComplexesToPlot = maxNbrComplexesToPlot, seed = seed,
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
    diannFile <- normalizePath(diannFile)
    diannLogFile <- normalizePath(diannLogFile)

    ## -------------------------------------------------------------------------
    ## Copy Rmd template and insert arguments
    ## -------------------------------------------------------------------------
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentInfo = experimentInfo, species = species,
             diannFile = diannFile, diannFileType = diannFileType,
             outLevel = outLevel, diannLogFile = diannLogFile, aName = aName,
             idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
             proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
             extraFeatureCols = extraFeatureCols,
             reportTitle = reportTitle, reportAuthor = reportAuthor,
             sampleAnnot = sampleAnnot,
             includeOnlySamples = includeOnlySamples,
             excludeSamples = excludeSamples, minScore = minScore,
             minPeptides = minPeptides, imputeMethod = imputeMethod,
             assaysForExport = assaysForExport, addHeatmaps = addHeatmaps,
             mergeGroups = mergeGroups,
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
             maxNbrComplexesToPlot = maxNbrComplexesToPlot, seed = seed,
             includeFeatureCollections = includeFeatureCollections,
             minSizeToKeepSet = minSizeToKeepSet,
             customComplexes = customComplexes, complexSpecies = complexSpecies,
             complexDbPath = complexDbPath, stringVersion = stringVersion,
             stringDir = stringDir, linkTableColumns = linkTableColumns,
             expType = "DIANN")
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
