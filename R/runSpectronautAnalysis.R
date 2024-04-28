#' Run analysis on Spectronaut data
#'
#' Launch an analysis workflow on quantifications obtained with
#' \code{Spectronaut}. Note that Spectronaut support in einprot is currently
#' experimental - please be aware that the interface may change, and interpret
#' results with caution.
#'
#' @inheritParams runMaxQuantAnalysis
#'
#' @param spectronautFile Character string pointing to the Spectronaut
#'     \code{Report.tsv} file. File paths will be expressed in canonical form
#'     (using \code{normalizePath()}) before they are processed.
#' @param spectronautFileType Character string indicating what type of input
#'     file \code{spectronautFile} represents. Either \code{"pg_pivot"}
#'     or \code{"long_format"}.
#' @param outLevel Character string indicating the desired output level.
#'     Currently only \code{"pg"} is supported.
#' @param spectronautLogFile Character string pointing to the Spectronaut
#'     log file. File paths will be expressed in canonical form (using
#'     \code{normalizePath()}) before they are processed.
#' @param iColPattern Character scalar defining a regular expression to
#'     identify sample columns (only used if \code{spectronautFileType} is
#'     \code{"pg_pivot"}. Typically one of \code{".PG.Quantity$"} or
#'     \code{".PG.IBAQ$"}.
#' @param aName Character scalar indicating the column to use for the main
#'     assay.
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
runSpectronautAnalysis <- function(
    templateRmd = system.file("extdata/process_basic_template.Rmd",
                              package = "einprot"),
    outputDir = ".", outputBaseName = "SpectronautAnalysis",
    reportTitle = "Spectronaut LFQ data processing", reportAuthor = "",
    forceOverwrite = FALSE,
    experimentInfo = list(), species, spectronautFile, spectronautFileType,
    outLevel, spectronautLogFile, aName,
    idCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                        "PG.ProteinGroups")),
    labelCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                           "PG.ProteinGroups")),
    geneIdCol = function(df) getFirstId(df, colName = "PG.Genes"),
    proteinIdCol = "PG.ProteinGroups",
    stringIdCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                              "PG.ProteinGroups"),
                                          combineWhen = "missing",
                                          makeUnique = FALSE),
    extraFeatureCols = NULL, iColPattern = ".PG.Quantity$",
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
    .checkArgumentsSpectronaut(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, reportTitle = reportTitle,
        reportAuthor = reportAuthor, forceOverwrite = forceOverwrite,
        experimentInfo = experimentInfo, species = species,
        spectronautFile = spectronautFile,
        spectronautFileType = spectronautFileType, outLevel = outLevel,
        spectronautLogFile = spectronautLogFile, aName = aName,
        idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
        proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
        extraFeatureCols = extraFeatureCols, iColPattern = iColPattern,
        sampleAnnot = sampleAnnot, includeOnlySamples = includeOnlySamples,
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
    spectronautFile <- normalizePath(spectronautFile)
    spectronautLogFile <- normalizePath(spectronautLogFile)

    ## -------------------------------------------------------------------------
    ## Copy Rmd template and insert arguments
    ## -------------------------------------------------------------------------
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentInfo = experimentInfo, species = species,
             spectronautFile = spectronautFile,
             spectronautFileType = spectronautFileType, outLevel = outLevel,
             spectronautLogFile = spectronautLogFile, aName = aName,
             idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
             proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
             extraFeatureCols = extraFeatureCols,
             reportTitle = reportTitle, reportAuthor = reportAuthor,
             iColPattern = iColPattern, sampleAnnot = sampleAnnot,
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
             expType = "Spectronaut")
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
