#' Run analysis on PD/TMT data
#'
#' Launch a workflow on quantifications obtained with
#' \code{Proteome Discoverer}.
#'
#' @inheritParams runMaxQuantAnalysis
#'
#' @param pdOutputFolder Character string pointing to the PD/TMT output folder.
#'     Should contain the files \code{pdResultName_InputFiles.txt},
#'     \code{pdResultName_StudyInformation.txt} and
#'     \code{pdResultName_Proteins.txt}. In order to generate the stand-alone
#'     pdf file with QC metrics, additionally the following files should
#'     be present:
#'     \code{pdResultName_PSMs.txt}, \code{pdResultName_PeptideGroups.txt},
#'     \code{pdResultName_MSMSSpectrumInfo.txt},
#'     \code{pdResultName_QuanSpectra.txt}. File paths will be expressed in
#'     canonical form (using \code{normalizePath()}) before they are processed.
#' @param pdResultName Character string providing the base name for the
#'     files in the \code{pdOutputFolder}.
#' @param inputLevel Character string specifying which of the PD files to use
#'     for the analysis. Currently only \code{"Proteins"} and
#'     \code{"PeptideGroups"} are supported.
#' @param pdAnalysisFile Character string pointing to the \code{pdAnalysis}
#'     file. File paths will be expressed in canonical form (using
#'     \code{normalizePath()}) before they are processed.
#' @param modificationsCol Character string pointing to a column containing
#'     modification details. \code{excludeUnmodifiedPeptides} and
#'     \code{keepModifications} will use information from this column. Only
#'     used if \code{inputLevel} is "PeptideGroups".
#' @param excludeUnmodifiedPeptides Logical scalar, whether to filter out
#'     peptides without modifications. Only used if \code{inputLevel} is
#'     "PeptideGroups".
#' @param keepModifications Character string (or \code{NULL}) indicating
#'     which modifications to retain in the analysis. Can be a regular
#'     expression, which will be matched against the \code{modificationsCol}.
#'     If \code{NULL} (the default), all rows are retained. Only used if
#'     \code{inputLevel} is "PeptideGroups".
#' @param iColPattern Regular expression identifying the columns of the PD
#'     \code{Proteins.txt} file to use for the analysis.
#' @param minScore,minDeltaScore Numeric, minimum score for a protein (or
#'     delta score for a peptide group) to be retained in the analysis.
#'     Set to \code{NULL} if no score filtering is desired.
#' @param minPeptides,minPSMs Numeric, minimum number of peptides for a protein
#'     (or PSMs for a peptide group) to be retained in the analysis.
#'     Set to \code{NULL} if no filtering on the number of peptides/PTMs is
#'     desired.
#' @param masterProteinsOnly Logical scalar indicating whether only master
#'     proteins (where the \code{Master} column value is
#'     \code{IsMasterProtein}) should be retained. Only used if
#'     \code{inputLevel} is \code{"Proteins"}.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' if (interactive()) {
#'     ## Basic analysis
#'     out <- runPDTMTAnalysis(
#'         outputDir = tempdir(),
#'         outputBaseName = "PD_TMT_basic",
#'         species = "baker's yeast",
#'         iColPattern = "^Abundance.F.+.Sample.",
#'         pdOutputFolder = system.file("extdata", "pdtmt_example",
#'                                      package = "einprot"),
#'         pdResultName = "Fig2_m23139_RTS_QC_varMods",
#'         inputLevel = "Proteins",
#'         pdAnalysisFile = system.file("extdata", "pdtmt_example",
#'                           "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
#'                           package = "einprot"),
#'         sampleAnnot = data.frame(
#'             sample = c("HIS4KO_S05", "HIS4KO_S06", "HIS4KO_S07", "HIS4KO_S08",
#'                        "MET6KO_S01", "MET6KO_S02", "MET6KO_S03", "MET6KO_S04",
#'                        "URA2KO_S09", "URA2KO_S10", "URA2KO_S11", "URA2KO_S12",
#'                        "WT_S13", "WT_S14", "WT_S15", "WT_S16"),
#'             group = c(rep("HIS4KO", 4), rep("MET6KO", 4), rep("URA2KO", 4),
#'                       rep("WT", 4))),
#'         stringIdCol = NULL
#'     )
#'     ## Output file
#'     out
#' }
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
#'     mutate relocate last_col
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
runPDTMTAnalysis <- function(
    templateRmd = system.file("extdata/process_basic_template.Rmd",
                              package = "einprot"),
    outputDir = ".", outputBaseName = "PDTMTAnalysis",
    reportTitle = "PD data processing", reportAuthor = "",
    forceOverwrite = FALSE,
    experimentInfo = list(), species, pdOutputFolder, pdResultName,
    inputLevel = "Proteins", pdAnalysisFile,
    idCol = function(df) combineIds(df, combineCols = c("Gene.Symbol",
                                                        "Accession")),
    labelCol = function(df) combineIds(df, combineCols = c("Gene.Symbol",
                                                           "Accession")),
    geneIdCol = function(df) getFirstId(df, colName = "Gene.Symbol"),
    proteinIdCol = "Accession",
    stringIdCol = function(df) combineIds(df, combineCols = c("Gene.Symbol",
                                                              "Accession"),
                                          combineWhen = "missing",
                                          makeUnique = FALSE),
    extraFeatureCols = NULL,
    modificationsCol = "Modifications", excludeUnmodifiedPeptides = FALSE,
    keepModifications = NULL, iColPattern, sampleAnnot,
    includeOnlySamples = "", excludeSamples = "",
    minScore = 2, minDeltaScore = 0.2, minPeptides = 2, minPSMs = 2,
    masterProteinsOnly = FALSE, imputeMethod = "MinProb",
    assaysForExport = NULL, addAbundanceValues = TRUE, 
    addHeatmaps = TRUE, mergeGroups = list(), comparisons = list(),
    ctrlGroup = "", allPairwiseComparisons = TRUE, singleFit = TRUE,
    subtractBaseline = FALSE, baselineGroup = "", normMethod = "none",
    spikeFeatures = NULL, stattest = "limma", minNbrValidValues = 2,
    minlFC = 0, samSignificance = FALSE, nperm = 250, volcanoAdjPvalThr = 0.05,
    volcanoLog2FCThr = 1, volcanoMaxFeatures = 25, volcanoLabelSign = "both",
    volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
    addInteractiveVolcanos = FALSE, interactiveDisplayColumns = NULL,
    interactiveGroupColumn = NULL, complexFDRThr = 0.1,
    maxNbrComplexesToPlot = 10, maxComplexSimilarity = 1, seed = 42,
    includeFeatureCollections = c(), minSizeToKeepSet = 2,
    customComplexes = list(),
    complexSpecies = "all", complexDbPath = NULL, stringVersion = "11.5",
    stringDir = NULL, linkTableColumns = c(),
    customYml = NULL, doRender = TRUE
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
    .checkArgumentsPDTMT(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, reportTitle = reportTitle,
        reportAuthor = reportAuthor, forceOverwrite = forceOverwrite,
        experimentInfo = experimentInfo, species = species,
        pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
        inputLevel = inputLevel,
        pdAnalysisFile = pdAnalysisFile, idCol = idCol, labelCol = labelCol,
        geneIdCol = geneIdCol, proteinIdCol = proteinIdCol,
        stringIdCol = stringIdCol, extraFeatureCols = extraFeatureCols,
        modificationsCol = modificationsCol,
        excludeUnmodifiedPeptides = excludeUnmodifiedPeptides,
        keepModifications = keepModifications,
        iColPattern = iColPattern, sampleAnnot = sampleAnnot,
        includeOnlySamples = includeOnlySamples,
        excludeSamples = excludeSamples,
        minScore = minScore, minDeltaScore = minDeltaScore,
        minPeptides = minPeptides, minPSMs = minPSMs,
        masterProteinsOnly = masterProteinsOnly,
        imputeMethod = imputeMethod, assaysForExport = assaysForExport,
        addAbundanceValues = addAbundanceValues, addHeatmaps = addHeatmaps,
        mergeGroups = mergeGroups, comparisons = comparisons,
        ctrlGroup = ctrlGroup, allPairwiseComparisons = allPairwiseComparisons,
        singleFit = singleFit,
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
    pdOutputFolder <- normalizePath(pdOutputFolder)
    pdAnalysisFile <- normalizePath(pdAnalysisFile)

    ## -------------------------------------------------------------------------
    ## Copy Rmd template and insert arguments
    ## -------------------------------------------------------------------------
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentInfo = experimentInfo, species = species,
             pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
             inputLevel = inputLevel, pdAnalysisFile = pdAnalysisFile,
             idCol = idCol, labelCol = labelCol, geneIdCol = geneIdCol,
             proteinIdCol = proteinIdCol, stringIdCol = stringIdCol,
             extraFeatureCols = extraFeatureCols,
             modificationsCol = modificationsCol,
             excludeUnmodifiedPeptides = excludeUnmodifiedPeptides,
             keepModifications = keepModifications,
             reportTitle = reportTitle, reportAuthor = reportAuthor,
             iColPattern = iColPattern, sampleAnnot = sampleAnnot,
             includeOnlySamples = includeOnlySamples,
             excludeSamples = excludeSamples,
             minScore = minScore, minDeltaScore = minDeltaScore,
             minPeptides = minPeptides, minPSMs = minPSMs,
             masterProteinsOnly = masterProteinsOnly,
             imputeMethod = imputeMethod, assaysForExport = assaysForExport,
             addAbundanceValues = addAbundanceValues, addHeatmaps = addHeatmaps,
             mergeGroups = mergeGroups, comparisons = comparisons,
             ctrlGroup = ctrlGroup,
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
             expType = "ProteomeDiscoverer")
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
