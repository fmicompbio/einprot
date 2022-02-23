#' Run analysis on MaxQuant data
#'
#' @param templateRmd Path to the template Rmd. Typically does not need to
#'     be modified.
#' @param outputDir Path to a directory where all output files will be
#'     written. Will be created if it doesn't exist.
#' @param outputBaseName Character string providing the 'base name' of the
#'     output files. All output files will start with this prefix.
#' @param reportTitle,reportAuthor Character scalars, giving the title and
#'     author for the result report.
#' @param forceOverwrite Logical, whether to force overwrite an existing
#'     Rmd file with the same \code{outputBaseName} in the \code{outputDir}.
#' @param experimentInfo Named list with information about the experiment.
#'     Each entry of the list must be a scalar value.
#' @param species Character scalar providing the species. Must be one of the
#'     supported species (see \code{getSupportedSpecies()}).
#' @param mqFile Character string pointing to the MaxQuant
#'     \code{proteinGroups.txt} file.
#' @param mqParameterFile Character string pointing to the MaxQuant
#'     parameter (xml) file.
#' @param aName Character string providing the desired name of the base assay
#'     in the output \code{SingleCellExperiment} object.
#' @param iColPattern Regular expression identifying the columns of the MaxQuant
#'     \code{proteinGroups.txt} file to use for the analysis. Typically either
#'     "^Intensity\\\\.", "^LFQ\\\\.intensity\\\\." or "^iBAQ\\\\."
#' @param sampleAnnot A \code{data.frame} with at least columns named
#'     \code{sample} and \code{group}, used to explicitly specify the group
#'     assignment for each sample. It can also contain a column named
#'     \code{batch}, in which case this will be used as a covariate in
#'     the \code{limma} tests.
#' @param includeOnlySamples,excludeSamples Character vectors defining specific
#'     samples to include or exclude from all analyses.
#' @param minScore Numeric, minimum score for a protein to be retained in the
#'     analysis.
#' @param minPeptides Numeric, minimum number of peptides for a protein to be
#'     retained in the analysis.
#' @param imputeMethod Character string defining the imputation method to use.
#' @param mergeGroups Named list of character vectors defining sample groups
#'     to merge to create new groups, that will be used for comparisons.
#'     Any specification of \code{comparisons} or \code{ctrlGroup} should
#'     be done in terms of the new (merged) group names.
#' @param comparisons List of character vectors defining comparisons to
#'     perform. The first element of each vector represents the
#'     denominator of the comparison. If not empty, \code{ctrlGroup} and
#'     \code{allPairwiseComparisons} are ignored.
#' @param ctrlGroup Character scalar defining the sample group to use as
#'     control group in comparisons.
#' @param allPairwiseComparisons Logical, should all pairwise comparisons be
#'     performed?
#' @param normMethod Character scalar indicating the normalization method to
#'     use.
#' @param stattest Either \code{"ttest"} or \code{"limma"}, the testing
#'     framework to use.
#' @param minNbrValidValues Numeric, the minimum number of valid values for a
#'     protein to be used for statistical testing.
#' @param minlFC Numeric, minimum log fold change to test against (only used
#'     if \code{stattest = "limma"}).
#' @param nperm Numeric, number of permutations to use in the statistical
#'     testing (only used if \code{stattest = "ttest"}).
#' @param volcanoAdjPvalThr Numeric, adjusted p-value threshold to determine
#'     which proteins to highlight in the volcano plots.
#' @param volcanoLog2FCThr Numeric, log-fold change threshold to determine
#'     which proteins to highlight in the volcano plots.
#' @param volcanoMaxFeatures Numeric, maximum number of significant features
#'     to label in the volcano plots.
#' @param volcanoS0 Numeric, S0 value to use to generate the significance
#'     curve in the volcano plots (only used if \code{stattest = "ttest"}).
#' @param volcanoFeaturesToLabel Character vector with features to always
#'     label in the volcano plots (regardless of significance).
#' @param addInteractiveVolcanos Logical scalar indicating whether to add
#'     interactive volcano plots to the html report. For experiments with
#'     many quantified features or many comparisons, setting this to
#'     \code{TRUE} can make the html report very large and difficult to
#'     interact with.
#' @param complexFDRThr Numeric, FDR threshold for significance in testing
#'     of complexes.
#' @param seed Numeric, random seed to use for any non-deterministic
#'     calculations.
#' @param includeFeatureCollections Character vector, a subset of
#'     c("complexes", "GO").
#' @param customComplexes List of character vectors providing custom complexes
#'     to test for significant differences between groups.
#' @param complexSpecies Either \code{"all"} or \code{"current"}, depending
#'     on whether complexes defined for all species, or only those defined
#'     for the current species, should be tested for significance.
#' @param complexDbPath Character string providing path to the complex DB
#'     file (generated with \code{makeComplexDB()}).
#' @param customYml Character string providing the path to a custom YAML file
#'     that can be used to overwrite default settings in the report. If set
#'     to \code{NULL} (default), no alterations are made.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return Invisibly, the path to the compiled html report.
#'
#' @importFrom xfun Rscript_call
#' @importFrom rmarkdown render
#' @importFrom readr read_file write_file
#'
runMaxQuantAnalysis <- function(
    templateRmd = system.file("extdata/process_MaxQuant_template.Rmd",
                              package = "einprot"),
    outputDir = ".", outputBaseName = "MaxQuantAnalysis",
    reportTitle = "MaxQuant LFQ data processing", reportAuthor = "",
    forceOverwrite = FALSE,
    experimentInfo, species, mqFile, mqParameterFile,
    aName, iColPattern, sampleAnnot,
    includeOnlySamples, excludeSamples,
    minScore = 10, minPeptides = 2, imputeMethod = "MinProb",
    mergeGroups = list(), comparisons = list(),
    ctrlGroup = "", allPairwiseComparisons = TRUE,
    normMethod = "none", stattest = "limma", minNbrValidValues = 2,
    minlFC = 0, nperm = 250, volcanoAdjPvalThr = 0.05,
    volcanoLog2FCThr = 1, volcanoMaxFeatures = 25,
    volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
    addInteractiveVolcanos = FALSE, complexFDRThr = 0.1, seed = 42,
    includeFeatureCollections, customComplexes = list(),
    complexSpecies = "all", complexDbPath, customYml = NULL
) {
    ## --------------------------------------------------------------------- ##
    ## Fix ctrlGroup/mergeGroups
    ## --------------------------------------------------------------------- ##
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

    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .checkArgumentsMaxQuant(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, reportTitle = reportTitle,
        reportAuthor = reportAuthor, forceOverwrite = forceOverwrite,
        experimentInfo = experimentInfo, species = species,
        mqFile = mqFile, mqParameterFile = mqParameterFile,
        aName = aName,
        iColPattern = iColPattern, sampleAnnot = sampleAnnot,
        includeOnlySamples = includeOnlySamples,
        excludeSamples = excludeSamples,
        minScore = minScore, minPeptides = minPeptides,
        imputeMethod = imputeMethod, mergeGroups = mergeGroups,
        comparisons = comparisons, ctrlGroup = ctrlGroup,
        allPairwiseComparisons = allPairwiseComparisons,
        normMethod = normMethod, stattest = stattest,
        minNbrValidValues = minNbrValidValues, minlFC = minlFC,
        nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
        volcanoLog2FCThr = volcanoLog2FCThr,
        volcanoMaxFeatures = volcanoMaxFeatures,
        volcanoS0 = volcanoS0, volcanoFeaturesToLabel = volcanoFeaturesToLabel,
        addInteractiveVolcanos = addInteractiveVolcanos,
        complexFDRThr = complexFDRThr, seed = seed,
        includeFeatureCollections = includeFeatureCollections,
        customComplexes = customComplexes, complexSpecies = complexSpecies,
        complexDbPath = complexDbPath, customYml = customYml)

    ## --------------------------------------------------------------------- ##
    ## Copy Rmd template and insert arguments
    ## --------------------------------------------------------------------- ##
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentInfo = experimentInfo, species = species,
             mqFile = mqFile, mqParameterFile = mqParameterFile,
             reportTitle = reportTitle, reportAuthor = reportAuthor,
             aName = aName,
             iColPattern = iColPattern, sampleAnnot = sampleAnnot,
             includeOnlySamples = includeOnlySamples,
             excludeSamples = excludeSamples,
             minScore = minScore, minPeptides = minPeptides,
             imputeMethod = imputeMethod, mergeGroups = mergeGroups,
             comparisons = comparisons, ctrlGroup = ctrlGroup,
             allPairwiseComparisons = allPairwiseComparisons,
             normMethod = normMethod, stattest = stattest,
             minNbrValidValues = minNbrValidValues, minlFC = minlFC,
             nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
             volcanoLog2FCThr = volcanoLog2FCThr,
             volcanoMaxFeatures = volcanoMaxFeatures,
             volcanoS0 = volcanoS0, volcanoFeaturesToLabel = volcanoFeaturesToLabel,
             addInteractiveVolcanos = addInteractiveVolcanos,
             complexFDRThr = complexFDRThr, seed = seed,
             includeFeatureCollections = includeFeatureCollections,
             customComplexes = customComplexes, complexSpecies = complexSpecies,
             complexDbPath = complexDbPath, customYml = customYml)
    )

    ## Read Rmd
    rmd <- readr::read_file(templateRmd)

    ## Determine where to insert the config chunk
    ## From https://community.rstudio.com/t/how-to-write-r-script-into-rmd-as-functioning-code-chunk/37453/2
    header_regex <- sprintf("\\{\\{%sStart\\}\\}(.*?)\\{\\{%sEnd\\}\\}",
                            confighook,
                            confighook)

    ## Replace hooks with config chunk
    output <- gsub(header_regex, configchunk, rmd)

    ## Similarly, add any custom yaml
    ymlhook <- "YmlParameters"
    header_regex_yml <- sprintf("\\{\\{%sStart\\}\\}(.*?)\\{\\{%sEnd\\}\\}",
                                ymlhook,
                                ymlhook)
    if (!is.null(customYml)) {
        customYml <- paste(readLines(customYml), collapse = "\n")
    } else {
        customYml <- ""
    }
    output <- gsub(header_regex_yml, customYml, output)

    ## Write output to file
    if (!dir.exists(outputDir)) {
        dir.create(outputDir, recursive = TRUE)
    }
    outputFile <- file.path(outputDir, paste0(outputBaseName, ".Rmd"))
    if (file.exists(outputFile) && !forceOverwrite) {
        stop(outputFile, " already exists and forceOverwrite = FALSE, stopping.")
    } else if (file.exists(outputFile) && forceOverwrite) {
        message(outputFile, " already exists but forceOverwrite = TRUE, overwriting.")
    }
    readr::write_file(output, file = outputFile)

    ## --------------------------------------------------------------------- ##
    ## Render the Rmd file
    ## --------------------------------------------------------------------- ##
    args <- list()
    args$input <- outputFile
    args$output_format <- "html_document"
    args$output_dir <- outputDir
    args$intermediates_dir <- outputDir
    args$quiet <- FALSE
    args$run_pandoc <- TRUE

    outputReport <- xfun::Rscript_call(
        rmarkdown::render,
        args
    )

    ## --------------------------------------------------------------------- ##
    ## Return (invisibly) the path to the rendered html file
    ## --------------------------------------------------------------------- ##
    invisible(outputReport)
}
