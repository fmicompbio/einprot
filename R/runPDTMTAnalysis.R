#' Run analysis on PD/TMT data
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @param templateRmd Path to the template Rmd. Typically does not need to
#'     be modified.
#' @param outputDir Path to a directory where all output files will be
#'     written. Will be created if it doesn't exist.
#' @param outputBaseName Character string providing the 'base name' of the
#'     output files. All output files will start with this prefix.
#' @param forceOverwrite Logical, whether to force overwrite an existing
#'     Rmd file with the same \code{outputBaseName} in the \code{outputDir}.
#' @param experimentId Numeric, the FMI experiment ID.
#' @param pdOutputFolder Character string pointing to the PD/TMT output folder.
#'     Should contain the files \code{pdResultName_InputFiles.txt},
#'     \code{pdResultName_StudyInformation.txt} and
#'     \code{pdResultName_Proteins.txt}
#' @param pdResultName Character string providing the base name for the
#'     files in the \code{pdOutputFolder}.
#' @param pdAnalysisFile Character string pointing to the \code{pdAnalysis}
#'     file
#' @param analysisDetails Character string, any specific details about the
#'     analysis to mention in the report.
#' @param cysAlkylation Character string.
#' @param sampleIs Character string with sample description (e.g., 'digested')
#' @param enzymes Character string indicating the enzymes that were used.
#' @param aName Character string providing the desired name of the base assay
#'     in the output \code{SingleCellExperiment} object.
#' @param iColPattern Regular expression identifying the columns of the PD
#'     \code{Proteins.txt} file to use for the analysis.
#' @param samplePattern Regular expression identifying the sample pattern, which
#'     will be removed from the sample ID to generate the group name.
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
#' @param ctrlGroup Character vector defining the sample group(s) to use as
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
#'
#' @return Invisibly, the path to the compiled html report.
#'
#' @importFrom xfun Rscript_call
#' @importFrom rmarkdown render
#' @importFrom readr read_file write_file
#'
runPDTMTAnalysis <- function(
    templateRmd = system.file("extdata/process_PD_TMT_template.Rmd",
                              package = "einprot"),
    outputDir = ".", outputBaseName = "PDTMTAnalysis",
    forceOverwrite = FALSE,
    experimentId, pdOutputFolder, pdResultName,
    pdAnalysisFile,
    analysisDetails, cysAlkylation, sampleIs, enzymes,
    aName, iColPattern, samplePattern,
    includeOnlySamples, excludeSamples,
    minScore = 10, minPeptides = 2, imputeMethod = "MinProb",
    mergeGroups = list(), comparisons = list(),
    ctrlGroup = "", allPairwiseComparisons = TRUE,
    normMethod = "none", stattest = "limma", minNbrValidValues = 2,
    minlFC = 0, nperm = 250, volcanoAdjPvalThr = 0.05,
    volcanoLog2FCThr = 1, volcanoMaxFeatures = 25,
    volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
    complexFDRThr = 0.1, seed = 42,
    includeFeatureCollections, customComplexes = list(),
    complexSpecies = "all", complexDbPath
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
    .checkArgumentsPDTMT(
        templateRmd = templateRmd, outputDir = outputDir,
        outputBaseName = outputBaseName, forceOverwrite = forceOverwrite,
        experimentId = experimentId, pdOutputFolder = pdOutputFolder,
        pdResultName = pdResultName, pdAnalysisFile = pdAnalysisFile,
        analysisDetails = analysisDetails,  cysAlkylation = cysAlkylation,
        sampleIs = sampleIs, enzymes = enzymes, aName = aName,
        iColPattern = iColPattern, samplePattern = samplePattern,
        includeOnlySamples = includeOnlySamples,
        excludeSamples = excludeSamples,
        minScore = minScore, minPeptides = minPeptides,
        imputeMethod = imputeMethod, mergeGroups = mergeGroups, comparisons = comparisons,
        ctrlGroup = ctrlGroup, allPairwiseComparisons = allPairwiseComparisons,
        normMethod = normMethod, stattest = stattest,
        minNbrValidValues = minNbrValidValues, minlFC = minlFC,
        nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
        volcanoLog2FCThr = volcanoLog2FCThr,
        volcanoMaxFeatures = volcanoMaxFeatures,
        volcanoS0 = volcanoS0, volcanoFeaturesToLabel = volcanoFeaturesToLabel,
        complexFDRThr = complexFDRThr, seed = seed,
        includeFeatureCollections = includeFeatureCollections,
        customComplexes = customComplexes, complexSpecies = complexSpecies,
        complexDbPath = complexDbPath)

    ## --------------------------------------------------------------------- ##
    ## Copy Rmd template and insert arguments
    ## --------------------------------------------------------------------- ##
    confighook <- "ConfigParameters"

    ## Concatenate Rmd chunk yml
    configchunk <- .generateConfigChunk(
        list(experimentId = experimentId, pdOutputFolder = pdOutputFolder,
             pdResultName = pdResultName, pdAnalysisFile = pdAnalysisFile,
             analysisDetails = analysisDetails,  cysAlkylation = cysAlkylation,
             sampleIs = sampleIs, enzymes = enzymes, aName = aName,
             iColPattern = iColPattern, samplePattern = samplePattern,
             includeOnlySamples = includeOnlySamples,
             excludeSamples = excludeSamples,
             minScore = minScore, minPeptides = minPeptides,
             imputeMethod = imputeMethod, mergeGroups = mergeGroups, comparisons = comparisons,
             ctrlGroup = ctrlGroup, allPairwiseComparisons = allPairwiseComparisons,
             normMethod = normMethod, stattest = stattest,
             minNbrValidValues = minNbrValidValues, minlFC = minlFC,
             nperm = nperm, volcanoAdjPvalThr = volcanoAdjPvalThr,
             volcanoLog2FCThr = volcanoLog2FCThr,
             volcanoMaxFeatures = volcanoMaxFeatures,
             volcanoS0 = volcanoS0, volcanoFeaturesToLabel = volcanoFeaturesToLabel,
             complexFDRThr = complexFDRThr, seed = seed,
             includeFeatureCollections = includeFeatureCollections,
             customComplexes = customComplexes, complexSpecies = complexSpecies,
             complexDbPath = complexDbPath)
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
