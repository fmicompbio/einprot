#' Check validity of arguments for PD TMT analysis
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
#'
#' @importFrom MsCoreUtils normalizeMethods
.checkArgumentsPDTMT <- function(
    templateRmd, outputDir, outputBaseName, forceOverwrite,
    experimentId, pdOutputFolder, pdResultName,
    pdAnalysisFile, analysisDetails, cysAlkylation,
    sampleIs, enzymes, aName, iColPattern, samplePattern, includeOnlySamples,
    excludeSamples, minScore, minPeptides, imputeMethod, mergeGroups,
    comparisons, ctrlGroup, allPairwiseComparisons, normMethod, stattest,
    minNbrValidValues, minlFC, nperm, volcanoAdjPvalThr, volcanoLog2FCThr,
    volcanoMaxFeatures, volcanoS0, volcanoFeaturesToLabel,
    addInteractiveVolcanos, complexFDRThr, seed,
    includeFeatureCollections, customComplexes, complexSpecies, complexDbPath
) {
    ## templateRmd
    .assertScalar(templateRmd, type = "character")
    if (!file.exists(templateRmd)) {
        stop("'templateRmd' must point to an existing file")
    }

    ## outputDir
    .assertScalar(outputDir, type = "character")

    ## outputBaseName
    .assertScalar(outputBaseName, type = "character")

    ## forceOverwrite
    .assertScalar(forceOverwrite, type = "logical")

    ## experimentId
    .assertScalar(experimentId, type = "numeric")

    ## PD files
    .assertScalar(pdOutputFolder, type = "character")
    .assertScalar(pdResultName, type = "character")
    if (!file.exists(file.path(pdOutputFolder, paste0(pdResultName, "_Proteins.txt")))) {
        stop("The file ",
             file.path(pdOutputFolder, paste0(pdResultName, "_Proteins.txt")),
             " doesn't exist")
    }
    if (!file.exists(file.path(pdOutputFolder, paste0(pdResultName, "_InputFiles.txt")))) {
        stop("The file ",
             file.path(pdOutputFolder, paste0(pdResultName, "_InputFiles.txt")),
             " doesn't exist")
    }
    if (!file.exists(file.path(pdOutputFolder, paste0(pdResultName, "_StudyInformation.txt")))) {
        stop("The file ",
             file.path(pdOutputFolder, paste0(pdResultName, "_StudyInformation.txt")),
             " doesn't exist")
    }

    ## More files
    .assertScalar(pdAnalysisFile, type = "character")
    if (!file.exists(pdAnalysisFile)) {
        stop("'pdAnalysisFile' must point to an existing file")
    }

    ## Analysis details
    .assertScalar(analysisDetails, type = "character")
    .assertScalar(cysAlkylation, type = "character")
    .assertScalar(sampleIs, type = "character")
    .assertScalar(enzymes, type = "character")

    ## Names and patterns
    .assertScalar(aName, type = "character")
    .assertScalar(iColPattern, type = "character",
                  validValues = c("^Abundance\\\\.F.+\\\\.Sample\\\\."))
    ## "^Abundances\\\\.Count\\\\.F.+\\\\.Sample\\\\."
    .assertScalar(samplePattern, type = "character")

    ## Samples to include or exclude
    .assertVector(includeOnlySamples, type = "character")
    .assertVector(excludeSamples, type = "character")
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }


    ## Score thresholds
    .assertScalar(minScore, type = "numeric")
    .assertScalar(minPeptides, type = "numeric")

    ## Method choices
    .assertScalar(imputeMethod, type = "character",
                  validValues = c("impSeqRob", "MinProb"))
    .assertScalar(normMethod, type = "character",
                  validValues = c(MsCoreUtils::normalizeMethods(), "none"))
    .assertScalar(stattest, type = "character",
                  validValues = c("limma", "ttest"))
    if (stattest == "ttest") {
        stop("'ttest' is currently not supported for PD/TMT data, due to the ",
             "extensive computational time")
    }

    ## Test parameters
    .assertScalar(minNbrValidValues, type = "numeric")
    .assertScalar(minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(nperm, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(volcanoMaxFeatures, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(volcanoS0, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertVector(volcanoFeaturesToLabel, type = "character")
    .assertVector(comparisons, type = "list")
    .assertVector(ctrlGroup, type = "character")
    .assertScalar(allPairwiseComparisons, type = "logical")
    .assertScalar(addInteractiveVolcanos, type = "logical")

    if (length(mergeGroups) > 0) {
        if (is.null(names(mergeGroups)) || any(names(mergeGroups) == "") ||
            any(duplicated(names(mergeGroups)))) {
            stop("'mergeGroups' must be a named list, without duplicated names")
        }

    }

    ## seed
    .assertScalar(seed, type = "numeric", rngIncl = c(1, Inf))

    ## Complexes
    .assertVector(includeFeatureCollections, type = "character",
                  validValues = c("complexes", "GO"))
    .assertVector(customComplexes, type = "list")
    .assertScalar(complexSpecies, type = "character",
                  validValues = c("current", "all"))
    .assertScalar(complexDbPath, type = "character")
    if (!file.exists(complexDbPath)) {
        stop("'complexDbPath' must point to an existing file")
    }

}

