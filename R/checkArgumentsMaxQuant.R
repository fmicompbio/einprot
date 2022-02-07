#' Check validity of arguments for MaxQuant analysis
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
#'
#' @importFrom MsCoreUtils normalizeMethods
.checkArgumentsMaxQuant <- function(
    templateRmd, outputDir, outputBaseName, reportTitle, reportAuthor, forceOverwrite,
    experimentInfo, species, mqFile, mqParameterFile,
    aName, iColPattern, samplePattern, includeOnlySamples,
    excludeSamples, minScore, minPeptides, imputeMethod, mergeGroups,
    comparisons, ctrlGroup, allPairwiseComparisons, normMethod, stattest,
    minNbrValidValues, minlFC, nperm, volcanoAdjPvalThr, volcanoLog2FCThr,
    volcanoMaxFeatures, volcanoS0, volcanoFeaturesToLabel,
    addInteractiveVolcanos, complexFDRThr, seed,
    includeFeatureCollections, customComplexes, complexSpecies, complexDbPath
) {
    ## templateRmd
    .assertScalar(x = templateRmd, type = "character")
    if (!file.exists(templateRmd)) {
        stop("'templateRmd' must point to an existing file")
    }

    ## Output specifications
    .assertScalar(x = outputDir, type = "character")
    .assertScalar(x = outputBaseName, type = "character")
    .assertScalar(x = reportTitle, type = "character")
    .assertScalar(x = reportAuthor, type = "character")
    .assertScalar(x = forceOverwrite, type = "logical")

    ## Experiment info
    .assertVector(x = experimentInfo, type = "list")
    .assertScalar(x = species, type = "character")
    tmp <- getSpeciesInfo(species) ## gives an error for unsupported species

    ## MQ files
    .assertScalar(x = mqFile, type = "character")
    if (!file.exists(mqFile)) {
        stop("'mqFile' must point to an existing file")
    }
    .assertScalar(x = mqParameterFile, type = "character")
    if (!file.exists(mqParameterFile)) {
        stop("'mqParameterFile' must point to an existing file")
    }

    ## Names and patterns
    .assertScalar(x = aName, type = "character")
    .assertScalar(x = iColPattern, type = "character",
                  validValues = c("^MS\\\\.MS\\\\.Count\\\\.",
                                  "^LFQ\\\\.intensity\\\\.",
                                  "^Intensity\\\\.",
                                  "^Sequence\\\\.coverage\\\\.",
                                  "^iBAQ\\\\."))
    .assertScalar(x = samplePattern, type = "character")

    ## Samples to include or exclude
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }


    ## Score thresholds
    .assertScalar(x = minScore, type = "numeric")
    .assertScalar(x = minPeptides, type = "numeric")

    ## Method choices
    .assertScalar(x = imputeMethod, type = "character",
                  validValues = c("impSeqRob", "MinProb"))
    .assertScalar(x = normMethod, type = "character",
                  validValues = c(MsCoreUtils::normalizeMethods(), "none"))
    .assertScalar(x = stattest, type = "character",
                  validValues = c("limma", "ttest"))

    ## Test parameters
    .assertScalar(x = minNbrValidValues, type = "numeric")
    .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = nperm, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoMaxFeatures, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoS0, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertVector(x = mergeGroups, type = "list")
    .assertVector(x = comparisons, type = "list")
    .assertScalar(x = ctrlGroup, type = "character")
    .assertScalar(x = allPairwiseComparisons, type = "logical")
    .assertScalar(x = addInteractiveVolcanos, type = "logical")

    if (length(mergeGroups) > 0) {
        if (is.null(names(mergeGroups)) || any(names(mergeGroups) == "") ||
            any(duplicated(names(mergeGroups)))) {
            stop("'mergeGroups' must be a named list, without duplicated names")
        }
    }

    ## seed
    .assertScalar(x = seed, type = "numeric", rngIncl = c(1, Inf))

    ## Complexes
    .assertVector(x = includeFeatureCollections, type = "character",
                  validValues = c("complexes", "GO"))
    .assertVector(x = customComplexes, type = "list")
    .assertScalar(x = complexSpecies, type = "character",
                  validValues = c("current", "all"))
    .assertScalar(x = complexDbPath, type = "character")
    if (!file.exists(complexDbPath)) {
        stop("'complexDbPath' must point to an existing file")
    }
}
