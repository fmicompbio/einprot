#' Check validity of arguments for MaxQuant analysis
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
#'
#' @importFrom MsCoreUtils normalizeMethods
.checkArgumentsMaxQuant <- function(
    templateRmd, outputDir, outputBaseName, forceOverwrite,
    experimentId, mqFile, mqParameterFile, analysisDetails, cysAlkylation,
    sampleIs, enzymes, aName, iColPattern, samplePattern, sampleAnnot, includeOnlySamples,
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

    ## mqFile
    .assertScalar(mqFile, type = "character")
    if (!file.exists(mqFile)) {
        stop("'mqFile' must point to an existing file")
    }

    ## mqParameterFile
    .assertScalar(mqParameterFile, type = "character")
    if (!file.exists(mqParameterFile)) {
        stop("'mqParameterFile' must point to an existing file")
    }

    ## Analysis details
    .assertScalar(analysisDetails, type = "character")
    .assertScalar(cysAlkylation, type = "character")
    .assertScalar(sampleIs, type = "character")
    .assertScalar(enzymes, type = "character")

    ## Names and patterns
    .assertScalar(aName, type = "character")
    .assertScalar(iColPattern, type = "character",
                  validValues = c("^MS\\\\.MS\\\\.Count\\\\.",
                                  "^LFQ\\\\.intensity\\\\.", "^Intensity\\\\.",
                                  "^Sequence\\\\.coverage\\\\.",
                                  "^Unique\\\\.peptides\\\\.",
                                  "^Razor\\\\.unique\\\\.peptides\\\\.",
                                  "^Peptides\\\\.", "^iBAQ\\\\."))
    .assertScalar(samplePattern, type = "character")
    .assertVector(x = sampleAnnot, type = "data.frame", allowNULL = TRUE)

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
    .assertVector(mergeGroups, type = "list")
    .assertVector(comparisons, type = "list")
    .assertScalar(ctrlGroup, type = "character")
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
