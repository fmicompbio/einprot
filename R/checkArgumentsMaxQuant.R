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
    sampleIs, enzymes, aName, iColPattern, samplePattern, includeOnlySamples,
    excludeSamples, minScore, minPeptides, imputeMethod,
    comparisons, ctrlGroup, allPairwiseComparisons, normMethod, stattest,
    minNbrValidValues, minlFC, nperm, volcanoAdjPvalThr, volcanoLog2FCThr,
    volcanoMaxFeatures, volcanoS0, volcanoFeaturesToLabel, complexFDRThr, seed,
    includeFeatureCollections, customComplexes, complexSpecies, complexDbPath
) {
    ## Check input values
    stopifnot(file.exists(mqFile))
    if (!(imputeMethod %in% c("MinProb", "impSeqRob"))) {
        stop("Misspecified imputeMethod")
    }
    if (!is.logical(allPairwiseComparisons)) {
        stop("Misspecified allPairwiseComparisons")
    }
    if (!(normMethod %in% c("none", MsCoreUtils::normalizeMethods()))) {
        stop("Misspecified normMethod")
    }
    if (!(stattest %in% c("limma", "ttest"))) {
        stop("Misspecified stattest")
    }
    if (!(is.numeric(minlFC) && minlFC >= 0)) {
        stop("Misspecified minlFC")
    }
    if (!(iColPattern %in% c("^MS\\\\.MS\\\\.Count\\\\.",
                             "^LFQ\\\\.intensity\\\\.", "^Intensity\\\\.",
                             "^Sequence\\\\.coverage\\\\.",
                             "^Unique\\\\.peptides\\\\.",
                             "^Razor\\\\.unique\\\\.peptides\\\\.",
                             "^Peptides\\\\.", "^iBAQ\\\\."))) {
        stop("Misspecified iColPattern")
    }
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }
}
