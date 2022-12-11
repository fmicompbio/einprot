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
    idCol, labelCol, geneIdCol, proteinIdCol,
    iColPattern, sampleAnnot, includeOnlySamples,
    excludeSamples, minScore, minPeptides, imputeMethod, mergeGroups,
    comparisons, ctrlGroup, allPairwiseComparisons, singleFit,
    subtractBaseline, baselineGroup, normMethod, spikeFeatures, stattest,
    minNbrValidValues, minlFC, samSignificance, nperm, volcanoAdjPvalThr,
    volcanoLog2FCThr, volcanoMaxFeatures, volcanoS0, volcanoFeaturesToLabel,
    addInteractiveVolcanos, complexFDRThr, maxNbrComplexesToPlot, seed,
    includeFeatureCollections, minSizeToKeepSet, customComplexes,
    complexSpecies, complexDbPath, customYml, doRender
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
    .assertScalar(x = doRender, type = "logical")

    ## Experiment info
    .assertVector(x = experimentInfo, type = "list")
    if (length(experimentInfo) > 0) {
        .assertVector(x = names(experimentInfo), type = "character")
    }
    tmp <- getSpeciesInfo(species) ## gives an error for unsupported species

    ## MQ files
    .assertScalar(x = mqFile, type = "character")
    if (!file.exists(mqFile)) {
        stop("'mqFile' must point to an existing file")
    }
    .assertScalar(x = mqParameterFile, type = "character", allowNULL = TRUE)
    if (!is.null(mqParameterFile) && !file.exists(mqParameterFile)) {
        stop("'mqParameterFile' must point to an existing file")
    }

    ## Samples to include or exclude
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    if ((length(includeOnlySamples) > 1 || includeOnlySamples != "") &&
        (length(excludeSamples) > 1 || excludeSamples != "")) {
        stop("Please specify max one of includeOnlySamples and excludeSamples")
    }

    ## Names and patterns
    .assertScalar(x = iColPattern, type = "character",
                  validValues = c("^MS\\\\.MS\\\\.Count\\\\.",
                                  "^LFQ\\\\.intensity\\\\.",
                                  "^Intensity\\\\.",
                                  "^Sequence\\\\.coverage\\\\.",
                                  "^iBAQ\\\\."))
    .assertVector(x = sampleAnnot, type = "data.frame")
    .assertVector(x = colnames(sampleAnnot), type = "character")
    stopifnot(all(c("sample", "group") %in% colnames(sampleAnnot)))
    .assertVector(x = sampleAnnot$group, type = "character")
    ics <- getIntensityColumns(inFile = mqFile,
                               iColPattern = gsub("\\\\", "\\", iColPattern,
                                                  fixed = TRUE),
                               includeOnlySamples = includeOnlySamples,
                               excludeSamples = excludeSamples,
                               stopIfEmpty = TRUE)
    ics <- gsub(gsub("\\\\", "\\", iColPattern,
                     fixed = TRUE), "", ics$iCols)
    msg <- setdiff(ics, sampleAnnot$sample)
    if (length(msg) > 0) {
        stop("Not all sample names are available in the sample annotation. ",
             "Missing samples: ", paste(msg, collapse = ","))
    }

    if (is(idCol, "function")) {
        stopifnot(length(formals(idCol)) == 1)
    } else {
        .assertVector(x = idCol, type = "character")
    }
    if (is(labelCol, "function")) {
        stopifnot(length(formals(labelCol)) == 1)
    } else {
        .assertVector(x = labelCol, type = "character")
    }
    if (is(geneIdCol, "function")) {
        stopifnot(length(formals(geneIdCol)) == 1)
    } else {
        .assertVector(x = geneIdCol, type = "character")
    }
    if (is(proteinIdCol, "function")) {
        stopifnot(length(formals(proteinIdCol)) == 1)
    } else {
        .assertVector(x = proteinIdCol, type = "character")
    }

    ## Score thresholds
    .assertScalar(x = minScore, type = "numeric")
    .assertScalar(x = minPeptides, type = "numeric")

    ## Method choices
    .assertScalar(x = imputeMethod, type = "character",
                  validValues = c("impSeqRob", "MinProb"))
    .assertScalar(x = normMethod, type = "character",
                  validValues = c(MsCoreUtils::normalizeMethods(), "none"))
    .assertVector(x = spikeFeatures, type = "character", allowNULL = TRUE)
    .assertScalar(x = stattest, type = "character",
                  validValues = c("limma", "ttest", "proDA"))

    ## Test parameters
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = samSignificance, type = "logical")
    .assertScalar(x = nperm, type = "numeric", rngIncl = c(1, Inf))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoMaxFeatures, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoS0, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = maxNbrComplexesToPlot, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minSizeToKeepSet, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertVector(x = mergeGroups, type = "list")
    .assertVector(x = comparisons, type = "list")
    .assertScalar(x = ctrlGroup, type = "character")
    .assertScalar(x = allPairwiseComparisons, type = "logical")
    .assertScalar(x = addInteractiveVolcanos, type = "logical")
    .assertScalar(x = singleFit, type = "logical")
    .assertScalar(x = subtractBaseline, type = "logical")
    .assertScalar(x = baselineGroup, type = "character")

    if (length(mergeGroups) > 0) {
        if (is.null(names(mergeGroups)) || any(names(mergeGroups) == "") ||
            any(duplicated(names(mergeGroups)))) {
            stop("'mergeGroups' must be a named list, without duplicated names")
        }
    }

    if (length(comparisons) > 0) {
        if (!all(vapply(comparisons, length, 0) == 2)) {
            stop("Each entry in 'comparisons' must have exactly two elements")
        }
    }

    ## seed
    .assertScalar(x = seed, type = "numeric", rngIncl = c(1, Inf))

    ## Complexes
    .assertVector(x = includeFeatureCollections, type = "character",
                  validValues = c("complexes", "GO"), allowNULL = TRUE)
    .assertVector(x = customComplexes, type = "list")
    if (length(customComplexes) > 0) {
        .assertVector(x = names(customComplexes), type = "character")
    }
    .assertScalar(x = complexSpecies, type = "character",
                  validValues = c("current", "all"), allowNULL = TRUE)
    .assertScalar(x = complexDbPath, type = "character", allowNULL = TRUE)
    if (!is.null(complexDbPath) && !file.exists(complexDbPath)) {
        stop("'complexDbPath' must point to an existing file")
    }

    .assertScalar(x = customYml, type = "character", allowNULL = TRUE)
    if (!is.null(customYml) && !file.exists(customYml)) {
        stop("'customYml' must point to an existing file")
    }
}
