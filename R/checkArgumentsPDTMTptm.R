#' Check validity of arguments for PD TMT PTM analysis
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
#'
#' @importFrom MsCoreUtils normalizeMethods
.checkArgumentsPDTMTptm <- function(
        templateRmd, outputDir, outputBaseName, reportTitle, reportAuthor, forceOverwrite,
        experimentInfo, species, sceProteins, scePeptides, assayForTests,
        assayImputation, proteinIdColProteins,
        proteinIdColPeptides, excludeUnmodifiedPeptides,
        comparisons, ctrlGroup, allPairwiseComparisons, singleFit,
        stattest, minNbrValidValues, minlFC, volcanoAdjPvalThr, volcanoLog2FCThr,
        volcanoMaxFeatures, volcanoFeaturesToLabel,
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

    ## Paths to SCE objects
    .assertScalar(x = sceProteins, type = "character")
    .assertScalar(x = scePeptides, type = "character")
    if (!file.exists(sceProteins)) {
        stop("'sceProteins' must point to an existing file")
    }
    if (!file.exists(scePeptides)) {
        stop("'scePeptides' must point to an existing file")
    }

    ## Names and patterns
    # stopifnot(all(c("sample", "group") %in% colnames(sampleAnnot)))
    if (is(proteinIdColProteins, "function")) {
        stopifnot(length(formals(proteinIdColProteins)) == 1)
    } else {
        .assertVector(x = proteinIdColProteins, type = "character", allowNULL = TRUE)
    }
    if (is(proteinIdColPeptides, "function")) {
        stopifnot(length(formals(proteinIdColPeptides)) == 1)
    } else {
        .assertVector(x = proteinIdColPeptides, type = "character", allowNULL = TRUE)
    }

    ## Method choices
    .assertScalar(x = stattest, type = "character",
                  validValues = c("limma"))

    ## Test parameters
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoMaxFeatures, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = maxNbrComplexesToPlot, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minSizeToKeepSet, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertVector(x = comparisons, type = "list")
    .assertScalar(x = ctrlGroup, type = "character")
    .assertScalar(x = allPairwiseComparisons, type = "logical")
    .assertScalar(x = addInteractiveVolcanos, type = "logical")
    .assertScalar(x = singleFit, type = "logical")

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

