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
        proteinIdColPeptides, modificationsCol, excludeUnmodifiedPeptides,
        keepModifications, comparisons, ctrlGroup, allPairwiseComparisons,
        singleFit, subtractBaseline, baselineGroup,
        testType, minNbrValidValues, minlFC, volcanoAdjPvalThr, volcanoLog2FCThr,
        volcanoMaxFeatures, volcanoFeaturesToLabel,
        addInteractiveVolcanos, interactiveDisplayColumns, seed, customYml, doRender
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

    .assertScalar(x = assayForTests, type = "character")
    .assertScalar(x = assayImputation, type = "character")

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

    .assertScalar(x = modificationsCol, type = "character")
    .assertScalar(x = excludeUnmodifiedPeptides, type = "logical")
    .assertScalar(x = keepModifications, type = "character", allowNULL = TRUE)

    ## Method choices
    .assertScalar(x = testType, type = "character",
                  validValues = c("interaction", "welch"))

    ## Test parameters
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoMaxFeatures, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertVector(x = comparisons, type = "list")
    .assertScalar(x = ctrlGroup, type = "character")
    .assertScalar(x = allPairwiseComparisons, type = "logical")
    .assertScalar(x = addInteractiveVolcanos, type = "logical")
    .assertVector(x = interactiveDisplayColumns, type = "character", allowNULL = TRUE)
    .assertScalar(x = singleFit, type = "logical")
    .assertScalar(x = subtractBaseline, type = "logical")
    .assertScalar(x = baselineGroup, type = "character")

    if (length(comparisons) > 0) {
        if (!all(vapply(comparisons, length, 0) == 2)) {
            stop("Each entry in 'comparisons' must have exactly two elements")
        }
    }

    ## seed
    .assertScalar(x = seed, type = "numeric", rngIncl = c(1, Inf))

    .assertScalar(x = customYml, type = "character", allowNULL = TRUE)
    if (!is.null(customYml) && !file.exists(customYml)) {
        stop("'customYml' must point to an existing file")
    }

}

