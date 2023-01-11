#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom mclust adjustedRandIndex
#' @importFrom genefilter rowttests
#'
.getThresholdCurve <- function(exprvals, fc, res, seed, nperm, volcanoS0,
                               volcanoAdjPvalThr) {
    ## Permute and calculate SAM statistics
    set.seed(seed)
    sampermL <- lapply(seq_len(nperm), function(i) {
        fc0 <- factor(sample(as.character(fc), size = length(fc)))
        ## Don't include the true grouping
        if (mclust::adjustedRandIndex(fc, fc0) != 1) {
            tts0 <- genefilter::rowttests(as.matrix(exprvals), fc0)
            tts0$statistic/(1 + tts0$statistic * volcanoS0/tts0$dm)
        }
    })
    sampermdf <- abs(do.call(cbind, sampermL))

    ## Define candidate cutoff values - use the observed SAM statistics plus the
    ## middle point between each consecutive pair
    cs <- sort(unique(abs(res$sam)))
    cs <- sort(unique(c(cs, cs + c(diff(cs)/2, 0))))

    ## Get mean number of FPs across permutations for each cutoff
    nbrFP <- sapply(cs, function(cval) {
        apply(sampermdf, 2, function(w) sum(w >= cval))
    })
    meannfp <- colMeans(nbrFP)

    ## Get number of significant hits for each cutoff
    signhits <- sapply(cs, function(cval) {
        sum(abs(res$sam) >= cval)
    })

    ## Get all cutoffs that pass the FDR threshold and extract the smallest one
    passthr <- which(meannfp/signhits < volcanoAdjPvalThr)
    if (length(passthr) > 0) {
        ta <- min(cs[passthr])
        ta2 <- cs[min(passthr)]
        lowerbound <- ceiling(100 * ta * volcanoS0)/100
        x <- seq(from = lowerbound,
                 to = max(lowerbound, 1.1 * max(abs(res$logFC))), by = 0.01)
    } else {
        ta <- Inf
        x <- NA
    }
    df <- sum(table(fc)) - 2
    list(x = x, ta = ta, s0 = volcanoS0, df = df)
}

#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr mutate across everything
#' @importFrom genefilter rowSds
#'
.addAbundanceValues <- function(res, sce, aName) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    stopifnot("pid" %in% colnames(res))
    stopifnot(all(rownames(sce) == res$pid))
    .assertScalar(x = aName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))

    abundance_values <- as.data.frame(SummarizedExperiment::assay(sce, aName))
    colnames(abundance_values) <- paste0(aName, ".",
                                         colnames(abundance_values))
    log_abundance_values <- log2(abundance_values) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    .fns = ~ ifelse(is.finite(.), ., NA)))
    res <- cbind(res, abundance_values)
    for (gr in unique(sce$group)) {
        res[[paste0(aName, ".", gr, ".avg")]] <-
            rowMeans(abundance_values[, sce$group == gr, drop = FALSE],
                     na.rm = TRUE)
        res[[paste0(aName, ".", gr, ".sd")]] <-
            genefilter::rowSds(abundance_values[, sce$group == gr,
                                                drop = FALSE],
                               na.rm = TRUE)
        res[[paste0("log2_", aName, ".", gr, ".avg")]] <-
            rowMeans(log_abundance_values[, sce$group == gr, drop = FALSE],
                     na.rm = TRUE)
        res[[paste0("log2_", aName, ".", gr, ".sd")]] <-
            genefilter::rowSds(log_abundance_values[, sce$group == gr,
                                                    drop = FALSE],
                               na.rm = TRUE)
    }
    return(res)
}

#' Run statistical test
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param comparisons A list of character vectors of length 2, each giving the
#'     two groups to be compared.
#' @param groupComposition A list providing the composition of each group
#'     used in any of the comparisons. If \code{NULL}, assumes that each
#'     group used in \code{comparisons} consists of a single group in the
#'     \code{group} column of \code{colData(sce)}.
#' @param testType Character scalar, either "limma", "ttest" or "proDA".
#' @param assayForTests Character scalar, the name of an assay of the
#'     \code{SummarizedExperiment} object with values that will be used to
#'     perform the test.
#' @param assayImputation Character scalar, the name of an assay of
#'     \code{sce} with logical values indicating whether an entry was imputed
#'     or not.
#' @param minNbrValidValues Numeric scalar, the minimum number of valid
#'     (non-imputed) values that must be present for a features to include it
#'     in the result table.
#' @param minlFC Non-negative numeric scalar, the logFC threshold to use for
#'     limma-treat. If \code{minlFC} = 0, \code{limma::eBayes} is used instead.
#' @param featureCollections List of CharacterLists with feature collections.
#' @param complexFDRThr Numeric scalar giving the significance (FDR) threshold
#'     below which a complex will be considered significant.
#' @param volcanoAdjPvalThr Numeric scalar giving the FDR threshold for
#'     significance (for later use in volcano plots).
#' @param volcanoLog2FCThr Numeric scalar giving the logFC threshold for
#'     significance (for later use in volcano plots).
#' @param baseFileName Character scalar or \code{NULL}, the base file name of
#'     the output text files. If \code{NULL}, no result files are generated.
#' @param seed Numeric scalar, the random seed to use for permutation (only
#'     used if \code{testType} is \code{"ttest"}).
#' @param samSignificance Logical scalar, indicating whether the SAM statistic
#'     should be used to determine significance (similar to the approach used by
#'     Perseus). Only used if \code{testType = "ttest"}. If \code{FALSE}, the
#'     p-values are adjusted using the Benjamini-Hochberg approach and used
#'     to determine significance.
#' @param nperm Numeric scalar, the number of permutations (only
#'     used if \code{testType} is \code{"ttest"}).
#' @param volcanoS0 Numeric scalar, the S0 value to use for creating
#'     significance curves (only used if \code{testType} is \code{"ttest"}).
#' @param addAbundanceValues Logical scalar, whether to extract abundance
#'     and add to the result table.
#' @param aName Character scalar, the name of the assay in the
#'     \code{SummarizedExperiment} object to get abundance values from (only
#'     required if \code{addAbundanceValues} is \code{TRUE}).
#' @param singleFit Logical scalar, whether to fit a single model to the full
#'     data set and extract relevant results using contrasts. If \code{FALSE},
#'     the data set will be subset for each comparison to only the relevant
#'     samples. Setting \code{singleFit} to \code{TRUE} is only supported
#'     for \code{testType = "limma"} or \code{"proDA"}.
#' @param subtractBaseline Logical scalar, whether to subtract the background/
#'     reference value for each feature in each batch before fitting the
#'     model. If \code{TRUE}, requires that a 'batch' column is available.
#' @param baselineGroup Character scalar representing the reference group.
#'     Only used if \code{subtractBaseline} is \code{TRUE}, in which case the
#'     abundance values for a given sample will be adjusted by subtracting the
#'     average value across all samples in the \code{baselineGroup} from the
#'     same batch as the original sample.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with seven components: \code{tests} (a list with test
#' results), \code{plotnotes} (the prior df used by limma), \code{plottitles}
#' (indicating the type of test), \code{plotsubtitles} (indicating the
#' significance thresholds), \code{featureCollections} (list of
#' feature sets, expanded with results from camera), \code{topsets}
#' (a list with the significant feature sets), \code{messages} (any
#' messages for the user), \code{design} (information about the
#' experimental design) and \code{curveparams} (information required to
#' create Perseus-like significance curves). In addition, if
#' \code{baseFileName} is not \code{NULL}, text files with test results
#' (including only features and feature sets passing the imposed significance
#' thresholds) are saved.
#'
#' @importFrom SummarizedExperiment assay rowData colData
#' @importFrom stats model.matrix p.adjust
#' @importFrom limma lmFit treat topTreat cameraPR ids2indices eBayes
#'     topTable
#' @importFrom proDA proDA test_diff
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate select contains arrange filter left_join
#'     across everything rename
#' @importFrom rlang .data
#' @importFrom S4Vectors mcols
#' @importFrom genefilter rowttests
#'
runTest <- function(sce, comparisons, groupComposition = NULL, testType,
                    assayForTests, assayImputation = NULL, minNbrValidValues = 2,
                    minlFC = 0, featureCollections = list(),
                    complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05,
                    volcanoLog2FCThr = 1, baseFileName = NULL, seed = 123,
                    samSignificance = TRUE, nperm = 250, volcanoS0 = 0.1,
                    addAbundanceValues = FALSE, aName = NULL, singleFit = TRUE,
                    subtractBaseline = FALSE, baselineGroup = "") {
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    .assertVector(x = sce, type = "SummarizedExperiment")
    stopifnot("group" %in% colnames(SummarizedExperiment::colData(sce)))
    .assertVector(x = sce$group, type = "character")
    .assertVector(x = comparisons, type = "list")
    .assertVector(x = groupComposition, type = "list", allowNULL = TRUE)
    for (comparison in comparisons) {
        .assertVector(x = comparison, type = "character", len = 2,
                      validValues = unique(c(sce$group, names(groupComposition))))
    }
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest", "proDA"))
    .assertScalar(x = assayForTests, type = "character",
                  validValues = assayNames(sce))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = assayNames(sce), allowNULL = TRUE)
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    if (testType == "limma") {
        .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    } else if (testType == "ttest") {
        .assertScalar(x = seed, type = "numeric", rngIncl = c(0, Inf))
        .assertScalar(x = nperm, type = "numeric", rngIncl = c(1, Inf))
        .assertScalar(x = volcanoS0, type = "numeric", rngIncl = c(0, Inf))
        .assertScalar(x = samSignificance, type = "logical")
    }
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertScalar(x = addAbundanceValues, type = "logical")
    if (addAbundanceValues) {
        .assertScalar(x = aName, type = "character")
    }
    .assertScalar(x = singleFit, type = "logical")
    .assertScalar(x = subtractBaseline, type = "logical")
    .assertScalar(x = baselineGroup, type = "character")
    if (subtractBaseline) {
        stopifnot(baselineGroup %in% sce$group)
        stopifnot("batch" %in% colnames(SummarizedExperiment::colData(sce)))
    }

    if (singleFit && testType == "ttest") {
        message("A single model fit is currently not supported for t-tests. ",
                "Changing to fit a separate model for each comparison.")
        singleFit <- FALSE
    }

    ## If comparisons is not named, add names
    comparisons <- .assignNamesToComparisons(comparisons)

    ## If there are entries in unlist(comparisons) that are not defined in
    ## groupComposition, add them to the latter
    stdf <- setdiff(unlist(comparisons), names(groupComposition))
    groupComposition <- c(groupComposition,
                          setNames(as.list(stdf), stdf))

    ## Check that all entries in groupComposition[comparisons] are in the
    ## group column
    stdf <- setdiff(unlist(groupComposition[unlist(comparisons)]),
                    sce$group)
    if (length(stdf) > 0) {
        stop("Missing group(s) in sce$groups: ", paste(stdf, collapse = ", "))
    }

    ## --------------------------------------------------------------------- ##
    ## Initialize result lists
    ## --------------------------------------------------------------------- ##
    plottitles <- list()
    plotsubtitles <- list()
    plotnotes <- list()
    tests <- list()
    curveparams <- list()
    topsets <- list()
    messages <- list()
    returndesign <- list()

    ## --------------------------------------------------------------------- ##
    ## Subset and define design
    ## --------------------------------------------------------------------- ##
    if (singleFit) {
        if ("batch" %in% colnames(SummarizedExperiment::colData(sce))) {
            if (subtractBaseline) {
                idxkeep <- which(sce$group != baselineGroup)
                exprvals <- getMatSubtractedBaseline(
                    sce[, idxkeep, drop = FALSE],
                    assayName = assayForTests,
                    baselineGroup = baselineGroup,
                    sceFull = sce)
                fc <- factor(structure(sce$group[idxkeep],
                                       names = colnames(sce)[idxkeep]))
                dfdes <- data.frame(fc = fc)
                design <- stats::model.matrix(~ fc, data = dfdes)
            } else {
                exprvals <- SummarizedExperiment::assay(sce, assayForTests,
                                                        withDimnames = TRUE)
                fc <- factor(structure(sce$group, names = colnames(sce)))
                bc <- structure(sce$batch, names = colnames(sce))
                dfdes <- data.frame(fc = fc, bc = bc)
                if (length(unique(bc)) == 1) {
                    messages <- paste0("Only one unique value for batch - ",
                                       "fitting a model without batch.")
                    design <- stats::model.matrix(~ fc, data = dfdes)
                } else {
                    design <- stats::model.matrix(~ bc + fc, data = dfdes)
                }
            }
        } else {
            exprvals <- SummarizedExperiment::assay(sce, assayForTests,
                                                    withDimnames = TRUE)
            fc <- factor(structure(sce$group, names = colnames(sce)))
            dfdes <- data.frame(fc = fc)
            design <- stats::model.matrix(~ fc, data = dfdes)
        }
        returndesign <- list(design = design, sampleData = dfdes,
                             contrasts = list())
        if (testType == "limma") {
            fit0 <- limma::lmFit(exprvals, design)
        } else if (testType == "proDA") {
            fit0 <- proDA::proDA(exprvals, design = design)
        }
    }
    for (comparisonName in names(comparisons)) {
        comparison <- comparisons[[comparisonName]]
        scesub <- sce[, sce$group %in% unlist(groupComposition[comparison])]
        ## Only consider features with at least a given number of valid values
        if (!is.null(assayImputation)) {
            imputedvals <- SummarizedExperiment::assay(scesub, assayImputation,
                                                       withDimnames = TRUE)
            keep <- rowSums(!imputedvals) >= minNbrValidValues
        } else {
            ## Need an object with the full set of rownames, to use for
            ## creating the final data frame
            imputedvals <- SummarizedExperiment::assay(scesub, assayForTests,
                                                       withDimnames = TRUE)
            keep <- rep(TRUE, nrow(scesub))
        }

        if (singleFit) {
            fit <- fit0[keep, ]
            contrast <-
                (1/length(groupComposition[[comparison[2]]])) *
                (colnames(design) %in%
                     paste0("fc", groupComposition[[comparison[2]]])) -
                (1/length(groupComposition[[comparison[1]]])) *
                (colnames(design) %in%
                     paste0("fc", groupComposition[[comparison[1]]]))
            returndesign$contrasts[[comparisonName]] <- contrast
            if (testType == "limma") {
                fit <- limma::contrasts.fit(fit, contrasts = contrast)
                if (minlFC == 0) {
                    fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                    res <- limma::topTable(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                } else {
                    fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                    res <- limma::topTreat(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                }
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::mutate(s2.prior = fit$s2.prior,
                                  weights = fit$weights,
                                  sigma = fit$sigma,
                                  se.logFC = sqrt(fit$s2.post) *
                                      fit$stdev.unscaled[, 1],
                                  df.total = fit$df.total)
                camerastat <- "t"
            } else if (testType == "proDA") {
                res <- proDA::test_diff(fit, contrast = contrast) %>%
                    dplyr::rename(pid = "name",
                                  t = "t_statistic",
                                  adj.P.Val = "adj_pval",
                                  logFC = "diff",
                                  P.Value = "pval")
                camerastat <- "t"
            }
        } else {
            ## Create a new vector with the "merged" group names
            ## Check that it has length 2
            ## Also check that no column is in both groups
            c1 <- which(scesub$group %in% groupComposition[[comparison[1]]])
            c2 <- which(scesub$group %in% groupComposition[[comparison[2]]])
            if (length(intersect(c1, c2)) > 0) {
                stop("The same original group is part of both groups ",
                     "to be compared")
            }
            if (!all(sort(union(c1, c2)) == seq_len(ncol(scesub)))) {
                stop("Subsetting error - not all samples seem to have ",
                     "an assigned group")
            }
            fc <- rep(NA_character_, ncol(scesub))
            fc[c1] <- comparison[1]
            fc[c2] <- comparison[2]

            if (testType %in% c("limma", "proDA")) {
                if ("batch" %in% colnames(SummarizedExperiment::colData(scesub))) {
                    fc <- factor(structure(fc, names = colnames(scesub)),
                                 levels = comparison)
                    bc <- structure(scesub$batch, names = colnames(scesub))
                    dfdes <- data.frame(fc = fc, bc = bc)
                    if (length(unique(bc)) == 1) {
                        messages[[comparisonName]] <-
                            paste0("Only one unique value for batch - ",
                                   "fitting a model without batch.")
                        design <- stats::model.matrix(~ fc, data = dfdes)
                    } else {
                        design <- stats::model.matrix(~ bc + fc, data = dfdes)
                    }
                } else {
                    fc <- factor(structure(fc, names = colnames(scesub)),
                                 levels = comparison)
                    dfdes <- data.frame(fc = fc)
                    design <- stats::model.matrix(~ fc, data = dfdes)
                }
                contrast <- (colnames(design) == paste0("fc", comparison[2])) -
                    (colnames(design) == paste0("fc", comparison[1]))
                returndesign[[comparisonName]] <-
                    list(design = design, sampleData = dfdes, contrast = contrast)
            } else if (testType == "ttest") {
                fc <- factor(fc, levels = rev(comparison))
            }
            if (subtractBaseline) {
                exprvals <- getMatSubtractedBaseline(scesub,
                                                     assayName = assayForTests,
                                                     baselineGroup = baselineGroup,
                                                     sceFull = sce)
            } else {
                exprvals <- SummarizedExperiment::assay(scesub, assayForTests,
                                                        withDimnames = TRUE)
            }
            exprvals <- exprvals[keep, , drop = FALSE]

            ## ------------------------------------------------------------- ##
            ## Run test
            ## ------------------------------------------------------------- ##
            if (testType == "limma") {
                fit <- limma::lmFit(exprvals, design)
                fit <- limma::contrasts.fit(fit, contrasts = contrast)
                if (minlFC == 0) {
                    fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                    res <- limma::topTable(fit, coef = 1, number = Inf,
                                           confint = TRUE, sort.by = "none")
                } else {
                    fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                    res <- limma::topTreat(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                }
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::mutate(s2.prior = fit$s2.prior,
                                  weights = fit$weights,
                                  sigma = fit$sigma,
                                  se.logFC = sqrt(fit$s2.post) *
                                      fit$stdev.unscaled[, 1],
                                  df.total = fit$df.total)
                camerastat <- "t"
            } else if (testType == "proDA") {
                fit <- proDA::proDA(exprvals, design = design)
                res <- proDA::test_diff(fit, contrast = contrast) %>%
                    dplyr::rename(pid = "name",
                                  t = "t_statistic",
                                  adj.P.Val = "adj_pval",
                                  logFC = "diff",
                                  P.Value = "pval")
                camerastat <- "t"
            } else if (testType == "ttest") {
                res <- genefilter::rowttests(exprvals, fac = fc)
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::rename(t = "statistic", logFC = "dm",
                                  P.Value = "p.value") %>%
                    dplyr::mutate(adj.P.Val = p.adjust(.data$P.Value,
                                                       method = "BH"),
                                  AveExpr = rowMeans(exprvals)) %>%
                    dplyr::mutate(sam = .data$t/(1 + .data$t * volcanoS0/.data$logFC))
                if (samSignificance) {
                    camerastat <- "sam"
                } else {
                    camerastat <- "t"
                }
            }
        }

        ## Calculate -log10(p). For features with p-value = 0, use half the
        ## smallest non-zero p-value as a proxy to be able to make volcano plots.
        res <- res %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::mutate(mlog10p = replace(
                .data$mlog10p, .data$P.Value == 0,
                -log10(min(.data$P.Value[which(.data$P.Value > 0)])/2))) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(scesub)) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select(tidyselect::any_of(
                        c("pid", "einprotGene", "einprotProtein",
                          "einprotLabel"))),
                by = "pid")

        ## ----------------------------------------------------------------- ##
        ## Test feature sets
        ## ----------------------------------------------------------------- ##
        featureCollections <- lapply(featureCollections, function(fcoll) {
            notna <- which(!is.na(res[[camerastat]]))
            camres <- limma::cameraPR(
                statistic = structure(res[[camerastat]][notna],
                                      names = res$pid[notna]),
                index = limma::ids2indices(as.list(fcoll), res$pid[notna],
                                           remove.empty = FALSE),
                sort = FALSE
            )

            if (!("FDR" %in% colnames(camres))) {
                camres$FDR <- stats::p.adjust(camres$PValue, method = "BH")
            }
            colnames(camres) <- paste0(comparisonName, "_", colnames(camres))
            stopifnot(all(rownames(camres) == names(fcoll)))
            S4Vectors::mcols(fcoll) <- cbind(S4Vectors::mcols(fcoll), camres)
            fcoll
        })

        ## Write test results for feature collections to text files
        topSets <- list()
        for (setname in names(featureCollections)) {
            tmpres <- as.data.frame(S4Vectors::mcols(featureCollections[[setname]]),
                                    optional = TRUE) %>%
                tibble::rownames_to_column("set") %>%
                dplyr::select(dplyr::any_of(c("set", "genes", "sharedGenes",
                                              "Source", "All.names", "PMID")),
                              dplyr::contains(comparisonName)) %>%
                dplyr::arrange(.data[[paste0(comparisonName, "_FDR")]]) %>%
                dplyr::filter(.data[[paste0(comparisonName, "_FDR")]] < complexFDRThr)
            topSets[[setname]] <- tmpres
            if (nrow(tmpres) > 0 && !is.null(baseFileName)) {
                write.table(tmpres,
                            file = paste0(baseFileName,
                                          paste0("_testres_", comparisonName,
                                                 "_camera_", setname, ".txt")),
                            row.names = FALSE, col.names = TRUE,
                            quote = FALSE, sep = "\t")
            }
        }

        ## ----------------------------------------------------------------- ##
        ## Get the threshold curve (replicating Perseus plots)
        ## ----------------------------------------------------------------- ##
        if (testType %in% c("limma", "proDA")) {
            curveparam <- list()
        } else if (testType == "ttest") {
            if (samSignificance) {
                curveparam <- .getThresholdCurve(
                    exprvals = exprvals, fc = fc, res = res, seed = seed,
                    nperm = nperm, volcanoS0 = volcanoS0,
                    volcanoAdjPvalThr = volcanoAdjPvalThr)
            } else {
                curveparam <- list()
            }
        }

        ## ----------------------------------------------------------------- ##
        ## Determine significant features
        ## ----------------------------------------------------------------- ##
        res <- data.frame(pid = rownames(imputedvals)) %>%
            dplyr::left_join(res, by = "pid")
        if (testType %in% c("limma", "proDA")) {
            res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
                abs(res$logFC) >= volcanoLog2FCThr
        } else if (testType == "ttest") {
            if (samSignificance) {
                res$showInVolcano <- abs(res$sam) >= curveparam$ta
            } else {
                res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
                    abs(res$logFC) >= volcanoLog2FCThr
            }
        }

        ## ----------------------------------------------------------------- ##
        ## Add abundance values and STRING IDs
        ## ----------------------------------------------------------------- ##
        if (addAbundanceValues) {
            res <- .addAbundanceValues(res = res, sce = scesub, aName = aName)
        }

        if ("IDsForSTRING" %in%
            colnames(SummarizedExperiment::rowData(scesub))) {
            res$IDsForSTRING <- SummarizedExperiment::rowData(
                scesub)$IDsForSTRING[match(res$pid, rownames(scesub))]
        }

        ## ----------------------------------------------------------------- ##
        ## Write results to file
        ## ----------------------------------------------------------------- ##
        if (!is.null(baseFileName)) {
            write.table(res %>%
                            dplyr::filter(.data$showInVolcano) %>%
                            dplyr::arrange(desc(.data$logFC)),
                        file = paste0(baseFileName, "_testres_", comparisonName,
                                      ".txt"),
                        row.names = FALSE, col.names = TRUE,
                        quote = FALSE, sep = "\t")
        }

        ## ----------------------------------------------------------------- ##
        ## Generate return values
        ## ----------------------------------------------------------------- ##
        if (testType == "limma") {
            if (minlFC == 0) {
                plottitle <- paste0(sub("_vs_", " vs ", comparisonName),
                                    ", limma")
            } else {
                plottitle <- paste0(sub("_vs_", " vs ", comparisonName),
                                    ", limma treat (H0: |log2FC| <= ", minlFC, ")")
            }
            plotnote <- paste0("df.prior = ", round(fit$df.prior, digits = 2))
            plotsubtitle <- paste0("Adj.p threshold = ", volcanoAdjPvalThr,
                                   ", |log2FC| threshold = ", volcanoLog2FCThr)
        } else if (testType == "proDA") {
            plottitle <- paste0(sub("_vs_", " vs ", comparisonName), ", proDA")
            plotnote <- ""
            plotsubtitle <- paste0("Adj.p threshold = ", volcanoAdjPvalThr,
                                   ", |log2FC| threshold = ", volcanoLog2FCThr)
        } else if (testType == "ttest") {
            plottitle <- paste0(sub("_vs_", " vs ", comparisonName), ", t-test")
            plotnote <- ""
            if (samSignificance) {
                plotsubtitle <- paste0("FDR threshold = ", volcanoAdjPvalThr,
                                       ", s0 = ", curveparam$s0)
            } else {
                plotsubtitle <- paste0("FDR threshold = ", volcanoAdjPvalThr)
            }
        }

        ## ----------------------------------------------------------------- ##
        ## Populate result lists
        ## ----------------------------------------------------------------- ##
        plottitles[[comparisonName]] <- plottitle
        plotsubtitles[[comparisonName]] <- plotsubtitle
        plotnotes[[comparisonName]] <- plotnote
        tests[[comparisonName]] <- res
        curveparams[[comparisonName]] <- curveparam
        topsets[[comparisonName]] <- topSets

    } ## end comparison

    return(list(plottitles = plottitles, plotsubtitles = plotsubtitles,
                plotnotes = plotnotes, tests = tests,
                curveparams = curveparams, topsets = topsets,
                messages = messages, design = returndesign,
                featureCollections = featureCollections))
}
