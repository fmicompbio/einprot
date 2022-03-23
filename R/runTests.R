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
#' @param testType Character scalar, either "limma" or "ttest".
#' @param assayForTests Character scalar, the name of an assay of the
#'     \code{SummarizedExperiment} object with values that will be used to
#'     perform the test.
#' @param assayImputation Character scalar, the name of an assay of the
#'     \code{SummarizedExperiment} object with logical values indicating
#'     whether an entry was imputed or not.
#' @param minNbrValidValues Numeric scalar, the minimum number of valid
#'     (non-imputed) values that must be present for a features to include it
#'     in the result table.
#' @param minlFC Non-negative numeric scalar, the logFC threshold to use for
#'     limma-treat.
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
#'     used if \code{stattest} is \code{"ttest"}).
#' @param nperm Numeric scalar, the number of permutations (only
#'     used if \code{stattest} is \code{"ttest"}).
#' @param volcanoS0 Numeric scalar, the S0 value to use for creating
#'     significance curves (only used if \code{stattest} is \code{"ttest"}).
#' @param addAbundanceValues Logical scalar, whether to extract abundance
#'     and add to the result table.
#' @param aName Character scalar, the name of the assay in the
#'     \code{SummarizedExperiment} object to get abundance values from (only
#'     required if \code{addAbundanceValues} is \code{TRUE}).
#' @param singleFit Logical scalar, whether to fit a single model to the full
#'     data set and extract relevant results using contrasts. If \code{FALSE},
#'     the data set will be subset for each comparison to only the relevant
#'     samples. Setting \code{singleFit} to \code{TRUE} is only supported
#'     for \code{testType = "limma"}.
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
#' @importFrom limma lmFit treat topTreat cameraPR ids2indices
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>% mutate select contains arrange filter left_join
#'     across everything rename
#' @importFrom rlang .data
#' @importFrom S4Vectors mcols
#' @importFrom genefilter rowttests
#'
runTest <- function(sce, comparisons, testType, assayForTests,
                    assayImputation, minNbrValidValues = 2,
                    minlFC = 0, featureCollections = list(),
                    complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05,
                    volcanoLog2FCThr = 1, baseFileName = NULL, seed = 123,
                    nperm = 250, volcanoS0 = 0.1, addAbundanceValues = FALSE,
                    aName = NULL, singleFit = TRUE) {
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    .assertVector(x = sce, type = "SummarizedExperiment")
    stopifnot("group" %in% colnames(SummarizedExperiment::colData(sce)))
    .assertVector(x = sce$group, type = "character")
    .assertVector(x = comparisons, type = "list")
    for (comparison in comparisons) {
        .assertVector(x = comparison, type = "character", len = 2,
                      validValues = unique(sce$group))
    }
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest"))
    .assertScalar(x = assayForTests, type = "character",
                  validValues = assayNames(sce))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = assayNames(sce))
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    if (testType == "limma") {
        .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    } else if (testType == "ttest") {
        .assertScalar(x = seed, type = "numeric", rngIncl = c(0, Inf))
        .assertScalar(x = nperm, type = "numeric", rngIncl = c(1, Inf))
        .assertScalar(x = volcanoS0, type = "numeric", rngIncl = c(0, Inf))
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

    if (singleFit && testType == "ttest") {
        message("A single model fit is currently not supported for t-tests. ",
                "Changing to fit a separate model for each comparison.")
        singleFit <- FALSE
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
        fc <- factor(structure(sce$group, names = colnames(sce)))
        if ("batch" %in% colnames(SummarizedExperiment::colData(sce))) {
            bc <- structure(sce$batch, names = colnames(sce))
            dfdes <- data.frame(fc = fc, bc = bc)
            if (length(unique(bc)) == 1) {
                messages <- paste0("Only one unique value for batch - ",
                                   "fitting a model without batch.")
                design <- stats::model.matrix(~ fc, data = dfdes)
            } else {
                design <- stats::model.matrix(~ bc + fc, data = dfdes)
            }
        } else {
            dfdes <- data.frame(fc = fc)
            design <- stats::model.matrix(~ fc, data = dfdes)
        }
        returndesign <- list(design = design, sampleData = dfdes,
                             contrasts = list())
        exprvals <- SummarizedExperiment::assay(sce, assayForTests,
                                                withDimnames = TRUE)
        fit0 <- limma::lmFit(exprvals, design)
    }
    for (comparison in comparisons) {
        scesub <- sce[, sce$group %in% comparison]
        ## Only consider features with at least a given number of valid values
        imputedvals <- SummarizedExperiment::assay(scesub, assayImputation,
                                                   withDimnames = TRUE)
        keep <- rowSums(!imputedvals) >= minNbrValidValues

        if (singleFit) {
            fit <- fit0[keep, ]
            contrast <- (colnames(design) == paste0("fc", comparison[2])) -
                (colnames(design) == paste0("fc", comparison[1]))
            returndesign$contrasts[[paste0(comparison[[2]], "_vs_",
                                           comparison[[1]])]] <- contrast
            fit <- limma::contrasts.fit(fit, contrasts = contrast)
            fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
            res <- limma::topTreat(fit, coef = 1,
                                   number = Inf, sort.by = "none") %>%
                tibble::rownames_to_column("pid") %>%
                dplyr::mutate(s2.prior = fit$s2.prior,
                              weights = fit$weights,
                              sigma = fit$sigma)
            camerastat <- "t"
        } else {
            if (testType == "limma") {
                fc <- factor(structure(scesub$group, names = colnames(scesub)),
                             levels = comparison)
                if ("batch" %in% colnames(SummarizedExperiment::colData(scesub))) {
                    bc <- structure(scesub$batch, names = colnames(scesub))
                    dfdes <- data.frame(fc = fc, bc = bc)
                    if (length(unique(bc)) == 1) {
                        messages[[paste0(comparison[[2]], "_vs_",
                                         comparison[[1]])]] <-
                            paste0("Only one unique value for batch - ",
                                   "fitting a model without batch.")
                        design <- stats::model.matrix(~ fc, data = dfdes)
                    } else {
                        design <- stats::model.matrix(~ bc + fc, data = dfdes)
                    }
                } else {
                    dfdes <- data.frame(fc = fc)
                    design <- stats::model.matrix(~ fc, data = dfdes)
                }
                contrast <- (colnames(design) == paste0("fc", comparison[2])) -
                    (colnames(design) == paste0("fc", comparison[1]))
                returndesign[[paste0(comparison[[2]], "_vs_",
                                     comparison[[1]])]] <-
                    list(design = design, sampleData = dfdes, contrast = contrast)
            } else if (testType == "ttest") {
                fc <- factor(scesub$group, levels = rev(comparison))
            }
            exprvals <- SummarizedExperiment::assay(scesub, assayForTests,
                                                    withDimnames = TRUE)
            exprvals <- exprvals[keep, , drop = FALSE]

            ## ------------------------------------------------------------- ##
            ## Run test
            ## ------------------------------------------------------------- ##
            if (testType == "limma") {
                fit <- limma::lmFit(exprvals, design)
                fit <- limma::contrasts.fit(fit, contrasts = contrast)
                fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                res <- limma::topTreat(fit, coef = 1,
                                       number = Inf, sort.by = "none") %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::mutate(s2.prior = fit$s2.prior,
                                  weights = fit$weights,
                                  sigma = fit$sigma)
                camerastat <- "t"
            } else if (testType == "ttest") {
                res <- genefilter::rowttests(exprvals, fac = fc)
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::rename(t = .data$statistic, logFC = .data$dm,
                                  P.Value = .data$p.value) %>%
                    dplyr::mutate(adj.P.Val = p.adjust(.data$P.Value,
                                                       method = "BH"),
                                  AveExpr = rowMeans(exprvals)) %>%
                    dplyr::mutate(sam = .data$t/(1 + .data$t * volcanoS0/.data$logFC))
                camerastat <- "sam"
            }
        }

        res <- res %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(scesub)) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select(.data$pid, .data$geneIdSingle,
                                  .data$proteinIdSingle),
                by = "pid")

        ## ----------------------------------------------------------------- ##
        ## Test feature sets
        ## ----------------------------------------------------------------- ##
        featureCollections <- lapply(featureCollections, function(fcoll) {
            camres <- limma::cameraPR(
                statistic = structure(res[[camerastat]], names = res$pid),
                index = limma::ids2indices(as.list(fcoll), res$pid,
                                           remove.empty = FALSE),
                sort = FALSE
            )

            if (!("FDR" %in% colnames(camres))) {
                camres$FDR <- stats::p.adjust(camres$PValue, method = "BH")
            }
            colnames(camres) <- paste0(comparison[2], "_vs_",
                                       comparison[1], "_", colnames(camres))
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
                              dplyr::contains(paste0(comparison[2], "_vs_",
                                                     comparison[1]))) %>%
                dplyr::arrange(.data[[paste0(comparison[2], "_vs_",
                                             comparison[1], "_FDR")]]) %>%
                dplyr::filter(.data[[paste0(comparison[2], "_vs_",
                                            comparison[1], "_FDR")]] < complexFDRThr)
            topSets[[setname]] <- tmpres
            if (nrow(tmpres) > 0 && !is.null(baseFileName)) {
                write.table(tmpres,
                            file = paste0(baseFileName,
                                          paste0("_testres_", comparison[2],
                                                 "_vs_", comparison[1],
                                                 "_camera_", setname, ".txt")),
                            row.names = FALSE, col.names = TRUE,
                            quote = FALSE, sep = "\t")
            }
        }

        ## ----------------------------------------------------------------- ##
        ## Get the threshold curve (replicating Perseus plots)
        ## ----------------------------------------------------------------- ##
        if (testType == "limma") {
            curveparam <- list()
        } else if (testType == "ttest") {
            curveparam <- .getThresholdCurve(
                exprvals = exprvals, fc = fc, res = res, seed = seed,
                nperm = nperm, volcanoS0 = volcanoS0,
                volcanoAdjPvalThr = volcanoAdjPvalThr)
        }

        ## ----------------------------------------------------------------- ##
        ## Determine significant features
        ## ----------------------------------------------------------------- ##
        res <- data.frame(pid = rownames(imputedvals)) %>%
            dplyr::left_join(res, by = "pid")
        if (testType == "limma") {
            res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
                abs(res$logFC) >= volcanoLog2FCThr
        } else if (testType == "ttest") {
            res$showInVolcano <- abs(res$sam) >= curveparam$ta
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
                        file = paste0(baseFileName, "_testres_", comparison[2],
                                      "_vs_", comparison[1], ".txt"),
                        row.names = FALSE, col.names = TRUE,
                        quote = FALSE, sep = "\t")
        }

        ## ----------------------------------------------------------------- ##
        ## Generate return values
        ## ----------------------------------------------------------------- ##
        if (testType == "limma") {
            plottitle <- paste0(comparison[2], " vs ", comparison[1],
                                ", limma treat (H0: |log2FC| <= ", minlFC, ")")
            plotnote <- paste0("df.prior = ", round(fit$df.prior, digits = 2))
            plotsubtitle <- paste0("Adj.p threshold = ", volcanoAdjPvalThr,
                                   ", |log2FC| threshold = ", volcanoLog2FCThr)
        } else if (testType == "ttest") {
            plottitle <- paste0(comparison[2], " vs ", comparison[1], ", t-test")
            plotnote <- ""
            plotsubtitle = paste0("FDR threshold = ", volcanoAdjPvalThr,
                                  ", s0 = ", curveparam$s0)
        }

        ## ----------------------------------------------------------------- ##
        ## Populate result lists
        ## ----------------------------------------------------------------- ##
        plottitles[[paste0(comparison[2], "_vs_", comparison[1])]] <- plottitle
        plotsubtitles[[paste0(comparison[2], "_vs_", comparison[1])]] <- plotsubtitle
        plotnotes[[paste0(comparison[2], "_vs_", comparison[1])]] <- plotnote
        tests[[paste0(comparison[2], "_vs_", comparison[1])]] <- res
        curveparams[[paste0(comparison[2], "_vs_", comparison[1])]] <- curveparam
        topsets[[paste0(comparison[2], "_vs_", comparison[1])]] <- topSets

    } ## end comparison

    return(list(plottitles = plottitles, plotsubtitles = plotsubtitles,
                plotnotes = plotnotes, tests = tests,
                curveparams = curveparams, topsets = topsets,
                messages = messages, design = returndesign,
                featureCollections = featureCollections))
}
