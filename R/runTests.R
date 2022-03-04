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
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom dplyr mutate across everything
#' @importFrom genefilter rowSds
#'
.addAbundanceValues <- function(res, qftsub, aName, iColPattern) {
    stopifnot(all(rownames(SummarizedExperiment::rowData(qftsub[[aName]])) == res$pid))

    ## If there are "iBAQ" columns in the row data, use those (it means that
    ## iBAQ is not the main assay, but if the values are available they should
    ## still be used)
    if (all(sub(iColPattern, "iBAQ.", colnames(qftsub[[aName]])) %in%
            colnames(rowData(qftsub[[aName]])))) {
        abundance_values <- as.data.frame(
            SummarizedExperiment::rowData(
                qftsub[[aName]])[, sub(iColPattern, "iBAQ.",
                                       colnames(qftsub[[aName]])), drop = FALSE]
        )
        colpat <- "iBAQ"
        finalName <- "iBAQ"
    } else {
        ## Otherwise, use the aName assay
        abundance_values <- as.data.frame(SummarizedExperiment::assay(qftsub[[aName]]))
        abundance_values[is.na(abundance_values)] <- 0
        colpat <- gsub("\\", "", sub("\\^", "", sub("\\.$", "", iColPattern)), fixed = TRUE)
        finalName <- aName
        colnames(abundance_values) <- gsub(iColPattern, paste0(finalName, "."),
                                           colnames(abundance_values))
    }

    log_abundance_values <- log2(abundance_values) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    .fns = ~ ifelse(is.finite(.), ., NA)))
    res <- cbind(res, abundance_values)
    for (gr in unique(SummarizedExperiment::colData(qftsub)$group)) {
        res[[paste0(finalName, ".", gr, ".avg")]] <-
            rowMeans(
                abundance_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                            drop = FALSE])
        res[[paste0(finalName, ".", gr, ".sd")]] <-
            genefilter::rowSds(
                abundance_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                            drop = FALSE])
        res[[paste0("log2_", finalName, ".", gr, ".avg")]] <-
            rowMeans(
                log_abundance_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                                drop = FALSE],
                na.rm = TRUE)
        res[[paste0("log2_", finalName, ".", gr, ".sd")]] <-
            genefilter::rowSds(
                log_abundance_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                                drop = FALSE],
                na.rm = TRUE)
    }
    return(list(res = res, colpat = finalName))
}

#' Run statistical test
#'
#' @param qft A \code{QFeatures} object.
#' @param comparison A character vector of length 2, giving the two groups to
#'     be compared.
#' @param testType Character scalar, either "limma" or "ttest".
#' @param assayForTests Character scalar, the name of an assay of the
#'     \code{QFeatures} object with values that will be used to perform the
#'     test.
#' @param assayImputation Character scalar, the name of an assay of the
#'     \code{QFeatures} object with logical values indicating whether an entry
#'     was imputed or not.
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
#'     significance curves.
#' @param addAbundanceValues Logical scalar, whether to extract abundance
#'     and add to the result table. For \code{MaxQuant}, iBAQ values will be
#'     used if they are available, otherwise the abundance values in the
#'     \code{aName} assay will be used. For other types of input, the
#'     abundance values in the \code{aName} assay will be used.
#' @param iColPattern Character scalar, a regular expression used to identify
#'     intensity columns in the input files (only required if
#'     \code{addAbundanceValues} is \code{TRUE}).
#' @param aName Character scalar, the name of the base assay in the
#'     \code{QFeatures} object (only required if \code{addAbundanceValues} is
#'     \code{TRUE}).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with six components: \code{res} (a data.frame with test
#' results), \code{plotnote} (the prior df used by limma), \code{plottitle}
#' (indicating the type of test), \code{plotsubtitle} (indicating the
#' significance thresholds), \code{featureCollections} (list of
#' feature sets, expanded with results from camera), \code{curveparam}
#' (information required to create Perseus-like significance curves). In
#' addition, if \code{baseFileName} is not \code{NULL}, text files with
#' test results (including only features and feature sets passing the
#' imposed significance thresholds) are saved.
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
runTest <- function(qft, comparison, testType, assayForTests,
                    assayImputation, minNbrValidValues = 2,
                    minlFC = 0, featureCollections = list(),
                    complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05,
                    volcanoLog2FCThr = 1, baseFileName = NULL, seed = 123,
                    nperm = 250, volcanoS0 = 0.1, addAbundanceValues = FALSE,
                    iColPattern = NULL, aName = NULL) {
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    .assertVector(x = qft, type = "QFeatures")
    stopifnot("group" %in% colnames(SummarizedExperiment::colData(qft)))
    .assertVector(x = qft$group, type = "character")
    .assertVector(x = comparison, type = "character", len = 2,
                  validValues = unique(qft$group))
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest"))
    .assertScalar(x = assayForTests, type = "character",
                  validValues = names(qft))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = names(qft))
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    if (testType == "limma") {
        .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    } else if (testType == "ttest") {
        .assertScalar(x = nperm, type = "numeric", rngIncl = c(1, Inf))
        .assertScalar(x = volcanoS0, type = "numeric", rngIncl = c(0, Inf))
    }
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertScalar(x = seed, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = addAbundanceValues, type = "logical")
    if (addAbundanceValues) {
        .assertScalar(x = iColPattern, type = "character")
        .assertScalar(x = aName, type = "character")
    }

    ## --------------------------------------------------------------------- ##
    ## Subset and define design
    ## --------------------------------------------------------------------- ##
    qftsub <- qft[, qft$group %in% comparison]
    if (testType == "limma") {
        fc <- factor(qftsub$group, levels = comparison)
        if ("batch" %in% colnames(SummarizedExperiment::colData(qftsub))) {
            bc <- qftsub$batch
            design <- stats::model.matrix(~ bc + fc)
        } else {
            design <- stats::model.matrix(~ fc)
        }
    } else if (testType == "ttest") {
        fc <- factor(qftsub$group, levels = rev(comparison))
    }

    exprvals <- SummarizedExperiment::assay(qftsub[[assayForTests]],
                                            withDimnames = TRUE)
    ## Only consider features with at least a given number of valid values
    imputedvals <- SummarizedExperiment::assay(qftsub[[assayImputation]],
                                               withDimnames = TRUE)
    keep <- rowSums(!imputedvals) >= minNbrValidValues
    exprvals <- exprvals[keep, , drop = FALSE]

    ## --------------------------------------------------------------------- ##
    ## Run test
    ## --------------------------------------------------------------------- ##
    if (testType == "limma") {
        fit <- limma::lmFit(exprvals, design)
        fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
        res <- limma::topTreat(fit, coef = paste0("fc", comparison[2]),
                               number = Inf, sort.by = "none") %>%
            tibble::rownames_to_column("pid") %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(qftsub[[assayForTests]])) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select(.data$pid, .data$geneIdSingle,
                                  .data$proteinIdSingle),
                by = "pid")
    } else if (testType == "ttest") {
        res <- genefilter::rowttests(exprvals, fac = fc)
        res <- res %>%
            tibble::rownames_to_column("pid") %>%
            dplyr::rename(t = .data$statistic, logFC = .data$dm,
                          P.Value = .data$p.value) %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value),
                          adj.P.Val = p.adjust(.data$P.Value, method = "BH"),
                          AveExpr = rowMeans(exprvals)) %>%
            dplyr::mutate(sam = t/(1 + t * volcanoS0/.data$logFC)) %>%
            dplyr::left_join(as.data.frame(rowData(qftsub[[assayForTests]])) %>%
                                 tibble::rownames_to_column("pid") %>%
                                 dplyr::select(.data$pid, .data$geneIdSingle,
                                               .data$proteinIdSingle),
                             by = "pid")
    }

    ## --------------------------------------------------------------------- ##
    ## Test feature sets
    ## --------------------------------------------------------------------- ##
    featureCollections <- lapply(featureCollections, function(fcoll) {
        if (testType == "limma") {
            camres <- limma::cameraPR(
                statistic = structure(res$t, names = res$pid),
                index = limma::ids2indices(as.list(fcoll), res$pid,
                                           remove.empty = FALSE),
                sort = FALSE
            )
        } else if (testType == "ttest") {
            camres <- limma::cameraPR(
                statistic = structure(res$sam, names = res$pid),
                index = limma::ids2indices(as.list(fcoll), res$pid,
                                           remove.empty = FALSE),
                sort = FALSE
            )
        }
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

    ## --------------------------------------------------------------------- ##
    ## Get the threshold curve (replicating Perseus plots)
    ## --------------------------------------------------------------------- ##
    if (testType == "limma") {
        curveparam <- list()
    } else if (testType == "ttest") {
        curveparam <- .getThresholdCurve(
            exprvals = exprvals, fc = fc, res = res, seed = seed,
            nperm = nperm, volcanoS0 = volcanoS0,
            volcanoAdjPvalThr = volcanoAdjPvalThr)
    }

    ## --------------------------------------------------------------------- ##
    ## Determine significant features
    ## --------------------------------------------------------------------- ##
    res <- data.frame(pid = rownames(imputedvals)) %>%
        dplyr::left_join(res, by = "pid")
    if (testType == "limma") {
        res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
            abs(res$logFC) >= volcanoLog2FCThr
    } else if (testType == "ttest") {
        res$showInVolcano <- abs(res$sam) >= curveparam$ta
    }

    ## --------------------------------------------------------------------- ##
    ## Add abundance values and STRING IDs
    ## --------------------------------------------------------------------- ##
    if (addAbundanceValues) {
        tmp <- .addAbundanceValues(res = res, qftsub = qftsub, aName = aName,
                                   iColPattern = iColPattern)
        res <- tmp$res
    } else {
        tmp <- NULL
    }

    if ("IDsForSTRING" %in%
        colnames(SummarizedExperiment::rowData(qftsub[[assayForTests]]))) {
        res$IDsForSTRING <- SummarizedExperiment::rowData(
            qftsub[[assayForTests]])$IDsForSTRING[match(
                rownames(res),
                rownames(qftsub[[assayForTests]])
            )]
    }

    ## --------------------------------------------------------------------- ##
    ## Write results to file
    ## --------------------------------------------------------------------- ##
    if (!is.null(baseFileName)) {
        write.table(res %>%
                        dplyr::filter(.data$showInVolcano) %>%
                        dplyr::arrange(desc(.data$logFC)),
                    file = paste0(baseFileName, "_testres_", comparison[2],
                                  "_vs_", comparison[1], ".txt"),
                    row.names = FALSE, col.names = TRUE,
                    quote = FALSE, sep = "\t")
    }

    ## --------------------------------------------------------------------- ##
    ## Generate return values
    ## --------------------------------------------------------------------- ##
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

    return(list(res = res, plotnote = plotnote, plottitle = plottitle,
                plotsubtitle = plotsubtitle, colpat = tmp$colpat, topSets = topSets,
                featureCollections = featureCollections, curveparam = curveparam))
}
