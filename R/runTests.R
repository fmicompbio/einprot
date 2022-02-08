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
.addiBAQvalues <- function(res, qftsub, aName, iColPattern) {
    stopifnot(all(rownames(SummarizedExperiment::rowData(qftsub[[aName]])) == res$pid))
    if (iColPattern == "^iBAQ\\.") {
        ## If we're using iBAQ values for the analysis, they will not be
        ## present in the rowData anymore.
        ibaq_values <- as.data.frame(SummarizedExperiment::assay(qftsub[[aName]]))
        ibaq_values[is.na(ibaq_values)] <- 0
    } else {
        ibaq_values <- as.data.frame(
            SummarizedExperiment::rowData(
                qftsub[[aName]])[, sub(iColPattern, "iBAQ.",
                                       colnames(qftsub[[aName]])), drop = FALSE]
        )
    }
    log_ibaq_values <- log2(ibaq_values) %>%
        dplyr::mutate(dplyr::across(dplyr::everything(),
                                    .fns = ~ ifelse(is.finite(.), ., NA)))
    res <- cbind(res, ibaq_values)
    for (gr in unique(SummarizedExperiment::colData(qftsub)$group)) {
        res[[paste0("iBAQ.", gr, ".avg")]] <-
            rowMeans(
                ibaq_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                            drop = FALSE])
        res[[paste0("iBAQ.", gr, ".sd")]] <-
            genefilter::rowSds(
                ibaq_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                            drop = FALSE])
        res[[paste0("log2_iBAQ.", gr, ".avg")]] <-
            rowMeans(
                log_ibaq_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                                drop = FALSE],
                na.rm = TRUE)
        res[[paste0("log2_iBAQ.", gr, ".sd")]] <-
            genefilter::rowSds(
                log_ibaq_values[, SummarizedExperiment::colData(qftsub)$group == gr,
                                drop = FALSE],
                na.rm = TRUE)
    }
    res
}

#' Run statistical test
#'
#' @param qft QFeatures object.
#' @param comparison Vector of length 2, giving the two groups to be compared.
#' @param testType Either "limma" or "ttest".
#' @param assayForTests A name of an assay of the QFeatures object.
#' @param assayImputation A name of an assay of the QFeatures object.
#' @param minNbrValidValues Numeric scalar, min number of valid values to run
#'     the test.
#' @param minlFC Min logFC threshold for treat.
#' @param featureCollections List of CharacterLists with feature collections.
#' @param complexFDRThr Numeric scalar giving the FDR threshold to plot
#'     significant complexes.
#' @param volcanoAdjPvalThr Numeric scalar giving the FDR threshold for
#'     significance.
#' @param volcanoLog2FCThr Numeric scalar giving the logFC threshold for
#'     significance.
#' @param baseFileName Base file name.
#' @param iColPattern iColPattern.
#' @param seed Random seed.
#' @param nperm Number of permutations.
#' @param volcanoS0 S0 value to use for creating significance curves.
#' @param aName The name of the base assay in the QFeatures object.
#'
#' @author Charlotte Soneson
#' @export
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
                    assayImputation, minNbrValidValues,
                    minlFC, featureCollections, complexFDRThr,
                    volcanoAdjPvalThr, volcanoLog2FCThr,
                    baseFileName, iColPattern, seed, nperm,
                    volcanoS0, aName) {
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest"))

    ## ----------------------------------------------------------------- ##
    ## Subset and define design
    ## ----------------------------------------------------------------- ##
    qftsub <- qft[, qft$group %in% comparison]
    if (testType == "limma") {
        fc <- factor(qftsub$group, levels = comparison)
        design <- stats::model.matrix(~ fc)
    } else if (testType == "ttest") {
        fc <- factor(qftsub$group, levels = rev(comparison))
    }

    exprvals <- SummarizedExperiment::assay(qftsub[[assayForTests]], withDimnames = TRUE)
    ## Only consider features with at least a given number of valid values
    imputedvals <- SummarizedExperiment::assay(qftsub[[assayImputation]],
                                               withDimnames = TRUE)
    keep <- rowSums(!imputedvals) >= minNbrValidValues
    exprvals <- exprvals[keep, , drop = FALSE]

    ## ----------------------------------------------------------------- ##
    ## Run test
    ## ----------------------------------------------------------------- ##
    if (testType == "limma") {
        fit <- limma::lmFit(exprvals, design)
        fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
        res <- limma::topTreat(fit, coef = paste0("fc", comparison[2]),
                               number = Inf, sort.by = "none") %>%
            tibble::rownames_to_column("pid") %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value))
    } else if (testType == "ttest") {
        res <- genefilter::rowttests(exprvals, fac = fc)
        res <- res %>%
            tibble::rownames_to_column("pid") %>%
            dplyr::rename(t = statistic, logFC = dm, P.Value = p.value) %>%
            dplyr::mutate(mlog10p = -log10(P.Value),
                          adj.P.Val = p.adjust(P.Value, method = "BH"),
                          AveExpr = rowMeans(exprvals)) %>%
            dplyr::mutate(sam = t/(1 + t * volcanoS0/logFC))
    }

    ## ----------------------------------------------------------------- ##
    ## Test feature sets
    ## ----------------------------------------------------------------- ##
    featureCollections <- lapply(featureCollections, function(fcoll) {
        if (testType == "limma") {
            camres <- limma::cameraPR(
                statistic = structure(res$t, names = res$pid),
                index = limma::ids2indices(as.list(fcoll), res$pid, remove.empty = FALSE),
                sort = FALSE
            )
        } else if (testType == "ttest") {
            camres <- cameraPR(
                statistic = structure(res$sam, names = res$pid),
                index = ids2indices(as.list(fcoll), res$pid, remove.empty = FALSE),
                sort = FALSE
            )
        }
        if (!("FDR" %in% colnames(camres))) {
            camres$FDR <- stats::p.adjust(camres$PValue, method = "BH")
        }
        colnames(camres) <- paste0(comparison[2], "_vs_", comparison[1], "_", colnames(camres))
        stopifnot(all(rownames(camres) == names(fcoll)))
        S4Vectors::mcols(fcoll) <- cbind(S4Vectors::mcols(fcoll), camres)
        fcoll
    })

    ## Write test results for complexes to text files
    if ("complexes" %in% names(featureCollections)) {
        tmpres <- as.data.frame(S4Vectors::mcols(featureCollections$complexes), optional = TRUE) %>%
            tibble::rownames_to_column("complex") %>%
            dplyr::select(.data$complex, .data$genes, .data$sharedGenes,
                          .data$Source, .data$All.names, .data$PMID,
                          dplyr::contains(paste0(comparison[2], "_vs_", comparison[1]))) %>%
            dplyr::arrange(.data[[paste0(comparison[2], "_vs_", comparison[1], "_FDR")]]) %>%
            dplyr::filter(.data[[paste0(comparison[2], "_vs_", comparison[1], "_FDR")]] < complexFDRThr)
        if (nrow(tmpres) > 0) {
            write.table(tmpres,
                        file = paste0(baseFileName,
                                      paste0("_testres_", comparison[2], "_vs_", comparison[1],
                                             "_camera_complexes.txt")),
                        row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
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
    ## Add iBAQ values and STRING IDs
    ## --------------------------------------------------------------------- ##
    res <- .addiBAQvalues(res = res, qftsub = qftsub, aName = aName,
                          iColPattern = iColPattern)

    res$IDsForSTRING <- SummarizedExperiment::rowData(
        qftsub[[assayForTests]])$IDsForSTRING[match(
            rownames(res),
            rownames(qftsub[[assayForTests]])
        )]

    ## --------------------------------------------------------------------- ##
    ## Generate return values
    ## --------------------------------------------------------------------- ##
    if (testType == "limma") {
        plottitle <- paste0(comparison[2], " vs ", comparison[1],
                            ", limma treat (H0: |log2FC| <= ", minlFC, ")")
        plotnote <- paste0("df.prior = ", round(fit$df.prior, digits = 2))
    } else if (testType == "ttest") {
        plottitle <- paste0(comparison[2], " vs ", comparison[1], ", t-test")
        plotnote <- ""
    }

    return(list(plotnote = plotnote, plottitle = plottitle, res = res,
                featureCollections = featureCollections, curveparam = curveparam))
}
