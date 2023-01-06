#' @export
#' @author Charlotte Soneson
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom dplyr select mutate left_join
#' @importFrom tibble rownames_to_column
#' @importFrom limma lmFit contrasts.fit eBayes topTable topTreat treat
#'
runInteractionPTMTest <- function(sceProteins, scePeptides, comparisons,
                                  groupComposition = NULL, assayForTests,
                                  assayImputation, minNbrValidValues, minlFC,
                                  volcanoAdjPvalThr, volcanoLog2FCThr,
                                  singleFit) {

    ## If comparisons is not named, add names
    for (i in seq_along(comparisons)) {
        if (is.null(names(comparisons)) ||
            (!is.null(names(comparisons)) &&
             (is.na(names(comparisons)[i]) || names(comparisons)[i] == ""))) {
            names(comparisons)[i] <- paste0(comparisons[[i]][2], "_vs_",
                                            comparisons[[i]][1])
        }
    }
    if (any(duplicated(names(comparisons)))) {
        stop("Duplicated comparison names not allowed: ",
             paste(names(comparisons)[duplicated(names(comparisons))],
                   collapse = ", "))
    }

    ## If there are entries in unlist(comparisons) that are not defined in
    ## groupComposition, add them to the latter
    stdf <- setdiff(unlist(comparisons), names(groupComposition))
    groupComposition <- c(groupComposition,
                          setNames(as.list(stdf), stdf))

    ## Check that all entries in groupComposition[comparisons] are in the
    ## group column
    stdf <- setdiff(unlist(groupComposition[unlist(comparisons)]),
                    intersect(sceProteins$group, scePeptides$group))
    if (length(stdf) > 0) {
        stop("Missing group(s) in sceProteins/scePeptides$groups: ",
             paste(stdf, collapse = ", "))
    }

    ## --------------------------------------------------------------------- ##
    ## Initialize result lists
    ## --------------------------------------------------------------------- ##
    returndesign <- list()
    tests <- list()
    plottitles <- list()
    plotsubtitles <- list()
    plotnotes <- list()
    curveparams <- list()
    topsets <- list()
    messages <- list()

    exprvals <- cbind(SummarizedExperiment::assay(scePeptides, assayForTests),
                      SummarizedExperiment::assay(sceProteins, assayForTests))
    df <- rbind(as.data.frame(SummarizedExperiment::colData(scePeptides)) %>%
                    dplyr::select(sample, group) %>%
                    dplyr::mutate(dataLevel = "peptide"),
                as.data.frame(SummarizedExperiment::colData(sceProteins)) %>%
                    dplyr::select(sample, group) %>%
                    dplyr::mutate(dataLevel = "protein"))

    if (singleFit) {
        design <- model.matrix(~ sample, data = df)
        for (gr in unique(df$group)) {
            design <- cbind(design, as.numeric(df$group == gr &
                                                   df$dataLevel == "peptide"))
            colnames(design)[ncol(design)] <- paste0(gr, "_pept")
        }
        returndesign <- list(design = design, sampleData = df, contrasts = list())
        fit0 <- limma::lmFit(exprvals, design)
    }

    for (comparisonName in names(comparisons)) {
        comparison <- comparisons[[comparisonName]]
        idx <- which(df$group %in% unlist(groupComposition[comparison]))
        ## Only consider features with at least a given number of valid values
        imputedvalsPeptides <- SummarizedExperiment::assay(
            scePeptides[, colnames(scePeptides) %in% df$sample[idx]],
            assayImputation, withDimnames = TRUE)
        imputedvalsProteins <- SummarizedExperiment::assay(
            sceProteins[, colnames(sceProteins) %in% df$sample[idx]],
            assayImputation, withDimnames = TRUE)
        keep <- (rowSums(!imputedvalsPeptides) >= minNbrValidValues) &
            (rowSums(!imputedvalsProteins) > minNbrValidValues)

        if (singleFit) {
            fit <- fit0[keep, ]
            contrast <- (colnames(design) == paste0(comparison[2], "_pept")) -
                (colnames(design) == paste0(comparison[1], "_pept"))
            returndesign$contrasts[[comparisonName]] <- contrast
            fit <- limma::contrasts.fit(fit, contrasts = contrast)
            if (minlFC == 0) {
                fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                res <- limma::topTable(fit, coef = 1,
                                       number = Inf, sort.by = "none")
            } else {
                fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                res <- limma::topTreat(fit, coef = 1,
                                       number = Inf, sort.by = "none")
            }
            res <- res %>%
                tibble::rownames_to_column("pid") %>%
                dplyr::mutate(s2.prior = fit$s2.prior,
                              weights = fit$weights,
                              sigma = fit$sigma)
            camerastat <- "t"
        } else {
            exprvalssub <- exprvals[keep, idx, drop = FALSE]
            dfsub <- df[idx, , drop = FALSE]
            design <- model.matrix(~ sample, data = dfsub)
            for (gr in unique(dfsub$group)) {
                design <- cbind(design, as.numeric(dfsub$group == gr &
                                                       dfsub$dataLevel == "peptide"))
                colnames(design)[ncol(design)] <- paste0(gr, "_pept")
            }
            contrast <- (colnames(design) == paste0(comparison[2], "_pept")) -
                (colnames(design) == paste0(comparison[1], "_pept"))
            returndesign[[comparisonName]] <-
                list(design = design, sampleData = dfsub, contrast = contrast)

            ## ------------------------------------------------------------- ##
            ## Run test
            ## ------------------------------------------------------------- ##
            fit <- limma::lmFit(exprvalssub, design)
            fit <- limma::contrasts.fit(fit, contrasts = contrast)
            if (minlFC == 0) {
                fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                res <- limma::topTable(fit, coef = 1, number = Inf, sort.by = "none")
            } else {
                fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                res <- limma::topTreat(fit, coef = 1,
                                       number = Inf, sort.by = "none")
            }
            res <- res %>%
                tibble::rownames_to_column("pid") %>%
                dplyr::mutate(s2.prior = fit$s2.prior,
                              weights = fit$weights,
                              sigma = fit$sigma)
            camerastat <- "t"
        }

        ## Calculate -log10(p). For features with p-value = 0, use half the
        ## smallest non-zero p-value as a proxy to be able to make volcano plots.
        res <- res %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::mutate(mlog10p = replace(
                .data$mlog10p, .data$P.Value == 0,
                -log10(min(.data$P.Value[which(.data$P.Value > 0)])/2))) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(scePeptides)) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select(.data$pid, .data$einprotGene,
                                  .data$einprotProtein, .data$einprotLabel),
                by = "pid")

        curveparam <- list()
        topSets <- list()

        res <- data.frame(pid = rownames(imputedvalsPeptides)) %>%
            dplyr::left_join(res, by = "pid")
        res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
            abs(res$logFC) >= volcanoLog2FCThr

        plottitle <- paste0(comparison[2], " vs ", comparison[1],
                            ", limma")
        plotnote <- paste0("df.prior = ", round(fit$df.prior, digits = 2))
        plotsubtitle <- paste0("Adj.p threshold = ", volcanoAdjPvalThr,
                               ", |log2FC| threshold = ", volcanoLog2FCThr)

        plottitles[[comparisonName]] <- plottitle
        plotsubtitles[[comparisonName]] <- plotsubtitle
        plotnotes[[comparisonName]] <- plotnote
        tests[[comparisonName]] <- res
        curveparams[[comparisonName]] <- curveparam
        topsets[[comparisonName]] <- topSets

    }

    return(list(plottitles = plottitles, plotsubtitles = plotsubtitles,
                plotnotes = plotnotes, tests = tests,
                curveparams = curveparams, topsets = topsets,
                messages = messages, design = returndesign))
}
