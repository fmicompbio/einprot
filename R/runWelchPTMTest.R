#' @export
#' @author Charlotte Soneson
#'
#' @importFrom SummarizedExperiment rowData
runWelchPTMTest <- function(sceProteins, scePeptides, comparisons,
                            volcanoAdjPvalThr, volcanoLog2FCThr) {

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

    tests <- list()
    plottitles <- list()
    plotsubtitles <- list()
    plotnotes <- list()
    curveparams <- list()
    topsets <- list()
    messages <- list()

    resProtein <- SummarizedExperiment::rowData(sceProteins)
    resPeptide <- SummarizedExperiment::rowData(scePeptides)

    for (comparisonName in names(comparisons)) {
        df <- data.frame(pid = rownames(resPeptide))
        df$logFC <- resPeptide[[paste0(comparisonName, ".logFC")]] -
            resProtein[[paste0(comparisonName, ".logFC")]]
        s2_1 <- resPeptide[[paste0(comparisonName, ".se.logFC")]] ^ 2
        s2_2 <- resProtein[[paste0(comparisonName, ".se.logFC")]] ^ 2
        df$se.logFC <- sqrt(s2_1 + s2_2)
        numer <- (s2_1 + s2_2) ^ 2
        denom <- (s2_1 ^ 2 / resPeptide[[paste0(comparisonName, ".df.total")]] +
                      s2_2 ^ 2 / resProtein[[paste0(comparisonName, ".df.total")]])
        df$df <- numer / denom
        df$t <- df$logFC / df$se.logFC
        df$P.Value <- 2 * stats::pt(abs(df$t), df$df, lower.tail = FALSE)
        df$adj.P.Val <- p.adjust(df$P.Value, method = "BH")

        ## Calculate -log10(p). For features with p-value = 0, use half the
        ## smallest non-zero p-value as a proxy to be able to make volcano plots.
        res <- df %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::mutate(mlog10p = replace(
                .data$mlog10p, .data$P.Value == 0,
                -log10(min(.data$P.Value[which(.data$P.Value > 0)])/2))) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(scePeptides)) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select("pid", "einprotGene",
                                  "einprotProtein", "einprotLabel"),
                by = "pid")

        curveparam <- list()
        topSets <- list()

        res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
            abs(res$logFC) >= volcanoLog2FCThr

        plottitle <- paste0(comparisonName, ", limma")
        plotnote <- ""
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
                messages = messages))
}
