#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom stats pt
.curvefun <- function(x, ta, s0, df) {
    -log10(2 * (1 - stats::pt(q = ta * (1 + s0/(abs(x)/ta - s0)), df = df)))
}

#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom dplyr filter select mutate left_join %>% matches
#'     group_by summarize
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment colData
#' @importFrom ggplot2 ggplot aes geom_bar position_dodge geom_errorbar
#'     theme_bw theme element_text labs scale_fill_manual geom_jitter
#' @importFrom stats sd
#'
.complexBarPlot <- function(res, prs, sce, cplx, colpat, groupmap) {
    bardata <- res %>%
        dplyr::filter(.data$pid %in% prs) %>%
        dplyr::select("pid", dplyr::matches(colpat)) %>%
        tidyr::gather(key = "sample", value = "Abundance", -"pid") %>%
        dplyr::mutate(sample = sub(paste0("^", colpat, "\\."),
                                   "", .data$sample)) %>%
        dplyr::left_join(as.data.frame(
            SummarizedExperiment::colData(sce)), by = "sample") %>%
        dplyr::filter(!is.na(.data$group))
    if (!is.null(groupmap)) {
        bardata <- bardata %>%
            dplyr::left_join(groupmap, by = "group")
    } else {
        bardata$mergegroup <- bardata$group
    }
    ggbar <- ggplot(
        bardata %>% dplyr::group_by(.data$pid, .data$mergegroup) %>%
            dplyr::summarize(
                mean_abundance = mean(.data$Abundance, na.rm = TRUE),
                sd_abundance = stats::sd(.data$Abundance, na.rm = TRUE),
                .groups = "drop"),
        aes(x = .data$pid, y = .data$mean_abundance,
            fill = .data$mergegroup)) +
        geom_bar(position = position_dodge(), stat = "identity",
                 colour = "black", linewidth = 0.3) +
        geom_errorbar(aes(ymin = .data$mean_abundance - .data$sd_abundance,
                          ymax = .data$mean_abundance + .data$sd_abundance),
                      linewidth = 0.3, width = 0.2,
                      position = position_dodge(width = 0.9)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12, angle = 90,
                                         hjust = 1, vjust = 0.5),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              title = element_text(size = 14)) +
        labs(x = "", y = paste0("Mean +/- SD ", colpat), title = cplx) +
        scale_fill_manual(name = "", values = c("steelblue", "firebrick2"))
    if (length(unique(bardata$sample)) <= 6) {
        ggbar <- ggbar +
            geom_jitter(data = bardata, aes(y = .data$Abundance,
                                            shape = .data$sample), size = 2,
                        position = position_dodge(width = 0.9))
    } else {
        ggbar <- ggbar +
            geom_jitter(data = bardata, aes(y = .data$Abundance), size = 2,
                        position = position_dodge(width = 0.9))
    }
    ggbar
}

#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom dplyr bind_rows filter arrange desc slice mutate
#' @importFrom forcats fct_reorder
#' @importFrom rlang .data
#' @importFrom ggplot2 ggplot aes geom_col coord_flip geom_text theme_minimal
#'     theme element_blank ggtitle element_line scale_y_continuous
#'     scale_fill_gradient2
#' @importFrom scales muted
.makeWaterfallPlot <- function(res, ntop, xv = "logFC",
                               volcind = "showInVolcano", title = "") {
    .assertVector(x = res, type = "data.frame")
    .assertScalar(x = ntop, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = xv, type = "character", validValues = colnames(res))
    .assertScalar(x = volcind, type = "character", validValues = colnames(res))
    .assertScalar(x = title, type = "character")

    a <- res %>%
        dplyr::filter(.data[[volcind]]) %>%
        dplyr::arrange(dplyr::desc(abs(.data[[xv]]))) %>%
        dplyr::slice(seq_len(ntop))
    rng <- c(-max(abs(a[[xv]])), max(abs(a[[xv]])))
    a <- a %>%
        dplyr::mutate(label_y = ifelse(.data[[xv]] < 0, rng[2]/20, rng[1]/20),
                      label_hjust = ifelse(.data[[xv]] < 0, 0, 1))
    ggplot2::ggplot(a, ggplot2::aes(x = forcats::fct_reorder(.data$pid, .data[[xv]]),
                                    y = .data[[xv]], fill = sign(.data[[xv]]))) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::geom_text(ggplot2::aes(label = .data$einprotLabel, y = .data$label_y,
                                        hjust = .data$label_hjust)) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.text.y = ggplot2::element_blank(),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title.y = ggplot2::element_blank(),
            legend.position = "none",
            panel.grid.major.y = ggplot2::element_blank(),
            panel.grid.minor.y = ggplot2::element_blank(),
            panel.grid.major.x = ggplot2::element_line(colour = "grey80",
                                                       linetype = "dashed"),
            panel.grid.minor.x = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(size = 12),
            axis.title.x = ggplot2::element_text(size = 14),
            title = ggplot2::element_text(size = 14)) +
        ggplot2::scale_y_continuous(limits = rng) +
        ggplot2::scale_fill_gradient2(low = scales::muted("red"),
                                      high = scales::muted("blue"),
                                      limits = c(-1, 1)) +
        ggplot2::labs(title = title, y = "log2(fold change)")
}

#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#' @importFrom ggplot2 ggplot aes theme_bw theme annotate labs coord_cartesian
#'     element_text geom_line
#' @importFrom rlang .data
#'
.makeBaseVolcano <- function(res, testType, xv, yv, plotnote, plottitle,
                             plotsubtitle, curveparam) {
    xr <- range(res[[xv]], na.rm = TRUE)
    xr <- c(-max(abs(xr), na.rm = TRUE), max(abs(xr), na.rm = TRUE))
    yr <- range(res[[yv]], na.rm = TRUE)

    ggbase <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[xv]],
                                                y = .data[[yv]])) +
        ggplot2::theme_bw() +
        ggplot2::coord_cartesian(xlim = xr, ylim = yr) +
        ggplot2::theme(axis.text = ggplot2::element_text(size = 12),
                       axis.title = ggplot2::element_text(size = 14),
                       title = ggplot2::element_text(size = 14)) +
        ggplot2::annotate("text", x = min(xr), y = 0, hjust = "left",
                          vjust = "bottom", label = plotnote) +
        ggplot2::labs(x = "log2(fold change)", y = "-log10(p-value)",
                      title = plottitle, subtitle = plotsubtitle)
    if (testType == "ttest") {
        if (length(curveparam) > 0 && is.finite(curveparam$ta)) {
            ggbase <- ggbase +
                ggplot2::geom_line(
                    data = data.frame(x = curveparam$x,
                                      y = .curvefun(x = curveparam$x,
                                                    ta = curveparam$ta,
                                                    s0 = curveparam$s0,
                                                    df = curveparam$df)),
                    ggplot2::aes(x = .data$x, y = .data$y), color = "red",
                    linetype = "dashed") +
                ggplot2::geom_line(
                    data = data.frame(x = -curveparam$x,
                                      y = .curvefun(x = curveparam$x,
                                                    ta = curveparam$ta,
                                                    s0 = curveparam$s0,
                                                    df = curveparam$df)),
                    ggplot2::aes(x = .data$x, y = .data$y), color = "red",
                    linetype = "dashed")
        }
    }

    ggbase
}

.getVolcanoColumns <- function(testType) {
    if (testType == "limma") {
        xv <- "logFC"
        yv <- "mlog10p"
        xvma = "AveExpr"
        apv <- "adj.P.Val"
        tv <- NULL
        volcind <- "showInVolcano"
    } else if (testType == "ttest") {
        xv <- "logFC"
        yv <- "mlog10p"
        xvma = NULL
        apv <- "adj.P.Val"
        tv <- "sam"
        volcind <- "showInVolcano"
    } else if (testType == "proDA") {
        xv <- "logFC"
        yv <- "mlog10p"
        xvma = NULL
        apv <- "adj.P.Val"
        tv <- NULL
        volcind <- "showInVolcano"
    }
    list(xv = xv, yv = yv, xvma = xvma, apv = apv, volcind = volcind, tv = tv)
}

#' Make volcano plots
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param res A \code{data.frame} object with test results (typically
#'     generated by \code{runTest}).
#' @param testType Character scalar indicating the type of test that was run,
#'     either \code{"ttest"}, \code{"limma"} or \code{"proDA"}.
#' @param xv,yv Character scalars indicating which columns of \code{res} that
#'     should be used as the x- and y-axis of the volcano plot, respectively.
#'     If \code{NULL}, will be determined based on \code{testType}.
#' @param xvma If not \code{NULL}, a character scalar indicating which
#'     column of \code{res} should be used as the x-axis for an MA plot. The
#'     y-axis column will be \code{xv}. If \code{NULL}, no MA plot is
#'     generated.
#' @param volcind Character scalar indicating which column in \code{res} that
#'     represents the "significance" column. This should be a logical
#'     column; rows with a value equal to \code{TRUE} will be colored
#'     in the plot. If \code{NULL}, will be determined based on \code{testType}.
#' @param plotnote Character scalar with a note to add to the plot.
#' @param plottitle Character scalar giving the title of the plot.
#' @param plotsubtitle Character scalar giving the subtitle of the plot.
#' @param volcanoFeaturesToLabel Character vector, features to label in the
#'     plot.
#' @param volcanoMaxFeatures Numeric scalar, the maximum number of features
#'     to color in the plot.
#' @param baseFileName Character scalar or \code{NULL}, the base file name of
#'     the output files. If \code{NULL}, no result files are generated.
#' @param comparisonString Character scalar giving the name of the comparison of
#'     interest. This is used to extract the appropriate column from the
#'     metadata of the feature collections to make the gene set plots.
#' @param groupComposition A list providing the composition of each group
#'     used in the comparisons indicated by \code{comparisonString}.
#'     If \code{NULL}, assumes that each group used in \code{comparisonString}
#'     consists of a single group in the \code{group} column of
#'     \code{colData(sce)}.
#' @param stringDb A STRINGdb object or \code{NULL}. If not \code{NULL},
#'     STRING network plots of up- and downregulated genes will be added to
#'     the output pdf file.
#' @param featureCollections A list of \code{CharacterList}s with feature
#'     collections. If there is a collection named \code{"complexes"},
#'     volcano plots and barplots for significant complexes are included in
#'     the output pdf files.
#' @param complexFDRThr Numeric scalar giving the FDR threshold for complexes
#'     to be considered significant.
#' @param maxNbrComplexesToPlot Numeric scalar, the largest number of
#'     significant complexes to generate separate volcano plots for.
#' @param curveparam List with curve parameters for creating the Perseus-like
#'     significance curves in the volcano plots.
#' @param abundanceColPat Character scalar providing the column pattern used
#'     to identify abundance columns in the result table, to make bar plots
#'     for significant complexes. Typically the name of an assay.
#'
#' @return A list with two plot objects; one ggplot object and one
#'     interactive ggiraph object. If \code{xvma} is not \code{NULL},
#'     the list will also contain a ggplot object for the MA plot.
#'     If \code{baseFileName} is not \code{NULL}, pdf files with volcano
#'     plots and bar plots for significant complexes will also be generated.
#'
#' @author Charlotte Soneson
#' @export
#' @importFrom ggplot2 geom_point aes labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggiraph geom_point_interactive
#' @importFrom dplyr filter arrange between row_number desc
#' @importFrom rlang .data
#' @importFrom utils stack
#'
plotVolcano <- function(sce, res, testType, xv = NULL, yv = NULL, xvma = NULL,
                        volcind = NULL, plotnote, plottitle, plotsubtitle,
                        volcanoFeaturesToLabel, volcanoMaxFeatures,
                        baseFileName, comparisonString, groupComposition = NULL,
                        stringDb, featureCollections, complexFDRThr,
                        maxNbrComplexesToPlot, curveparam,
                        abundanceColPat = "") {

    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertVector(x = res, type = "data.frame")
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest", "proDA"))

    cols <- .getVolcanoColumns(testType = testType)
    if (!is.null(xv)) {
        cols$xv <- xv
    }
    if (!is.null(yv)) {
        cols$yv <- yv
    }
    if (!is.null(xvma)) {
        cols$xvma <- xvma
    }
    if (!is.null(volcind)) {
        cols$volcind <- volcind
    }

    .assertScalar(x = cols$xv, type = "character", validValues = colnames(res))
    .assertScalar(x = cols$yv, type = "character", validValues = colnames(res))
    .assertScalar(x = cols$xvma, type = "character", allowNULL = TRUE,
                  validValues = colnames(res))
    .assertScalar(x = cols$volcind, type = "character",
                  validValues = colnames(res))
    .assertScalar(x = plotnote, type = "character")
    .assertScalar(x = plottitle, type = "character")
    .assertScalar(x = plotsubtitle, type = "character")
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertScalar(x = volcanoMaxFeatures, type = "numeric",
                  rngIncl = c(0, Inf))
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertScalar(x = comparisonString, type = "character")
    .assertVector(x = groupComposition, type = "list", allowNULL = TRUE)
    .assertVector(x = stringDb, type = "STRINGdb", allowNULL = TRUE)
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = maxNbrComplexesToPlot, type = "numeric", rngIncl = c(0, Inf))
    .assertVector(x = curveparam, type = "list")
    .assertScalar(x = abundanceColPat, type = "character")

    ## If the 'einprotLabel' column is not available, create it using the 'pid' column
    if (!("einprotLabel" %in% colnames(res))) {
        res$einprotLabel <- res$pid
    }

    ## Make "base" volcano plot
    ggbase <- .makeBaseVolcano(res = res, testType = testType, xv = cols$xv, yv = cols$yv,
                               plotnote = plotnote, plottitle = plottitle,
                               plotsubtitle = plotsubtitle, curveparam = curveparam)

    ## Add labels for significant features
    ggtest <- ggbase +
        ggplot2::geom_point(fill = "lightgrey", color = "grey",
                            pch = 21, size = 1.5) +
        ggplot2::geom_point(data = res %>%
                       dplyr::filter(.data[[cols$volcind]]),
                   fill = "red", color = "grey", pch = 21, size = 1.5) +
        ggrepel::geom_text_repel(
            data = res %>%
                dplyr::filter(
                    .data[[cols$volcind]] |
                        .data$pid %in% volcanoFeaturesToLabel
                ) %>%
                dplyr::arrange(
                    dplyr::desc(abs(.data[[cols$xv]]) + abs(.data[[cols$yv]]))
                ) %>%
                dplyr::filter(dplyr::between(dplyr::row_number(), 0,
                                             volcanoMaxFeatures) |
                                  .data$pid %in% volcanoFeaturesToLabel),
            aes(label = .data$einprotLabel), max.overlaps = Inf, size = 4,
            min.segment.length = 0.1)

    ## Interactive version
    ggint <- ggbase +
        ggiraph::geom_point_interactive(
            aes(tooltip = .data$einprotLabel), fill = "lightgrey", color = "grey",
            pch = 21, size = 1.5) +
        ggiraph::geom_point_interactive(
            data = res %>% dplyr::filter(.data[[cols$volcind]]),
            aes(tooltip = .data$einprotLabel), fill = "red", color = "grey",
            pch = 21, size = 1.5)

    ## MA plot
    if (!is.null(cols$xvma)) {
        yrma <- range(res[[cols$xv]], na.rm = TRUE)
        yrma <- c(-max(abs(res[[cols$xv]]), na.rm = TRUE), max(abs(yrma), na.rm = TRUE))
        ggma <- ggplot(res, aes(x = .data[[cols$xvma]], y = .data[[cols$xv]])) +
            geom_hline(yintercept = 0, color = "black", linetype = "dashed") +
            geom_point(fill = "lightgrey", color = "grey", pch = 21, size = 1.5) +
            theme_bw() + coord_cartesian(ylim = yrma) +
            theme(axis.text = element_text(size = 12),
                  axis.title = element_text(size = 14),
                  title = element_text(size = 14)) +
            labs(x = "Average abundance", y = "log2(fold change)",
                 title = plottitle, subtitle = plotsubtitle) +
            geom_point(data = res %>%
                           dplyr::filter(.data[[cols$volcind]]),
                       fill = "red", color = "grey", pch = 21, size = 1.5) +
            geom_text_repel(
                data = res %>%
                    dplyr::filter(.data[[cols$volcind]] |
                                      .data$pid %in% volcanoFeaturesToLabel) %>%
                    dplyr::arrange(desc(abs(.data[[cols$xv]]) + abs(.data[[cols$yv]]))) %>%
                    dplyr::filter(dplyr::between(row_number(), 0, volcanoMaxFeatures) |
                                      .data$pid %in% volcanoFeaturesToLabel),
                aes(label = .data$einprotLabel), max.overlaps = Inf, size = 4,
                min.segment.length = 0.1)
    } else {
        ggma <- NULL
    }

    ## Waterfall plot
    if (any(!is.na(res[[cols$volcind]]) & res[[cols$volcind]])) {
        ggwf <- .makeWaterfallPlot(res = res, ntop = volcanoMaxFeatures,
                                   xv = cols$xv,
                                   volcind = cols$volcind, title = plottitle)
    } else {
        ggwf <- NULL
    }

    ## --------------------------------------------------------------------- ##
    ## Write to file, including STRING plots
    ## --------------------------------------------------------------------- ##
    if (!is.null(baseFileName)) {
        pdf(paste0(baseFileName, paste0("_volcano_", comparisonString, ".pdf")),
            width = 10.5, height = 7.5)
        print(ggtest)

        ## STRING plots of up- and downregulated proteins
        if (!is.null(stringDb) && "IDsForSTRING" %in% colnames(res)) {
            res0 <- res %>%
                dplyr::filter(.data[[cols$volcind]]) %>%
                dplyr::arrange(dplyr::desc(abs(.data[[cols$xv]]) + abs(.data[[cols$yv]]))) %>%
                dplyr::filter(dplyr::between(dplyr::row_number(), 0, volcanoMaxFeatures))
            res0 <- stringDb$map(res0, "IDsForSTRING", removeUnmappedRows = TRUE)
            if (any(res0[[cols$xv]] > 0)) {
                stringDb$plot_network(res0 %>% dplyr::filter(.data[[cols$xv]] > 0) %>%
                                          dplyr::pull("STRING_id"))
            }
            if (any(res0[[cols$xv]] < 0)) {
                stringDb$plot_network(res0 %>% dplyr::filter(.data[[cols$xv]] < 0) %>%
                                          dplyr::pull("STRING_id"))
            }
        }

        if (!is.null(ggwf)) {
            print(ggwf)
        }
        dev.off()
    }

    ## --------------------------------------------------------------------- ##
    ## Create a volcano plot for each significantly enriched complex
    ## --------------------------------------------------------------------- ##
    if ("complexes" %in% names(featureCollections)) {
        ## Find significant complexes
        idx <- which(
            mcols(featureCollections$complexes)[, paste0(comparisonString,
                                                         "_FDR")] < complexFDRThr &
                mcols(featureCollections$complexes)[, paste0(comparisonString,
                                                             "_NGenes")] > 1
        )
        tmpcomplx <- mcols(featureCollections$complexes)[idx, , drop = FALSE]
        tmpcomplx <- tmpcomplx[order(tmpcomplx[paste0(comparisonString,
                                                      "_PValue")]), , drop = FALSE]
        cplxs <- rownames(tmpcomplx)
        cplxs <- cplxs[seq_len(min(length(cplxs), maxNbrComplexesToPlot))]

        if (length(cplxs) > 0 && !is.null(baseFileName)) {
            pdf(paste0(baseFileName, "_volcano_", comparisonString,
                       "_complexes.pdf"), width = 10.5, height = 7.5)
            for (cplx in cplxs) {
                prs <- featureCollections$complexes[[cplx]]
                cplxpval <- signif(mcols(
                    featureCollections$complexes)[cplx, paste0(comparisonString,
                                                               "_PValue")], digits = 3)
                cplxfdr <- signif(mcols(
                    featureCollections$complexes)[cplx, paste0(comparisonString,
                                                               "_FDR")], digits = 3)
                if (length(intersect(prs, res$pid)) > 1) {
                    gg <- ggbase +
                        ggplot2::geom_point(
                            fill = "lightgrey", color = "grey",
                            pch = 21, size = 1.5) +
                        ggplot2::geom_point(
                            data = res %>%
                                dplyr::filter(.data$pid %in% prs),
                            fill = "red", color = "grey", pch = 21, size = 1.5) +
                        ggrepel::geom_text_repel(
                            data = res %>%
                                dplyr::filter(.data$pid %in% prs),
                            aes(label = .data$pid), max.overlaps = Inf, size = 4,
                            min.segment.length = 0.1) +
                        ggplot2::labs(caption = paste0(cplx, ", PValue = ", cplxpval,
                                                       ", FDR = ", cplxfdr))
                    print(gg)

                    ## Bar plot
                    if (is.null(groupComposition)) {
                        groupmap <- NULL
                    } else {
                        groupmap <- utils::stack(groupComposition) %>%
                            setNames(c("group", "mergegroup"))
                    }
                    print(.complexBarPlot(
                        res = res, prs = prs, sce = sce, cplx = cplx,
                        colpat = abundanceColPat, groupmap = groupmap))
                }
            }
            dev.off()
        }
    }
    return(list(gg = ggtest, ggint = ggiraph::girafe(ggobj = ggint),
                ggma = ggma, ggwf = ggwf))
}
