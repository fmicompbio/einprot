#' @noRd
#' @keywords internal
#' @importFrom stats pt
.curvefun <- function(x, ta, s0, df) {
    -log10(2 * (1 - stats::pt(q = ta * (1 + s0/(abs(x)/ta - s0)), df = df)))
}

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
.complexBarPlot <- function(res, prs, qft, cplx) {
    bardata <- res %>%
        dplyr::filter(.data$pid %in% prs) %>%
        dplyr::select(.data$pid, dplyr::matches("^iBAQ")) %>%
        tidyr::gather(key = "sample", value = "iBAQ", -.data$pid) %>%
        dplyr::mutate(sample = sub("^iBAQ.", "", .data$sample)) %>%
        dplyr::left_join(as.data.frame(colData(qft)), by = "sample") %>%
        dplyr::filter(!is.na(.data$group))
    ggbar <- ggplot(bardata %>% dplyr::group_by(.data$pid, .data$group) %>%
                        dplyr::summarize(mean_iBAQ = mean(.data$iBAQ, na.rm = TRUE),
                                         sd_iBAQ = stats::sd(.data$iBAQ, na.rm = TRUE),
                                         .groups = "drop"),
                    aes(x = .data$pid, y = .data$mean_iBAQ,
                        fill = .data$group)) +
        geom_bar(position = position_dodge(), stat = "identity",
                 colour = "black", size = 0.3) +
        geom_errorbar(aes(ymin = .data$mean_iBAQ - .data$sd_iBAQ,
                          ymax = .data$mean_iBAQ + .data$sd_iBAQ),
                      size = 0.3, width = 0.2,
                      position = position_dodge(width = 0.9)) +
        theme_bw() +
        theme(axis.text.x = element_text(size = 12, angle = 90,
                                         hjust = 1,
                                         vjust = 0.5),
              axis.text.y = element_text(size = 12),
              axis.title = element_text(size = 14),
              title = element_text(size = 14)) +
        labs(x = "", y = "Mean +/- SD iBAQ", title = cplx) +
        scale_fill_manual(name = "",
                          values = c("steelblue", "firebrick2"))
    if (length(unique(bardata$sample)) <= 6) {
        ggbar <- ggbar +
            geom_jitter(data = bardata, aes(y = .data$iBAQ,
                                            shape = .data$sample), size = 2,
                        position = position_dodge(width = 0.9))
    } else {
        ggbar <- ggbar +
            geom_jitter(data = bardata, aes(y = .data$iBAQ), size = 2,
                        position = position_dodge(width = 0.9))
    }
    ggbar
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

    ggbase <- ggplot2::ggplot(res, ggplot2::aes(x = .data[[xv]], y = .data[[yv]])) +
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
        if (is.finite(curveparam$ta)) {
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

#' Make volcano plots
#'
#' @param qft A \code{QFeatures} object.
#' @param res A \code{data.frame} object with test results (typically
#'     generated by \code{runTest}).
#' @param testType Character scalar indicating the type of test that was run,
#'     either \code{"ttest"} or \code{"limma"}.
#' @param xv,yv Character scalars indicating which columns of \code{res} that
#'     should be used as the x- and y-axis, respectively.
#' @param volcind Character scalar indicating which column in \code{res} that
#'     represents the "significance" column. This should be a logical
#'     column; rows with a value equal to \code{TRUE} will be colored
#'     in the plot.
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
#' @param stringDb A STRINGdb object or \code{NULL}. If not \code{NULL},
#'     STRING network plots of up- and downregulated genes will be added to
#'     the output pdf file.
#' @param featureCollections A list of \code{CharacterList}s with feature
#'     collections. If there is a collection named \code{"complexes"},
#'     volcano plots and barplots for significant complexes are included in
#'     the output pdf files.
#' @param complexFDRThr Numeric scalar giving the FDR threshold for complexes
#'     to be considered significant.
#' @param curveparam List with curve parameters for creating the Perseus-like
#'     significance curves in the volcano plots.
#'
#' @return A list with two plot objects; one ggplot object and one
#'     interactive ggiraph object. If \code{baseFileName} is not \code{NULL},
#'     pdf files with volcano plots and bar plots for significant
#'     complexes will also be generated.
#'
#' @export
#' @importFrom ggplot2 geom_point aes labs
#' @importFrom ggrepel geom_text_repel
#' @importFrom ggiraph geom_point_interactive
#' @importFrom dplyr filter arrange between row_number desc
#' @importFrom rlang .data
#'
plotVolcano <- function(qft, res, testType, xv, yv, volcind,
                        plotnote, plottitle, plotsubtitle,
                        volcanoFeaturesToLabel, volcanoMaxFeatures,
                        baseFileName, comparisonString, stringDb,
                        featureCollections, complexFDRThr, curveparam) {

    .assertVector(x = qft, type = "QFeatures")
    .assertVector(x = res, type = "data.frame")
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest"))
    .assertScalar(x = xv, type = "character", validValues = colnames(res))
    .assertScalar(x = yv, type = "character", validValues = colnames(res))
    .assertScalar(x = volcind, type = "character",
                  validValues = colnames(res))
    .assertScalar(x = plotnote, type = "character")
    .assertScalar(x = plottitle, type = "character")
    .assertScalar(x = plotsubtitle, type = "character")
    .assertVector(x = volcanoFeaturesToLabel, type = "character")
    .assertScalar(x = volcanoMaxFeatures, type = "numeric",
                  rngIncl = c(0, Inf))
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertScalar(x = comparisonString, type = "character")
    .assertVector(x = stringDb, type = "STRINGdb", allowNULL = TRUE)
    .assertVector(x = featureCollections, type = "list")
    .assertScalar(x = complexFDRThr, type = "numeric", rngIncl = c(0, 1))
    .assertVector(x = curveparam, type = "list")

    ## Make "base" volcano plot
    ggbase <- .makeBaseVolcano(res = res, testType = testType, xv = xv, yv = yv,
                               plotnote = plotnote, plottitle = plottitle,
                               plotsubtitle = plotsubtitle, curveparam = curveparam)

    ## Add labels for significant features
    ggtest <- ggbase +
        ggplot2::geom_point(fill = "lightgrey", color = "grey", pch = 21, size = 1.5) +
        ggplot2::geom_point(data = res %>%
                       dplyr::filter(.data[[volcind]]),
                   fill = "red", color = "grey", pch = 21, size = 1.5) +
        ggrepel::geom_text_repel(
            data = res %>%
                dplyr::filter(
                    .data[[volcind]] |
                        .data$pid %in% volcanoFeaturesToLabel
                ) %>%
                dplyr::arrange(
                    dplyr::desc(abs(.data[[xv]]) + abs(.data[[yv]]))
                ) %>%
                dplyr::filter(dplyr::between(dplyr::row_number(), 0,
                                             volcanoMaxFeatures) |
                                  .data$pid %in% volcanoFeaturesToLabel),
            aes(label = .data$pid), max.overlaps = Inf, size = 4,
            min.segment.length = 0.1)

    ## Interactive version
    ggint <- ggbase +
        ggiraph::geom_point_interactive(
            aes(tooltip = .data$pid), fill = "lightgrey", color = "grey",
            pch = 21, size = 1.5) +
        ggiraph::geom_point_interactive(
            data = res %>% dplyr::filter(.data[[volcind]]),
            aes(tooltip = .data$pid), fill = "red", color = "grey",
            pch = 21, size = 1.5)

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
                dplyr::filter(.data[[volcind]]) %>%
                dplyr::arrange(dplyr::desc(abs(.data[[xv]]) + abs(.data[[yv]]))) %>%
                dplyr::filter(dplyr::between(dplyr::row_number(), 0, volcanoMaxFeatures))
            res0 <- stringDb$map(res0, "IDsForSTRING", removeUnmappedRows = TRUE)
            if (any(res0[[xv]] > 0)) {
                stringDb$plot_network(res0 %>% dplyr::filter(.data[[xv]] > 0) %>%
                                          dplyr::pull("STRING_id"))
            }
            if (any(res0[[xv]] < 0)) {
                stringDb$plot_network(res0 %>% dplyr::filter(.data[[xv]] < 0) %>%
                                          dplyr::pull("STRING_id"))
            }
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
                    print(.complexBarPlot(res = res, prs = prs, qft = qft,
                                          cplx = cplx))
                }
            }
            dev.off()
        }
    }
    return(list(gg = ggtest, ggint = ggiraph::girafe(ggobj = ggint)))
}
