#' Plot heatmap of missing values
#'
#' @param sce A \code{SummarizedExperiment} object.
#' @param assayMissing Character scalar indicating the name of a
#'     logical assay of \code{sce} representing the missingness pattern.
#'     "FALSE" entries should represent present values, while
#'     "TRUE" entries represent missing values.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' sce <- importExperiment(system.file("extdata", "mq_example",
#'                                     "1356_proteinGroups.txt",
#'                                     package = "einprot"),
#'                         iColPattern = "^iBAQ\\.")$sce
#' SummarizedExperiment::assay(sce, "iBAQ")[
#'     SummarizedExperiment::assay(sce, "iBAQ") == 0] <- NA
#' SummarizedExperiment::assay(sce, "missing") <-
#'     is.na(SummarizedExperiment::assay(sce, "iBAQ"))
#' plotMissingValuesHeatmap(sce, "missing")
#'
#' @return A ComplexHeatmap object.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom SummarizedExperiment assay assayNames
#'
plotMissingValuesHeatmap <- function(sce, assayMissing) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayMissing, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    if (any(is.na(SummarizedExperiment::assay(sce, assayMissing)))) {
        stop("Assay contains missing values")
    }

    col_fun = circlize::colorRamp2(c(0, 1), c("grey50", "white"))
    ComplexHeatmap::Heatmap(
        SummarizedExperiment::assay(sce, assayMissing) + 0,
        col = col_fun, name = "imputed",
        column_title = "Missing value pattern (white = missing)",
        cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE,
        show_heatmap_legend = FALSE)
}

#' Plot the fraction of detected features in each sample
#'
#' @param dfNA A \code{DFrame} with at least columns named "name" and
#'     "pNA" representing the sample name and the fraction of missing values
#'     (e.g., returned by \code{QFeatures::nNA}).
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes theme geom_bar theme_bw labs expand_limits
#'     geom_text
#' @importFrom dplyr filter %>%
#' @importFrom rlang .data
#'
plotFractionDetectedPerSample <- function(dfNA) {
    .assertVector(x = dfNA, type = "data.frame")
    stopifnot(all(c("sample", "pNA") %in% colnames(dfNA)))

    ggplot2::ggplot(dfNA,
        ggplot2::aes(x = .data$sample, y = 100 - .data$pNA,
                     label = paste0(round(100 - .data$pNA, 1), "%"))) +
        ggplot2::geom_bar(stat = "identity") + ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                                  vjust = 0.5)) +
        ggplot2::labs(x = "", y = "Detected features (%)") +
        ggplot2::expand_limits(y = 100) +
        ggplot2::geom_text(vjust = 1.5, color = "white", size = 3)
}

#' Plot the distribution of the number of samples where features are detected
#'
#' @param dfNA A \code{DFrame} with at least columns "name", "nNA" and
#'     "pNA", representing the feature name and the number and fraction of
#'     missing values (e.g., returned by \code{QFeatures::nNA}).
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A ggplot object.
#'
#' @importFrom ggplot2 ggplot aes geom_bar labs
#' @importFrom rlang .data
#' @importFrom dplyr filter %>% group_by tally pull mutate
#'
plotDetectedInSamples <- function(dfNA) {
    .assertVector(x = dfNA, type = "data.frame")
    stopifnot(all(c("pNA", "nNA") %in% colnames(dfNA)))

    ## Get the total number of samples
    totN <- dfNA %>%
        dplyr::mutate(totN = round(.data$nNA/.data$pNA * 100)) %>%
        dplyr::filter(!is.na(.data$totN)) %>%
        dplyr::pull(totN) %>%
        unique()
    stopifnot(length(totN) == 1)

    ## Count number of features observed in a given number of samples
    plotdf <- dfNA %>%
        dplyr::group_by(.data$nNA) %>%
        dplyr::tally() %>%
        dplyr::mutate(nObs = totN - .data$nNA)

    ## Expand with zeros if necessary
    missingN <- setdiff(c(0, seq_len(totN)), plotdf$nObs)
    if (length(missingN) > 0) {
        plotdf <- plotdf %>%
            dplyr::bind_rows(
                data.frame(nObs = missingN,
                           n = 0) %>%
                    dplyr::mutate(nNA = totN - .data$nObs)
            )
    }
    plotdf <- plotdf %>%
        dplyr::mutate(nObs = factor(.data$nObs, levels = c(0, seq_len(totN))))
    ggplot2::ggplot(plotdf,
           ggplot2::aes(x = .data$nObs, y = .data$n)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Number of samples in which a feature is detected",
                      y = "Number of features")
}
