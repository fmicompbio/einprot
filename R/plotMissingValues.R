#' Plot heatmap of missing values
#'
#' @param qft A \code{QFeatures} object.
#' @param assayMissing Character scalar indicating the name of a
#'     logical assay of \code{qft} representing the missingness pattern.
#'     "FALSE" entries should represent present values, while
#'     "TRUE" entries represent missing values.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A ComplexHeatmap object.
#'
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom SummarizedExperiment assay
#'
plotMissingValuesHeatmap <- function(qft, assayMissing) {
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = assayMissing, type = "character",
                  validValues = names(qft))

    col_fun = circlize::colorRamp2(c(0, 1), c("grey50", "white"))
    ComplexHeatmap::Heatmap(
        SummarizedExperiment::assay(qft[[paste0("imputed_", aName)]]) + 0,
        col = col_fun, name = "imputed",
        column_title = "Missing value pattern (white = missing)",
        cluster_rows = TRUE, cluster_columns = TRUE, show_row_names = FALSE,
        show_heatmap_legend = FALSE)
}

#' Plot the fraction of detected features in each sample
#'
#' @param dfNA A \code{DFrame} with columns "assay", "name", "nNA" and
#'     "pNA" (such as those returned by \code{QFeatures::nNA}).
#' @param aName Character scalar indicating the assay for which to plot
#'     the fraction of detected features.
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
plotFractionDetectedPerSample <- function(dfNA, aName) {
    .assertVector(x = dfNA, type = "DFrame")
    .assertScalar(x = aName, type = "character",
                  validValues = unique(dfNA$assay))
    stopifnot(all(c("name", "pNA") %in% colnames(dfNA)))

    ggplot2::ggplot(
        as.data.frame(dfNA) %>%
            dplyr::filter(.data$assay == aName),
        ggplot2::aes(x = .data$name, y = 100 - .data$pNA,
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
#' @param dfNA A \code{DFrame} with columns "assay", "name", "nNA" and
#'     "pNA" (such as those returned by \code{QFeatures::nNA}).
#' @param aName Character scalar indicating the assay for which to plot
#'     the the detection patterns.
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
plotDetectedInSamples <- function(dfNA, aName) {
    .assertVector(x = dfNA, type = "DFrame")
    .assertScalar(x = aName, type = "character",
                  validValues = unique(dfNA$assay))
    stopifnot(all(c("name", "pNA", "nNA") %in% colnames(dfNA)))

    ## Get the total number of samples
    totN <- as.data.frame(dfNA) %>%
        dplyr::filter(.data$assay == aName) %>%
        dplyr::mutate(totN = round(.data$nNA/.data$pNA * 100)) %>%
        dplyr::filter(!is.na(.data$totN)) %>%
        dplyr::pull(totN) %>%
        unique()
    stopifnot(length(totN) == 1)
    ggplot2::ggplot(as.data.frame(dfNA) %>%
                        dplyr::filter(.data$assay == aName) %>%
                        dplyr::group_by(.data$nNA) %>%
                        dplyr::tally(),
           ggplot2::aes(x = totN - .data$nNA, y = .data$n)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::theme_bw() +
        ggplot2::labs(x = "Number of samples in which a feature is detected",
                      y = "Number of features")
}
