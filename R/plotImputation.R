#' Plot distribution of imputed and unimputed values
#'
#' @param qft A \code{QFeatures} object.
#' @param assayToPlot Character scalar indicating the name of a numeric
#'     assay of \code{qft} to use for plotting.
#' @param assayImputation Character scalar indicating the name of a
#'     logical assay of \code{qft} to use for filling the distribution plots.
#' @param iColPattern Character scalar indicating a regular expression to
#'     remove from the sample names in the plot labels.
#' @param xlab Character scalar providing the x-axis label for the plot.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A ggplot object.
#'
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr left_join mutate
#' @importFrom ggplot2 ggplot aes geom_histogram theme_bw facet_wrap labs
#'     scale_fill_manual
#' @importFrom rlang .data
#'
plotImputationDistribution <- function(qft, assayToPlot, assayImputation,
                                       iColPattern = "", xlab = "") {
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = assayToPlot, type = "character",
                  validValues = names(qft))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = names(qft))
    .assertScalar(x = iColPattern, type = "character")
    .assertScalar(x = xlab, type = "character")

    plotdf <- as.data.frame(
        SummarizedExperiment::assay(qft[[assayToPlot]])) %>%
        tibble::rownames_to_column("pid") %>%
        tidyr::gather(key = "sample", value = "log2intensity", -pid) %>%
        dplyr::left_join(
            as.data.frame(
                SummarizedExperiment::assay(qft[[assayImputation]])) %>%
                tibble::rownames_to_column("pid") %>%
                tidyr::gather(key = "sample", value = "imputed", -pid),
            by = c("pid", "sample")
        ) %>%
        dplyr::mutate(sample = sub(iColPattern, "", sample))
    ggplot2::ggplot(plotdf, ggplot2::aes(x = .data$log2intensity,
                                         fill = .data$imputed)) +
        ggplot2::geom_histogram(bins = 50) +
        ggplot2::facet_wrap(~ sample) +
        ggplot2::theme_bw() + ggplot2::labs(x = xlab) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "grey",
                                              `FALSE` = "firebrick1"))
}
