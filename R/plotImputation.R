#' Plot distribution of imputed and unimputed values
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param assayToPlot Character scalar indicating the name of a numeric
#'     assay of \code{sce} to use for plotting.
#' @param assayImputation Character scalar indicating the name of a
#'     logical assay of \code{sce} to use for filling the distribution plots.
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
plotImputationDistribution <- function(sce, assayToPlot, assayImputation,
                                       xlab = "") {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayToPlot, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = xlab, type = "character")

    plotdf <- as.data.frame(
        SummarizedExperiment::assay(sce, assayToPlot)) %>%
        tibble::rownames_to_column("pid") %>%
        tidyr::gather(key = "sample", value = "log2intensity", -.data$pid) %>%
        dplyr::left_join(
            as.data.frame(
                SummarizedExperiment::assay(sce, assayImputation)) %>%
                tibble::rownames_to_column("pid") %>%
                tidyr::gather(key = "sample", value = "imputed", -.data$pid),
            by = c("pid", "sample")
        )
    ggplot2::ggplot(plotdf, ggplot2::aes(x = .data$log2intensity,
                                         fill = .data$imputed)) +
        ggplot2::geom_histogram(bins = 50) +
        ggplot2::facet_wrap(~ sample) +
        ggplot2::theme_bw() + ggplot2::labs(x = xlab) +
        ggplot2::scale_fill_manual(values = c(`TRUE` = "grey",
                                              `FALSE` = "firebrick1"))
}
