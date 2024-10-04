#' Plot distribution of imputed and unimputed values
#'
#' Create a plot showing the distribution of both imputed and observed
#' abundance values contained in a \code{SummarizedExperiment} object.
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param assayToPlot Character scalar indicating the name of a numeric
#'     assay of \code{sce} to use for plotting.
#' @param assayImputation Character scalar indicating the name of a
#'     logical assay of \code{sce} to use for filling the distribution plots.
#' @param xlab Character scalar providing the x-axis label for the plot.
#' @param plotType Character scalar indicating the type of plot to make
#'     (either "histogram" or "density").
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A ggplot object.
#'
#' @examples
#' sce <- readRDS(system.file("extdata", "mq_example", "1356_sce.rds",
#'                            package = "einprot"))
#' plotImputationDistribution(sce, assayToPlot = "log2_LFQ.intensity",
#'                            assayImputation = "imputed_LFQ.intensity")
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
                                       xlab = "", plotType = "histogram") {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayToPlot, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = plotType, type = "character",
                  validValues = c("histogram", "density"))

    plotdf <- as.data.frame(
        SummarizedExperiment::assay(sce, assayToPlot)) %>%
        tibble::rownames_to_column("pid") %>%
        tidyr::gather(key = "sample", value = "log2intensity", -"pid") %>%
        dplyr::left_join(
            as.data.frame(
                SummarizedExperiment::assay(sce, assayImputation)) %>%
                tibble::rownames_to_column("pid") %>%
                tidyr::gather(key = "sample", value = "imputed", -"pid"),
            by = c("pid", "sample")
        )
    if (plotType == "histogram") {
        ggplot2::ggplot(plotdf, ggplot2::aes(x = .data$log2intensity,
                                             fill = .data$imputed)) +
            ggplot2::geom_histogram(bins = 50) +
            ggplot2::facet_wrap(~ sample) +
            ggplot2::theme_bw() + ggplot2::labs(x = xlab) +
            ggplot2::scale_fill_manual(values = c(`TRUE` = "grey",
                                                  `FALSE` = "firebrick1"))
    } else if (plotType == "density") {
        ggplot2::ggplot(plotdf, ggplot2::aes(x = .data$log2intensity,
                                             color = .data$imputed)) +
            ggplot2::geom_density(linewidth = 1.5) +
            ggplot2::facet_wrap(~ sample) +
            ggplot2::theme_bw() + ggplot2::labs(x = xlab) +
            ggplot2::scale_color_manual(values = c(`TRUE` = "grey",
                                                   `FALSE` = "firebrick1"))
    } else {
        ## Should never end up here as the parameter is checked above
        #nocov start
        stop("Unknown value of the plotType parameter")
        #nocov end
    }
}
