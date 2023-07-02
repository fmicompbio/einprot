#' Make boxplot of intensities
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative). The
#'     \code{colData} must have columns named \code{"sample"} and
#'     \code{"group"}, for grouping and coloring the values, respectively.
#' @param assayName Character scalar, the name of the assay of \code{sce}
#'     to use for the plots.
#' @param doLog Logical scalar, whether to log-transform the y-axis.
#' @param ylab Character scalar, the label to use for the y-axis.
#'
#' @returns A ggplot2 object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
#'                       package = "einprot")
#' samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
#'              "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
#'              "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
#' out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
#'                         includeOnlySamples = samples)
#' sampleAnnot <- data.frame(sample = samples,
#'                           group = gsub("_IP.*", "", samples))
#' sce <- addSampleAnnots(out$sce, sampleAnnot = sampleAnnot)
#' makeIntensityBoxplots(sce, assayName = "iBAQ", doLog = TRUE,
#'                       ylab = "log intensity")
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw scale_y_log10
#'     theme_bw theme labs element_text rel
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment assay colData
#' @importFrom rlang .data
#'
makeIntensityBoxplots <- function(sce, assayName, doLog, ylab) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = doLog, type = "logical")
    .assertScalar(x = ylab, type = "character")
    stopifnot(all(c("sample", "group") %in%
                      colnames(SummarizedExperiment::colData(sce))))

    gg <- ggplot2::ggplot(as.data.frame(
        SummarizedExperiment::assay(sce, assayName)) %>%
            tidyr::gather(key = "col_id", value = "intensity") %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::colData(sce)) %>%
                    tibble::rownames_to_column("col_id"),
                by = "col_id"),
           ggplot2::aes(x = .data$sample, y = .data$intensity,
                        fill = .data$group)) +
        ggplot2::geom_boxplot(alpha = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(angle = 90,
                                                hjust = 1, vjust = 0.5)) +
        ggplot2::labs(x = "", y = ylab)
    if (doLog) {
        gg <- gg + ggplot2::scale_y_log10()
    }
    if (length(unique(sce$group)) > 15) {
        gg <- gg +
            ggplot2::theme(
                legend.text = ggplot2::element_text(size = ggplot2::rel(0.75)))
    }
    gg
}

#' Make mean-vs-SD plot
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param assayName Character scalar, the name of the assay of \code{sce}
#'     to use for the plots.
#' @param xlab,ylab Character scalars, the labels to use for the x/y-axis,
#'     respectively.
#'
#' @returns A ggplot2 object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @examples
#' mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
#'                       package = "einprot")
#' out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.")
#' makeMeanSDPlot(out$sce, assayName = "iBAQ", xlab = "Mean", ylab = "SD")
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw labs
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate ungroup
#' @importFrom rlang .data
#' @importFrom stats sd
#'
makeMeanSDPlot <- function(sce, assayName, xlab = "Mean", ylab = "SD") {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = ylab, type = "character")

    gg <- ggplot2::ggplot(as.data.frame(
        SummarizedExperiment::assay(sce, assayName)) %>%
            tibble::rownames_to_column("pid") %>%
            tidyr::gather(key = "col_id", value = "intensity",
                          -"pid") %>%
            dplyr::group_by(.data$pid) %>%
            dplyr::mutate(mean_intensity = mean(.data$intensity),
                          sd_intensity = stats::sd(.data$intensity)) %>%
            dplyr::ungroup(),
        ggplot2::aes(x = .data$mean_intensity, y = .data$sd_intensity)) +
        ggplot2::geom_point(alpha = 0.05) +
        ggplot2::geom_smooth() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = xlab, y = ylab)

    gg
}

#' Construct SA plot from limma results
#'
#' Given a list of data frames with limma test results, create an SA plot for
#' each contrast.
#'
#' @param testList List of test results, typically generated using
#'     \code{runTest()}.
#'
#' @returns A cowplot object
#'
#' @export
#' @author Charlotte Soneson
#'
#' @importFrom ggplot2 ggplot geom_point geom_line aes labs coord_cartesian
#'     theme_bw
#' @importFrom cowplot plot_grid
#'
makeSAPlot <- function(testList) {
    .assertVector(x = testList, type = "list")

    xrng <- range(unlist(lapply(testList, function(df) df$AveExpr)),
                  na.rm = TRUE)
    yrng <- range(unlist(lapply(testList, function(df) sqrt(df$sigma))),
                  na.rm = TRUE)
    saplots <- lapply(names(testList), function(nm) {
        df <- testList[[nm]]
        ggplot(df) +
            geom_point(aes(x = .data$AveExpr, y = sqrt(.data$sigma)),
                       alpha = 0.5, size = 0.5) +
            geom_line(aes(x = .data$AveExpr, y = sqrt(sqrt(.data$s2.prior))),
                      color = "blue", linewidth = 1) +
            labs(x = "Average log abundance",
                 y = "sqrt(residual standard deviation)",
                 title = nm) +
            coord_cartesian(xlim = xrng, ylim = yrng) +
            theme_bw()
    })
    cowplot::plot_grid(plotlist = saplots, ncol = 3)
}
