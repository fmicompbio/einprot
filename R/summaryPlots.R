#' Make boxplot of intensities
#'
#' @param qft QFeatures object.
#' @param assayName Character scalar, the name of the assay of \code{qft}
#'     to use for the plots.
#' @param doLog Logical scalar, whether to log-transform the y-axis.
#' @param ylab Character scalar, the label of the y-axis.
#'
#' @return A ggplot2 object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw scale_y_log10
#'     theme_bw theme labs
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment assay colData
#' @importFrom rlang .data
#'
makeIntensityBoxplots <- function(qft, assayName, doLog, ylab) {
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = assayName, type = "character")
    .assertScalar(x = doLog, type = "logical")
    .assertScalar(x = ylab, type = "character")

    gg <- ggplot2::ggplot(as.data.frame(
        SummarizedExperiment::assay(qft[[assayName]])) %>%
            tidyr::gather(key = "col_id", value = "intensity") %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::colData(qft)) %>%
                    tibble::rownames_to_column("col_id"),
                by = "col_id"),
           ggplot2::aes(x = .data$sample, y = .data$intensity,
                        fill = .data$group)) +
        ggplot2::geom_boxplot(alpha = 0.5) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = element_text(angle = 90,
                                                  hjust = 1, vjust = 0.5)) +
        ggplot2::labs(x = "", y = ylab)
    if (doLog) {
        gg <- gg + ggplot2::scale_y_log10()
    }
    gg
}

#' Make mean-vs-SD plot
#'
#' @param qft QFeatures object.
#' @param assayName Character scalar, the name of the assay of \code{qft}
#'     to use for the plots.
#' @param xlab,ylab Character scalars, the labels of the x/y-axis,
#'     respectively.
#'
#' @return A ggplot2 object.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw labs
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate
#' @importFrom rlang .data
#'
makeMeanSDPlot <- function(qft, assayName, xlab, ylab) {
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = assayName, type = "character")
    .assertScalar(x = xlab, type = "character")
    .assertScalar(x = ylab, type = "character")

    gg <- ggplot2::ggplot(as.data.frame(
        SummarizedExperiment::assay(qft[[assayName]])) %>%
            tibble::rownames_to_column("pid") %>%
            tidyr::gather(key = "col_id", value = "intensity",
                          -.data$pid) %>%
            dplyr::group_by(.data$pid) %>%
            dplyr::mutate(mean_intensity = mean(.data$intensity),
                          sd_intensity = sd(.data$intensity)),
        ggplot2::aes(x = .data$mean_intensity, y = .data$sd_intensity)) +
        ggplot2::geom_point(alpha = 0.05) +
        ggplot2::geom_smooth() +
        ggplot2::theme_bw() +
        ggplot2::labs(x = xlab, y = ylab)

    gg
}

