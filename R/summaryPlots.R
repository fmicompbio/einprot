#' @export
#' @author Charlotte Soneson
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_bw scale_y_log10
#'     theme_bw theme labs
#' @importFrom tidyr gather
#' @importFrom dplyr left_join
#' @importFrom tibble rownames_to_column
#' @importFrom SummarizedExperiment assay colData
#'
makeIntensityBoxplots <- function(qft, aName, doLog, ylab) {
    gg <- ggplot(as.data.frame(assay(qft[[aName]])) %>%
               tidyr::gather(key = "col_id", value = "intensity") %>%
               dplyr::left_join(as.data.frame(colData(qft)) %>%
                                    tibble::rownames_to_column("col_id"),
                                by = "col_id"),
           aes(x = sample, y = intensity, fill = group)) +
        geom_boxplot(alpha = 0.5) + theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(x = "", y = ylab)
    if (doLog) {
        gg <- gg + scale_y_log10()
    }
    gg
}

#' @export
#' @author Charlotte Soneson
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_smooth theme_bw labs
#' @importFrom SummarizedExperiment assay
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr gather
#' @importFrom dplyr group_by mutate
makeMeanSDPlot <- function(qft, aName, xlab, ylab) {
    ggplot(as.data.frame(assay(qft[[aName]])) %>%
               tibble::rownames_to_column("pid") %>%
               tidyr::gather(key = "col_id", value = "intensity", -pid) %>%
               dplyr::group_by(pid) %>%
               dplyr::mutate(mean_intensity = mean(intensity),
                             sd_intensity = sd(intensity)),
           aes(x = mean_intensity, y = sd_intensity)) +
        geom_point(alpha = 0.05) + geom_smooth() +
        theme_bw() +
        labs(x = xlab, y = ylab)
}

