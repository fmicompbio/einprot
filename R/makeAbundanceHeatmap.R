#' @noRd
#' @keywords internal
.gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[seq_len(n)]
}

#' Generate heatmap
#'
#' Generate a heatmap from a defined assay.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @param sce A `SummarizedExperiment` object (or a derivative).
#' @param assayToPlot Character scalar giving the name of the assay in
#'     \code{sce} to plot.
#' @param doCenter Logical scalar, whether to center the abundance values by
#'     row before creating the heatmap.
#' @param settings Character scalar or \code{NULL}. Setting this to either
#'     `"report"` or `"export"` creates heatmaps with specific settings used
#'     in `einprot` reports and when exporting the heatmap to a pdf. Setting
#'     it to \code{NULL} allows any argument to be passed to
#'     \code{ComplexHeatmap::Heatmap} via the `...` argument.
#' @param ... If \code{settings} is \code{NULL}, additional arguments passed to
#'     \code{ComplexHeatmap::Heatmap}.
#'
#' @return A ComplexHeatmap object.
#'
#' @importFrom ComplexHeatmap Heatmap columnAnnotation
#' @importFrom SummarizedExperiment assay assayNames
#' @importFrom grid gpar
#' @importFrom circlize colorRamp2
#'
makeAbundanceHeatmap <- function(sce, assayToPlot, doCenter,
                                 settings = "report", ...) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayToPlot, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = doCenter, type = "logical")
    .assertScalar(x = settings, type = "character", allowNULL = TRUE,
                  validValues = c("report", "export"))

    if (!"pNA" %in% colnames(SummarizedExperiment::rowData(sce))) {
        SummarizedExperiment::rowData(sce)$pNA <- NA_real_
    }
    pNAcol <- circlize::colorRamp2(c(0, 100), c("white", "darkviolet"))

    groupcols = .gg_color_hue(length(unique(sce$group)))
    names(groupcols) <- levels(factor(sce$group))

    if (doCenter) {
        mat <- t(scale(t(SummarizedExperiment::assay(sce, assayToPlot)),
                       center = TRUE, scale = FALSE))
        nm <- paste0(assayToPlot, "\ncentered")
    } else {
        mat <- SummarizedExperiment::assay(sce, assayToPlot)
        nm <- assayToPlot
    }
    if (!is.null(settings) && settings == "report") {
        ComplexHeatmap::Heatmap(
            mat,
            name = nm,
            show_row_names = FALSE,
            use_raster = TRUE,
            top_annotation = ComplexHeatmap::columnAnnotation(
                group = sce$group,
                col = list(group = groupcols))
        )
    } else if (!is.null(settings) && settings == "export") {
        ComplexHeatmap::Heatmap(
            mat,
            name = nm,
            show_row_names = TRUE, show_row_dend = FALSE,
            cluster_columns = TRUE, column_split = sce$group,
            use_raster = TRUE,
            top_annotation = ComplexHeatmap::columnAnnotation(
                group = sce$group,
                annotation_name_side = "left",
                col = list(group = groupcols)),
            bottom_annotation = ComplexHeatmap::columnAnnotation(
                group = sce$group,
                annotation_name_side = "left",
                col = list(group = groupcols)),
            right_annotation = ComplexHeatmap::rowAnnotation(
                pNA = SummarizedExperiment::rowData(sce)$pNA,
                col = list(pNA = pNAcol)),
            heatmap_legend_param = list(direction = "horizontal"),
            column_title_gp = grid::gpar(fontsize = 6)
        )
    } else if (is.null(settings)) {
        ComplexHeatmap::Heatmap(
            mat, name = nm, ...
        )
    } else {
        stop("Unknown value of the settings parameter")
    }

}

