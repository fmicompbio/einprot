#' Compile an R script for launching an adapted iSEE instance
#'
#' The function is intended to be used within the einprot workflows, and
#' assumes that the analysis has been performed as in these.
#'
#' @param iSEEScript Character scalar providing the path where the generated
#'     R script will be saved.
#' @param sceFile Character scalar providing the path to the SCE object that
#'     will be used by the script.
#' @param aName Character scalar providing the base assay name.
#' @param tests Named list with results from statistical tests.
#' @param assayForPlots Character scalar, the assay that should be used for
#'     feature and sample assay plot panels.
#'
#' @return The path to the generated script.
#'
#' @author Charlotte Soneson
#' @export
#'
makeiSEEScript <- function(iSEEScript, sceFile, aName, tests, assayForPlots) {
    .assertScalar(x = iSEEScript, type = "character")
    .assertScalar(x = sceFile, type = "character")
    .assertScalar(x = tools::file_ext(sceFile), type = "character",
                  validValues = "rds")
    .assertScalar(x = aName, type = "character")
    .assertVector(x = tests, type = "list")
    .assertVector(x = names(tests), type = "character")
    .assertScalar(x = assayForPlots, type = "character")

    ## Assemble a script that can be sourced to run iSEE
    ## Load packages, read SCE object and define ECM
    iSEECode <- c(
        "library(iSEE)",
        "library(iSEEu)",
        "library(shiny)",
        paste0("sce <- readRDS('", sceFile, "')"),
        "imp_color_fun <- function(n) {",
        "    structure(c('grey', 'firebrick1'), names = c('TRUE', 'FALSE'))",
        "}",
        "ecm <- ExperimentColorMap(",
        "    assays = list(",
        paste0("        imputed_", aName, " = imp_color_fun"),
        "    )",
        ")",
        "app <- iSEE(sce, colormap = ecm, initial = list("
    )

    ## Add volcano plots
    for (nm in names(tests)) {
        iSEECode <- c(
            iSEECode,
            c(paste0("    VolcanoPlot(XAxis = 'Row data', YAxis = '", nm, ".P.Value',"),
              paste0("                XAxisRowData = '", nm, ".logFC', "),
              "                RowSelectionDynamicSource = TRUE, PanelWidth = 4L),")
        )
    }

    ## Add MA plots
    for (nm in names(tests)) {
        iSEECode <- c(
            iSEECode,
            c(paste0("    MAPlot(XAxis = 'Row data', YAxis = '", nm, ".logFC',"),
              paste0("           XAxisRowData = '", nm, ".AveExpr', "),
              paste0("           PValueField = '", nm, ".P.Value', "),
              "           RowSelectionDynamicSource = TRUE, PanelWidth = 4L),")
        )
    }

    ## Add other panels
    iSEECode <- c(
        iSEECode,
        "    RowDataTable(PanelWidth = 6L, RowSelectionDynamicSource = TRUE),",
        "    FeatureAssayPlot(PanelWidth = 3L, ",
        paste0("                     Assay = '", assayForPlots, "',"),
        "                     XAxis = 'Column data', XAxisColumnData = 'group',",
        "                     ColorBy = 'Feature name',",
        paste0("                     ColorByFeatureNameAssay = 'imputed_", aName, "',"),
        "                     ColorByFeatureSource = 'RowDataTable1',",
        "                     YAxisFeatureSource = 'RowDataTable1', PointSize = 5),",
        "    ReducedDimensionPlot(PanelWidth = 3L, PointSize = 5, ColorBy = 'Column data',",
        "                         ColorByColumnData = 'group'),",
        "    FeatureSetTable(PanelWidth = 8L), ",
        "    ComplexHeatmapPlot(PanelWidth = 4L, ",
        paste0("                       Assay = 'iBAQ', "),
        "                       RowSelectionDynamicSource = TRUE, CustomRows = FALSE,",
        "                       ColumnData = 'group', ShowColumnSelection = FALSE,",
        "                       OrderColumnSelection = FALSE), ",
        "    SampleAssayPlot(PanelWidth = 4L, ",
        paste0("                    Assay = '", assayForPlots, "',"),
        "                    XAxis = 'Sample name'),",
        "    ColumnDataPlot(PanelWidth = 4L, YAxis = 'pNA', XAxis = 'Column data',",
        "                   XAxisColumnData = 'group', ColorBy = 'Column data',",
        "                   PointSize = 5, ColorByColumnData = 'group'),",
        "    RowDataPlot(PanelWidth = 4L, YAxis = 'Score')",
        "))",
        "shiny::runApp(app)"
    )

    writeLines(iSEECode, con = iSEEScript, sep = "\n")
    invisible(iSEEScript)
}
