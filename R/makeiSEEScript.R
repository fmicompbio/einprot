#' Compile an R script for launching an adapted iSEE instance
#'
#' @param iSEEScript Name of R script
#' @param sceFile Path to SCE object
#' @param aname Base assay name
#' @param tests Test results
#' @param assayForTests Assay used for tests
#'
#' @author Charlotte Soneson
#' @export
#'
makeiSEEScript <- function(iSEEScript, sceFile, aName, tests, assayForTests) {
    ## Assemble a script that can be sourced to run iSEE
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
    for (nm in names(tests)) {
        iSEECode <- c(
            iSEECode,
            c(paste0("    VolcanoPlot(XAxis = 'Row data', YAxis = '", nm, ".P.Value',"),
              paste0("                XAxisRowData = '", nm, ".logFC', "),
              "                RowSelectionDynamicSource = TRUE, PanelWidth = 4L),")
        )
    }
    for (nm in names(tests)) {
        iSEECode <- c(
            iSEECode,
            c(paste0("    MAPlot(XAxis = 'Row data', YAxis = '", nm, ".logFC',"),
              paste0("           XAxisRowData = '", nm, ".AveExpr', "),
              paste0("           PValueField = '", nm, ".P.Value', "),
              "           RowSelectionDynamicSource = TRUE, PanelWidth = 4L),")
        )
    }
    iSEECode <- c(
        iSEECode,
        "    RowDataTable(PanelWidth = 6L, RowSelectionDynamicSource = TRUE),",
        "    FeatureAssayPlot(PanelWidth = 3L, ",
        paste0("                     Assay = '", assayForTests, "',"),
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
        paste0("                    Assay = '", assayForTests, "',"),
        "                    XAxis = 'Sample name'),",
        "    ColumnDataPlot(PanelWidth = 4L, YAxis = 'pNA', XAxis = 'Column data',",
        "                   XAxisColumnData = 'group', ColorBy = 'Column data',",
        "                   PointSize = 5, ColorByColumnData = 'group'),",
        "    RowDataPlot(PanelWidth = 4L, YAxis = 'Score')",
        "))",
        "shiny::runApp(app)"
    )

    writeLines(iSEECode, con = iSEEScript, sep = "\n")
}
