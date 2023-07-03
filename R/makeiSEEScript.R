#' Compile an R script for launching an adapted \code{iSEE} instance
#'
#' The function is intended to be used within the \code{einprot} workflows, and
#' assumes that the analysis has been performed as in these.
#'
#' @param iSEEScript Character scalar providing the path where the generated
#'     R script will be saved.
#' @param sceFile Character scalar providing the path to the
#'     \code{SingleCellExperiment} object that will be used by the script.
#' @param aName Character scalar providing the base assay name.
#' @param tests Named list with results from statistical tests.
#' @param assayForPlots Character scalar, the assay that should be used for
#'     feature and sample assay plot panels.
#' @param assayForHeatmaps Character scalar, the assay that should be used for
#'     heatmap panels.
#' @param includeFeatureSetTable Logical scalar, whether to include a
#'     feature set table panel.
#'
#' @returns The path to the generated script.
#'
#' @examples
#' sceFile <- system.file("extdata", "mq_example", "1356_sce.rds",
#'                        package = "einprot")
#' sce <- readRDS(sceFile)
#' isee <- makeiSEEScript(tempfile(fileext = "_isee.R"), sceFile = sceFile,
#'                        aName = "LFQ.intensity",
#'                        tests = S4Vectors::metadata(sce)$testres$tests,
#'                        assayForPlots = "log2_LFQ.intensity",
#'                        assayForHeatmaps = "iBAQ",
#'                        includeFeatureSetTable = FALSE)
#' file.exists(isee)
#' head(readLines(isee), 10)
#'
#' @author Charlotte Soneson
#' @export
#'
#' @importFrom utils read.csv
#'
makeiSEEScript <- function(iSEEScript, sceFile, aName, tests, assayForPlots,
                           assayForHeatmaps, includeFeatureSetTable) {
    .assertScalar(x = iSEEScript, type = "character")
    .assertScalar(x = sceFile, type = "character")
    .assertScalar(x = tools::file_ext(sceFile), type = "character",
                  validValues = "rds")
    .assertScalar(x = aName, type = "character")
    .assertVector(x = tests, type = "list")
    if (length(tests) > 0) {
        .assertVector(x = names(tests), type = "character")
    }
    .assertScalar(x = assayForPlots, type = "character")
    .assertScalar(x = assayForHeatmaps, type = "character")
    .assertScalar(x = includeFeatureSetTable, type = "logical")

    tour <- utils::read.csv(system.file("extdata", "iSEEtour.csv",
                                        package = "einprot"))
    if (!includeFeatureSetTable) {
        tour <- tour[!tour$element %in% c("#FeatureSetTable1",
                                          "#ComplexHeatmapPlot1"), ]
    }
    if (length(tests) == 0) {
        tour <- tour[!tour$element %in% c("#VolcanoPlot1",
                                          "#VolcanoPlot1_VisualBoxOpen",
                                          "#MAPlot1"), ]
    }
    tourFile <- sub("\\.R$", "_tour.csv", iSEEScript)
    write.table(tour, file = tourFile,
                sep = ",", row.names = FALSE, col.names = TRUE,
                quote = TRUE)

    snames <- tryCatch({
        sce <- readRDS(sceFile)
        c(xn = colnames(sce)[1], yn = colnames(sce)[2])
    }, error = function(e) {
        c(xn = NA, yn = NA)
    })

    ## Assemble a script that can be sourced to run iSEE
    ## Load packages, read SCE object and define ECM
    iSEECode <- c(
        "library(iSEE)",
        "library(iSEEu)",
        "library(shiny)",
        paste0("sce <- readRDS('", sceFile, "')"),
        paste0("tour <- read.csv('", tourFile, "')"),
        "panelDefaults(TooltipRowData = c('einprotLabel'))",
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
            c(paste0("    VolcanoPlot(XAxis = 'Row data', YAxis = '", nm,
                     ".P.Value',"),
              paste0("                XAxisRowData = '", nm, ".logFC', "),
              paste0("                RowSelectionDynamicSource = TRUE, ",
                     "PanelWidth = 4L),"))
        )
    }

    ## Add MA plots
    for (nm in names(tests)) {
        iSEECode <- c(
            iSEECode,
            c(paste0("    MAPlot(XAxis = 'Row data', YAxis = '", nm,
                     ".logFC',"),
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
        if (includeFeatureSetTable) "    FeatureSetTable(PanelWidth = 8L), ",
        "    ComplexHeatmapPlot(PanelWidth = 4L, ",
        paste0("                       Assay = '", assayForHeatmaps, "', "),
        "                       RowSelectionDynamicSource = TRUE, CustomRows = FALSE,",
        "                       ColumnData = 'group', ShowColumnSelection = FALSE,",
        "                       OrderColumnSelection = FALSE), ",
        "    SampleAssayPlot(PanelWidth = 4L, ",
        paste0("                    Assay = '", assayForPlots, "',"),
        "                    XAxis = 'Sample name',",
        paste0("                    XAxisSampleName = '", snames["xn"], "',"),
        paste0("                    YAxisSampleName = '", snames["yn"], "'),"),
        "    ColumnDataPlot(PanelWidth = 4L, YAxis = 'pNA', XAxis = 'Column data',",
        "                   XAxisColumnData = 'group', ColorBy = 'Column data',",
        "                   PointSize = 5, ColorByColumnData = 'group'),",
        "    RowDataPlot(PanelWidth = 4L, YAxis = 'Score')",
        "), tour = tour)",
        "shiny::runApp(app)"
    )

    writeLines(iSEECode, con = iSEEScript, sep = "\n")
    invisible(iSEEScript)
}
