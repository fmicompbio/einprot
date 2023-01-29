test_that("creating an iSEE script works", {
    expect_error(makeiSEEScript(iSEEScript = 1, sceFile = "file.rds",
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'iSEEScript' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = 1,
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'sceFile' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = c("file1.rds", "file2.rds"),
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'sceFile' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file.R",
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "All values in 'tools::file_extsceFile' must be one of: rds",
                 fixed = TRUE)
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = 1, tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'aName' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = c("Intensity", "iBAQ"),
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'aName' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = 1,
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'tests' must be of class 'list'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(1, 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'namestests' must not be NULL")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = 1,
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'assayForPlots' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = c("Intensity", "iBAQ"),
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = TRUE),
                 "'assayForPlots' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = 1,
                                includeFeatureSetTable = TRUE),
                 "'assayForHeatmaps' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = c("Intensity", "iBAQ"),
                                includeFeatureSetTable = TRUE),
                 "'assayForHeatmaps' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = 1),
                 "'includeFeatureSetTable' must be of class 'logical'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity",
                                assayForHeatmaps = "Intensity",
                                includeFeatureSetTable = c(TRUE, FALSE)),
                 "'includeFeatureSetTable' must have length 1")

    expect_warning({
        iss <- makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                              sceFile = "file1.rds",
                              aName = "Intensity",
                              tests = list(t1 = 1, t2 = 2),
                              assayForPlots = "Intensity",
                              assayForHeatmaps = "Intensity",
                              includeFeatureSetTable = TRUE)
    }, "cannot open compressed file")
    expect_type(iss, "character")
    expect_true(file.exists(iss))
    rl <- readLines(iss)
    expect_true(any(grepl("FeatureSetTable", rl)))

    expect_warning({
        iss <- makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                              sceFile = "file1.rds",
                              aName = "Intensity",
                              tests = list(t1 = 1, t2 = 2),
                              assayForPlots = "Intensity",
                              assayForHeatmaps = "Intensity",
                              includeFeatureSetTable = FALSE)
    }, "cannot open compressed file")
    expect_type(iss, "character")
    expect_true(file.exists(iss))
    rl <- readLines(iss)
    expect_false(any(grepl("FeatureSetTable", rl)))
    expect_true(any(grepl("MAPlot", rl)))
    expect_true(any(grepl("VolcanoPlot", rl)))

    expect_warning({
        iss <- makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                              sceFile = "file1.rds",
                              aName = "Intensity",
                              tests = list(),
                              assayForPlots = "Intensity",
                              assayForHeatmaps = "Intensity",
                              includeFeatureSetTable = FALSE)
    }, "cannot open compressed file")
    expect_type(iss, "character")
    expect_true(file.exists(iss))
    rl <- readLines(iss)
    expect_false(any(grepl("FeatureSetTable", rl)))
    expect_false(any(grepl("MAPlot", rl)))
    expect_false(any(grepl("VolcanoPlot", rl)))

})
