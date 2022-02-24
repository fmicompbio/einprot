test_that("creating an iSEE script works", {
    expect_error(makeiSEEScript(iSEEScript = 1, sceFile = "file.rds",
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "'iSEEScript' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = 1,
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "'sceFile' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = c("file1.rds", "file2.rds"),
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "'sceFile' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file.R",
                                aName = "Intensity", tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "All values in 'tools::file_extsceFile' must be one of: rds",
                 fixed = TRUE)
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = 1, tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "'aName' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = c("Intensity", "iBAQ"),
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = "Intensity"),
                 "'aName' must have length 1")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = 1,
                                assayForPlots = "Intensity"),
                 "'tests' must be of class 'list'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(1, 2),
                                assayForPlots = "Intensity"),
                 "'namestests' must not be NULL")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = 1),
                 "'assayForPlots' must be of class 'character'")
    expect_error(makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                                sceFile = "file1.rds",
                                aName = "Intensity",
                                tests = list(t1 = 1, t2 = 2),
                                assayForPlots = c("Intensity", "iBAQ")),
                 "'assayForPlots' must have length 1")

    iss <- makeiSEEScript(iSEEScript = tempfile(fileext = ".R"),
                          sceFile = "file1.rds",
                          aName = "Intensity",
                          tests = list(t1 = 1, t2 = 2),
                          assayForPlots = "Intensity")
    expect_type(iss, "character")
    expect_true(file.exists(iss))
})
