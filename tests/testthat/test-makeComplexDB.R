test_that("makeComplexDB works", {
    expect_error(makeComplexDB(dbDir = 1, customComplexTxt = NULL),
                 "'dbDir' must be of class 'character'")
    expect_error(makeComplexDB(dbDir = c("dir1", "dir2"),
                               customComplexTxt = NULL),
                 "'dbDir' must have length 1")
    expect_error(makeComplexDB(dbDir = "subdir", customComplexTxt = 1),
                 "'customComplexTxt' must be of class 'character'")
    expect_error(makeComplexDB(dbDir = "subdir",
                               customComplexTxt = c("txt1", "txt2")),
                 "'customComplexTxt' must have length 1")
    expect_error(makeComplexDB(dbDir = "subdir", customComplexTxt = "missing"),
                 "file.exists(customComplexTxt) is not TRUE", fixed = TRUE)
})
