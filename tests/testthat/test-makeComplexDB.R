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

    ## Check that species exist in babelgene
    ## As of May 2022, neither C elegans nor S pombe has a 'common name'
    sps <- babelgene::species()
    expect_true(any(grepl("mouse", sps$common_name)))
    expect_true(any(grepl("baker's yeast", sps$common_name)))
    expect_true(any(grepl("Caenorhabditis elegans", sps$scientific_name)))
    expect_true(any(grepl("Schizosaccharomyces pombe 972h-", sps$scientific_name)))
})
