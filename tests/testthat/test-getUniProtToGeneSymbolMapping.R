test_that("getUniProtToGeneSymbolMapping works", {
    expect_error(getUniProtToGeneSymbolMapping(1), "Unknown species")
    expect_error(getUniProtToGeneSymbolMapping("wrongSpecies"), "Unknown species")

    df <- getUniProtToGeneSymbolMapping("fruitfly")
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 2)
    expect_named(df, c("UniProtID", "GeneSymbol"))
})
