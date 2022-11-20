test_that("getUniProtToGeneSymbolMapping works", {
    expect_error(getUniProtToGeneSymbolMapping(1), "Unknown species")
    expect_error(getUniProtToGeneSymbolMapping("wrongSpecies"), "Unknown species")

    expect_error(.getUniProtToGeneSymbolMappingFile("mouse"),
                 "'spi' must be of class 'list'")
    expect_error(.getUniProtToGeneSymbolMappingFile(list(a = "mouse")),
                 "%in% names(spi) is not TRUE", fixed = TRUE)
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("mouse")),
                 "MOUSE_10090_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("human")),
                 "HUMAN_9606_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("roundworm")),
                 "CAEEL_6239_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("zebrafish")),
                 "DANRE_7955_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("fruitfly")),
                 "DROME_7227_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("baker's yeast")),
                 "YEAST_559292_idmapping.dat.gz")
    expect_equal(.getUniProtToGeneSymbolMappingFile(getSpeciesInfo("fission yeast")),
                 "SCHPO_284812_idmapping.dat.gz")
    expect_error(.getUniProtToGeneSymbolMappingFile(list(species = "missing")),
                 "Unsupported species")

    df <- getUniProtToGeneSymbolMapping("fruitfly")
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 2)
    expect_named(df, c("UniProtID", "GeneSymbol"))
})
