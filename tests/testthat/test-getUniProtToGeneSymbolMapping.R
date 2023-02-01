test_that("getUniProtToIDMapping works", {
    expect_error(getUniProtToIDMapping(1), "Unknown species")
    expect_error(getUniProtToIDMapping("wrongSpecies"), "Unknown species")
    expect_error(getUniProtToIDMapping("fruitfly", targetId = 1),
                 "'targetId' must be of class 'character'")

    expect_error(.getUniProtToIDMappingFile("mouse"),
                 "'spi' must be of class 'list'")
    expect_error(.getUniProtToIDMappingFile(list(a = "mouse")),
                 "%in% names(spi) is not TRUE", fixed = TRUE)
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("mouse")),
                 "MOUSE_10090_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("human")),
                 "HUMAN_9606_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("roundworm")),
                 "CAEEL_6239_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("zebrafish")),
                 "DANRE_7955_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("fruitfly")),
                 "DROME_7227_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("baker's yeast")),
                 "YEAST_559292_idmapping.dat.gz")
    expect_equal(.getUniProtToIDMappingFile(getSpeciesInfo("fission yeast")),
                 "SCHPO_284812_idmapping.dat.gz")
    expect_error(.getUniProtToIDMappingFile(list(species = "missing")),
                 "Unsupported species")

    df <- getUniProtToIDMapping("fruitfly")
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 2)
    expect_named(df, c("UniProtID", "Gene_Name"))

    df <- getUniProtToIDMapping("fruitfly", targetId = "RefSeq_NT")
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 2)
    expect_named(df, c("UniProtID", "RefSeq_NT"))
})
