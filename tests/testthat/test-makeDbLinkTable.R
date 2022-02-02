test_that("making the link table works", {
    expect_error(.makeLinkFromId(1))
    expect_error(.makeLinkFromId("Q7YTG1", "WrongType"))
    expect_error(.makeLinkFromId(c("Q7YTG1", "SPBC460.01c")))

    expect_equal(.makeLinkFromId(""), "")
    expect_equal(.makeLinkFromId("Q7YTG1", "UniProt"),
                 '<a href="https://www.uniprot.org/uniprot/Q7YTG1" target="_blank"> Q7YTG1</a>')
    expect_equal(.makeLinkFromId("Q7YTG1", "AlphaFold"),
                 '<a href="https://alphafold.ebi.ac.uk/entry/Q7YTG1" target="_blank"> Q7YTG1</a>')
    expect_equal(.makeLinkFromId("SPBC460.01c", "PomBase"),
                 '<a href="https://www.pombase.org/gene/SPBC460.01c" target="_blank"> SPBC460.01c</a>')
    expect_equal(.makeLinkFromId("WBGene00001330", "WormBase"),
                 '<a href="https://wormbase.org/species/c_elegans/gene/WBGene00001330" target="_blank"> WBGene00001330</a>')

    expect_error(getConvTable(type = "WrongType"))
    expect_error(getConvTable(type = c("Pombase", "Wormbase")))
    expect_error(getConvTable(type = 1))

    ## Without species-specific columns
    dblt <- makeDbLinkTable(data.frame(id = c("B5BP45", "O13282")), idCol = "id",
                            speciesCommon = "fission yeast")
    expect_s3_class(dblt, "data.frame")
    expect_equal(ncol(dblt), 3)
    expect_equal(nrow(dblt), 2)
    expect_named(dblt, c("id", "UniProt", "AlphaFold"))
    expect_equal(dblt$id, c("B5BP45", "O13282"))

    ## With Pombase column
    dblt2 <- makeDbLinkTable(
        data.frame(id = c("B5BP45", "O13282")), idCol = "id",
        speciesCommon = "fission yeast",
        convTablePombase = data.frame(PomBaseID = c("SPBC460.01c", "SPCC5E4.03c"),
                                      UniProtID = c("B5BP45", "O13282")))
    expect_s3_class(dblt2, "data.frame")
    expect_equal(ncol(dblt2), 4)
    expect_equal(nrow(dblt2), 2)
    expect_named(dblt2, c("id", "UniProt", "AlphaFold", "PomBase"))

    ## With Pombase conversion table, but ignore
    dblt3 <- makeDbLinkTable(
        data.frame(id = c("B5BP45", "O13282")), idCol = "id",
        speciesCommon = "fission yeast", addSpeciesSpecificColumns = FALSE,
        convTablePombase = data.frame(PomBaseID = c("SPBC460.01c", "SPCC5E4.03c"),
                                      UniProtID = c("B5BP45", "O13282")))
    expect_s3_class(dblt3, "data.frame")
    expect_equal(ncol(dblt3), 3)
    expect_equal(nrow(dblt3), 2)
    expect_named(dblt3, c("id", "UniProt", "AlphaFold"))

    ## With Wormbase column
    dblt4 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm",
        convTableWormbase = data.frame(UniProtID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
                                       UniProtKB.ID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
                                       WormBaseID = c("WBGene00001330", "WBGene00001328",
                                                      "WBGene00001330", "WBGene00001328"),
                                       check.names = FALSE))
    expect_s3_class(dblt4, "data.frame")
    expect_equal(ncol(dblt4), 5)
    expect_equal(nrow(dblt4), 2)
    expect_named(dblt4, c("gid", "pid", "UniProt", "AlphaFold", "WormBase"))
    expect_equal(grep(";", dblt4$UniProt), c(1, 2))
    expect_equal(grep(";", dblt4$AlphaFold), c(1, 2))
    expect_equal(grep(";", dblt4$WormBase), integer(0))

    ## With Wormbase conversion table, but ignore
    dblt5 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm", addSpeciesSpecificColumns = FALSE,
        convTableWormbase = data.frame(UniProtID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
                                       UniProtKB.ID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
                                       WormBaseID = c("WBGene00001330", "WBGene00001328",
                                                      "WBGene00001330", "WBGene00001328"),
                                       check.names = FALSE))
    expect_s3_class(dblt5, "data.frame")
    expect_equal(ncol(dblt5), 4)
    expect_equal(nrow(dblt5), 2)
    expect_named(dblt5, c("gid", "pid", "UniProt", "AlphaFold"))
    expect_equal(grep(";", dblt5$UniProt), c(1, 2))
    expect_equal(grep(";", dblt5$AlphaFold), c(1, 2))

    ## Skip the rest of the tests if no internet connection is available
    skip_if_offline()
    dfp <- getConvTable(type = "Pombase")
    expect_equal(ncol(dfp), 2)
    expect_named(dfp, c("PomBaseID", "UniProtID"))

    dfw <- getConvTable(type = "Wormbase")
    expect_equal(ncol(dfw), 3)
    expect_named(dfw, c("UniProtID", "UniProtKB.ID", "WormBaseID"))

    ## Don't add more tests below here unless they should be skipped if offline
})
