test_that("making the link table works", {
    ## .makeLinkFromId
    ## --------------------------------------------------------------------- ##
    expect_error(.makeLinkFromId(1),
                 "'id' must be of class 'character'")
    expect_error(.makeLinkFromId(c("Q7YTG1", "SPBC460.01c")),
                 "'id' must have length 1")
    expect_error(.makeLinkFromId("Q7YTG1", 1),
                 "'linktype' must be of class 'character'")
    expect_error(.makeLinkFromId("Q7YTG1", c("AlphaFold", "UniProt")),
                 "'linktype' must have length 1")
    expect_error(.makeLinkFromId("Q7YTG1", "WrongType"),
                 "All values in 'linktype' must be one of")

    expect_equal(.makeLinkFromId(""), "")
    expect_equal(.makeLinkFromId("Q7YTG1", "UniProt"),
                 '<a href="https://www.uniprot.org/uniprot/Q7YTG1" target="_blank"> Q7YTG1</a>')
    expect_equal(.makeLinkFromId("Q7YTG1", "AlphaFold"),
                 '<a href="https://alphafold.ebi.ac.uk/entry/Q7YTG1" target="_blank"> Q7YTG1</a>')
    expect_equal(.makeLinkFromId("SPBC460.01c", "PomBase"),
                 '<a href="https://www.pombase.org/gene/SPBC460.01c" target="_blank"> SPBC460.01c</a>')
    expect_equal(.makeLinkFromId("WBGene00001330", "WormBase"),
                 '<a href="https://wormbase.org/species/c_elegans/gene/WBGene00001330" target="_blank"> WBGene00001330</a>')
    expect_equal(.makeLinkFromId("Q7YTG1-1", "AlphaFold", removeSuffix = TRUE),
                 '<a href="https://alphafold.ebi.ac.uk/entry/Q7YTG1" target="_blank"> Q7YTG1</a>')
    expect_equal(.makeLinkFromId("Q7YTG1-1", "AlphaFold", removeSuffix = FALSE),
                 '<a href="https://alphafold.ebi.ac.uk/entry/Q7YTG1-1" target="_blank"> Q7YTG1-1</a>')

    ## getConvTable
    ## --------------------------------------------------------------------- ##
    expect_error(getConvTable(type = "WrongType"),
                 "All values in 'type' must be one of")
    expect_error(getConvTable(type = c("PomBase", "WormBase")),
                 "'type' must have length 1")
    expect_error(getConvTable(type = 1),
                 "'type' must be of class 'character'")

    ## makeDbLinkTable
    ## --------------------------------------------------------------------- ##
    expect_error(makeDbLinkTable(df = 1, idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'df' must be of class 'data.frame'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = 1,
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'idCol' must be of class 'character'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = c("id", "id"),
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'idCol' must have length 1")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "missing",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "All values in 'idCol' must be one of")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = 1,
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'speciesCommon' must be of class 'character'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = c("mouse", "human"),
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'speciesCommon' must have length 1")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "missing",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "All values in 'speciesCommon' must be one of")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = 1,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'addSpeciesSpecificColumns' must be of class 'logical'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = c(TRUE, FALSE),
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'addSpeciesSpecificColumns' must have length 1")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = 1,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'convTablePomBase' must be of class 'data.frame'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = 1,
                                 removeSuffix = TRUE,
                                 signifDigits = 3),
                 "'convTableWormBase' must be of class 'data.frame'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = 1,
                                 signifDigits = 3),
                 "'removeSuffix' must be of class 'logical'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = c(TRUE, FALSE),
                                 signifDigits = 3),
                 "'removeSuffix' must have length 1")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = "3"),
                 "'signifDigits' must be of class 'numeric'")
    expect_error(makeDbLinkTable(df = data.frame(id = c("B5BP45", "O13282")),
                                 idCol = "id",
                                 speciesCommon = "fission yeast",
                                 addSpeciesSpecificColumns = TRUE,
                                 convTablePomBase = NULL,
                                 convTableWormBase = NULL,
                                 removeSuffix = TRUE,
                                 signifDigits = c(3, 4)),
                 "'signifDigits' must have length 1")


    ## Without species-specific columns
    dblt <- makeDbLinkTable(data.frame(id = c("B5BP45", "O13282", "B5BP45"),
                                       numcol = c(1.23456, 0.00034561, 7625.23)),
                            idCol = "id", speciesCommon = "fission yeast",
                            signifDigits = 3)
    expect_s3_class(dblt, "data.frame")
    expect_equal(ncol(dblt), 4)
    expect_equal(nrow(dblt), 3)
    expect_named(dblt, c("id", "numcol", "UniProt", "AlphaFold"))
    expect_equal(dblt$id, c("B5BP45", "O13282", "B5BP45"))
    expect_equal(dblt$numcol, c(1.23, 0.000346, 7630))
    expect_equal(dblt$UniProt, c(
        '<a href="https://www.uniprot.org/uniprot/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://www.uniprot.org/uniprot/O13282" target="_blank"> O13282</a>',
        '<a href="https://www.uniprot.org/uniprot/B5BP45" target="_blank"> B5BP45</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt$AlphaFold, c(
        '<a href="https://alphafold.ebi.ac.uk/entry/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/O13282" target="_blank"> O13282</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/B5BP45" target="_blank"> B5BP45</a>'),
        ignore_attr = TRUE)

    ## As above, but different number of significant digits
    dblt <- makeDbLinkTable(data.frame(id = c("B5BP45", "O13282", "O13282"),
                                       numcol = c(1.23456, 0.00034561, 7625.23)),
                            idCol = "id", speciesCommon = "fission yeast",
                            signifDigits = 1)
    expect_s3_class(dblt, "data.frame")
    expect_equal(ncol(dblt), 4)
    expect_equal(nrow(dblt), 3)
    expect_named(dblt, c("id", "numcol", "UniProt", "AlphaFold"))
    expect_equal(dblt$id, c("B5BP45", "O13282", "O13282"))
    expect_equal(dblt$numcol, c(1, 0.0003, 8000))
    expect_equal(dblt$UniProt, c(
        '<a href="https://www.uniprot.org/uniprot/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://www.uniprot.org/uniprot/O13282" target="_blank"> O13282</a>',
        '<a href="https://www.uniprot.org/uniprot/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt$AlphaFold, c(
        '<a href="https://alphafold.ebi.ac.uk/entry/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/O13282" target="_blank"> O13282</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)

    ## With Pombase column
    dblt2 <- makeDbLinkTable(
        data.frame(id = c("B5BP45", "O13282")), idCol = "id",
        speciesCommon = "fission yeast",
        convTablePomBase = data.frame(
            PomBaseID = c("SPCC5E4.03c", "SPBC460.01c"),
            UniProtID = c("O13282", "B5BP45")))
    expect_s3_class(dblt2, "data.frame")
    expect_equal(ncol(dblt2), 4)
    expect_equal(nrow(dblt2), 2)
    expect_named(dblt2, c("id", "UniProt", "AlphaFold", "PomBase"))
    expect_equal(dblt2$id, c("B5BP45", "O13282"))
    expect_equal(dblt2$UniProt, c(
        '<a href="https://www.uniprot.org/uniprot/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://www.uniprot.org/uniprot/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt2$AlphaFold, c(
        '<a href="https://alphafold.ebi.ac.uk/entry/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt2$PomBase, c(
        '<a href=\"https://www.pombase.org/gene/SPBC460.01c\" target=\"_blank\"> SPBC460.01c</a>',
        '<a href=\"https://www.pombase.org/gene/SPCC5E4.03c\" target=\"_blank\"> SPCC5E4.03c</a>'),
        ignore_attr = TRUE)

    ## With Pombase conversion table, but ignore
    dblt3 <- makeDbLinkTable(
        data.frame(id = c("B5BP45", "O13282")), idCol = "id",
        speciesCommon = "fission yeast", addSpeciesSpecificColumns = FALSE,
        convTablePomBase = data.frame(
            PomBaseID = c("SPBC460.01c", "SPCC5E4.03c"),
            UniProtID = c("B5BP45", "O13282")))
    expect_s3_class(dblt3, "data.frame")
    expect_equal(ncol(dblt3), 3)
    expect_equal(nrow(dblt3), 2)
    expect_named(dblt3, c("id", "UniProt", "AlphaFold"))
    expect_equal(dblt3$id, c("B5BP45", "O13282"))
    expect_equal(dblt3$UniProt, c(
        '<a href="https://www.uniprot.org/uniprot/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://www.uniprot.org/uniprot/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt3$AlphaFold, c(
        '<a href="https://alphafold.ebi.ac.uk/entry/B5BP45" target="_blank"> B5BP45</a>',
        '<a href="https://alphafold.ebi.ac.uk/entry/O13282" target="_blank"> O13282</a>'),
        ignore_attr = TRUE)

    ## With Wormbase column
    dblt4 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm",
        convTableWormBase = data.frame(
            UniProtID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
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
    expect_equal(dblt4$gid, c("eps-8", "epi-1"))
    expect_equal(dblt4$pid, c("Q7YTG1;O18250", "C1P641;C1P640"))
    expect_equal(dblt4$UniProt, c(
        '<a href=\"https://www.uniprot.org/uniprot/Q7YTG1\" target=\"_blank\"> Q7YTG1</a>;<a href=\"https://www.uniprot.org/uniprot/O18250\" target=\"_blank\"> O18250</a>',
        '<a href=\"https://www.uniprot.org/uniprot/C1P641\" target=\"_blank\"> C1P641</a>;<a href=\"https://www.uniprot.org/uniprot/C1P640\" target=\"_blank\"> C1P640</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt4$AlphaFold, c(
        '<a href=\"https://alphafold.ebi.ac.uk/entry/Q7YTG1\" target=\"_blank\"> Q7YTG1</a>;<a href=\"https://alphafold.ebi.ac.uk/entry/O18250\" target=\"_blank\"> O18250</a>',
        '<a href=\"https://alphafold.ebi.ac.uk/entry/C1P641\" target=\"_blank\"> C1P641</a>;<a href=\"https://alphafold.ebi.ac.uk/entry/C1P640\" target=\"_blank\"> C1P640</a>'),
        ignore_attr = TRUE)
    expect_equal(dblt4$WormBase, c(
        '<a href="https://wormbase.org/species/c_elegans/gene/WBGene00001330" target="_blank"> WBGene00001330</a>',
        '<a href="https://wormbase.org/species/c_elegans/gene/WBGene00001328" target="_blank"> WBGene00001328</a>'),
        ignore_attr = TRUE)

    ## With Wormbase conversion table, but ignore
    dblt5 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm", addSpeciesSpecificColumns = FALSE,
        convTableWormBase = data.frame(
            UniProtID = c("O18250", "C1P641", "Q7YTG1", "C1P640"),
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
    expect_equal(dblt5$AlphaFold, dblt4$AlphaFold)
    expect_equal(dblt5$UniProt, dblt4$UniProt)

    ## With Wormbase column, but missing conversion
    dblt6 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm",
        convTableWormBase = data.frame(
            UniProtID = c("O18250", "C1P641"),
            UniProtKB.ID = c("O18250", "C1P641"),
            WormBaseID = c("WBGene00001330", "WBGene00001328"),
            check.names = FALSE))
    expect_s3_class(dblt6, "data.frame")
    expect_equal(ncol(dblt6), 5)
    expect_equal(nrow(dblt6), 2)
    expect_named(dblt6, c("gid", "pid", "UniProt", "AlphaFold", "WormBase"))
    expect_equal(grep(";", dblt6$UniProt), c(1, 2))
    expect_equal(grep(";", dblt6$AlphaFold), c(1, 2))
    expect_equal(grep(";", dblt6$WormBase), integer(0))
    expect_equal(dblt6$gid, c("eps-8", "epi-1"))
    expect_equal(dblt6$pid, c("Q7YTG1;O18250", "C1P641;C1P640"))
    expect_equal(dblt6$UniProt, dblt4$UniProt)
    expect_equal(dblt6$AlphaFold, dblt4$AlphaFold)
    expect_equal(dblt6$WormBase, dblt4$WormBase)

    ## With Wormbase column, but missing conversion (2)
    dblt7 <- makeDbLinkTable(
        data.frame(gid = c("eps-8", "epi-1"),
                   pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
        idCol = "pid", speciesCommon = "roundworm",
        convTableWormBase = data.frame(UniProtID = c("O18250"),
                                       UniProtKB.ID = c("O18250"),
                                       WormBaseID = c("WBGene00001330"),
                                       check.names = FALSE))
    expect_s3_class(dblt7, "data.frame")
    expect_equal(ncol(dblt7), 5)
    expect_equal(nrow(dblt7), 2)
    expect_named(dblt7, c("gid", "pid", "UniProt", "AlphaFold", "WormBase"))
    expect_equal(grep(";", dblt7$UniProt), c(1, 2))
    expect_equal(grep(";", dblt7$AlphaFold), c(1, 2))
    expect_equal(grep(";", dblt7$WormBase), integer(0))
    expect_equal(dblt7$gid, c("eps-8", "epi-1"))
    expect_equal(dblt7$pid, c("Q7YTG1;O18250", "C1P641;C1P640"))
    expect_equal(dblt7$UniProt, dblt4$UniProt)
    expect_equal(dblt7$AlphaFold, dblt4$AlphaFold)
    expect_equal(dblt7$WormBase, c(dblt4$WormBase[1], ""),
                 ignore_attr = TRUE)

    ## Skip the rest of the tests if no internet connection is available
    # skip_if_offline()
    dfp <- getConvTable(type = "PomBase")
    expect_equal(ncol(dfp), 2)
    expect_named(dfp, c("PomBaseID", "UniProtID"))

    dfw <- getConvTable(type = "WormBase")
    expect_equal(ncol(dfw), 3)
    expect_named(dfw, c("UniProtID", "UniProtKB.ID", "WormBaseID"))

    ## Don't add more tests below here unless they should be skipped if offline
})
