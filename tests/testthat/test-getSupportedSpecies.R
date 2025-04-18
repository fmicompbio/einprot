test_that("geting species information works", {
    df <- getSupportedSpecies()
    expect_s3_class(df, "data.frame")
    expect_equal(ncol(df), 3)
    expect_named(df, c("taxId", "species", "speciesCommon"))
    expect_equal(nrow(df), 20)

    expect_warning(getSpeciesInfo("missing"),
                   "Unknown species missing")
    expect_error(getSpeciesInfo(list(x = 1)),
                 "is.character(species) || is.numeric(species) is not TRUE",
                 fixed = TRUE)
    expect_error(getSpeciesInfo(c("mouse", "human")),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    expect_type(getSpeciesInfo("mouse"), "list")
    expect_equal(length(getSpeciesInfo("mouse")), 3)
    expect_named(getSpeciesInfo("mouse"),
                 c("species", "speciesCommon", "taxId"))

    expect_equal(getSpeciesInfo("Mouse")$species, "Mus musculus")
    expect_equal(getSpeciesInfo("MoUsE")$speciesCommon, "mouse")
    expect_equal(getSpeciesInfo("mouse")$taxId, 10090)
    expect_equal(getSpeciesInfo("Mus musculus")$species, "Mus musculus")
    expect_equal(getSpeciesInfo("Mus musculus")$speciesCommon, "mouse")
    expect_equal(getSpeciesInfo("Mus musculus")$taxId, 10090)
    expect_equal(getSpeciesInfo(10090)$species, "Mus musculus")
    expect_equal(getSpeciesInfo(10090)$speciesCommon, "mouse")
    expect_equal(getSpeciesInfo(10090)$taxId, 10090)

    expect_equal(getSpeciesInfo("Human")$species, "Homo sapiens")
    expect_equal(getSpeciesInfo("HuMaN")$speciesCommon, "human")
    expect_equal(getSpeciesInfo("human")$taxId, 9606)
    expect_equal(getSpeciesInfo("Homo sapiens")$species, "Homo sapiens")
    expect_equal(getSpeciesInfo("Homo sapiens")$speciesCommon, "human")
    expect_equal(getSpeciesInfo("Homo sapiens")$taxId, 9606)
    expect_equal(getSpeciesInfo(9606)$species, "Homo sapiens")
    expect_equal(getSpeciesInfo(9606)$speciesCommon, "human")
    expect_equal(getSpeciesInfo(9606)$taxId, 9606)

    expect_equal(getSpeciesInfo("Roundworm")$species, "Caenorhabditis elegans")
    expect_equal(getSpeciesInfo("RounDWorm")$speciesCommon, "roundworm")
    expect_equal(getSpeciesInfo("roundworm")$taxId, 6239)
    expect_equal(getSpeciesInfo("Caenorhabditis elegans")$species,
                 "Caenorhabditis elegans")
    expect_equal(getSpeciesInfo("Caenorhabditis elegans")$speciesCommon,
                 "roundworm")
    expect_equal(getSpeciesInfo("Caenorhabditis elegans")$taxId, 6239)
    expect_equal(getSpeciesInfo(6239)$species, "Caenorhabditis elegans")
    expect_equal(getSpeciesInfo(6239)$speciesCommon, "roundworm")
    expect_equal(getSpeciesInfo(6239)$taxId, 6239)

    expect_equal(getSpeciesInfo("Zebrafish")$species, "Danio rerio")
    expect_equal(getSpeciesInfo("ZebraFish")$speciesCommon, "zebrafish")
    expect_equal(getSpeciesInfo("zebrafish")$taxId, 7955)
    expect_equal(getSpeciesInfo("Danio rerio")$species, "Danio rerio")
    expect_equal(getSpeciesInfo("Danio rerio")$speciesCommon, "zebrafish")
    expect_equal(getSpeciesInfo("Danio rerio")$taxId, 7955)
    expect_equal(getSpeciesInfo(7955)$species, "Danio rerio")
    expect_equal(getSpeciesInfo(7955)$speciesCommon, "zebrafish")
    expect_equal(getSpeciesInfo(7955)$taxId, 7955)

    expect_equal(getSpeciesInfo("Fruitfly")$species, "Drosophila melanogaster")
    expect_equal(getSpeciesInfo("FruitFly")$speciesCommon, "fruitfly")
    expect_equal(getSpeciesInfo("fruitfly")$taxId, 7227)
    expect_equal(getSpeciesInfo("Drosophila melanogaster")$species,
                 "Drosophila melanogaster")
    expect_equal(getSpeciesInfo("Drosophila melanogaster")$speciesCommon,
                 "fruitfly")
    expect_equal(getSpeciesInfo("Drosophila melanogaster")$taxId, 7227)
    expect_equal(getSpeciesInfo(7227)$species, "Drosophila melanogaster")
    expect_equal(getSpeciesInfo(7227)$speciesCommon, "fruitfly")
    expect_equal(getSpeciesInfo(7227)$taxId, 7227)

    expect_equal(getSpeciesInfo("Baker's yeast")$species,
                 "Saccharomyces cerevisiae")
    expect_equal(getSpeciesInfo("Baker's Yeast")$speciesCommon,
                 "baker's yeast")
    expect_equal(getSpeciesInfo("baker's yeast")$taxId, 4932)
    expect_equal(getSpeciesInfo("Saccharomyces cerevisiae")$species,
                 "Saccharomyces cerevisiae")
    expect_equal(getSpeciesInfo("Saccharomyces cerevisiae")$speciesCommon,
                 "baker's yeast")
    expect_equal(getSpeciesInfo("Saccharomyces cerevisiae")$taxId, 4932)
    expect_equal(getSpeciesInfo(4932)$species, "Saccharomyces cerevisiae")
    expect_equal(getSpeciesInfo(4932)$speciesCommon, "baker's yeast")
    expect_equal(getSpeciesInfo(4932)$taxId, 4932)

    expect_equal(getSpeciesInfo("Fission Yeast")$species,
                 "Schizosaccharomyces pombe 972h-")
    expect_equal(getSpeciesInfo("Fission yeast")$speciesCommon,
                 "fission yeast")
    expect_equal(getSpeciesInfo("fission yeast")$taxId, 284812)
    expect_equal(getSpeciesInfo("Schizosaccharomyces pombe 972h-")$species,
                 "Schizosaccharomyces pombe 972h-")
    expect_equal(
        getSpeciesInfo("Schizosaccharomyces pombe 972h-")$speciesCommon,
        "fission yeast"
    )
    expect_equal(getSpeciesInfo("Schizosaccharomyces pombe 972h-")$taxId,
                 284812)
    expect_equal(getSpeciesInfo(284812)$species,
                 "Schizosaccharomyces pombe 972h-")
    expect_equal(getSpeciesInfo(284812)$speciesCommon, "fission yeast")
    expect_equal(getSpeciesInfo(284812)$taxId, 284812)

})
