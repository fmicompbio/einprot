test_that("preparing feature collections works", {
    ## .createPattern
    ## -----------------------------------------------------------------------
    expect_equal(.createPattern(2), "([^;]+;[^;]+);")
    expect_equal(.createPattern(3), "([^;]+;[^;]+;[^;]+);")

    ## .replaceIdsInList
    ## -----------------------------------------------------------------------
    chl <- IRanges::CharacterList(a = c("g3", "g2", "g1", "g7", "g4"),
                                  b = c("g5", "g6", "g1"))
    dfConv <- data.frame(orig = paste0("g", 1:6),
                         new = paste0("h", 1:6))
    pat <- .createPattern(2)

    expect_error(.replaceIdsInList(chl = 1, dfConv = dfConv,
                                   currentIdCol = "orig", newIdCol = "new",
                                   pat = pat),
                 "'chl' must be of class 'CharacterList'")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = 1, newIdCol = "new",
                                   pat = pat),
                 "'currentIdCol' must be of class 'character'")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = c("orig", "new"),
                                   newIdCol = "new", pat = pat),
                 "'currentIdCol' must have length 1")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = "missing",
                                   newIdCol = "new", pat = pat),
                 "All values in 'currentIdCol' must be one of")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = "orig", newIdCol = 1,
                                   pat = pat),
                 "'newIdCol' must be of class 'character'")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = "orig",
                                   newIdCol = c("orig", "new"), pat = pat),
                 "'newIdCol' must have length 1")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = "orig",
                                   newIdCol = "missing", pat = pat),
                 "All values in 'newIdCol' must be one of")
    expect_error(.replaceIdsInList(chl = chl, dfConv = dfConv,
                                   currentIdCol = "orig", newIdCol = "new",
                                   pat = 1),
                 "'pat' must be of class 'character'")

    chl1 <- .replaceIdsInList(chl = chl, dfConv = dfConv,
                              currentIdCol = "orig", newIdCol = "new",
                              pat = pat)
    expect_s4_class(chl1, "CharacterList")
    expect_length(chl1, 2)
    expect_named(chl1, c("a", "b"))
    expect_equal(chl1[["a"]], c("h1", "h2", "h3", "h4"))
    expect_equal(chl1[["b"]], c("h1", "h5", "h6"))
    expect_s4_class(S4Vectors::mcols(chl1), "DFrame")
    expect_equal(S4Vectors::mcols(chl1)$sharedGenes,
                 c("h1;h2; h3;h4", "h1;h5; h6"), ignore_attr = TRUE)
    expect_equal(S4Vectors::mcols(chl1)$nSharedGenes, c(4, 3),
                 ignore_attr = TRUE)

    ## many-to-one mappings
    chl <- IRanges::CharacterList(a = c("g3", "g2", "g1", "g7", "g4"),
                                  b = c("g5", "g6", "g1"))
    dfConv <- data.frame(orig = paste0("g", 1:6),
                         new = paste0("h", c(1, 1, 3, 3, 5, 5)))
    pat <- .createPattern(2)
    chl1 <- .replaceIdsInList(chl = chl, dfConv = dfConv,
                              currentIdCol = "orig", newIdCol = "new",
                              pat = pat)
    expect_s4_class(chl1, "CharacterList")
    expect_length(chl1, 2)
    expect_named(chl1, c("a", "b"))
    expect_equal(chl1[["a"]], c("h1", "h3"))
    expect_equal(chl1[["b"]], c("h1", "h5"))
    expect_s4_class(S4Vectors::mcols(chl1), "DFrame")
    expect_equal(S4Vectors::mcols(chl1)$sharedGenes,
                 c("h1;h3", "h1;h5"), ignore_attr = TRUE)
    expect_equal(S4Vectors::mcols(chl1)$nSharedGenes, c(2, 2),
                 ignore_attr = TRUE)

    ## one-to-many mappings
    chl <- IRanges::CharacterList(a = c("g3", "g2", "g1", "g7", "g4"),
                                  b = c("g5", "g6", "g1"))
    dfConv <- data.frame(orig = paste0("g", c(1:6, 1:3)),
                         new = c(paste0("h", 1:6), paste0("i", 1:3)))
    pat <- .createPattern(2)
    chl1 <- .replaceIdsInList(chl = chl, dfConv = dfConv,
                              currentIdCol = "orig", newIdCol = "new",
                              pat = pat)
    expect_s4_class(chl1, "CharacterList")
    expect_length(chl1, 2)
    expect_named(chl1, c("a", "b"))
    expect_equal(chl1[["a"]], c("h1", "h2", "h3", "h4", "i1", "i2", "i3"))
    expect_equal(chl1[["b"]], c("h1", "h5", "h6", "i1"))
    expect_s4_class(S4Vectors::mcols(chl1), "DFrame")
    expect_equal(S4Vectors::mcols(chl1)$sharedGenes,
                 c("h1;h2; h3;h4; i1;i2; i3", "h1;h5; h6;i1"), ignore_attr = TRUE)
    expect_equal(S4Vectors::mcols(chl1)$nSharedGenes, c(7, 4),
                 ignore_attr = TRUE)

    ## no matches
    chl <- IRanges::CharacterList(a = c("g3", "g2", "g1", "g7", "g4"),
                                  b = c("g8", "g9"))
    dfConv <- data.frame(orig = paste0("g", c(1:6, 1:3)),
                         new = c(paste0("h", 1:6), paste0("i", 1:3)))
    pat <- .createPattern(2)
    chl1 <- .replaceIdsInList(chl = chl, dfConv = dfConv,
                              currentIdCol = "orig", newIdCol = "new",
                              pat = pat)
    expect_s4_class(chl1, "CharacterList")
    expect_length(chl1, 2)
    expect_named(chl1, c("a", "b"))
    expect_equal(chl1[["a"]], c("h1", "h2", "h3", "h4", "i1", "i2", "i3"))
    expect_equal(chl1[["b"]], character(0))
    expect_s4_class(S4Vectors::mcols(chl1), "DFrame")
    expect_equal(S4Vectors::mcols(chl1)$sharedGenes,
                 c("h1;h2; h3;h4; i1;i2; i3", ""), ignore_attr = TRUE)
    expect_equal(S4Vectors::mcols(chl1)$nSharedGenes, c(7, 0),
                 ignore_attr = TRUE)

    ## prepareFeatureCollections
    ## -----------------------------------------------------------------------
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 100)$sce
    sce <- fixFeatureIds(sce)

    args0 <- list(
        sce = sce,
        idCol = "Gene.names",
        includeFeatureCollections = "complexes",
        complexDbPath = system.file("extdata", "complexes",
                                    "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                    package = "einprot"),
        speciesInfo = getSpeciesInfo("mouse"),
        complexSpecies = "current",
        customComplexes = list(),
        minSizeToKeep = 2
    )

    ## Fails with wrong arguments
    args <- args0
    args$sce <- 1
    expect_error(do.call(prepareFeatureCollections, args),
                 "'sce' must be of class 'SummarizedExperiment'")

    args <- args0
    args$idCol <- 1
    expect_error(do.call(prepareFeatureCollections, args),
                 "'idCol' must be of class 'character'")
    args$idCol <- c("Gene.names", "Majority.protein.IDs")
    expect_error(do.call(prepareFeatureCollections, args),
                 "'idCol' must have length 1")
    args$idCol <- "missing"
    expect_error(do.call(prepareFeatureCollections, args),
                 "All values in 'idCol' must be one of")

    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(prepareFeatureCollections, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "missing"
    expect_error(do.call(prepareFeatureCollections, args),
                 "All values in 'includeFeatureCollections' must be one of")

    args <- args0
    args$speciesInfo <- "mouse"
    expect_error(do.call(prepareFeatureCollections, args),
                 "'speciesInfo' must be of class 'list'")
    args$speciesInfo <- list("mouse", "Mus musculus")
    expect_error(do.call(prepareFeatureCollections, args),
                 "'namesspeciesInfo' must not be NULL")
    args$speciesInfo <- list(one = "mouse", two = "Mus musculus")
    expect_error(do.call(prepareFeatureCollections, args),
                 "All values in 'namesspeciesInfo' must be one of")
    args$speciesInfo <- list(speciesCommon = "missing", species = "missing")
    expect_error(do.call(prepareFeatureCollections, args),
                 "No complex database available for the current species")

    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(prepareFeatureCollections, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(prepareFeatureCollections, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "missing"
    expect_error(do.call(prepareFeatureCollections, args),
                 "All values in 'complexSpecies' must be one of")

    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(prepareFeatureCollections, args),
                 "'customComplexes' must be of class 'list'")
    args <- args0
    args$customComplexes <- list(c("gene1", "gene2"))
    expect_error(do.call(prepareFeatureCollections, args),
                 "'namescustomComplexes' must not be NULL")

    args <- args0
    args$minSizeToKeep <- "2"
    expect_error(do.call(prepareFeatureCollections, args),
                 "'minSizeToKeep' must be of class 'numeric'")
    args$minSizeToKeep <- c(1, 2)
    expect_error(do.call(prepareFeatureCollections, args),
                 "'minSizeToKeep' must have length 1")

    ## Works with correct arguments
    fcoll <- do.call(prepareFeatureCollections, args0)
    expect_type(fcoll, "list")
    expect_named(fcoll, "complexes")
    expect_s4_class(fcoll$complexes, "CharacterList")
    expect_length(fcoll$complexes, 1)
    expect_true(all(lengths(fcoll$complexes) >= 2))
    expect_equal(lengths(fcoll$complexes), 2, ignore_attr = TRUE)
    expect_true(all(unlist(fcoll$complexes) %in% rownames(sce)))
    mcc <- S4Vectors::mcols(fcoll$complexes)
    expect_s4_class(mcc, "DFrame")
    expect_named(mcc, c("Species.common", "Source", "PMID", "All.names",
                        "genes", "nGenes", "sharedGenes", "nSharedGenes"))
    expect_equal(lengths(fcoll$complexes), mcc$nSharedGenes)
    expect_equal(mcc$Species.common, rep("mouse", 1))

    ## Wrong column for IDs - empty output
    args <- args0
    args$idCol <- "Majority.protein.IDs"
    expect_length(do.call(prepareFeatureCollections, args)$complexes, 0)

    ## complexSpecies
    args <- args0
    args$complexSpecies <- "all"
    fcoll <- do.call(prepareFeatureCollections, args)
    expect_length(fcoll$complexes, 45)
    expect_true(all(unlist(fcoll$complexes) %in% rownames(args$sce)))
    expect_equal(substr(names(fcoll$complexes)[1], 1, 5), "mouse")

    ## customComplexes
    args <- args0
    args$customComplexes <- list(my_complex = c("Dhx9", "Krt10", "missing"))
    fcoll <- do.call(prepareFeatureCollections, args)
    expect_type(fcoll, "list")
    expect_s4_class(fcoll$complexes, "CharacterList")
    expect_length(fcoll$complexes, 2)
    expect_true(all(unlist(fcoll$complexes) %in% rownames(args$sce)))
    expect_equal(names(fcoll$complexes)[1], "my_complex")
    expect_equal(fcoll$complexes[[1]], c("Dhx9", "Krt10"))
    mcc <- S4Vectors::mcols(fcoll$complexes)
    expect_s4_class(mcc, "DFrame")
    expect_named(mcc, c("Species.common", "Source", "PMID", "All.names",
                        "genes", "nGenes", "sharedGenes", "nSharedGenes"))
    expect_equal(lengths(fcoll$complexes), mcc$nSharedGenes)
    expect_equal(mcc$Species.common, rep("mouse", 2))
    expect_equal(mcc$Source[1], "custom")
    expect_equal(mcc$All.names[1], "my_complex")
    expect_equal(mcc$genes[1], "Dhx9;Krt10;missing", ignore_attr = TRUE)
    expect_equal(mcc$nGenes[1], 3, ignore_attr = TRUE)
    expect_equal(mcc$sharedGenes[1], "Dhx9;Krt10", ignore_attr = TRUE)
    expect_equal(mcc$nSharedGenes[1], 2, ignore_attr = TRUE)

    ## minSizeToKeep
    args <- args0
    args$minSizeToKeep <- 3
    args$speciesInfo$species <- "mouse"
    args$speciesInfo$speciesCommon <- "missing"  ## should use species instead
    expect_length(do.call(prepareFeatureCollections, args)$complexes, 0)

    ## Don't extract anything
    args <- args0
    fcoll <- prepareFeatureCollections(sce = args$sce, idCol = args$idCol,
                                       includeFeatureCollections = NULL,
                                       complexDbPath = args$complexDbPath,
                                       speciesInfo = args$speciesInfo,
                                       complexSpecies = args$complexSpecies,
                                       customComplexes = args$customComplexes,
                                       minSizeToKeep = args$minSizeToKeep)
    expect_length(fcoll, 0)
    expect_type(fcoll, "list")

    ## Only a custom complex
    args <- args0
    fcoll <- prepareFeatureCollections(
        sce = args$sce, idCol = args$idCol,
        includeFeatureCollections = NULL,
        complexDbPath = args$complexDbPath,
        speciesInfo = args$speciesInfo,
        complexSpecies = args$complexSpecies,
        customComplexes = list(my_complex = c("Dhx9", "Krt10", "missing")),
        minSizeToKeep = args$minSizeToKeep)
    expect_length(fcoll, 1)
    expect_named(fcoll, "complexes")
    expect_type(fcoll, "list")
    expect_s4_class(fcoll$complexes, "CharacterList")
    expect_named(fcoll$complexes, "my_complex")
    expect_equal(fcoll$complexes[[1]], c("Dhx9", "Krt10"))
    expect_length(fcoll$complexes, 1)
    expect_true(all(unlist(fcoll$complexes) %in% rownames(args$sce)))
    mcc <- S4Vectors::mcols(fcoll$complexes)
    expect_s4_class(mcc, "DFrame")
    expect_named(mcc, c("Species.common", "Source", "PMID", "All.names",
                        "genes", "nGenes", "sharedGenes", "nSharedGenes"))
    expect_equal(lengths(fcoll$complexes), mcc$nSharedGenes)
    expect_equal(mcc$Species.common, rep("mouse", 1))
    expect_equal(mcc$Source[1], "custom")
    expect_equal(mcc$All.names[1], "my_complex")
    expect_equal(mcc$genes[1], "Dhx9;Krt10;missing", ignore_attr = TRUE)
    expect_equal(mcc$nGenes[1], 3, ignore_attr = TRUE)
    expect_equal(mcc$sharedGenes[1], "Dhx9;Krt10", ignore_attr = TRUE)
    expect_equal(mcc$nSharedGenes[1], 2, ignore_attr = TRUE)

    ## Only a custom complex, rows with multiple genes
    args <- args0
    fcoll <- prepareFeatureCollections(
        sce = args$sce, idCol = args$idCol,
        includeFeatureCollections = NULL,
        complexDbPath = args$complexDbPath,
        speciesInfo = args$speciesInfo,
        complexSpecies = args$complexSpecies,
        customComplexes = list(my_complex = c("Gm13145", "Zfp518b", "Gnai1",
                                              "Gnai2", "Gm3424", "Zfp462"),
                               my_complex2 = c("Zfp462")),
        minSizeToKeep = args$minSizeToKeep)
    expect_length(fcoll, 1)
    expect_named(fcoll, "complexes")
    expect_type(fcoll, "list")
    expect_s4_class(fcoll$complexes, "CharacterList")
    expect_named(fcoll$complexes, c("my_complex", "my_complex2"))
    expect_equal(sort(fcoll$complexes[[1]]),
                 sort(c("Gm13145", "Znf518b", "Gnai1",
                        "1700001E04Rik", "Zfp462.B1AWL2",
                        "Zfp462.B1AWL4", "Zfp462.B1AWL5")))
    expect_equal(sort(fcoll$complexes[[2]]),
                 sort(c("Zfp462.B1AWL2", "Zfp462.B1AWL4", "Zfp462.B1AWL5")))
    expect_length(fcoll$complexes, 2)
    expect_true(all(unlist(fcoll$complexes) %in% rownames(args$sce)))
    mcc <- S4Vectors::mcols(fcoll$complexes)
    expect_s4_class(mcc, "DFrame")
    expect_named(mcc, c("Species.common", "Source", "PMID", "All.names",
                        "genes", "nGenes", "sharedGenes", "nSharedGenes"))
    expect_equal(lengths(fcoll$complexes), mcc$nSharedGenes)
    expect_equal(mcc$Species.common, rep("mouse", 2))
    expect_equal(mcc$Source, c("custom", "custom"))
    expect_equal(mcc$All.names, c("my_complex", "my_complex2"))
    expect_equal(mcc$genes, c("Gm13145;Zfp518b;Gnai1;Gnai2;Gm3424;Zfp462",
                              "Zfp462"), ignore_attr = TRUE)
    expect_equal(mcc$nGenes, c(6, 1), ignore_attr = TRUE)
    expect_equal(mcc$sharedGenes,
                 c("Gm13145;Zfp462.B1AWL2;Zfp462.B1AWL4;Zfp462.B1AWL5;Znf518b;Gnai1;1700001E04Rik",
                   "Zfp462.B1AWL2;Zfp462.B1AWL4;Zfp462.B1AWL5"),
                 ignore_attr = TRUE)
    expect_equal(mcc$nSharedGenes, c(7, 3), ignore_attr = TRUE)
})
