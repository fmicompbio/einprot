test_that("fixing feature IDs works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce

    ## Fail with wrong arguments
    expect_error(fixFeatureIds(sce = 1, primaryIdCol = "Gene.names",
                               secondaryIdCol = "Majority.protein.IDs"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = 1,
                               secondaryIdCol = "Majority.protein.IDs"),
                 "'primaryIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = c("Gene.names", "Gene.names"),
                               secondaryIdCol = "Majority.protein.IDs"),
                 "'primaryIdCol' must have length 1")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = "missing",
                               secondaryIdCol = "Majority.protein.IDs"),
                 "All values in 'primaryIdCol' must be one of")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                               secondaryIdCol = 1),
                 "'secondaryIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                               secondaryIdCol = c("Majority.protein.IDs",
                                                "Majority.protein.IDs")),
                 "'secondaryIdCol' must have length 1")
    expect_error(fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                               secondaryIdCol = "missing"),
                 "All values in 'secondaryIdCol' must be one of")

    ## Test that it does the right thing
    ## All gene and protein names are unique
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")

    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                          secondaryIdCol = "Majority.protein.IDs")
    expect_equal(rownames(sce1), gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Majority.protein.IDs",
                          secondaryIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)

    ## Missing gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce[36:45, ]
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs

    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                          secondaryIdCol = "Majority.protein.IDs")
    midx <- which(gns == "")  ## indices for missing gene IDs
    eidx <- which(gns != "")  ## indices for present gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    expect_equal(rownames(sce1)[eidx], gns[eidx])
    expect_equal(rownames(sce1)[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[eidx], gns[eidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Majority.protein.IDs",
                          secondaryIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)

    ## Duplicated gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Gene.names",
                          secondaryIdCol = "Majority.protein.IDs")
    midx <- which(gns == "")  ## indices for missing gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    didx <- setdiff(which(gns %in% gns[duplicated(gns)]), midx)  ## indices for duplicated gene IDs
    uidx <- setdiff(seq_along(gns), c(midx, didx))  ## indices for unique gene IDs
    expect_length(uidx, 38)
    expect_length(midx, 4)
    expect_length(didx, 3)
    expect_equal(rownames(sce1)[uidx], gns[uidx])
    expect_equal(rownames(sce1)[midx], pns[midx])
    expect_equal(rownames(sce1)[didx], paste0(gns[didx], ".", pns[didx]))
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Majority.protein.IDs",
                          secondaryIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)

    ## --------------------------------------------------------------------- ##
    ## PD data
    sce <- importExperiment(
        inFile = system.file("extdata", "pdtmt_example",
                             "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                             package = "einprot"),
        iColPattern = "^Abundance\\.F.+\\.Sample\\.")$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.Symbol
    pns <- SummarizedExperiment::rowData(sce)$Accession
    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Gene.Symbol",
                          secondaryIdCol = "Accession")
    midx <- which(gns == "")  ## indices for missing gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    didx <- setdiff(which(gns %in% gns[duplicated(gns)]), midx)  ## indices for duplicated gene IDs
    uidx <- setdiff(seq_along(gns), c(midx, didx))  ## indices for unique gene IDs
    expect_length(uidx, 1620)
    expect_length(midx, 75)
    expect_length(didx, 28)
    expect_equal(rownames(sce1)[uidx], gns[uidx])
    expect_equal(rownames(sce1)[midx], pns[midx])
    expect_equal(rownames(sce1)[didx], paste0(gns[didx], ".", pns[didx]))
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, primaryIdCol = "Accession",
                          secondaryIdCol = "Gene.Symbol")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$primaryIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$secondaryIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)
})
