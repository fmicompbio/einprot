test_that("fixing feature IDs works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce

    ## Fail with wrong arguments
    expect_error(fixFeatureIds(sce = 1, geneIdCol = "Gene.names",
                               proteinIdCol = "Majority.protein.IDs"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = 1,
                               proteinIdCol = "Majority.protein.IDs"),
                 "'geneIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = c("Gene.names", "Gene.names"),
                               proteinIdCol = "Majority.protein.IDs"),
                 "'geneIdCol' must have length 1")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = "missing",
                               proteinIdCol = "Majority.protein.IDs"),
                 "All values in 'geneIdCol' must be one of")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                               proteinIdCol = 1),
                 "'proteinIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                               proteinIdCol = c("Majority.protein.IDs",
                                                "Majority.protein.IDs")),
                 "'proteinIdCol' must have length 1")
    expect_error(fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                               proteinIdCol = "missing"),
                 "All values in 'proteinIdCol' must be one of")

    ## Test that it does the right thing
    ## All gene and protein names are unique
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")

    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
    expect_equal(rownames(sce1), gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)

    ## Missing gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce[36:45, ]
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs

    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
    midx <- which(gns == "")  ## indices for missing gene IDs
    eidx <- which(gns != "")  ## indices for present gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    expect_equal(rownames(sce1)[eidx], gns[eidx])
    expect_equal(rownames(sce1)[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[eidx], gns[eidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)

    ## Duplicated gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
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
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, gns)
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
    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Gene.Symbol",
                          proteinIdCol = "Accession")
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
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce, geneIdCol = "Accession",
                          proteinIdCol = "Gene.Symbol")
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$geneIdSingle, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$proteinIdSingle, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, pns)
})
