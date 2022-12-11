test_that("fixing feature IDs works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce

    ## Fail with wrong arguments
    expect_error(fixFeatureIds(sce = 1, idCol = "Gene.names",
                               labelCol = "Gene.names",
                               geneIdCol = "Gene.names",
                               proteinIdCol = "Majority.protein.IDs"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(fixFeatureIds(sce = sce, idCol = 1,
                               labelCol = "Gene.names",
                               geneIdCol = "Gene.names",
                               proteinIdCol = "Majority.protein.IDs"),
                 "'idCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, idCol = "Gene.names",
                               labelCol = 1,
                               geneIdCol = "Gene.names",
                               proteinIdCol = "Majority.protein.IDs"),
                 "'labelCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, idCol = "Gene.names",
                               labelCol = "Gene.names",
                               geneIdCol = 1,
                               proteinIdCol = "Majority.protein.IDs"),
                 "'geneIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce, idCol = "Gene.names",
                               labelCol = "Gene.names",
                               geneIdCol = "Gene.names",
                               proteinIdCol = 1),
                 "'proteinIdCol' must be of class 'character'")

    ## Test that it does the right thing
    ## All gene and protein names are unique
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")

    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                          labelCol = "Majority.protein.IDs",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
    expect_equal(rownames(sce1), gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    ## Missing gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce[36:45, ]
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs

    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
    midx <- which(gns == "")  ## indices for missing gene IDs
    eidx <- which(gns != "")  ## indices for present gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    expect_equal(rownames(sce1)[eidx], gns[eidx])
    expect_equal(rownames(sce1)[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[eidx], gns[eidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## Duplicated gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
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
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs"))
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## Change the separator
    scesep <- sce
    rowData(scesep)$Majority.protein.IDs <- gsub(";", ",",
                                                rowData(scesep)$Majority.protein.IDs)
    rowData(scesep)$Gene.names <- gsub(";", ",",
                                       rowData(scesep)$Gene.names)
    sce1 <- fixFeatureIds(sce = scesep,
                          idCol = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names"),
                                                          splitSeparator = ","),
                          labelCol = "Gene.names",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.names",
                                                              separator = ","),
                          proteinIdCol = function(df) getFirstId(df, colName = "Majority.protein.IDs",
                                                                 separator = ","))
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## --------------------------------------------------------------------- ##
    ## PD data
    sce <- importExperiment(
        inFile = system.file("extdata", "pdtmt_example",
                             "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                             package = "einprot"),
        iColPattern = "^Abundance\\.F.+\\.Sample\\.")$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.Symbol
    pns <- SummarizedExperiment::rowData(sce)$Accession
    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Gene.Symbol", "Accession")),
                          labelCol = "Gene.Symbol",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.Symbol"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Accession"))
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
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(sce = sce,
                          idCol = function(df) combineIds(df, c("Accession", "Gene.Symbol")),
                          labelCol = "Gene.Symbol",
                          geneIdCol = function(df) getFirstId(df, colName = "Gene.Symbol"),
                          proteinIdCol = function(df) getFirstId(df, colName = "Accession"))
    expect_equal(rownames(sce1), pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
})
