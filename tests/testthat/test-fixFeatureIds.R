test_that("fixing feature IDs works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
                 "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
                 "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
    ecol <- paste0("iBAQ.", samples)
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "Intensity",
                                    sep = "\t", nrows = 25)

    ## Fail with wrong arguments
    expect_error(fixFeatureIds(qft = 1, geneIdCol = "Gene.names",
                               proteinIdCol = "Majority.protein.IDs"),
                 "'qft' must be of class 'QFeatures'")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = 1,
                               proteinIdCol = "Majority.protein.IDs"),
                 "'geneIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = c("Gene.names", "Gene.names"),
                               proteinIdCol = "Majority.protein.IDs"),
                 "'geneIdCol' must have length 1")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = "missing",
                               proteinIdCol = "Majority.protein.IDs"),
                 "All values in 'geneIdCol' must be one of")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                               proteinIdCol = 1),
                 "'proteinIdCol' must be of class 'character'")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                               proteinIdCol = c("Majority.protein.IDs",
                                                "Majority.protein.IDs")),
                 "'proteinIdCol' must have length 1")
    expect_error(fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                               proteinIdCol = "missing"),
                 "All values in 'proteinIdCol' must be one of")

    ## Test that it does the right thing
    ## All gene and protein names are unique
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "Intensity",
                                    sep = "\t", nrows = 25)
    gns <- SummarizedExperiment::rowData(qft[[1]])$Gene.names
    pns <- SummarizedExperiment::rowData(qft[[1]])$Majority.protein.IDs

    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
    expect_equal(rownames(qft1[[1]]), vapply(strsplit(gns, ";"),
                                             .subset, 1, FUN.VALUE = ""))
    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(qft1[[1]]), vapply(strsplit(pns, ";"),
                                             .subset, 1, FUN.VALUE = ""))

    ## Missing gene names
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "Intensity",
                                    sep = "\t", nrows = 45)[36:45, ]
    gns <- SummarizedExperiment::rowData(qft[[1]])$Gene.names
    pns <- SummarizedExperiment::rowData(qft[[1]])$Majority.protein.IDs

    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
    midx <- which(gns == "")  ## indices for missing gene IDs
    eidx <- which(gns != "")  ## indices for present gene IDs
    expect_equal(rownames(qft1[[1]])[eidx], vapply(strsplit(gns[eidx], ";"),
                                                   .subset, 1, FUN.VALUE = ""))
    expect_equal(rownames(qft1[[1]])[midx], vapply(strsplit(pns[midx], ";"),
                                                   .subset, 1, FUN.VALUE = ""))
    ## Protein names are still there and unique
    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(qft1[[1]]), vapply(strsplit(pns, ";"),
                                             .subset, 1, FUN.VALUE = ""))

    ## Duplicated gene names
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "Intensity",
                                    sep = "\t", nrows = 45)
    gns <- SummarizedExperiment::rowData(qft[[1]])$Gene.names
    pns <- SummarizedExperiment::rowData(qft[[1]])$Majority.protein.IDs
    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Gene.names",
                          proteinIdCol = "Majority.protein.IDs")
    midx <- which(gns == "")  ## indices for missing gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    didx <- setdiff(which(gns %in% gns[duplicated(gns)]), midx)  ## indices for duplicated gene IDs
    uidx <- setdiff(seq_along(gns), c(midx, didx))  ## indices for unique gene IDs
    expect_length(uidx, 38)
    expect_length(midx, 4)
    expect_length(didx, 3)
    expect_equal(rownames(qft1[[1]])[uidx], gns[uidx])
    expect_equal(rownames(qft1[[1]])[midx], pns[midx])
    expect_equal(rownames(qft1[[1]])[didx], paste0(gns[didx], ".", pns[didx]))
    ## Protein names are still there and unique
    qft1 <- fixFeatureIds(qft = qft, geneIdCol = "Majority.protein.IDs",
                          proteinIdCol = "Gene.names")
    expect_equal(rownames(qft1[[1]]), pns)
})
