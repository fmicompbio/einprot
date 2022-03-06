test_that("importExperiment works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
                 "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
                 "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")

    ## Fails with wrong arguments
    expect_error(importExperiment(inFile = 1, iColPattern = "^iBAQ\\."),
                 "'inFile' must be of class 'character'")
    expect_error(importExperiment(inFile = c(mqFile, mqFile),
                                  iColPattern = "^iBAQ\\."),
                 "'inFile' must have length 1")
    expect_error(importExperiment(inFile = "missing",
                                  iColPattern = "^iBAQ\\."),
                 "file.exists(inFile) is not TRUE", fixed = TRUE)
    expect_error(importExperiment(inFile = mqFile, iColPattern = 1),
                 "'iColPattern' must be of class 'character'")
    expect_error(importExperiment(inFile = mqFile,
                                  iColPattern = c("^iBAQ\\.", "^iBAQ\\.")),
                 "'iColPattern' must have length 1")
    expect_error(importExperiment(inFile = mqFile, iColPattern = "iBAQ"),
                 "All values in 'iColPattern' must be one of")
    ## Valid iColPattern, but not for MQ
    expect_error(importExperiment(
        inFile = mqFile, iColPattern = "^Abundances\\.Count\\.F.+\\.Sample\\."),
        "unable to find an inherited method for function")
    expect_error(importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                                  includeOnlySamples = 1),
                 "'includeOnlySamples' must be of class 'character'")
    expect_error(importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                                  excludeSamples = 1),
                 "'excludeSamples' must be of class 'character'")
    expect_error(importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                                  includeOnlySamples = "Adnp_IP04",
                                  excludeSamples = "Adnp_IP05"),
                 "Please specify max one of includeOnlySamples and exclude")

    ## Extract some values to compare to later
    tmp <- read.delim(mqFile, sep = "\t", nrow = 20)
    ibaq_adnp_ip06_3 <- tmp$iBAQ.Adnp_IP06[3]
    lfq_chd4bf_ip08_8 <- tmp$LFQ.intensity.Chd4BF_IP08[8]
    razor_up_chd4bf_ip09_4 <- tmp$Razor...unique.peptides.Chd4BF_IP09[4]

    ## Without specifying samples to include/exclude
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "iBAQ")[3, "Adnp_IP06"],
                 ibaq_adnp_ip06_3)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[8, "Chd4BF_IP08"],
                 lfq_chd4bf_ip08_8)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[4, "Chd4BF_IP09"],
                 razor_up_chd4bf_ip09_4)
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Specifying samples to include
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "iBAQ")[3, "Adnp_IP06"],
        ibaq_adnp_ip06_3)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[8, "Chd4BF_IP08"],
        lfq_chd4bf_ip08_8)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[4, "Chd4BF_IP09"],
        razor_up_chd4bf_ip09_4)
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Specifying samples to exclude
    out <- importExperiment(
        inFile = mqFile, iColPattern = "^iBAQ\\.",
        excludeSamples = setdiff(samples, c("Chd4BF_IP08", "Adnp_IP06",
                                            "Chd4BF_IP09")), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "iBAQ")[3, "Adnp_IP06"],
        ibaq_adnp_ip06_3)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[8, "Chd4BF_IP08"],
        lfq_chd4bf_ip08_8)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[4, "Chd4BF_IP09"],
        razor_up_chd4bf_ip09_4)
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Different iColPattern
    out <- importExperiment(inFile = mqFile, iColPattern = "^LFQ\\.intensity\\.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "LFQ.intensity")
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("LFQ.intensity", "MS.MS.Count",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "iBAQ",
                   "Identification.type"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "iBAQ")[3, "Adnp_IP06"],
        ibaq_adnp_ip06_3)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[8, "Chd4BF_IP08"],
        lfq_chd4bf_ip08_8)
    expect_equal(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[4, "Chd4BF_IP09"],
        razor_up_chd4bf_ip09_4)
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

})
