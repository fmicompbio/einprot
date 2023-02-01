test_that("importExperiment works", {
    ## --------------------------------------------------------------------- ##
    ## Import MQ data
    ## --------------------------------------------------------------------- ##
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    mqSamples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
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

    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 9)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP08"), colnames(tmp))))
    }

    ## Without escaping the periods
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ.",
                             nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 9)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP08"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to include
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## Without escaping period
    out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to exclude
    ## -------------------------------------------------------------------------
    out <- importExperiment(
        inFile = mqFile, iColPattern = "^iBAQ\\.",
        excludeSamples = setdiff(mqSamples, c("Chd4BF_IP08", "Adnp_IP06",
                                              "Chd4BF_IP09")), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## Without escaping period
    out <- importExperiment(
        inFile = mqFile, iColPattern = "^iBAQ.",
        excludeSamples = setdiff(mqSamples, c("Chd4BF_IP08", "Adnp_IP06",
                                              "Chd4BF_IP09")), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "iBAQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Different iColPattern
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = mqFile,
                            iColPattern = "^LFQ\\.intensity\\.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "LFQ.intensity")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("LFQ.intensity", "MS.MS.Count",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "iBAQ",
                   "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## Without escaping period
    out <- importExperiment(inFile = mqFile,
                            iColPattern = "^LFQ.intensity.",
                            includeOnlySamples = c("Chd4BF_IP08", "Adnp_IP06",
                                                   "Chd4BF_IP09"), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "LFQ.intensity")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("Adnp_IP06", "Chd4BF_IP08",
                                      "Chd4BF_IP09"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("LFQ.intensity", "MS.MS.Count",
                   "Intensity", "Sequence.coverage", "Unique.peptides",
                   "Razor.unique.peptides", "Peptides", "iBAQ",
                   "Identification.type"))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "iBAQ")[, "Adnp_IP06"] == tmp$iBAQ.Adnp_IP06))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "LFQ.intensity")[, "Chd4BF_IP08"] ==
            tmp$LFQ.intensity.Chd4BF_IP08))
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Razor.unique.peptides")[, "Chd4BF_IP09"] ==
            tmp$Razor...unique.peptides.Chd4BF_IP09))
    expect_true(all(c("Peptides", "Unique.peptides", "Majority.protein.IDs",
                      "Gene.names", "Score", "Potential.contaminant",
                      "Reverse", "Only.identified.by.site") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("iBAQ", "MS.MS.Count", "LFQ.intensity",
                  "Intensity", "Sequence.coverage", "Unique.peptides",
                  "Razor.unique.peptides", "Peptides", "Identification.type")) {
        if (nms == "Razor.unique.peptides") {
            nmstmp <- "Razor...unique.peptides"
        } else {
            nmstmp <- nms
        }
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP07"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".Chd4BF_IP08"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".Chd4BF_IP07"), colnames(tmp))))
    }

    ## --------------------------------------------------------------------- ##
    ## Import PD data
    ## --------------------------------------------------------------------- ##
    pdFile <- system.file("extdata", "pdtmt_example",
                          "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                          package = "einprot")
    pdSamples <- c("HIS4KO_S05", "HIS4KO_S06", "HIS4KO_S07", "HIS4KO_S08",
                   "MET6KO_S01", "MET6KO_S02", "MET6KO_S03", "MET6KO_S04",
                   "URA2KO_S09", "URA2KO_S10", "URA2KO_S11", "URA2KO_S12",
                   "WT_S13", "WT_S14", "WT_S15", "WT_S16")

    ## Read plan text file
    tmp <- read.delim(pdFile, sep = "\t", nrow = 20)

    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = pdFile,
                            iColPattern = "^Abundance\\.F.+\\.Sample\\.",
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundance")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 16)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundance", "Abundances.count", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S01"), colnames(tmp))))
    }

    ## Without escaping periods
    out <- importExperiment(inFile = pdFile,
                            iColPattern = "^Abundance.F.+.Sample.",
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundance")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 16)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundance", "Abundances.count", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S01"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to include
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = pdFile,
                            iColPattern = "^Abundance\\.F.+\\.Sample\\.",
                            includeOnlySamples = c("HIS4KO_S06", "WT_S16",
                                                   "MET6KO_S01", "HIS4KO_S07",
                                                   "URA2KO_S10"),
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundance")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 5)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("HIS4KO_S06", "HIS4KO_S07", "MET6KO_S01",
                                      "URA2KO_S10", "WT_S16"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundance", "Abundances.count", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S02"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S02"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S02"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to exclude
    ## -------------------------------------------------------------------------
    out <- importExperiment(inFile = pdFile,
                            iColPattern = "^Abundance\\.F.+\\.Sample\\.",
                            excludeSamples = setdiff(
                                pdSamples, c("HIS4KO_S06", "WT_S16",
                                  "MET6KO_S01", "HIS4KO_S07", "URA2KO_S10")),
                            nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundance")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 5)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(colnames(out$sce), c("HIS4KO_S06", "HIS4KO_S07", "MET6KO_S01",
                                      "URA2KO_S10", "WT_S16"))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundance", "Abundances.count", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S02"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S02"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S02"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Different iColPattern
    ## -------------------------------------------------------------------------
    out <- importExperiment(
        inFile = pdFile, iColPattern = "^Abundances\\.Count\\.F.+\\.Sample\\.",
        nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundances.count")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 16)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundances.count", "Abundance", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S01"), colnames(tmp))))
    }

    ## Without escaping periods
    out <- importExperiment(
        inFile = pdFile, iColPattern = "^Abundances.Count.F.+.Sample.",
        nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundances.count")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 16)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundances.count", "Abundance", "Abundances.normalized",
                   "Abundances.grouped.count", "Abundances.grouped.CV",
                   "Abundances.grouped"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S01"), colnames(tmp))))
    }

    ## -------------------------------------------------------------------------
    ## Another different iColPattern (currently not supported)
    ## -------------------------------------------------------------------------
    expect_warning(out <- importExperiment(
        inFile = pdFile, iColPattern = "^Abundances\\.Grouped\\.", nrows = 20),
        "Note that the specified iColPattern may match different")
    expect_error(out <- importExperiment(
        inFile = pdFile, iColPattern = "\\.Spectral\\.Count$", nrows = 20))

    expect_warning(out <- importExperiment(
        inFile = pdFile, iColPattern = "^Abundances.Grouped.", nrows = 20),
        "Note that the specified iColPattern may match different")
    expect_error(out <- importExperiment(
        inFile = pdFile, iColPattern = ".Spectral.Count$", nrows = 20))

    expect_warning(out <- importExperiment(
        inFile = pdFile, iColPattern = "^Abundances.Grouped.",
        nrows = 20),
        "Note that the specified iColPattern may match different")
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "Abundances.grouped")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 16)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("Abundances.grouped", "Abundances.count", "Abundance",
                   "Abundances.normalized", "Abundances.grouped.count",
                   "Abundances.grouped.CV"))
    idx_not_na <- c(1, 3, 5, 8, 9, 10, 12, 13, 14, 15, 16, 18, 20)
    idx_na <- setdiff(seq_len(20), idx_not_na)
    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_not_na, "HIS4KO_S06"] ==
            tmp$Abundances.Grouped.HIS4KO_S06[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped")[idx_na, "HIS4KO_S06"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_not_na, "WT_S16"] ==
            tmp$Abundances.Grouped.CV.WT_S16[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.CV")[idx_na, "WT_S16"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Grouped.Count.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.grouped.count")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_not_na, "HIS4KO_S07"] ==
            tmp$Abundance.F12.129C.Sample.HIS4KO_S07[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundance")[idx_na, "HIS4KO_S07"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_not_na, "MET6KO_S01"] ==
            tmp$Abundances.Normalized.MET6KO_S01[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.normalized")[idx_na, "MET6KO_S01"])))

    expect_true(all(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_not_na, "URA2KO_S10"] ==
            tmp$Abundances.Count.F12.131N.Sample.URA2KO_S10[idx_not_na]))
    expect_true(all(is.na(SummarizedExperiment::assay(
        out$sce, "Abundances.count")[idx_na, "URA2KO_S10"])))

    expect_true(all(SummarizedExperiment::rowData(out$sce)$Accession ==
                        tmp$Accession))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Number.of.Peptides ==
                        tmp$Number.of.Peptides))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Gene.Symbol ==
                        tmp$Gene.Symbol))
    expect_true(all(SummarizedExperiment::rowData(out$sce)$Modifications ==
                        tmp$Modifications))
    expect_true(all(c("Accession", "Number.of.Peptides",
                      "Score.Sequest.HT.Sequest.HT", "Gene.Symbol") %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))

    ## Check that no sample-specific columns remain in the rowData
    ## (but they should be there in the temp data loaded above)
    for (nms in c("Abundance", "Abundances.count", "Abundances.normalized",
                  "Abundances.grouped.count", "Abundances.grouped.CV",
                  "Abundances.grouped")) {
        nmstmp <- dplyr::case_when(
            nms == "Abundance" ~ "Abundance.F.+.Sample",
            nms == "Abundances.count" ~ "Abundances.Count.F.+.Sample",
            nms == "Abundances.normalized" ~ "Abundances.Normalized.F.+.Sample",
            nms == "Abundances.grouped.count" ~ "Abundances.Grouped.Count",
            nms == "Abundances.grouped.CV" ~ "Abundances.Grouped.CV.in.Percent",
            nms == "Abundances.grouped" ~ "Abundances.Grouped"
        )
        expect_false(any(grepl(paste0(nms, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nms, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".MET6KO_S01"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_false(any(grepl(paste0(nmstmp, ".URA2KO_S10"),
                               colnames(SummarizedExperiment::rowData(out$sce)))))
        expect_true(any(grepl(paste0(nmstmp, ".MET6KO_S01"), colnames(tmp))))
    }
})
