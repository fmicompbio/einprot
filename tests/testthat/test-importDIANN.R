test_that("multiplication works", {
    pgmat <- system.file("extdata", "diann_example", "PXD028735.pg_matrix.tsv",
                         package = "einprot")
    prmat <- system.file("extdata", "diann_example", "PXD028735.pr_matrix.tsv",
                         package = "einprot")
    mainrep <- system.file("extdata", "diann_example", "PXD028735.report.tsv",
                           package = "einprot")

    ## Fails with wrong arguments
    ## -------------------------------------------------------------------------
    expect_error(importDIANN(inFile = 1),
                 "'inFile' must be of class 'character'")
    expect_error(importExperiment(inFile = c(pgmat, prmat)),
                 "'inFile' must have length 1")
    expect_error(importExperiment(inFile = "missing"),
                 "file.exists(inFile) is not TRUE", fixed = TRUE)

    expect_error(importDIANN(inFile = pgmat, fileType = 1),
                 "'fileType' must be of class 'character'")
    expect_error(importDIANN(inFile = pgmat, fileType = c("pg_matrix", "pr_matrix")),
                 "'fileType' must have length 1")
    expect_error(importDIANN(inFile = pgmat, fileType = "missing"),
                 "values in 'fileType' must be one of")

    expect_error(importDIANN(inFile = pgmat, outLevel = 1),
                 "'outLevel' must be of class 'character'")
    expect_error(importDIANN(inFile = pgmat, outLevel = c("pg", "pr")),
                 "'outLevel' must have length 1")
    expect_error(importDIANN(inFile = pgmat, outLevel = "missing"),
                 "values in 'outLevel' must be one of")

    expect_error(importDIANN(inFile = pgmat, fileType = "pg_matrix",
                             outLevel = "pr"),
                 "To obtain precursor output, please provide")
    expect_error(importDIANN(inFile = prmat, fileType = "pr_matrix",
                             outLevel = "pg"),
                 "To obtain protein group output, please provide")

    expect_error(importDIANN(inFile = pgmat, includeOnlySamples = 1),
                 "'includeOnlySamples' must be of class 'character'")
    expect_error(importDIANN(inFile = pgmat, excludeSamples = 1),
                 "'excludeSamples' must be of class 'character'")
    expect_error(importDIANN(inFile = pgmat, includeOnlySamples = "s1",
                             excludeSamples = "s2"),
                 "Please specify max one of includeOnlySamples and exclude")

    expect_error(importDIANN(inFile = pgmat, stopIfEmpty = 1),
                 "'stopIfEmpty' must be of class 'logical'")
    expect_error(importDIANN(inFile = pgmat, stopIfEmpty = c(TRUE, FALSE)),
                 "'stopIfEmpty' must have length 1")

    expect_error(importDIANN(inFile = pgmat, aName = 1),
                 "'aName' must be of class 'character'")
    expect_error(importDIANN(inFile = pgmat, aName = c("N1", "N2")),
                 "'aName' must have length 1")
    expect_error(importDIANN(inFile = mainrep, fileType = "main_report",
                             outLevel = "pg", aName = "MaxLFQ"),
                 "aName %in% colnames(tmp) is not TRUE", fixed = TRUE)

    ## Extract some values to compare to later
    pgtmp <- read.delim(pgmat, sep = "\t", nrow = 20)
    prtmp <- read.delim(prmat, sep = "\t", nrow = 20)
    maintmp <- read.delim(mainrep, sep = "\t")

    ## Works with correct arguments
    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude - pg_matrix
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = pgmat, fileType = "pg_matrix", outLevel = "pg",
                       includeOnlySamples = "", excludeSamples = "",
                       stopIfEmpty = FALSE, aName = "MaxLFQ", nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "MaxLFQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 6L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce), "MaxLFQ")
    expect_equal(SummarizedExperiment::assay(
        out$sce, "MaxLFQ")[, "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
        pgtmp[, "X.scratch.cpanse.PXD028735.dia.LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01.mzML"], ignore_attr = TRUE)
    infoCols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                  "First.Protein.Description")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude - pr_matrix
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = prmat, fileType = "pr_matrix", outLevel = "pr",
                       includeOnlySamples = "", excludeSamples = "",
                       stopIfEmpty = FALSE, aName = "MaxLFQ", nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "MaxLFQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 6L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce), "MaxLFQ")
    expect_equal(SummarizedExperiment::assay(
        out$sce, "MaxLFQ")[, "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
        prtmp[, "X.scratch.cpanse.PXD028735.dia.LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01.mzML"], ignore_attr = TRUE)
    infoCols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                  "First.Protein.Description", "Proteotypic", "Stripped.Sequence",
                  "Modified.Sequence", "Precursor.Charge", "Precursor.Id")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     prtmp[, ic])
    }
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude - main_report/pg
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = mainrep, fileType = "main_report",
                       outLevel = "pg",
                       includeOnlySamples = "", excludeSamples = "",
                       stopIfEmpty = FALSE, aName = "PG.MaxLFQ")
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "PG.MaxLFQ")
    expect_equal(nrow(out$sce), 72L) ## length(unique(maintmp$Protein.Group))
    expect_equal(ncol(out$sce), 6L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(sort(rownames(out$sce)), sort(unique(maintmp$Protein.Group)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("PG.MaxLFQ", "PG.Quantity", "PG.Normalised", "PG.Q.Value",
                   "Global.PG.Q.Value", "Lib.PG.Q.Value", "Protein.Ids"))

    ## Manual checks
    aNames <- SummarizedExperiment::assayNames(out$sce)
    sub1 <- subset(maintmp, Run == "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01")
    for (an in aNames) {
        sub2 <- sub1[match(c("P0A7J3", "O94826", "G5EA36;I3L328"),
                           sub1$Protein.Group), an]
        sub2[is.na(sub2)] <- 0
        expect_equal(SummarizedExperiment::assay(
            out$sce, an)[c("P0A7J3", "O94826", "G5EA36;I3L328"),
                         "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
            sub2, ignore_attr = TRUE)
    }

    infoCols <- c("Protein.Group")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude - main_report/pr
    ## -------------------------------------------------------------------------
    expect_error(
        out <- importDIANN(inFile = mainrep, fileType = "main_report",
                           outLevel = "pr",
                           includeOnlySamples = "", excludeSamples = "",
                           stopIfEmpty = FALSE, aName = "Precursor.Quantity"),
        "Not yet implemented")

    ## -------------------------------------------------------------------------
    ## Specifying samples to include - pg_matrix
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = pgmat, fileType = "pg_matrix", outLevel = "pg",
                       includeOnlySamples = "Condition_A", excludeSamples = "",
                       stopIfEmpty = FALSE, aName = "MaxLFQ", nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "MaxLFQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce), "MaxLFQ")
    expect_equal(SummarizedExperiment::assay(
        out$sce, "MaxLFQ")[, "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
        pgtmp[, "X.scratch.cpanse.PXD028735.dia.LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01.mzML"], ignore_attr = TRUE)
    infoCols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                  "First.Protein.Description")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

    ## -------------------------------------------------------------------------
    ## Specifying samples to exclude - pg_matrix
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = pgmat, fileType = "pg_matrix", outLevel = "pg",
                       includeOnlySamples = "", excludeSamples = "Condition_B",
                       stopIfEmpty = FALSE, aName = "MaxLFQ", nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "MaxLFQ")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce), "MaxLFQ")
    expect_equal(SummarizedExperiment::assay(
        out$sce, "MaxLFQ")[, "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
        pgtmp[, "X.scratch.cpanse.PXD028735.dia.LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01.mzML"], ignore_attr = TRUE)
    infoCols <- c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                  "First.Protein.Description")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

    ## -------------------------------------------------------------------------
    ## Specifying samples to include - main_report/pg
    ## -------------------------------------------------------------------------
    out <- importDIANN(inFile = mainrep, fileType = "main_report",
                       outLevel = "pg",
                       includeOnlySamples = "Condition_A", excludeSamples = "",
                       stopIfEmpty = FALSE, aName = "PG.MaxLFQ")
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "PG.MaxLFQ")
    expect_equal(nrow(out$sce), 62L) ## length(unique(maintmp$Protein.Group[grep("Condition_A", maintmp$Run)]))
    expect_equal(ncol(out$sce), 3L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(sort(rownames(out$sce)), sort(unique(maintmp$Protein.Group[grep("Condition_A", maintmp$Run)])))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("PG.MaxLFQ", "PG.Quantity", "PG.Normalised", "PG.Q.Value",
                   "Global.PG.Q.Value", "Lib.PG.Q.Value", "Protein.Ids"))

    ## Manual checks
    aNames <- SummarizedExperiment::assayNames(out$sce)
    sub1 <- subset(maintmp, Run == "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01")
    for (an in aNames) {
        sub2 <- sub1[match(c("P0A7J3", "O94826", "G5EA36;I3L328"),
                           sub1$Protein.Group), an]
        sub2[is.na(sub2)] <- 0
        expect_equal(SummarizedExperiment::assay(
            out$sce, an)[c("P0A7J3", "O94826", "G5EA36;I3L328"),
                         "LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01"],
            sub2, ignore_attr = TRUE)
    }

    infoCols <- c("Protein.Group")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    expect_false(any(grepl("scratch", colnames(out$sce))))
    expect_false(any(grepl("mzML", colnames(out$sce))))

})
