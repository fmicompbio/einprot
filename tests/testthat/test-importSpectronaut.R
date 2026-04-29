test_that("Spectronaut import works", {
    pgpivot <- system.file("extdata", "spectronaut_example",
                           "SN_unit_test_subset.tsv",
                           package = "einprot")

    ## Fails with wrong arguments
    ## -------------------------------------------------------------------------
    expect_error(importSpectronaut(inFile = 1),
                 "'inFile' must be of class 'character'")
    expect_error(importSpectronaut(inFile = c(pgpivot, pgpivot)),
                 "'inFile' must have length 1")
    expect_error(importSpectronaut(inFile = "missing"),
                 "file.exists(inFile) is not TRUE", fixed = TRUE)

    expect_error(importSpectronaut(inFile = pgpivot, fileType = 1),
                 "'fileType' must be of class 'character'")
    expect_error(importSpectronaut(inFile = pgpivot, fileType = c("long_format", "pg_pivot")),
                 "'fileType' must have length 1")
    expect_error(importSpectronaut(inFile = pgpivot, fileType = "missing"),
                 "values in 'fileType' must be one of")

    expect_error(importSpectronaut(inFile = pgpivot, outLevel = 1),
                 "'outLevel' must be of class 'character'")
    expect_error(importSpectronaut(inFile = pgpivot, outLevel = c("pg", "pg")),
                 "'outLevel' must have length 1")
    expect_error(importSpectronaut(inFile = pgpivot, outLevel = "missing"),
                 "values in 'outLevel' must be one of")

    expect_error(importSpectronaut(inFile = pgpivot, outLevel = "pg",
                                   iColPattern = 1),
                 "'iColPattern' must be of class 'character'")
    expect_error(importSpectronaut(inFile = pgpivot, outLevel = "pg",
                                   iColPattern = c(".PG.Quantity$", ".PG.Quantity$")),
                 "'iColPattern' must have length 1")
    expect_error(importSpectronaut(inFile = pgpivot, outLevel = "pg",
                                   iColPattern = "missing"),
                 "Invalid iColPattern")

    expect_error(importSpectronaut(inFile = pgpivot, includeOnlySamples = 1),
                 "'includeOnlySamples' must be of class 'character'")
    expect_error(importSpectronaut(inFile = pgpivot, excludeSamples = 1),
                 "'excludeSamples' must be of class 'character'")
    expect_error(importSpectronaut(inFile = pgpivot, includeOnlySamples = "s1",
                             excludeSamples = "s2"),
                 "Please specify max one of includeOnlySamples and exclude")

    expect_error(importSpectronaut(inFile = pgpivot, stopIfEmpty = 1),
                 "'stopIfEmpty' must be of class 'logical'")
    expect_error(importSpectronaut(inFile = pgpivot, stopIfEmpty = c(TRUE, FALSE)),
                 "'stopIfEmpty' must have length 1")

    expect_error(importSpectronaut(inFile = pgpivot, fileType = "pg_pivot",
                                   outLevel = "pg", filtersDF = 1),
                 "'filtersDF' must be of class 'list'")

    ## Extract some values to compare to later
    pgtmp <- read.delim(pgpivot, sep = "\t", nrow = 20)

    ## Works with correct arguments
    ## -------------------------------------------------------------------------
    ## Without specifying samples to include/exclude - pg_pivot
    ## -------------------------------------------------------------------------
    out <- importSpectronaut(inFile = pgpivot, fileType = "pg_pivot",
                             outLevel = "pg", iColPattern = ".PG.Quantity$",
                             includeOnlySamples = "", excludeSamples = "",
                             stopIfEmpty = FALSE, aName = "PG.Quantity",
                             filtersDF = list(), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "PG.Quantity")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 6L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("PG.Quantity", "PG.IsSingleHit",
                   "PG.NrOfPrecursorsUsedForQuantification",
                   "PG.MS1Quantity", "PG.IBAQ", "PG.Cscore..Run.Wise.",
                   "PG.QValue..Run.Wise.", "PG.RunEvidenceCount"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "PG.Quantity")[, "LFQ_24min_500ng_2Th_3ms_A_02"],
        pgtmp[, "X.2..LFQ_24min_500ng_2Th_3ms_A_02.raw.PG.Quantity"], ignore_attr = TRUE)
    infoCols <- c("PG.Qvalue", "PG.Cscore", "PG.ProteinGroups", "PG.Genes",
                  "PG.ProteinDescriptions",
                  "PG.Completeness", "PG.CV")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to include - pg_pivot
    ## -------------------------------------------------------------------------
    out <- importSpectronaut(inFile = pgpivot, fileType = "pg_pivot",
                             outLevel = "pg", iColPattern = ".PG.Quantity$",
                             includeOnlySamples = "_A_", excludeSamples = "",
                             stopIfEmpty = FALSE, aName = "PG.Quantity",
                             filtersDF = list(), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "PG.Quantity")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("PG.Quantity", "PG.IsSingleHit",
                   "PG.NrOfPrecursorsUsedForQuantification",
                   "PG.MS1Quantity", "PG.IBAQ", "PG.Cscore..Run.Wise.",
                   "PG.QValue..Run.Wise.", "PG.RunEvidenceCount"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "PG.Quantity")[, "LFQ_24min_500ng_2Th_3ms_A_02"],
        pgtmp[, "X.2..LFQ_24min_500ng_2Th_3ms_A_02.raw.PG.Quantity"], ignore_attr = TRUE)
    infoCols <- c("PG.Qvalue", "PG.Cscore", "PG.ProteinGroups", "PG.Genes",
                  "PG.ProteinDescriptions",
                  "PG.Completeness", "PG.CV")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }

    ## -------------------------------------------------------------------------
    ## Specifying samples to exclude - pg_pivot
    ## -------------------------------------------------------------------------
    out <- importSpectronaut(inFile = pgpivot, fileType = "pg_pivot",
                             outLevel = "pg", iColPattern = ".PG.Quantity$",
                             includeOnlySamples = "", excludeSamples = "_B_",
                             stopIfEmpty = FALSE, aName = "PG.Quantity",
                             filtersDF = list(), nrows = 20)
    expect_type(out, "list")
    expect_named(out, c("sce", "aName"))
    expect_equal(out$aName, "PG.Quantity")
    expect_equal(nrow(out$sce), 20)
    expect_equal(ncol(out$sce), 3L)
    expect_s4_class(out$sce, "SingleCellExperiment")
    expect_equal(rownames(out$sce), as.character(seq_len(20)))
    expect_equal(SummarizedExperiment::assayNames(out$sce),
                 c("PG.Quantity", "PG.IsSingleHit",
                   "PG.NrOfPrecursorsUsedForQuantification",
                   "PG.MS1Quantity", "PG.IBAQ", "PG.Cscore..Run.Wise.",
                   "PG.QValue..Run.Wise.", "PG.RunEvidenceCount"))
    expect_equal(SummarizedExperiment::assay(
        out$sce, "PG.Quantity")[, "LFQ_24min_500ng_2Th_3ms_A_02"],
        pgtmp[, "X.2..LFQ_24min_500ng_2Th_3ms_A_02.raw.PG.Quantity"], ignore_attr = TRUE)
    infoCols <- c("PG.Qvalue", "PG.Cscore", "PG.ProteinGroups", "PG.Genes",
                  "PG.ProteinDescriptions",
                  "PG.Completeness", "PG.CV")
    expect_true(all(infoCols %in%
                        colnames(SummarizedExperiment::rowData(out$sce))))
    for (ic in infoCols) {
        expect_equal(SummarizedExperiment::rowData(out$sce)[, ic],
                     pgtmp[, ic])
    }
})
