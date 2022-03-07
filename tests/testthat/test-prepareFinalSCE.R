test_that("assembling the SCE works", {
    out <- runTest(
        sce = sce_mq_final, comparison = c("Adnp", "RBC_ctrl"), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )

    args0 <- list(
        sce = sce_mq_final,
        baseFileName = NULL,
        featureCollections = out$featureCollections,
        expType = "MaxQuant"
    )

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## sce
    args <- args0
    args$sce <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'sce' must be of class 'SingleCellExperiment'")

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'baseFileName' must be of class 'character'")

    ## featureCollections
    args <- args0
    args$featureCollections <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'featureCollections' must be of class 'list'")

    ## expType
    args <- args0
    args$expType <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'expType' must be of class 'character'")
    args$expType <- c("MaxQuant", "ProteomeDiscoverer")
    expect_error(do.call(prepareFinalSCE, args),
                 "'expType' must have length 1")
    args$expType <- "missing"
    expect_error(do.call(prepareFinalSCE, args),
                 "All values in 'expType' must be one of")


    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    sce <- do.call(prepareFinalSCE, args0)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(nrow(sce), 70)
    expect_equal(ncol(sce), 9)
    expect_true(all(c("iBAQ", "log2_iBAQ", "log2_iBAQ_withNA", "imputed_iBAQ",
                      "MS.MS.Count", "LFQ.intensity", "Intensity",
                      "Sequence.coverage", "Unique.peptides", "Peptides") %in%
                        SummarizedExperiment::assayNames(sce)))
    expect_true(all(c("sample", "group") %in%
                        colnames(SummarizedExperiment::colData(sce))))
    expect_true(all(c("Gene.names", "Majority.protein.IDs") %in%
                        colnames(SummarizedExperiment::rowData(sce))))
    expect_false(any(c("Peptide.counts.all", "Peptide.counts.razor.unique",
                       "Peptide.counts.unique", "Fasta.headers",
                       "Peptide.IDs", "Peptide.is.razor", "Mod.peptide.IDs",
                       "Evidence.IDs", "MS.MS.IDs", "Best.MS.MS",
                       "Sequence.lengths", "Oxidation.M.site.IDs",
                       "Oxidation.M.site.positions") %in%
                         colnames(SummarizedExperiment::rowData(sce))))
})
