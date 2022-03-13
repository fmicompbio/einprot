test_that("assembling the SCE works", {
    out <- runTest(
        sce = sce_mq_final, comparison = c("Adnp", "RBC_ctrl"), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )

    args0_mq <- list(
        sce = sce_mq_final,
        baseFileName = tempfile(),
        featureCollections = out$featureCollections,
        expType = "MaxQuant"
    )

    args0_pd <- list(
        sce = sce_pd_final,
        baseFileName = tempfile(),
        featureCollections = fcoll_pd_final,
        expType = "ProteomeDiscoverer"
    )

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## sce
    args <- args0_mq
    args$sce <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'sce' must be of class 'SingleCellExperiment'")

    ## baseFileName
    args <- args0_mq
    args$baseFileName <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'baseFileName' must be of class 'character'")

    ## featureCollections
    args <- args0_mq
    args$featureCollections <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'featureCollections' must be of class 'list'")

    ## expType
    args <- args0_mq
    args$expType <- 1
    expect_error(do.call(prepareFinalSCE, args),
                 "'expType' must be of class 'character'")
    args$expType <- c("MaxQuant", "ProteomeDiscoverer")
    expect_error(do.call(prepareFinalSCE, args),
                 "'expType' must have length 1")
    args$expType <- "missing"
    expect_error(do.call(prepareFinalSCE, args),
                 "All values in 'expType' must be one of")


    ## Works with correct arguments - MaxQuant
    ## --------------------------------------------------------------------- ##
    sce <- do.call(prepareFinalSCE, args0_mq)
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
    expect_true(file.exists(paste0(args0_mq$baseFileName,
                                   "_sce_extra_annots.tsv")))

    ## Works with correct arguments - MaxQuant, add log2_iBAQ_withNA
    ## --------------------------------------------------------------------- ##
    args <- args0_mq
    SummarizedExperiment::assays(args$sce)[["log2_iBAQ_withNA"]] <- NULL
    sce <- do.call(prepareFinalSCE, args)
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

    ## Works with correct arguments - ProteomeDiscoverer
    ## --------------------------------------------------------------------- ##
    sce <- do.call(prepareFinalSCE, args0_pd)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(nrow(sce), 70)
    expect_equal(ncol(sce), 16)
    expect_true(all(c("Abundance", "Abundances.count", "Abundances.normalized",
                      "Abundances.grouped.count", "Abundances.grouped.CV",
                      "Abundances.grouped", "log2_Abundance",
                      "log2_Abundance_withNA", "imputed_Abundance") %in%
                        SummarizedExperiment::assayNames(sce)))
    expect_true(all(c("sample", "group") %in%
                        colnames(SummarizedExperiment::colData(sce))))
    expect_true(all(c("Gene.Symbol", "Accession") %in%
                        colnames(SummarizedExperiment::rowData(sce))))
    expect_false(any(c("Sequence", "GO.Accessions",
                       "Proteins.Unique.Sequence.ID") %in%
                         colnames(SummarizedExperiment::rowData(sce))))
    expect_true(file.exists(paste0(args0_pd$baseFileName,
                                   "_sce_extra_annots.tsv")))
})
