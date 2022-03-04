test_that("assembling the SCE works", {
    out <- runTest(
        qft = qft_mq_final, comparison = c("Adnp", "RBC_ctrl"), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )

    args0 <- list(
        qft = qft_mq_final,
        aName = "iBAQ",
        testResults = out$res,
        iColPattern = "^iBAQ\\.",
        iColsAll = getIntensityColumns(mqFile, "^iBAQ\\.")$iColsAll,
        baseFileName = NULL,
        nbrNA = nbr_na_mq,
        featureCollections = out$featureCollections,
        expType = "MaxQuant"
    )

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## qft
    args <- args0
    args$qft <- 1
    expect_error(do.call(assembleSCE, args),
                 "'qft' must be of class 'QFeatures'")

    ## aName
    args <- args0
    args$aName <- 1
    expect_error(do.call(assembleSCE, args),
                 "'aName' must be of class 'character'")
    args$aName <- c("iBAQ", "log2_iBAQ")
    expect_error(do.call(assembleSCE, args),
                 "'aName' must have length 1")
    args$aName <- "missing"
    expect_error(do.call(assembleSCE, args),
                 "All values in 'aName' must be one of")

    ## testResults
    args <- args0
    args$testResults <- 1
    expect_error(do.call(assembleSCE, args),
                 "'testResults' must be of class 'data.frame'")
    args$testResults <- as.matrix(args0$testResults)
    expect_error(do.call(assembleSCE, args),
                 "'testResults' must be of class 'data.frame'")

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(assembleSCE, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("pat1", "pat2")
    expect_error(do.call(assembleSCE, args),
                 "'iColPattern' must have length 1")

    ## iColsAll
    args <- args0
    args$iColsAll <- 1
    expect_error(do.call(assembleSCE, args),
                 "'iColsAll' must be of class 'character'")

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(assembleSCE, args),
                 "'baseFileName' must be of class 'character'")

    ## nbrNA
    args <- args0
    args$nbrNA <- 1
    expect_error(do.call(assembleSCE, args),
                 "'nbrNA' must be of class 'list'")
    args <- args0
    names(args$nbrNA) <- paste0("a", names(args$nbrNA))
    expect_error(do.call(assembleSCE, args),
                 "All values in 'namesnbrNA' must be one of")

    ## featureCollections
    args <- args0
    args$featureCollections <- 1
    expect_error(do.call(assembleSCE, args),
                 "'featureCollections' must be of class 'list'")

    ## expType
    args <- args0
    args$expType <- 1
    expect_error(do.call(assembleSCE, args),
                 "'expType' must be of class 'character'")
    args$expType <- c("MaxQuant", "ProteomeDiscoverer")
    expect_error(do.call(assembleSCE, args),
                 "'expType' must have length 1")
    args$expType <- "missing"
    expect_error(do.call(assembleSCE, args),
                 "All values in 'expType' must be one of")


    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    sce <- do.call(assembleSCE, args0)
    expect_s4_class(sce, "SingleCellExperiment")
    expect_equal(nrow(sce), 70)
    expect_equal(ncol(sce), 9)
    expect_true(all(c("iBAQ", "log2_iBAQ", "log2_iBAQ_withNA", "imputed_iBAQ",
                      "MS.MS.Count", "LFQ.intensity", "Intensity",
                      "Sequence.coverage", "Unique.peptides", "Peptides") %in%
                        SummarizedExperiment::assayNames(sce)))
    expect_true(all(c("sample", "group", "nNA", "pNA") %in%
                        colnames(SummarizedExperiment::colData(sce))))
    expect_true(all(c("Gene.names", "Majority.protein.IDs", "pid",
                      "adj.P.Val", "nNA", "pNA") %in%
                        colnames(SummarizedExperiment::rowData(sce))))
})
