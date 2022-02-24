test_that("volcano plots work", {
    ## Preparation
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
                 "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
                 "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
    ecol <- paste0("iBAQ.", samples)
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "iBAQ",
                                    sep = "\t", nrows = 70)
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    qft <- addSampleAnnots(qft, iColPattern = "^iBAQ\\.",
                           sampleAnnot = sampleAnnot, mergeGroups = list())
    qft <- fixFeatureIds(qft)
    qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ")
    qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ_withNA")
    qft <- QFeatures::zeroIsNA(qft, "iBAQ")
    qft <- QFeatures::infIsNA(qft, "log2_iBAQ")
    qft <- QFeatures::infIsNA(qft, "log2_iBAQ_withNA")
    nbr_na <- QFeatures::nNA(qft, i = seq_along(qft))
    set.seed(123)
    qft <- QFeatures::impute(qft, method = "MinProb", i = "log2_iBAQ")
    tmp <- qft[["log2_iBAQ"]]
    SummarizedExperiment::assay(tmp) <- !is.finite(SummarizedExperiment::assay(tmp))
    qft <- QFeatures::addAssay(qft, tmp, name = "imputed_iBAQ")
    fcoll <- prepareFeatureCollections(
        qft = qft, idCol = "Gene.names",
        includeFeatureCollections = "complexes",
        complexDbPath = system.file("extdata", "complexes",
                                    "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                    package = "einprot"),
        speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
        customComplexes = list(), minSizeToKeep = 2)
    out_limma <- runTest(
        qft = qft, comparison = c("Adnp", "RBC_ctrl"), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addiBAQvalues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )
    out_ttest <- runTest(
        qft = qft, comparison = c("Adnp", "RBC_ctrl"), testType = "ttest",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addiBAQvalues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )
    string_db <- STRINGdb::STRINGdb$new(
        version = "11.5", species = getSpeciesInfo("mouse")$taxId,
        score_threshold = 400, input_directory = "")

    ## .curvefun
    expect_equal(.curvefun(x = 2, ta = 1, s0 = 0.2, df = 4),
                 0.4830618, tolerance = 1e-5)

    ## .makeBaseVolcano
    expect_warning(
        expect_warning(
            expect_s3_class(.makeBaseVolcano(
                res = out_limma$res, testType = "limma",
                xv = "log2FC", yv = "mlog10p",
                plotnote = out_limma$plotnote,
                plottitle = out_limma$plottitle,
                plotsubtitle = out_limma$plotsubtitle,
                curveparam = out_limma$curveparam), "ggplot"),
        "no non-missing arguments"),
    "no non-missing arguments")

    expect_warning(
        expect_warning(
            expect_s3_class(.makeBaseVolcano(
                res = out_ttest$res, testType = "ttest",
                xv = "log2FC", yv = "mlog10p",
                plotnote = out_ttest$plotnote,
                plottitle = out_ttest$plottitle,
                plotsubtitle = out_ttest$plotsubtitle,
                curveparam = out_ttest$curveparam), "ggplot"),
            "no non-missing arguments"),
        "no non-missing arguments")

    ## .complexBarPlot
    expect_s3_class(.complexBarPlot(
        res = out_limma$res,
        prs = fcoll$complexes[[1]],
        qft = qft,
        cplx = names(fcoll$complexes)[1]
    ), "ggplot")

    expect_s3_class(.complexBarPlot(
        res = out_ttest$res,
        prs = fcoll$complexes[[1]],
        qft = qft,
        cplx = names(fcoll$complexes)[1]
    ), "ggplot")

    ## plotVolcano
    ## Fails with wrong arguments
    args0 <- list(
        qft = qft, res = out_ttest$res, testType = "ttest",
        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
        plotnote = out_ttest$plotnote,
        plottitle = out_ttest$plottitle,
        plotsubtitle = out_ttest$plotsubtitle,
        volcanoFeaturesToLabel = c("Chd3"),
        volcanoMaxFeatures = 10,
        baseFileName = NULL,
        comparisonString = "RBC_ctrl_vs_Adnp",
        stringDb = string_db,
        featureCollections = out_ttest$featureCollections,
        complexFDRThr = 0.1, curveparam = out_ttest$curveparam
    )

    ## qft
    args <- args0
    args$qft <- 1
    expect_error(do.call(plotVolcano, args),
                 "'qft' must be of class 'QFeatures'")

    ## res
    args <- args0
    args$res <- 1
    expect_error(do.call(plotVolcano, args),
                 "'res' must be of class 'data.frame'")
    args$res <- as.matrix(args0$res)
    expect_error(do.call(plotVolcano, args),
                 "'res' must be of class 'data.frame'")

    ## testType
    args <- args0
    args$testType <- 1
    expect_error(do.call(plotVolcano, args),
                 "'testType' must be of class 'character'")
    args$testType <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in 'testType' must be one of")
    args$testType <- c("ttest", "limma")
    expect_error(do.call(plotVolcano, args),
                 "'testType' must have length 1")

    ## xv
    args <- args0
    args$xv <- 1
    expect_error(do.call(plotVolcano, args),
                 "'xv' must be of class 'character'")
    args$xv <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'xv' must have length 1")
    args$xv <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in 'xv' must be one of")

    ## yv
    args <- args0
    args$yv <- 1
    expect_error(do.call(plotVolcano, args),
                 "'yv' must be of class 'character'")
    args$yv <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'yv' must have length 1")
    args$yv <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in 'yv' must be one of")

    ## volcind
    args <- args0
    args$volcind <- 1
    expect_error(do.call(plotVolcano, args),
                 "'volcind' must be of class 'character'")
    args$volcind <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'volcind' must have length 1")
    args$volcind <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in 'volcind' must be one of")

    ## plotnote
    args <- args0
    args$plotnote <- 1
    expect_error(do.call(plotVolcano, args),
                 "'plotnote' must be of class 'character'")
    args$plotnote <- c("note1", "note2")
    expect_error(do.call(plotVolcano, args),
                 "'plotnote' must have length 1")

    ## plottitle
    args <- args0
    args$plottitle <- 1
    expect_error(do.call(plotVolcano, args),
                 "'plottitle' must be of class 'character'")
    args$plottitle <- c("note1", "note2")
    expect_error(do.call(plotVolcano, args),
                 "'plottitle' must have length 1")

    ## plotsubtitle
    args <- args0
    args$plotsubtitle <- 1
    expect_error(do.call(plotVolcano, args),
                 "'plotsubtitle' must be of class 'character'")
    args$plotsubtitle <- c("note1", "note2")
    expect_error(do.call(plotVolcano, args),
                 "'plotsubtitle' must have length 1")

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(plotVolcano, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "2"
    expect_error(do.call(plotVolcano, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(1, 2)
    expect_error(do.call(plotVolcano, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(plotVolcano, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(plotVolcano, args),
                 "'baseFileName' must be of class 'character'")
    args$baseFileName <- c("pat1", "pat2")
    expect_error(do.call(plotVolcano, args),
                 "'baseFileName' must have length 1")

    ## comparisonString
    args <- args0
    args$comparisonString <- 1
    expect_error(do.call(plotVolcano, args),
                 "'comparisonString' must be of class 'character'")
    args$comparisonString <- c("note1", "note2")
    expect_error(do.call(plotVolcano, args),
                 "'comparisonString' must have length 1")

    ## stringDb
    args <- args0
    args$stringDb <- 1
    expect_error(do.call(plotVolcano, args),
                 "'stringDb' must be of class 'STRINGdb'")

    ## featureCollections
    args <- args0
    args$featureCollections <- "2"
    expect_error(do.call(plotVolcano, args),
                 "'featureCollections' must be of class 'list'")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "2"
    expect_error(do.call(plotVolcano, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(plotVolcano, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(plotVolcano, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$complexFDRThr <- 2
    expect_error(do.call(plotVolcano, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## curveparam
    args <- args0
    args$curveparam <- 1
    expect_error(do.call(plotVolcano, args),
                 "'curveparam' must be of class 'list'")


    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    outl <- plotVolcano(qft = qft, res = out_limma$res, testType = "limma",
                        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
                        plotnote = out_limma$plotnote,
                        plottitle = out_limma$plottitle,
                        plotsubtitle = out_limma$plotsubtitle,
                        volcanoFeaturesToLabel = c("Chd3"),
                        volcanoMaxFeatures = 10,
                        baseFileName = NULL,
                        comparisonString = "RBC_ctrl_vs_Adnp",
                        stringDb = string_db,
                        featureCollections = out_limma$featureCollections,
                        complexFDRThr = 0.1, curveparam = out_limma$curveparam)
    expect_type(outl, "list")
    expect_length(outl, 2)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")

    outl <- plotVolcano(qft = qft, res = out_ttest$res, testType = "ttest",
                        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
                        plotnote = out_ttest$plotnote,
                        plottitle = out_ttest$plottitle,
                        plotsubtitle = out_ttest$plotsubtitle,
                        volcanoFeaturesToLabel = c("Chd3"),
                        volcanoMaxFeatures = 10,
                        baseFileName = NULL,
                        comparisonString = "RBC_ctrl_vs_Adnp",
                        stringDb = string_db,
                        featureCollections = out_ttest$featureCollections,
                        complexFDRThr = 0.1, curveparam = out_ttest$curveparam)
    expect_type(outl, "list")
    expect_length(outl, 2)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")

    ## Save to file, no STRINGdb object
    bfn <- tempfile()
    outl <- plotVolcano(qft = qft, res = out_ttest$res, testType = "ttest",
                        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
                        plotnote = out_ttest$plotnote,
                        plottitle = out_ttest$plottitle,
                        plotsubtitle = out_ttest$plotsubtitle,
                        volcanoFeaturesToLabel = c("Chd3"),
                        volcanoMaxFeatures = 10,
                        baseFileName = bfn,
                        comparisonString = "RBC_ctrl_vs_Adnp",
                        stringDb = NULL,
                        featureCollections = out_ttest$featureCollections,
                        complexFDRThr = 0.1, curveparam = out_ttest$curveparam)
    expect_type(outl, "list")
    expect_length(outl, 2)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp.pdf")))

    ## Save to file, download string db object
    skip_if_offline()
    bfn <- tempfile()
    outl <- plotVolcano(qft = qft, res = out_ttest$res, testType = "ttest",
                        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
                        plotnote = out_ttest$plotnote,
                        plottitle = out_ttest$plottitle,
                        plotsubtitle = out_ttest$plotsubtitle,
                        volcanoFeaturesToLabel = c("Chd3"),
                        volcanoMaxFeatures = 10,
                        baseFileName = bfn,
                        comparisonString = "RBC_ctrl_vs_Adnp",
                        stringDb = string_db,
                        featureCollections = out_ttest$featureCollections,
                        complexFDRThr = 0.1, curveparam = out_ttest$curveparam)
    expect_type(outl, "list")
    expect_length(outl, 2)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp.pdf")))
})
