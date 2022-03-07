test_that("volcano plots work", {

    out_limma <- runTest(
        sce = sce_mq_final, comparison = c("Adnp", "RBC_ctrl"), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
    )
    out_ttest <- runTest(
        sce = sce_mq_final, comparison = c("Adnp", "RBC_ctrl"), testType = "ttest",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, iColPattern = "^iBAQ\\.", aName = "iBAQ"
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
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ"
    ), "ggplot")

    expect_s3_class(.complexBarPlot(
        res = out_ttest$res,
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ"
    ), "ggplot")

    ## plotVolcano
    ## Fails with wrong arguments
    args0 <- list(
        sce = sce_mq_final, res = out_ttest$res, testType = "ttest",
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
        complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
        curveparam = out_ttest$curveparam,
        abundanceColPat = "iBAQ"
    )

    ## sce
    args <- args0
    args$sce <- 1
    expect_error(do.call(plotVolcano, args),
                 "'sce' must be of class 'SummarizedExperiment'")

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

    ## xvma
    args <- args0
    args$xvma <- 1
    expect_error(do.call(plotVolcano, args),
                 "'xvma' must be of class 'character'")
    args$xvma <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'xvma' must have length 1")
    args$xvma <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in 'xvma' must be one of")

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

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "2"
    expect_error(do.call(plotVolcano, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(plotVolcano, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(plotVolcano, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## curveparam
    args <- args0
    args$curveparam <- 1
    expect_error(do.call(plotVolcano, args),
                 "'curveparam' must be of class 'list'")

    ## abundanceColPat
    args <- args0
    args$abundanceColPat <- 1
    expect_error(do.call(plotVolcano, args),
                 "'abundanceColPat' must be of class 'character'")
    args$abundanceColPat <- c("note1", "note2")
    expect_error(do.call(plotVolcano, args),
                 "'abundanceColPat' must have length 1")

    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    expect_warning(
        outl <- plotVolcano(sce = sce_mq_final, res = out_limma$res, testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma$plotnote,
                            plottitle = out_limma$plottitle,
                            plotsubtitle = out_limma$plotsubtitle,
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_limma$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma$curveparam, abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 3)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggma, "ggplot")

    expect_warning(
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$res, testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnote,
                            plottitle = out_ttest$plottitle,
                            plotsubtitle = out_ttest$plotsubtitle,
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparam, abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 3)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)

    ## Save to file, no STRINGdb object
    bfn <- tempfile()
    wns <- capture_warnings({
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$res, testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnote,
                            plottitle = out_ttest$plottitle,
                            plotsubtitle = out_ttest$plotsubtitle,
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = NULL,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparam, abundanceColPat = "iBAQ")})
    for (wn in wns) {
        expect_match(wn, "rows containing missing values")
    }
    expect_type(outl, "list")
    expect_length(outl, 3)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp.pdf")))

    ## Save to file, download string db object
    skip_if_offline()
    bfn <- tempfile()
    wns <- capture_warnings({
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$res, testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnote,
                            plottitle = out_ttest$plottitle,
                            plotsubtitle = out_ttest$plotsubtitle,
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparam, abundanceColPat = "iBAQ")})
    for (wn in wns) {
        expect_match(wn, "rows containing missing values")
    }
    expect_type(outl, "list")
    expect_length(outl, 3)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, "_volcano_RBC_ctrl_vs_Adnp.pdf")))
})
