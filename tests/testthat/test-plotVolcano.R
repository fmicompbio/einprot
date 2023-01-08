test_that("volcano plots work", {

    out_limma <- runTest(
        sce = sce_mq_final, comparisons = list(c("Adnp", "RBC_ctrl")), testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "iBAQ", singleFit = FALSE
    )
    out_limma_merged <- runTest(
        sce = sce_mq_final, comparisons = list(c("adnp_rbc_complement", "adnp_rbc")),
        groupComposition = list(adnp_rbc = c("Adnp", "RBC_ctrl"),
                                adnp_rbc_complement = "Chd4BF"),
        testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "iBAQ", singleFit = FALSE
    )
    out_ttest <- runTest(
        sce = sce_mq_final, comparisons = list(c("Adnp", "RBC_ctrl")), testType = "ttest",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "iBAQ", singleFit = FALSE
    )
    out_proda <- runTest(
        sce = sce_mq_final, comparisons = list(c("Adnp", "RBC_ctrl")), testType = "proDA",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "iBAQ", singleFit = FALSE
    )
    string_db <- STRINGdb::STRINGdb$new(
        version = "11.5", species = getSpeciesInfo("mouse")$taxId,
        score_threshold = 400, input_directory = "")

    ## ---------------------------------------------------------------------- ##
    ## .curvefun
    ## ---------------------------------------------------------------------- ##
    expect_equal(.curvefun(x = 2, ta = 1, s0 = 0.2, df = 4),
                 0.4830618, tolerance = 1e-5)

    ## ---------------------------------------------------------------------- ##
    ## .makeWaterfallPlot
    ## ---------------------------------------------------------------------- ##
    expect_error(.makeWaterfallPlot(res = 1, ntop = 10, xv = "logFC",
                                    volcind = "showInVolcano", title = ""),
                 "'res' must be of class 'data.frame'")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = "1",
                                    xv = "logFC",
                                    volcind = "showInVolcano", title = ""),
                 "'ntop' must be of class 'numeric'")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = c(1, 10),
                                    xv = "logFC",
                                    volcind = "showInVolcano", title = ""),
                 "'ntop' must have length 1")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = 1,
                                    volcind = "showInVolcano", title = ""),
                 "'xv' must be of class 'character'")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = c("logFC", "P.Value"),
                                    volcind = "showInVolcano", title = ""),
                 "'xv' must have length 1")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "missing",
                                    volcind = "showInVolcano", title = ""),
                 "All values in 'xv' must be one of")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "logFC",
                                    volcind = 1, title = ""),
                 "'volcind' must be of class 'character'")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "logFC",
                                    volcind = c("showInVolcano", "P.Value"),
                                    title = ""),
                 "'volcind' must have length 1")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "logFC",
                                    volcind = "missing", title = ""),
                 "All values in 'volcind' must be one of")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "logFC",
                                    volcind = "showInVolcano", title = 1),
                 "'title' must be of class 'character'")
    expect_error(.makeWaterfallPlot(res = out_limma$tests[[1]], ntop = 10,
                                    xv = "logFC",
                                    volcind = "showInVolcano",
                                    title = c("a", "b")),
                 "'title' must have length 1")

    out <- expect_s3_class(.makeWaterfallPlot(
        res = out_limma$tests[[1]], ntop = 10, xv = "logFC",
        volcind = "showInVolcano", title = ""), "ggplot")
    expect_s3_class(out$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(out$data)))
    expect_equal(rownames(out$data)[which.min(out$data$P.Value)], "Adnp")

    out <- expect_s3_class(.makeWaterfallPlot(
        res = out_ttest$tests[[1]], ntop = 10, xv = "logFC",
        volcind = "showInVolcano", title = ""), "ggplot")
    expect_s3_class(out$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(out$data)))
    expect_equal(rownames(out$data)[which.min(out$data$P.Value)], "Adnp")

    ## ---------------------------------------------------------------------- ##
    ## .makeBaseVolcano
    ## ---------------------------------------------------------------------- ##
    out <- expect_s3_class(.makeBaseVolcano(
        res = out_limma$tests[[1]], testType = "limma",
        xv = "logFC", yv = "mlog10p",
        plotnote = out_limma$plotnotes[[1]],
        plottitle = out_limma$plottitles[[1]],
        plotsubtitle = out_limma$plotsubtitles[[1]],
        curveparam = out_limma$curveparams[[1]]), "ggplot")
    expect_s3_class(out$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(out$data)))
    expect_equal(rownames(out$data)[which.min(out$data$P.Value)], "Adnp")
    expect_lt(out$data["Adnp", "logFC"], 0)
    expect_true(all(out$data$mlog10p[!is.na(out$data$logFC)] >= 0))
    expect_equal(out$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(out$data["Adnp", "showInVolcano"])
    expect_true(all(out$data$showInVolcano[which(abs(out$data$logFC) >= 1 &
                                                     out$data$adj.P.Val <= 0.05)]))

    out <- expect_s3_class(.makeBaseVolcano(
        res = out_ttest$tests[[1]], testType = "ttest",
        xv = "logFC", yv = "mlog10p",
        plotnote = out_ttest$plotnotes[[1]],
        plottitle = out_ttest$plottitles[[1]],
        plotsubtitle = out_ttest$plotsubtitles[[1]],
        curveparam = out_ttest$curveparams[[1]]), "ggplot")
    expect_s3_class(out$data, "data.frame")
    expect_true(all(c("pid", "logFC", "sam", "AveExpr", "mlog10p") %in%
                        colnames(out$data)))
    expect_equal(rownames(out$data)[which.min(out$data$P.Value)], "Adnp")
    expect_lt(out$data["Adnp", "logFC"], 0)
    expect_true(all(out$data$mlog10p[!is.na(out$data$logFC)] >= 0))
    expect_equal(out$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(out$data["Adnp", "showInVolcano"])
    expect_true(min(abs(out$data$sam[which(out$data$showInVolcano)])) >
                    max(abs(out$data$sam[which(!out$data$showInVolcano)])))

    out <- expect_s3_class(.makeBaseVolcano(
        res = out_proda$tests[[1]], testType = "proDA",
        xv = "logFC", yv = "mlog10p",
        plotnote = out_proda$plotnotes[[1]],
        plottitle = out_proda$plottitles[[1]],
        plotsubtitle = out_proda$plotsubtitles[[1]],
        curveparam = out_proda$curveparams[[1]]), "ggplot")
    expect_s3_class(out$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "avg_abundance", "mlog10p") %in%
                        colnames(out$data)))
    expect_equal(rownames(out$data)[which.min(out$data$P.Value)], "Chd4")
    expect_lt(out$data["Adnp", "logFC"], 0)
    expect_true(all(out$data$mlog10p[!is.na(out$data$logFC)] >= 0))
    expect_equal(out$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(out$data["Adnp", "showInVolcano"])
    expect_true(all(out$data$showInVolcano[which(abs(out$data$logFC) >= 1 &
                                                     out$data$adj.P.Val <= 0.05)]))

    ## ---------------------------------------------------------------------- ##
    ## .complexBarPlot
    ## ---------------------------------------------------------------------- ##
    ## limma
    out <- .complexBarPlot(
        res = out_limma$tests[[1]],
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ",
        groupmap = NULL
    )
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("pid", "mergegroup", "mean_abundance", "sd_abundance"))
    expect_equal(out$data$mean_abundance[out$data$pid == "Arnt" &
                                             out$data$mergegroup == "Adnp"],
                 mean(SummarizedExperiment::assay(
                     sce_mq_final, "iBAQ")["Arnt",
                                           sce_mq_final$group == "Adnp"],
                     na.rm = TRUE))
    expect_s3_class(out$layers[[3]]$data, "data.frame")
    expect_equal(nrow(out$layers[[3]]$data), 6)
    expect_named(out$layers[[3]]$data, c("pid", "sample", "Abundance",
                                         "group", "mergegroup"))
    expect_equal(out$layers[[3]]$data$Abundance[
        out$layers[[3]]$data$pid == "Arnt" &
            out$layers[[3]]$data$sample == "Adnp_IP05"],
                 SummarizedExperiment::assay(
                     sce_mq_final, "iBAQ")["Arnt", "Adnp_IP05"])

    ## t-test
    out <- .complexBarPlot(
        res = out_ttest$tests[[1]],
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ",
        groupmap = NULL
    )
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("pid", "mergegroup", "mean_abundance", "sd_abundance"))
    expect_equal(out$data$mean_abundance[out$data$pid == "Arnt" &
                                             out$data$mergegroup == "Adnp"],
                 mean(SummarizedExperiment::assay(
                     sce_mq_final, "iBAQ")["Arnt",
                                           sce_mq_final$group == "Adnp"],
                     na.rm = TRUE))
    expect_s3_class(out$layers[[3]]$data, "data.frame")
    expect_equal(nrow(out$layers[[3]]$data), 6)
    expect_named(out$layers[[3]]$data, c("pid", "sample", "Abundance",
                                         "group", "mergegroup"))
    expect_equal(out$layers[[3]]$data$Abundance[
        out$layers[[3]]$data$pid == "Arnt" &
            out$layers[[3]]$data$sample == "Adnp_IP05"],
        SummarizedExperiment::assay(
            sce_mq_final, "iBAQ")["Arnt", "Adnp_IP05"])

    ## Merge groups without accounting for the grouping
    out <- .complexBarPlot(
        res = out_limma_merged$tests[[1]],
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ",
        groupmap = NULL
    )
    expect_error(print(out), "Insufficient values in manual scale")

    ## Merge groups after accounting for the grouping
    out <- .complexBarPlot(
        res = out_limma_merged$tests[[1]],
        prs = fcoll_mq_final$complexes[[1]],
        sce = sce_mq_final,
        cplx = names(fcoll_mq_final$complexes)[1],
        colpat = "iBAQ",
        groupmap = data.frame(group = c("Adnp", "Chd4BF", "RBC_ctrl"),
                              mergegroup = c("adnp_rbc", "adnp_rbc_complement", "adnp_rbc"))
    )
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 4)
    expect_equal(nrow(out$data), 2)
    expect_named(out$data, c("pid", "mergegroup", "mean_abundance", "sd_abundance"))
    expect_equal(out$data$mean_abundance[out$data$pid == "Arnt" &
                                             out$data$mergegroup == "adnp_rbc"],
                 mean(SummarizedExperiment::assay(
                     sce_mq_final, "iBAQ")["Arnt",
                                           sce_mq_final$group %in% c("Adnp", "RBC_ctrl")],
                     na.rm = TRUE))
    expect_s3_class(out$layers[[3]]$data, "data.frame")
    expect_equal(nrow(out$layers[[3]]$data), 9)
    expect_named(out$layers[[3]]$data, c("pid", "sample", "Abundance",
                                         "group", "mergegroup"))
    expect_equal(out$layers[[3]]$data$Abundance[
        out$layers[[3]]$data$pid == "Arnt" &
            out$layers[[3]]$data$sample == "Adnp_IP05"],
        SummarizedExperiment::assay(
            sce_mq_final, "iBAQ")["Arnt", "Adnp_IP05"])


    ## ---------------------------------------------------------------------- ##
    ## plotVolcano
    ## ---------------------------------------------------------------------- ##
    ## Fails with wrong arguments
    args0 <- list(
        sce = sce_mq_final, res = out_ttest$tests[[1]], testType = "ttest",
        xv = "logFC", yv = "mlog10p", volcind = "showInVolcano",
        plotnote = out_ttest$plotnotes[[1]],
        plottitle = out_ttest$plottitles[[1]],
        plotsubtitle = out_ttest$plotsubtitles[[1]],
        volcanoFeaturesToLabel = c("Chd3"),
        volcanoMaxFeatures = 10,
        baseFileName = NULL,
        comparisonString = "RBC_ctrl_vs_Adnp",
        stringDb = string_db,
        featureCollections = out_ttest$featureCollections,
        complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
        curveparam = out_ttest$curveparams[[1]],
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
                 "'$colsxv' must be of class 'character'", fixed = TRUE)
    args$xv <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'$colsxv' must have length 1", fixed = TRUE)
    args$xv <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in '$colsxv' must be one of", fixed = TRUE)

    ## yv
    args <- args0
    args$yv <- 1
    expect_error(do.call(plotVolcano, args),
                 "'$colsyv' must be of class 'character'", fixed = TRUE)
    args$yv <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'$colsyv' must have length 1", fixed = TRUE)
    args$yv <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in '$colsyv' must be one of", fixed = TRUE)

    ## xvma
    args <- args0
    args$xvma <- 1
    expect_error(do.call(plotVolcano, args),
                 "'$colsxvma' must be of class 'character'", fixed = TRUE)
    args$xvma <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'$colsxvma' must have length 1", fixed = TRUE)
    args$xvma <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in '$colsxvma' must be one of", fixed = TRUE)

    ## volcind
    args <- args0
    args$volcind <- 1
    expect_error(do.call(plotVolcano, args),
                 "'$colsvolcind' must be of class 'character'", fixed = TRUE)
    args$volcind <- c("logFC", "P.Value")
    expect_error(do.call(plotVolcano, args),
                 "'$colsvolcind' must have length 1", fixed = TRUE)
    args$volcind <- "missing"
    expect_error(do.call(plotVolcano, args),
                 "All values in '$colsvolcind' must be one of", fixed = TRUE)

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
    ## limma
    expect_warning(
        outl <- plotVolcano(sce = sce_mq_final, res = out_limma$tests[[1]],
                            testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma$plotnotes[[1]],
                            plottitle = out_limma$plottitles[[1]],
                            plotsubtitle = out_limma$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_limma$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma$curveparams[[1]],
                            abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggma, "ggplot")
    expect_s3_class(outl$ggwf, "ggplot")
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "Adnp")
    expect_lt(outl$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(outl$gg$data["Adnp", "showInVolcano"])
    expect_true(all(outl$gg$data$showInVolcano[which(abs(outl$gg$data$logFC) >= 1 &
                                                         outl$gg$data$adj.P.Val <= 0.05)]))
    expect_equal(outl$gg$data, outl$ggma$data)

    ## limma with einprotLabel column missing (should use pid column)
    ## also all values in volcind column FALSE -> no ggwf plot
    restemp <- out_limma$tests[[1]]
    restemp$einprotLabel <- NULL
    restemp$showInVolcano <- FALSE
    expect_warning(
        outl <- plotVolcano(sce = sce_mq_final, res = restemp,
                            testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma$plotnotes[[1]],
                            plottitle = out_limma$plottitles[[1]],
                            plotsubtitle = out_limma$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_limma$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma$curveparams[[1]],
                            abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggma, "ggplot")
    expect_null(outl$ggwf)
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "Adnp")
    expect_lt(outl$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_equal(outl$gg$data, outl$ggma$data)

    ## limma with merging, write to file
    bfn <- tempfile()
    wns <- capture_warnings(
        outl <- plotVolcano(sce = sce_mq_final, res = out_limma_merged$tests[[1]],
                            testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma_merged$plotnotes[[1]],
                            plottitle = out_limma_merged$plottitles[[1]],
                            plotsubtitle = out_limma_merged$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = names(out_limma_merged$tests)[1],
                            groupComposition = list(adnp_rbc = c("Adnp", "RBC_ctrl"),
                                                    adnp_rbc_complement = "Chd4BF"),
                            stringDb = NULL,
                            featureCollections = out_limma_merged$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma_merged$curveparams[[1]],
                            abundanceColPat = "iBAQ"))
    expect_true(length(wns) > 0)
    expect_match(wns[1], ".*rows containing missing values.*")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggma, "ggplot")
    expect_s3_class(outl$ggwf, "ggplot")
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "Mbd3")
    expect_lt(outl$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(outl$gg$data["Mbd3", "showInVolcano"])
    expect_true(all(outl$gg$data$showInVolcano[which(abs(outl$gg$data$logFC) >= 1 &
                                                         outl$gg$data$adj.P.Val <= 0.05)]))
    expect_equal(outl$gg$data, outl$ggma$data)

    ## t-test
    expect_warning(
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$tests[[1]],
                            testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnotes[[1]],
                            plottitle = out_ttest$plottitles[[1]],
                            plotsubtitle = out_ttest$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparams[[1]],
                            abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_s3_class(outl$ggwf, "ggplot")
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "Adnp")
    expect_lt(outl$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(outl$gg$data["Adnp", "showInVolcano"])
    expect_true(min(abs(outl$gg$data$sam[which(outl$gg$data$showInVolcano)])) >
                    max(abs(outl$gg$data$sam[which(!outl$gg$data$showInVolcano)])))

    ## Don't specify arguments
    expect_warning(
        outl2 <- plotVolcano(sce = sce_mq_final, res = out_ttest$tests[[1]],
                             testType = "ttest",
                             xv = NULL, yv = NULL, xvma = NULL,
                             volcind = NULL,
                             plotnote = out_ttest$plotnotes[[1]],
                             plottitle = out_ttest$plottitles[[1]],
                             plotsubtitle = out_ttest$plotsubtitles[[1]],
                             volcanoFeaturesToLabel = c("Chd3"),
                             volcanoMaxFeatures = 10,
                             baseFileName = NULL,
                             comparisonString = "RBC_ctrl_vs_Adnp",
                             stringDb = string_db,
                             featureCollections = out_ttest$featureCollections,
                             complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                             curveparam = out_ttest$curveparams[[1]],
                             abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl2, "list")
    expect_length(outl2, 4)
    expect_s3_class(outl2$gg, "ggplot")
    expect_s3_class(outl2$ggint, "girafe")
    expect_null(outl2$ggma)
    expect_s3_class(outl2$ggwf, "ggplot")
    expect_s3_class(outl2$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl2$gg$data)))
    expect_equal(rownames(outl2$gg$data)[which.min(outl2$gg$data$P.Value)], "Adnp")
    expect_lt(outl2$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl2$gg$data$mlog10p[!is.na(outl2$gg$data$logFC)] >= 0))
    expect_equal(outl2$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(outl2$gg$data["Adnp", "showInVolcano"])
    expect_true(min(abs(outl2$gg$data$sam[which(outl2$gg$data$showInVolcano)])) >
                    max(abs(outl2$gg$data$sam[which(!outl2$gg$data$showInVolcano)])))

    ## Don't specify arguments, proDA
    expect_warning(
        outl2pr <- plotVolcano(sce = sce_mq_final, res = out_proda$tests[[1]],
                               testType = "proDA",
                               xv = NULL, yv = NULL, xvma = NULL,
                               volcind = NULL,
                               plotnote = out_proda$plotnotes[[1]],
                               plottitle = out_proda$plottitles[[1]],
                               plotsubtitle = out_proda$plotsubtitles[[1]],
                               volcanoFeaturesToLabel = c("Chd3"),
                               volcanoMaxFeatures = 10,
                               baseFileName = NULL,
                               comparisonString = "RBC_ctrl_vs_Adnp",
                               stringDb = string_db,
                               featureCollections = out_ttest$featureCollections,
                               complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                               curveparam = out_proda$curveparams[[1]],
                               abundanceColPat = "iBAQ"),
        "rows containing missing values")
    expect_type(outl2pr, "list")
    expect_length(outl2pr, 4)
    expect_s3_class(outl2pr$gg, "ggplot")
    expect_s3_class(outl2pr$ggint, "girafe")
    expect_null(outl2pr$ggma)
    expect_s3_class(outl2pr$ggwf, "ggplot")
    expect_s3_class(outl2pr$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "avg_abundance", "mlog10p") %in%
                        colnames(outl2pr$gg$data)))
    expect_equal(rownames(outl2pr$gg$data)[which.min(outl2pr$gg$data$P.Value)], "Chd4")
    expect_lt(outl2pr$gg$data["Adnp", "logFC"], 0)
    expect_true(all(outl2pr$gg$data$mlog10p[!is.na(outl2pr$gg$data$logFC)] >= 0))
    expect_equal(outl2pr$gg$data["Adnp", "IDsForSTRING"], "Adnp")
    expect_true(outl2pr$gg$data["Adnp", "showInVolcano"])

    ## Save to file, no STRINGdb object
    bfn <- tempfile()
    wns <- capture_warnings({
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$tests[[1]],
                            testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnotes[[1]],
                            plottitle = out_ttest$plottitles[[1]],
                            plotsubtitle = out_ttest$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = NULL,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.001, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparams[[1]],
                            abundanceColPat = "iBAQ")})
    for (wn in wns) {
        expect_match(wn, "rows containing missing values")
    }
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_s3_class(outl$ggwf, "ggplot")
    expect_false(file.exists(paste0(bfn, "_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, ".pdf")))

    ## Save to file, no STRINGdb object (limma)
    bfn <- tempfile()
    wns <- capture_warnings({
        outl <- plotVolcano(sce = sce_mq_final, res = out_limma$tests[[1]],
                            testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma$plotnotes[[1]],
                            plottitle = out_limma$plottitles[[1]],
                            plotsubtitle = out_limma$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = NULL,
                            featureCollections = out_limma$featureCollections,
                            complexFDRThr = 0.8, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma$curveparams[[1]],
                            abundanceColPat = "iBAQ")})
    for (wn in wns) {
        expect_match(wn, "rows containing missing values")
    }
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggwf, "ggplot")
    expect_true(file.exists(paste0(bfn, "_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, ".pdf")))

    ## Save to file, download string db object
    skip_if_offline()
    bfn <- tempfile()
    wns <- capture_warnings({
        outl <- plotVolcano(sce = sce_mq_final, res = out_ttest$tests[[1]],
                            testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnotes[[1]],
                            plottitle = out_ttest$plottitles[[1]],
                            plotsubtitle = out_ttest$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = c("Chd3"),
                            volcanoMaxFeatures = 10,
                            baseFileName = bfn,
                            comparisonString = "RBC_ctrl_vs_Adnp",
                            stringDb = string_db,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparams[[1]],
                            abundanceColPat = "iBAQ")})
    for (wn in wns) {
        expect_match(wn, "rows containing missing values")
    }
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_s3_class(outl$ggwf, "ggplot")
    expect_true(file.exists(paste0(bfn, "_complexes.pdf")))
    expect_true(file.exists(paste0(bfn, ".pdf")))

    ## ---------------------------------------------------------------------- ##
    ## PD data
    ## ---------------------------------------------------------------------- ##
    out_limma <- runTest(
        sce = sce_pd_final, comparisons = list(c("HIS4KO", "WT")), testType = "limma",
        assayForTests = "log2_Abundance", assayImputation = "imputed_Abundance",
        minNbrValidValues = 2, minlFC = 0.5, featureCollections = fcoll_pd_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "Abundance", singleFit = FALSE
    )
    out_ttest <- runTest(
        sce = sce_pd_final, comparisons = list(c("HIS4KO", "WT")), testType = "ttest",
        assayForTests = "log2_Abundance", assayImputation = "imputed_Abundance",
        minNbrValidValues = 2, minlFC = 0.5, featureCollections = fcoll_pd_final,
        complexFDRThr = 0.1, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "Abundance", singleFit = FALSE
    )

    expect_warning(
        outl <- plotVolcano(sce = sce_pd_final, res = out_limma$tests[[1]],
                            testType = "limma",
                            xv = "logFC", yv = "mlog10p", xvma = "AveExpr",
                            volcind = "showInVolcano",
                            plotnote = out_limma$plotnotes[[1]],
                            plottitle = out_limma$plottitles[[1]],
                            plotsubtitle = out_limma$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = "",
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "WT_vs_HIS4KO",
                            stringDb = NULL,
                            featureCollections = out_limma$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_limma$curveparams[[1]],
                            abundanceColPat = "Abundance"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_s3_class(outl$ggma, "ggplot")
    expect_s3_class(outl$ggwf, "ggplot")
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "ACH1")
    expect_lt(outl$gg$data["ACH1", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["ACH1", "IDsForSTRING"], "ACH1")
    expect_true(outl$gg$data["ACH1", "showInVolcano"])
    expect_true(all(outl$gg$data$showInVolcano[which(abs(outl$gg$data$logFC) >= 1 &
                                                         outl$gg$data$adj.P.Val <= 0.05)]))
    expect_equal(outl$gg$data, outl$ggma$data)
    ## Text labels
    tmp <- outl$gg$layers[[4]]$data
    expect_s3_class(tmp, "data.frame")
    expect_true(all(tmp$showInVolcano))

    expect_warning(
        outl <- plotVolcano(sce = sce_pd_final, res = out_ttest$tests[[1]],
                            testType = "ttest",
                            xv = "logFC", yv = "mlog10p", xvma = NULL,
                            volcind = "showInVolcano",
                            plotnote = out_ttest$plotnotes[[1]],
                            plottitle = out_ttest$plottitles[[1]],
                            plotsubtitle = out_ttest$plotsubtitles[[1]],
                            volcanoFeaturesToLabel = "TCP1",
                            volcanoMaxFeatures = 10,
                            baseFileName = NULL,
                            comparisonString = "WT_vs_HIS4KO",
                            stringDb = NULL,
                            featureCollections = out_ttest$featureCollections,
                            complexFDRThr = 0.1, maxNbrComplexesToPlot = 10,
                            curveparam = out_ttest$curveparams[[1]],
                            abundanceColPat = "Abundance"),
        "rows containing missing values")
    expect_type(outl, "list")
    expect_length(outl, 4)
    expect_s3_class(outl$gg, "ggplot")
    expect_s3_class(outl$ggint, "girafe")
    expect_null(outl$ggma)
    expect_s3_class(outl$ggwf, "ggplot")
    expect_s3_class(outl$gg$data, "data.frame")
    expect_true(all(c("pid", "logFC", "t", "AveExpr", "mlog10p") %in%
                        colnames(outl$gg$data)))
    expect_equal(rownames(outl$gg$data)[which.min(outl$gg$data$P.Value)], "ACH1")
    expect_lt(outl$gg$data["ACH1", "logFC"], 0)
    expect_true(all(outl$gg$data$mlog10p[!is.na(outl$gg$data$logFC)] >= 0))
    expect_equal(outl$gg$data["ACH1", "IDsForSTRING"], "ACH1")
    expect_true(outl$gg$data["ACH1", "showInVolcano"])
    expect_true(min(abs(outl$gg$data$sam[which(outl$gg$data$showInVolcano)])) >
                    max(abs(outl$gg$data$sam[which(!outl$gg$data$showInVolcano)])))
    ## Text labels
    tmp <- outl$gg$layers[[6]]$data
    expect_s3_class(tmp, "data.frame")
    expect_true("TCP1" %in% rownames(tmp))
    expect_true(all(tmp$showInVolcano[tmp$pid != "TCP1"]))
    expect_false(tmp$showInVolcano[tmp$pid == "TCP1"])

})
