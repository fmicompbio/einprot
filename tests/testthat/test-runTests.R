test_that("testing works", {

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    args0 <- list(
        sce = sce_mq_final,
        comparisons = list(c("Adnp", "RBC_ctrl")),
        testType = "limma",
        assayForTests = "log2_iBAQ",
        assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2,
        minlFC = 0,
        featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        baseFileName = NULL,
        seed = 123,
        nperm = 25,
        volcanoS0 = 0.1,
        addAbundanceValues = TRUE,
        aName = "iBAQ",
        singleFit = FALSE,
        subtractBaseline = FALSE,
        baselineGroup = ""
    )

    ## sce
    args <- args0
    args$sce <- 1
    expect_error(do.call(runTest, args),
                 "'sce' must be of class 'SummarizedExperiment'")
    args <- args0
    args$sce$group <- as.numeric(as.factor(args$sce$group))
    expect_error(do.call(runTest, args),
                 "'$scegroup' must be of class 'character'", fixed = TRUE)

    ## comparison
    args <- args0
    args$comparisons <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(c("Adnp"))
    expect_error(do.call(runTest, args),
                 "'comparison' must have length 2")
    args$comparisons <- list(c("missing", "Adnp"))
    expect_error(do.call(runTest, args),
                 "All values in 'comparison' must be one of")

    ## testType
    args <- args0
    args$testType <- 1
    expect_error(do.call(runTest, args),
                 "'testType' must be of class 'character'")
    args$testType <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'testType' must be one of")
    args$testType <- c("ttest", "limma")
    expect_error(do.call(runTest, args),
                 "'testType' must have length 1")

    ## assayForTests
    args <- args0
    args$assayForTests <- 1
    expect_error(do.call(runTest, args),
                 "'assayForTests' must be of class 'character'")
    args$assayForTests <- c("log2_iBAQ", "imputed_iBAQ")
    expect_error(do.call(runTest, args),
                 "'assayForTests' must have length 1")
    args$assayForTests <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'assayForTests' must be one of")

    ## assayImputation
    args <- args0
    args$assayImputation <- 1
    expect_error(do.call(runTest, args),
                 "'assayImputation' must be of class 'character'")
    args$assayImputation <- c("log2_iBAQ", "imputed_iBAQ")
    expect_error(do.call(runTest, args),
                 "'assayImputation' must have length 1")
    args$assayImputation <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'assayImputation' must be one of")

    ## minnbrValidValues
    args <- args0
    args$minNbrValidValues <- "2"
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "2"
    expect_error(do.call(runTest, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runTest, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## featureCollections
    args <- args0
    args$featureCollections <- "2"
    expect_error(do.call(runTest, args),
                 "'featureCollections' must be of class 'list'")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "2"
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$complexFDRThr <- 2
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$volcanoAdjPvalThr <- 2
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(runTest, args),
                 "'baseFileName' must be of class 'character'")
    args$baseFileName <- c("pat1", "pat2")
    expect_error(do.call(runTest, args),
                 "'baseFileName' must have length 1")

    ## seed
    args <- args0
    args$testType <- "ttest"
    args$seed <- "2"
    expect_error(do.call(runTest, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runTest, args),
                 "'seed' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## nperm
    args <- args0
    args$testType <- "ttest"
    args$nperm <- "2"
    expect_error(do.call(runTest, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runTest, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoS0
    args <- args0
    args$testType <- "ttest"
    args$volcanoS0 <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## addAbundanceValues
    args <- args0
    args$addAbundanceValues <- "2"
    expect_error(do.call(runTest, args),
                 "'addAbundanceValues' must be of class 'logical'")
    args$addAbundanceValues <- c(TRUE, FALSE)
    expect_error(do.call(runTest, args),
                 "'addAbundanceValues' must have length 1")

    ## aName
    args <- args0
    args$addAbundanceValues <- TRUE
    args$aName <- 1
    expect_error(do.call(runTest, args),
                 "'aName' must be of class 'character'")
    args$aName <- c("pat1", "pat2")
    expect_error(do.call(runTest, args),
                 "'aName' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- "2"
    expect_error(do.call(runTest, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runTest, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- "1"
    expect_error(do.call(runTest, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runTest, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$subtractBaseline <- TRUE
    args$baselineGroup <- 1
    expect_error(do.call(runTest, args),
                 "'baselineGroup' must be of class 'character'")
    args$baselineGroup <- "missing"
    expect_error(do.call(runTest, args),
                 "baselineGroup %in% sce$group is not TRUE", fixed = TRUE)
    args$baselineGroup <- "Adnp"
    expect_error(do.call(runTest, args),
                 "%in% colnames(SummarizedExperiment::colData(sce)) is not TRUE", fixed = TRUE)

    ## Works with correct arguments (a single test)
    ## --------------------------------------------------------------------- ##
    out <- do.call(runTest, args0)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$topsets[[1]], "list")
    expect_s3_class(out$topsets[[1]]$complexes, "data.frame")
    expect_type(out$design, "list")
    expect_named(out$design, "RBC_ctrl_vs_Adnp")
    expect_type(out$design$RBC_ctrl_vs_Adnp, "list")
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast"))
    expect_named(out$design$RBC_ctrl_vs_Adnp$sampleData, "fc")
    expect_equal(out$design$RBC_ctrl_vs_Adnp$contrast, c(0, 1))
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])

    ## t-test, write results to file
    args <- args0
    args$testType <- "ttest"
    args$singleFit <- TRUE
    args$baseFileName <- tempfile()
    args$addAbundanceValues <- FALSE
    expect_message(out1 <- do.call(runTest, args), "A single model fit")
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt")))
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp_camera_complexes.txt")))
    expect_type(out1, "list")
    expect_length(out1, 9)
    expect_named(out1, c("plottitles", "plotsubtitles", "plotnotes",
                         "tests", "curveparams", "topsets", "messages",
                         "design", "featureCollections"))
    expect_s3_class(out1$tests[[1]], "data.frame")
    expect_type(out1$plotnotes[[1]], "character")
    expect_type(out1$plottitles[[1]], "character")
    expect_type(out1$plotsubtitles[[1]], "character")
    expect_type(out1$topsets[[1]], "list")
    expect_s3_class(out1$topsets[[1]]$complexes, "data.frame")
    expect_type(out1$design, "list")
    expect_equal(length(out1$design), 0)
    expect_type(out1$featureCollections, "list")
    expect_type(out1$curveparams, "list")
    expect_named(out1$curveparams[[1]], c("x", "ta", "s0", "df"))
    expect_equal(nrow(out1$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out1$tests[[1]])))
    expect_false("iBAQ.Adnp_IP04" %in% colnames(out1$tests[[1]]))
    expect_equal(out1$tests[[1]]$pid, rownames(args$sce))
    expect_equal(out1$plotnotes[[1]], "")
    expect_equal(out1$plottitles[[1]], "RBC_ctrl vs Adnp, t-test")
    expect_s4_class(out1$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out1$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out1$featureCollections$complexes)))
    expect_equal(out1$tests[[1]]$pid[1:5], out1$tests[[1]]$IDsForSTRING[1:5])
    expect_false(any(grepl("iBAQ", colnames(out1$tests[[1]]))))
    expect_equal(out1$tests[[1]]$logFC[1],
                 mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Should correlate with limma results
    expect_equal(out$tests[[1]]$pid, out1$tests[[1]]$pid)
    idx <- which(!is.na(out$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_equal(out$tests[[1]]$AveExpr, out1$tests[[1]]$AveExpr, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC, out1$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(out$tests[[1]]$t[idx], out1$tests[[1]]$sam[idx]), 0.9)
    expect_lt(cor(out$tests[[1]]$t[idx], out1$tests[[1]]$sam[idx]), 0.99)

    ## Works also with singleFit = TRUE
    args <- args0
    args$singleFit <- TRUE
    out2 <- do.call(runTest, args)
    expect_type(out2, "list")
    expect_length(out2, 9)
    expect_named(out2, c("plottitles", "plotsubtitles", "plotnotes",
                         "tests", "curveparams", "topsets", "messages",
                         "design", "featureCollections"))
    expect_equal(length(out2$tests), 1)
    expect_s3_class(out2$tests[[1]], "data.frame")
    expect_type(out2$plotnotes[[1]], "character")
    expect_type(out2$plottitles[[1]], "character")
    expect_type(out2$plotsubtitles[[1]], "character")
    expect_type(out2$topsets[[1]], "list")
    expect_s3_class(out2$topsets[[1]]$complexes, "data.frame")
    expect_type(out2$design, "list")
    expect_named(out2$design, c("design", "sampleData", "contrasts"))
    expect_true(is.matrix(out2$design$design))
    expect_equal(colnames(out2$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out2$design$sampleData, "fc")
    expect_equal(out2$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_type(out2$featureCollections, "list")
    expect_type(out2$curveparams, "list")
    expect_equal(nrow(out2$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out2$tests[[1]])))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out2$tests[[1]]))
    expect_equal(out2$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out2$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out2$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out2$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out2$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out2$featureCollections$complexes)))
    expect_equal(out2$tests[[1]]$pid[1:5], out2$tests[[1]]$IDsForSTRING[1:5])
    expect_true(any(grepl("iBAQ", colnames(out2$tests[[1]]))))
    expect_equal(out2$tests[[1]]$logFC[1],
                 mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Should correlate with limma results with singleFit = FALSE
    expect_equal(out$tests[[1]]$pid, out2$tests[[1]]$pid)
    idx <- which(!is.na(out$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_false(all(out$tests[[1]]$AveExpr == out2$tests[[1]]$AveExpr))
    expect_equal(out$tests[[1]]$logFC, out2$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(out$tests[[1]]$t[idx], out2$tests[[1]]$t[idx]), 0.9)
    expect_lt(cor(out$tests[[1]]$t[idx], out2$tests[[1]]$t[idx]), 0.99)


    ## Identical results with singleFit = TRUE/FALSE when there are only two groups
    args <- args0
    args$sce <- args$sce[, args$sce$group %in% c("RBC_ctrl", "Adnp")]
    args$singleFit <- TRUE
    outtrue <- do.call(runTest, args)
    args$singleFit <- FALSE
    outfalse <- do.call(runTest, args)
    expect_equal(names(outtrue), names(outfalse))
    expect_equal(outtrue$tests[[1]], outfalse$tests[[1]])
    for (p in c("curveparams", "plottitles", "plotnotes", "plotsubtitles",
                "topsets", "messages")) {
        expect_equal(outtrue[[p]], outfalse[[p]])
    }
    expect_equal(outtrue$design$design, outfalse$design$RBC_ctrl_vs_Adnp$design)
    expect_equal(outtrue$design$sampleData, outfalse$design$RBC_ctrl_vs_Adnp$sampleData)
    expect_equal(outtrue$design$contrasts$RBC_ctrl_vs_Adnp,
                 outfalse$design$RBC_ctrl_vs_Adnp$contrast)

    ## Different assay
    args <- args0
    args$assayForTests <- "Intensity"
    args$aName <- "LFQ.intensity"
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "LFQ.intensity.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$LFQ.intensity.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "LFQ.intensity")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])

    ## Assay with missing values
    args <- args0
    args$assayForTests <- "log2_iBAQ_withNA"
    wn <- capture_warnings({
        out <- do.call(runTest, args)
    })
    expect_match(wn, "Partial NA coefficients")
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])

    ## proDA
    args <- args0
    args$assayForTests <- "log2_iBAQ_withNA"
    args$testType <- "proDA"
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, proDA")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])

    ## With batch effect
    args <- args0
    args$sce$batch <- c("B1", "B2", "B3", "B1", "B2", "B3", "B1", "B2", "B3")
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## With batch effect, subtract baseline
    args <- args0
    args$sce$batch <- c("B1", "B2", "B3", "B1", "B2", "B3", "B1", "B2", "B3")
    args$subtractBaseline <- TRUE
    args$baselineGroup <- "Chd4BF"
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## With batch effect, single batch
    args <- args0
    args$sce$batch <- "B1"
    out <- do.call(runTest, args)
    expect_true(grepl("Only one unique", out$messages[[1]]))
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## With batch effect, single fit
    args <- args0
    args$sce$batch <- c("B1", "B2", "B1", "B2", "B1", "B2", "B1", "B2", "B1")
    args$singleFit <- TRUE
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## With batch effect, single fit, subtract baseline
    args <- args0
    args$sce$batch <- c("B1", "B2", "B3", "B1", "B2", "B3", "B1", "B2", "B3")
    args$singleFit <- TRUE
    args$subtractBaseline <- TRUE
    args$baselineGroup <- "Chd4BF"
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## With batch effect, single fit, single batch
    args <- args0
    args$sce$batch <- "B1"
    args$singleFit <- TRUE
    out <- do.call(runTest, args)
    expect_true(grepl("Only one unique", out$messages))
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## -----------------------------------------------------------------------
    ## PD data
    args0_pd <- list(
        sce = sce_pd_final,
        comparisons = list(c("HIS4KO", "WT"), c("WT", "MET6KO")),
        testType = "limma",
        assayForTests = "log2_Abundance",
        assayImputation = "imputed_Abundance",
        minNbrValidValues = 2,
        minlFC = 0.5,
        featureCollections = fcoll_pd_final,
        complexFDRThr = 0.1,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        baseFileName = NULL,
        seed = 123,
        nperm = 25,
        volcanoS0 = 0.1,
        addAbundanceValues = TRUE,
        aName = "Abundance",
        singleFit = FALSE
    )
    out <- do.call(runTest, args0_pd)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_equal(length(out$plottitles), 2)
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[2]]), 70)
    expect_true(all(c("adj.P.Val", "Abundance.HIS4KO_S07",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_equal(out$tests[[2]]$pid, rownames(sce_pd_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.5)")
    expect_equal(out$plottitles[[2]], "MET6KO vs WT, limma treat (H0: |log2FC| <= 0.5)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("WT_vs_HIS4KO_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$Abundance.HIS4KO_S05,
                 SummarizedExperiment::assay(args0_pd$sce, "Abundance")[, "HIS4KO_S05"],
                 ignore_attr = TRUE)

    ## With batch column
    tmp <- sce_pd_final
    tmp$batch <- rep(c("b1", "b2"), 8)
    args <- args0_pd
    args$sce <- tmp
    outb <- do.call(runTest, args)
    expect_type(outb, "list")
    expect_length(outb, 9)
    expect_named(outb, c("plottitles", "plotsubtitles", "plotnotes",
                         "tests", "curveparams", "topsets", "messages",
                         "design", "featureCollections"))
    expect_s3_class(outb$tests[[1]], "data.frame")
    expect_type(outb$plotnotes[[2]], "character")
    expect_type(outb$plottitles[[1]], "character")
    expect_type(outb$featureCollections, "list")
    expect_type(outb$curveparams, "list")
    expect_equal(nrow(outb$tests[[1]]), 70)
    expect_true(all(c("adj.P.Val", "Abundance.HIS4KO_S07",
                      "showInVolcano", "IDsForSTRING") %in% colnames(outb$tests[[1]])))
    expect_equal(outb$tests[[2]]$pid, rownames(sce_pd_final))
    expect_equal(substr(outb$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(outb$plottitles[[1]], "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.5)")
    expect_s4_class(outb$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(outb$featureCollections$complexes), "DFrame")
    expect_true("WT_vs_HIS4KO_FDR" %in%
                    colnames(S4Vectors::mcols(outb$featureCollections$complexes)))
    expect_equal(outb$tests[[1]]$Abundance.HIS4KO_S05,
                 SummarizedExperiment::assay(args0_pd$sce, "Abundance")[, "HIS4KO_S05"],
                 ignore_attr = TRUE)
    ## Check that test results are different compared to without batch
    expect_equal(out$tests[[2]]$pid, outb$tests[[1]]$pid)
    expect_equal(which(!is.na(out$tests[[2]]$logFC)),
                 which(!is.na(outb$tests[[2]]$logFC)))
    nonaidx <- which(!is.na(outb$tests[[2]]$logFC))
    expect_gt(length(nonaidx), 10)
    expect_false(all(round(outb$tests[[2]]$P.Value[nonaidx], 13) ==
                         round(out$tests[[2]]$P.Value[nonaidx], 13)))
    expect_true(all(round(outb$tests[[2]]$logFC[nonaidx], 13) ==
                        round(out$tests[[2]]$logFC[nonaidx], 13)))
    expect_true(all(round(outb$tests[[2]]$AveExpr[nonaidx], 13) ==
                        round(out$tests[[2]]$AveExpr[nonaidx], 13)))
    expect_true(all(round(outb$tests[[2]]$Abundance.MET6KO_S02[nonaidx], 13) ==
                        round(out$tests[[2]]$Abundance.MET6KO_S02[nonaidx], 13)))
    expect_true(all(round(outb$tests[[2]]$log2_Abundance.WT.sd[nonaidx], 13) ==
                        round(out$tests[[2]]$log2_Abundance.WT.sd[nonaidx], 13)))
})
