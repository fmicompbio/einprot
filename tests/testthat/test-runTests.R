test_that("testing works", {

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    args0 <- list(
        sce = sce_mq_final,
        comparisons = list(c("Adnp", "RBC_ctrl")),
        groupComposition = NULL,
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
        samSignificance = TRUE,
        nperm = 25,
        volcanoS0 = 0.1,
        addAbundanceValues = TRUE,
        aName = "iBAQ",
        singleFit = FALSE,
        subtractBaseline = FALSE,
        baselineGroup = "",
        extraColumns = NULL
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

    args <- args0
    rownames(args$sce)[2] <- rownames(args$sce)[1]
    expect_error(do.call(runTest, args),
                 "The row names of sce cannot contain duplicated entries.")

    args <- args0
    args$sce$sampleweights <- LETTERS[seq_len(ncol(args$sce))]
    expect_error(do.call(runTest, args),
                 "'$scesampleweights' must be of class 'numeric'",
                 fixed = TRUE)
    args$sce$sampleweights <- -seq_len(ncol(args$sce))
    expect_error(do.call(runTest, args),
                 "'$scesampleweights' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## comparisons
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

    ## groupComposition
    args <- args0
    args$groupComposition <- 1
    expect_error(do.call(runTest, args),
                 "'groupComposition' must be of class 'list'")

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

    ## samSignificance
    args <- args0
    args$testType <- "ttest"
    args$samSignificance <- "2"
    expect_error(do.call(runTest, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(runTest, args),
                 "'samSignificance' must have length 1")

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

    ## extraColumns
    args <- args0
    args$extraColumns <- 1
    expect_error(do.call(runTest, args),
                 "'extraColumns' must be of class 'character'")

    ## Fails if duplicated comparison names
    args <- args0
    args$comparisons <- list(c("Adnp", "RBC_ctrl"),
                             RBC_ctrl_vs_Adnp = c("Adnp", "Chd4BF"))
    expect_error(do.call(runTest, args),
                 "Duplicated comparison names not allowed")

    ## Fails if groupComposition contains nonexistent group
    args <- args0
    args$groupComposition <- list(tmpgroup = c("Adnp", "missing"))
    args$comparisons <- list(c("Chd4BF", "tmpgroup"))
    expect_error(do.call(runTest, args), "Missing group(s) in sce$groups", fixed = TRUE)

    ## Fails if the groups overlap
    args <- args0
    args$groupComposition <- list(tmpgroup = c("Adnp", "Chd4BF"),
                                  Adnp = "Adnp")
    args$comparisons <- list(c("tmpgroup", "Adnp"))
    expect_error(do.call(runTest, args),
                 "The same original group is part of both groups")

    ## Works with correct arguments
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
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
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
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "Zmym4", "Zmym3"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.359509, -12.540088),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "Zmym4", "Zmym3"), "t"],
                 c(-19.610901, -19.118253, -8.672101, -7.880256, -7.508872),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)
    ## save for later comparison with t-test
    outsave <- out

    ## -------------------------------------------------------------------------
    ## Works also with singleFit = TRUE (don't put anything just before here,
    ## we're comparing to the results from the previous test)
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
    expect_named(out2$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(out2$design$design))
    expect_equal(colnames(out2$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out2$design$sampleData, "fc")
    expect_equal(out2$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_null(out2$design$sampleWeights)
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
    ## Compare to values calculated manually
    expect_equal(out2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.212135, -6.804194),
                 tolerance = 0.001)
    expect_equal(out2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-19.287015, -17.572886, -8.083965, -8.052675, -6.657392),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out2$tests[[1]]$logFC), na.rm = TRUE), 110)
    ## Should correlate with limma results with singleFit = FALSE
    expect_equal(out$tests[[1]]$pid, out2$tests[[1]]$pid)
    idx <- which(!is.na(out$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_false(all(out$tests[[1]]$AveExpr == out2$tests[[1]]$AveExpr))
    expect_equal(out$tests[[1]]$logFC, out2$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(out$tests[[1]]$t[idx], out2$tests[[1]]$t[idx]), 0.9)
    expect_lt(cor(out$tests[[1]]$t[idx], out2$tests[[1]]$t[idx]), 0.99)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out2$tests[[1]]$logFC + qt(p = 0.975, df = out2$tests[[1]]$df.total) *
                     out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$logFC - qt(p = 0.975, df = out2$tests[[1]]$df.total) *
                     out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out2$tests[[1]]$t),
                               out2$tests[[1]]$df.total, lower.tail = FALSE),
                 out2$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out2$tests[[1]]$logFC / out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## singleFit = TRUE, with sampleweights all equal to 2 - should be the
    ## same as without weights (=out2)
    args <- args0
    args$singleFit <- TRUE
    args$sce$sampleweight <-
        c(Adnp_IP04 = 2, Adnp_IP05 = 2, Adnp_IP06 = 2,
          Chd4BF_IP07 = 2, Chd4BF_IP08 = 2, Chd4BF_IP09 = 2,
          RBC_ctrl_IP01 = 2, RBC_ctrl_IP02 = 2, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
    out3 <- do.call(runTest, args)
    expect_type(out3, "list")
    expect_length(out3, 9)
    expect_named(out3, c("plottitles", "plotsubtitles", "plotnotes",
                         "tests", "curveparams", "topsets", "messages",
                         "design", "featureCollections"))
    expect_equal(length(out3$tests), 1)
    expect_s3_class(out3$tests[[1]], "data.frame")
    expect_type(out3$plotnotes[[1]], "character")
    expect_type(out3$plottitles[[1]], "character")
    expect_type(out3$plotsubtitles[[1]], "character")
    expect_type(out3$topsets[[1]], "list")
    expect_s3_class(out3$topsets[[1]]$complexes, "data.frame")
    expect_type(out3$design, "list")
    expect_named(out3$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(out3$design$design))
    expect_equal(colnames(out3$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out3$design$sampleData, "fc")
    expect_equal(out3$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_equal(out3$design$sampleWeights, rep(2, 9), ignore_attr = TRUE)
    expect_type(out3$featureCollections, "list")
    expect_type(out3$curveparams, "list")
    expect_equal(nrow(out3$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out3$tests[[1]])))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out3$tests[[1]]))
    expect_equal(out3$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out3$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out3$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out3$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out3$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out3$featureCollections$complexes)))
    expect_equal(out3$tests[[1]]$pid[1:5], out3$tests[[1]]$IDsForSTRING[1:5])
    expect_true(any(grepl("iBAQ", colnames(out3$tests[[1]]))))
    expect_equal(out3$tests[[1]]$logFC[1],
                 mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Compare to values calculated manually
    expect_equal(out3$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.212135, -6.804194),
                 tolerance = 0.001)
    expect_equal(out3$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-19.287015, -17.572886, -8.083965, -8.052675, -6.657392),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out3$tests[[1]]$logFC), na.rm = TRUE), 110)
    ## Should correlate with limma results with singleFit = FALSE
    expect_equal(out2$tests[[1]]$pid, out3$tests[[1]]$pid)
    idx <- which(!is.na(out2$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_true(all(out2$tests[[1]]$AveExpr[idx] == out3$tests[[1]]$AveExpr[idx]))
    expect_equal(out2$tests[[1]]$logFC, out3$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(out2$tests[[1]]$t[idx], out3$tests[[1]]$t[idx]), 0.9999)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out3$tests[[1]]$logFC + qt(p = 0.975, df = out3$tests[[1]]$df.total) *
                     out3$tests[[1]]$se.logFC,
                 out3$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out3$tests[[1]]$logFC - qt(p = 0.975, df = out3$tests[[1]]$df.total) *
                     out3$tests[[1]]$se.logFC,
                 out3$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out3$tests[[1]]$t),
                               out3$tests[[1]]$df.total, lower.tail = FALSE),
                 out3$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out3$tests[[1]]$logFC / out3$tests[[1]]$se.logFC,
                 out3$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## singleFit = TRUE, with sampleweights, do both RBC_ctrl vs Adnp and the opposite
    args <- args0
    args$singleFit <- TRUE
    args$comparisons <- list(c("Adnp", "RBC_ctrl"),
                             c("RBC_ctrl", "Adnp"))
    args$sce$sampleweight <-
        c(Adnp_IP04 = 1, Adnp_IP05 = 6, Adnp_IP06 = 2,
          Chd4BF_IP07 = 6, Chd4BF_IP08 = 1, Chd4BF_IP09 = 5,
          RBC_ctrl_IP01 = 7, RBC_ctrl_IP02 = 1, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_equal(length(out$tests), 2)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$topsets[[1]], "list")
    expect_s3_class(out$topsets[[1]]$complexes, "data.frame")
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$plotsubtitles[[2]], "character")
    expect_type(out$topsets[[2]], "list")
    expect_s3_class(out$topsets[[2]]$complexes, "data.frame")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(out$design$design))
    expect_equal(colnames(out$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out$design$sampleData, "fc")
    expect_equal(out$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_equal(out$design$contrasts$Adnp_vs_RBC_ctrl, c(0, 0, -1))
    expect_equal(out$design$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$sampleData)])
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_equal(nrow(out$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[2]])))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out$tests[[1]]))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out$tests[[2]]))
    expect_equal(out$tests[[1]]$pid, rownames(args$sce))
    expect_equal(out$tests[[2]]$pid, rownames(args$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(out$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(out$plottitles[[2]], "Adnp vs RBC_ctrl, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_true("Adnp_vs_RBC_ctrl_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])
    expect_true(any(grepl("iBAQ", colnames(out$tests[[1]]))))
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.928299, -13.854497, -9.686006, -9.275760, -7.401809),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 -c(-7.928299, -13.854497, -9.686006, -9.275760, -7.401809),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-23.289784, -24.947386, -9.138002, -10.096597, -7.640067),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 -c(-23.289784, -24.947386, -9.138002, -10.096597, -7.640067),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$logFC), na.rm = TRUE), 110)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## singleFit = FALSE, with sampleweights
    args <- args0
    args$singleFit <- FALSE
    args$comparisons <- list(c("Adnp", "RBC_ctrl"),
                             c("RBC_ctrl", "Adnp"))
    args$sce$sampleweight <-
        c(Adnp_IP04 = 1, Adnp_IP05 = 6, Adnp_IP06 = 2,
          Chd4BF_IP07 = 6, Chd4BF_IP08 = 1, Chd4BF_IP09 = 5,
          RBC_ctrl_IP01 = 7, RBC_ctrl_IP02 = 1, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_equal(length(out$tests), 2)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$topsets[[1]], "list")
    expect_s3_class(out$topsets[[1]]$complexes, "data.frame")
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$plotsubtitles[[2]], "character")
    expect_type(out$topsets[[2]], "list")
    expect_s3_class(out$topsets[[2]]$complexes, "data.frame")
    expect_type(out$design, "list")
    expect_named(out$design, c("RBC_ctrl_vs_Adnp", "Adnp_vs_RBC_ctrl"))
    expect_type(out$design$RBC_ctrl_vs_Adnp, "list")
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_named(out$design$RBC_ctrl_vs_Adnp$sampleData, "fc")
    expect_equal(out$design$RBC_ctrl_vs_Adnp$contrast, c(0, 1))
    expect_equal(out$design$RBC_ctrl_vs_Adnp$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$RBC_ctrl_vs_Adnp$sampleData)])
    expect_type(out$design$Adnp_vs_RBC_ctrl, "list")
    expect_named(out$design$Adnp_vs_RBC_ctrl, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_named(out$design$Adnp_vs_RBC_ctrl$sampleData, "fc")
    expect_equal(out$design$Adnp_vs_RBC_ctrl$contrast, c(0, 1))
    expect_equal(out$design$Adnp_vs_RBC_ctrl$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$Adnp_vs_RBC_ctrl$sampleData)])
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparams, "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out$tests[[1]]))
    expect_equal(out$tests[[1]]$pid, rownames(args$sce))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$pid[1:5], out$tests[[1]]$IDsForSTRING[1:5])
    expect_true(any(grepl("iBAQ", colnames(out$tests[[1]]))))
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.928299, -13.854497, -9.686006, -9.275760, -7.401809),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-23.746663, -25.691437, -11.242890, -9.452125, -8.781027),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$logFC), na.rm = TRUE), 110)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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

    ## -------------------------------------------------------------------------
    ## singleFit = TRUE, no imputation filtering, add extra column,
    ## add two abundances, write to file
    args <- args0
    args$baseFileName <- tempfile()
    args$singleFit <- TRUE
    args$assayImputation <- NULL
    args$extraColumns <- c("Gene.names", "Peptides")
    args$aName <- c("iBAQ", "LFQ.intensity")
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
    expect_named(out2$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(out2$design$design))
    expect_equal(colnames(out2$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out2$design$sampleData, "fc")
    expect_equal(out2$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_null(out2$design$sampleWeights)
    expect_type(out2$featureCollections, "list")
    expect_type(out2$curveparams, "list")
    expect_equal(nrow(out2$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val", "Gene.names", "Peptides",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out2$tests[[1]])))
    expect_true("iBAQ.Adnp_IP04" %in% colnames(out2$tests[[1]]))
    expect_true("LFQ.intensity.Adnp_IP04" %in% colnames(out2$tests[[1]]))
    expect_equal(out2$tests[[1]]$iBAQ.Adnp_IP04,
                 assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$LFQ.intensity.Adnp_IP04,
                 assay(args$sce, "LFQ.intensity")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$iBAQ.RBC_ctrl_IP01,
                 assay(args$sce, "iBAQ")[, "RBC_ctrl_IP01"],
                 ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$LFQ.intensity.RBC_ctrl_IP01,
                 assay(args$sce, "LFQ.intensity")[, "RBC_ctrl_IP01"],
                 ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$iBAQ.Adnp.avg,
                 rowMeans(out2$tests[[1]][, c("iBAQ.Adnp_IP04", "iBAQ.Adnp_IP05",
                                              "iBAQ.Adnp_IP06")], na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$LFQ.intensity.Adnp.avg,
                 rowMeans(out2$tests[[1]][, c("LFQ.intensity.Adnp_IP04",
                                              "LFQ.intensity.Adnp_IP05",
                                              "LFQ.intensity.Adnp_IP06")], na.rm = TRUE),
                 ignore_attr = TRUE)
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
    ## Compare to values calculated manually
    expect_equal(out2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.212135, -6.804194),
                 tolerance = 0.001)
    expect_equal(out2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-18.506681, -17.525983, -7.876344, -8.305873, -6.299199),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out2$tests[[1]]$logFC), na.rm = TRUE), 150)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out2$tests[[1]]$logFC + qt(p = 0.975, df = out2$tests[[1]]$df.total) *
                     out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out2$tests[[1]]$logFC - qt(p = 0.975, df = out2$tests[[1]]$df.total) *
                     out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out2$tests[[1]]$t),
                               out2$tests[[1]]$df.total, lower.tail = FALSE),
                 out2$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out2$tests[[1]]$logFC / out2$tests[[1]]$se.logFC,
                 out2$tests[[1]]$t, ignore_attr = TRUE)
    ## Check the saved file
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt")))
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp_camera_complexes.txt")))
    fl <- read.delim(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt"))
    expect_equal(nrow(fl), sum(out2$tests$RBC_ctrl_vs_Adnp$showInVolcano))
    expect_true(all(fl$pid %in% rownames(args$sce)))
    tmpsce <- args$sce[fl$pid, ]
    tmpres <- out2$tests$RBC_ctrl_vs_Adnp[match(fl$pid, out2$tests$RBC_ctrl_vs_Adnp$pid), ]
    expect_equal(fl$pid, tmpres$pid)
    expect_equal(fl$logFC, tmpres$logFC)
    expect_equal(fl$P.Value, tmpres$P.Value)
    expect_equal(fl$iBAQ.Adnp_IP06, assay(tmpsce, "iBAQ")[, "Adnp_IP06"],
                 ignore_attr = TRUE)
    expect_equal(fl$LFQ.intensity.Adnp_IP06, assay(tmpsce, "LFQ.intensity")[, "Adnp_IP06"],
                 ignore_attr = TRUE)
    expect_equal(fl$iBAQ.RBC_ctrl.avg,
                 rowMeans(assay(tmpsce, "iBAQ")[, c("RBC_ctrl_IP01", "RBC_ctrl_IP02",
                                                    "RBC_ctrl_IP03")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$LFQ.intensity.RBC_ctrl.avg,
                 rowMeans(assay(tmpsce, "LFQ.intensity")[, c("RBC_ctrl_IP01",
                                                             "RBC_ctrl_IP02",
                                                             "RBC_ctrl_IP03")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$iBAQ.RBC_ctrl.avg,
                 rowMeans(fl[, c("iBAQ.RBC_ctrl_IP01", "iBAQ.RBC_ctrl_IP02",
                                 "iBAQ.RBC_ctrl_IP03")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$LFQ.intensity.RBC_ctrl.avg,
                 rowMeans(fl[, c("LFQ.intensity.RBC_ctrl_IP01",
                                 "LFQ.intensity.RBC_ctrl_IP02",
                                 "LFQ.intensity.RBC_ctrl_IP03")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$iBAQ.RBC_ctrl.sd,
                 sqrt(matrixStats::rowVars(assay(
                     tmpsce, "iBAQ")[, c("RBC_ctrl_IP01", "RBC_ctrl_IP02",
                                         "RBC_ctrl_IP03")],
                     na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$LFQ.intensity.RBC_ctrl.sd,
                 sqrt(matrixStats::rowVars(assay(
                     tmpsce, "LFQ.intensity")[, c("RBC_ctrl_IP01", "RBC_ctrl_IP02",
                                                  "RBC_ctrl_IP03")],
                     na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_iBAQ.Adnp.avg,
                 rowMeans(log2(assay(tmpsce, "iBAQ"))[, c("Adnp_IP04", "Adnp_IP05",
                                                          "Adnp_IP06")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    tmplfq <- assay(tmpsce, "LFQ.intensity")  ## need to replace 0 by NA
    tmplfq[tmplfq == 0] <- NA
    expect_equal(fl$log2_LFQ.intensity.Adnp.avg,
                 rowMeans(log2(tmplfq)[, c("Adnp_IP04", "Adnp_IP05",
                                           "Adnp_IP06")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_iBAQ.Adnp.sd,
                 sqrt(matrixStats::rowVars(log2(assay(
                     tmpsce, "iBAQ"))[, c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06")],
                     na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_LFQ.intensity.Adnp.sd,
                 sqrt(matrixStats::rowVars(log2(tmplfq)[, c("Adnp_IP04", "Adnp_IP05",
                                                            "Adnp_IP06")],
                     na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    expect_equal(outsave$tests[[1]]$pid, out1$tests[[1]]$pid)
    idx <- which(!is.na(outsave$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_equal(outsave$tests[[1]]$AveExpr, out1$tests[[1]]$AveExpr, ignore_attr = TRUE)
    expect_equal(outsave$tests[[1]]$logFC, out1$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(outsave$tests[[1]]$t[idx], out1$tests[[1]]$sam[idx]), 0.9)
    expect_lt(cor(outsave$tests[[1]]$t[idx], out1$tests[[1]]$sam[idx]), 0.99)

    ## Check the saved file
    fl <- read.delim(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt"))
    expect_equal(nrow(fl), sum(out1$tests$RBC_ctrl_vs_Adnp$showInVolcano, na.rm = TRUE))
    expect_true(all(fl$pid %in% rownames(args$sce)))
    tmpsce <- args$sce[fl$pid, ]
    tmpres <- out1$tests$RBC_ctrl_vs_Adnp[match(fl$pid, out1$tests$RBC_ctrl_vs_Adnp$pid), ]
    expect_equal(fl$pid, tmpres$pid)
    expect_equal(fl$logFC, tmpres$logFC)
    expect_equal(fl$P.Value, tmpres$P.Value)
    expect_null(fl$iBAQ.Adnp_IP06) ## addAbundanceValues = FALSE

    ## -------------------------------------------------------------------------
    ## t-test, don't use SAM statistic for significance
    args <- args0
    args$testType <- "ttest"
    args$samSignificance <- FALSE
    args$singleFit <- TRUE
    args$addAbundanceValues <- FALSE
    expect_message(out1b <- do.call(runTest, args), "A single model fit")
    expect_type(out1b, "list")
    expect_length(out1b, 9)
    expect_named(out1b, c("plottitles", "plotsubtitles", "plotnotes",
                          "tests", "curveparams", "topsets", "messages",
                          "design", "featureCollections"))
    expect_s3_class(out1b$tests[[1]], "data.frame")
    expect_type(out1b$plotnotes[[1]], "character")
    expect_type(out1b$plottitles[[1]], "character")
    expect_type(out1b$plotsubtitles[[1]], "character")
    expect_type(out1b$topsets[[1]], "list")
    expect_s3_class(out1b$topsets[[1]]$complexes, "data.frame")
    expect_type(out1b$design, "list")
    expect_equal(length(out1b$design), 0)
    expect_type(out1b$featureCollections, "list")
    expect_type(out1b$curveparams, "list")
    expect_equal(out1b$curveparams[[1]], list())
    expect_equal(nrow(out1b$tests[[1]]), 150)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out1b$tests[[1]])))
    expect_false("iBAQ.Adnp_IP04" %in% colnames(out1b$tests[[1]]))
    expect_equal(out1b$tests[[1]]$pid, rownames(args$sce))
    expect_equal(out1b$plotnotes[[1]], "")
    expect_equal(out1b$plottitles[[1]], "RBC_ctrl vs Adnp, t-test")
    expect_s4_class(out1b$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out1b$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out1b$featureCollections$complexes)))
    expect_equal(out1b$tests[[1]]$pid[1:5], out1b$tests[[1]]$IDsForSTRING[1:5])
    expect_false(any(grepl("iBAQ", colnames(out1b$tests[[1]]))))
    expect_equal(out1b$tests[[1]]$logFC[1],
                 mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Compare to results with SAM significance
    expect_equal(out1$tests[[1]]$pid, out1b$tests[[1]]$pid)
    expect_equal(out1$tests[[1]]$t, out1b$tests[[1]]$t)
    expect_equal(out1$tests[[1]]$sam, out1b$tests[[1]]$sam)
    expect_equal(out1$tests[[1]]$logFC, out1b$tests[[1]]$logFC)
    expect_equal(out1$tests[[1]]$P.Value, out1b$tests[[1]]$P.Value)
    expect_equal(out1$tests[[1]]$adj.P.Val, out1b$tests[[1]]$adj.P.Val)
    idx <- which(!is.na(out1$tests[[1]]$showInVolcano))
    expect_false(all(out1$tests[[1]]$showInVolcano[idx] ==
                         out1b$tests[[1]]$showInVolcano[idx]))

    ## t-test, very small adj p-value threshold
    args <- args0
    args$testType <- "ttest"
    args$singleFit <- TRUE
    args$volcanoAdjPvalThr <- 0
    expect_message(outsm <- do.call(runTest, args), "A single model fit")
    expect_type(outsm, "list")
    expect_length(outsm, 9)
    expect_true(all(!outsm$tests[[1]]$showInVolcano[!is.na(outsm$tests[[1]]$showInVolcano)]))

    ## -------------------------------------------------------------------------
    ## Merged groups, save to file
    args <- args0
    args$baseFileName <- tempfile()
    args$groupComposition <- list(rbc_adnp = c("RBC_ctrl", "Adnp"))
    args$comparisons <- list(c("Adnp", "RBC_ctrl"), c("rbc_adnp", "Chd4BF"))
    outm <- do.call(runTest, args)
    expect_type(outm, "list")
    expect_length(outm, 9)
    expect_named(outm, c("plottitles", "plotsubtitles", "plotnotes",
                         "tests", "curveparams", "topsets", "messages",
                         "design", "featureCollections"))
    expect_length(outm$tests, 2)
    expect_named(outm$tests, c("RBC_ctrl_vs_Adnp", "Chd4BF_vs_rbc_adnp"))
    expect_s3_class(outm$tests[[1]], "data.frame")
    expect_type(outm$plotnotes[[1]], "character")
    expect_type(outm$plottitles[[1]], "character")
    expect_type(outm$plotsubtitles[[1]], "character")
    expect_type(outm$topsets[[1]], "list")
    expect_s3_class(outm$topsets[[1]]$complexes, "data.frame")
    expect_s3_class(outm$tests[[2]], "data.frame")
    expect_type(outm$plotnotes[[2]], "character")
    expect_type(outm$plottitles[[2]], "character")
    expect_type(outm$plotsubtitles[[2]], "character")
    expect_type(outm$topsets[[2]], "list")
    expect_s3_class(outm$topsets[[2]]$complexes, "data.frame")
    expect_type(outm$design, "list")
    expect_named(outm$design, c("RBC_ctrl_vs_Adnp", "Chd4BF_vs_rbc_adnp"))
    expect_type(outm$design$RBC_ctrl_vs_Adnp, "list")
    expect_named(outm$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                 "sampleWeights"))
    expect_named(outm$design$RBC_ctrl_vs_Adnp$sampleData, "fc")
    expect_equal(outm$design$RBC_ctrl_vs_Adnp$contrast, c(0, 1))
    expect_null(outm$design$RBC_ctrl_vs_Adnp$sampleWeights)
    expect_type(outm$design$Chd4BF_vs_rbc_adnp, "list")
    expect_named(outm$design$Chd4BF_vs_rbc_adnp, c("design", "sampleData", "contrast",
                                                   "sampleWeights"))
    expect_named(outm$design$Chd4BF_vs_rbc_adnp$sampleData, "fc")
    expect_equal(outm$design$Chd4BF_vs_rbc_adnp$contrast, c(0, 1))
    expect_null(outm$design$Chd4BF_vs_rbc_adnp$sampleWeights)
    expect_type(outm$featureCollections, "list")
    expect_type(outm$curveparams, "list")
    expect_equal(nrow(outm$tests[[1]]), 150)
    expect_equal(nrow(outm$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in%
                        colnames(outm$tests[[1]])))
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in%
                        colnames(outm$tests[[2]])))
    expect_equal(outm$tests[[1]]$pid, rownames(args0$sce))
    expect_equal(outm$tests[[2]]$pid, rownames(args0$sce))
    expect_equal(substr(outm$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(outm$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(outm$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(outm$plottitles[[2]], "Chd4BF vs rbc_adnp, limma")
    expect_s4_class(outm$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(outm$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(outm$featureCollections$complexes)))
    expect_true("Chd4BF_vs_rbc_adnp_FDR" %in%
                    colnames(S4Vectors::mcols(outm$featureCollections$complexes)))
    expect_equal(outm$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(outm$tests[[2]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_true(all(paste0("iBAQ.", c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
                                      "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
                                      "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")) %in%
                        colnames(outm$tests[[2]])))
    expect_equal(outm$tests[[1]]$pid[1:5], outm$tests[[1]]$IDsForSTRING[1:5])
    expect_equal(outm$tests[[2]]$pid[1:5], outm$tests[[2]]$IDsForSTRING[1:5])
    ## Compare to values calculated manually
    expect_equal(outm$tests[[1]][c("Adnp", "Chd4", "Dhx9", "Zmym4", "Zmym3"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.359509, -12.540088),
                 tolerance = 0.001)
    expect_equal(outm$tests[[1]][c("Adnp", "Chd4", "Dhx9", "Zmym4", "Zmym3"), "t"],
                 c(-19.610901, -19.118253, -8.672101, -7.880256, -7.508872),
                 tolerance = 0.001)
    expect_equal(outm$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Mta1.E9PX23", "Atp5c1"), "logFC"],
                 c(12.704305, 14.285579, 10.107046, 8.235868),
                 tolerance = 0.001)
    expect_equal(outm$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Mta1.E9PX23", "Atp5c1"), "t"],
                 c(9.621916, 9.245356, 7.524309, 7.062959),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(outm$tests[[1]]$logFC + qt(p = 0.975, df = outm$tests[[1]]$df.total) *
                     outm$tests[[1]]$se.logFC,
                 outm$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(outm$tests[[1]]$logFC - qt(p = 0.975, df = outm$tests[[1]]$df.total) *
                     outm$tests[[1]]$se.logFC,
                 outm$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(outm$tests[[1]]$t),
                               outm$tests[[1]]$df.total, lower.tail = FALSE),
                 outm$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(outm$tests[[1]]$logFC / outm$tests[[1]]$se.logFC,
                 outm$tests[[1]]$t, ignore_attr = TRUE)
    ## Check the saved file
    expect_true(file.exists(paste0(args$baseFileName, "_testres_Chd4BF_vs_rbc_adnp.txt")))
    expect_true(file.exists(paste0(args$baseFileName, "_testres_Chd4BF_vs_rbc_adnp_camera_complexes.txt")))
    fl <- read.delim(paste0(args$baseFileName, "_testres_Chd4BF_vs_rbc_adnp.txt"))
    expect_equal(nrow(fl), sum(outm$tests$Chd4BF_vs_rbc_adnp$showInVolcano, na.rm = TRUE))
    expect_true(all(fl$pid %in% rownames(args$sce)))
    tmpsce <- args$sce[fl$pid, ]
    tmpres <- outm$tests$Chd4BF_vs_rbc_adnp[match(fl$pid, outm$tests$Chd4BF_vs_rbc_adnp$pid), ]
    expect_equal(fl$pid, tmpres$pid)
    expect_equal(fl$logFC, tmpres$logFC)
    expect_equal(fl$P.Value, tmpres$P.Value)
    expect_equal(fl$iBAQ.Adnp_IP06, assay(tmpsce, "iBAQ")[, "Adnp_IP06"],
                 ignore_attr = TRUE)
    # expect_equal(fl$iBAQ.rbc_adnp.avg,
    #              rowMeans(assay(tmpsce, "iBAQ")[, c("RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03",
    #                                                 "Adnp_IP04", "Adnp_IP05", "Adnp_IP06")],
    #                       na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$iBAQ.RBC_ctrl.avg,
                 rowMeans(fl[, c("iBAQ.RBC_ctrl_IP01", "iBAQ.RBC_ctrl_IP02", "iBAQ.RBC_ctrl_IP03")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$iBAQ.RBC_ctrl.sd,
                 sqrt(matrixStats::rowVars(assay(tmpsce, "iBAQ")[, c("RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")],
                                           na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_iBAQ.Adnp.avg,
                 rowMeans(log2(assay(tmpsce, "iBAQ"))[, c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_iBAQ.Adnp.sd,
                 sqrt(matrixStats::rowVars(log2(assay(tmpsce, "iBAQ"))[, c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06")],
                                           na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## Single fit, merged groups
    args <- args0
    args$groupComposition <- list(rbc_adnp = c("RBC_ctrl", "Adnp"))
    args$comparisons <- list(c("Adnp", "RBC_ctrl"), c("rbc_adnp", "Chd4BF"))
    args$singleFit <- TRUE
    outm2 <- do.call(runTest, args)
    expect_type(outm2, "list")
    expect_length(outm2, 9)
    expect_named(outm2, c("plottitles", "plotsubtitles", "plotnotes",
                          "tests", "curveparams", "topsets", "messages",
                          "design", "featureCollections"))
    expect_length(outm2$tests, 2)
    expect_named(outm2$tests, c("RBC_ctrl_vs_Adnp", "Chd4BF_vs_rbc_adnp"))
    expect_s3_class(outm2$tests[[1]], "data.frame")
    expect_type(outm2$plotnotes[[1]], "character")
    expect_type(outm2$plottitles[[1]], "character")
    expect_type(outm2$plotsubtitles[[1]], "character")
    expect_type(outm2$topsets[[1]], "list")
    expect_s3_class(outm2$topsets[[1]]$complexes, "data.frame")
    expect_s3_class(outm2$tests[[2]], "data.frame")
    expect_type(outm2$plotnotes[[2]], "character")
    expect_type(outm2$plottitles[[2]], "character")
    expect_type(outm2$plotsubtitles[[2]], "character")
    expect_type(outm2$topsets[[2]], "list")
    expect_s3_class(outm2$topsets[[2]]$complexes, "data.frame")
    expect_type(outm2$design, "list")
    expect_named(outm2$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(outm2$design$design))
    expect_equal(colnames(outm2$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(outm2$design$sampleData, "fc")
    expect_equal(outm2$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_equal(outm2$design$contrasts$Chd4BF_vs_rbc_adnp, c(0, 1, -0.5))
    expect_null(outm2$design$sampleWeights)
    expect_type(outm2$featureCollections, "list")
    expect_type(outm2$curveparams, "list")
    expect_equal(nrow(outm2$tests[[1]]), 150)
    expect_equal(nrow(outm2$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in%
                        colnames(outm2$tests[[1]])))
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in%
                        colnames(outm2$tests[[2]])))
    expect_equal(outm2$tests[[1]]$pid, rownames(args0$sce))
    expect_equal(outm2$tests[[2]]$pid, rownames(args0$sce))
    expect_equal(substr(outm2$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(outm2$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(outm2$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(outm2$plottitles[[2]], "Chd4BF vs rbc_adnp, limma")
    expect_s4_class(outm2$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(outm2$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(outm2$featureCollections$complexes)))
    expect_true("Chd4BF_vs_rbc_adnp_FDR" %in%
                    colnames(S4Vectors::mcols(outm2$featureCollections$complexes)))
    expect_equal(outm2$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(outm2$tests[[2]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(outm2$tests[[1]]$pid[1:5], outm2$tests[[1]]$IDsForSTRING[1:5])
    expect_equal(outm2$tests[[2]]$pid[1:5], outm2$tests[[2]]$IDsForSTRING[1:5])
    expect_true(any(grepl("iBAQ", colnames(outm2$tests[[1]]))))
    expect_equal(outm2$tests[[1]]$logFC[1],
                 mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Compare to values calculated manually
    expect_equal(outm2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "logFC"],
                 c(-7.840171, -13.171775, -8.960746, -8.212135, -6.804194),
                 tolerance = 0.001)
    expect_equal(outm2$tests[[1]][c("Adnp", "Chd4", "Dhx9", "RBM8", "Ssb"), "t"],
                 c(-19.287015, -17.572886, -8.083965, -8.052675, -6.657392),
                 tolerance = 0.001)
    ## Compare to values calculated manually
    expect_equal(outm2$tests[[2]][c("Chd4", "Mta1.F8WHY8", "Mbd3", "Adnp"), "logFC"],
                 c(8.607417, 14.285579, 12.704305, 3.298679),
                 tolerance = 0.001)
    expect_equal(outm2$tests[[2]][c("Chd4", "Mta1.F8WHY8", "Mbd3", "Adnp"), "t"],
                 c(13.267899, 10.619452, 10.367206, 9.170229),
                 tolerance = 0.001)
    ## Should correlate with limma results with singleFit = FALSE
    expect_equal(outm$tests[[1]]$pid, outm2$tests[[1]]$pid)
    idx <- which(!is.na(outm$tests[[1]]$logFC))
    expect_gt(length(idx), 10)
    expect_false(all(outm$tests[[1]]$AveExpr == outm2$tests[[1]]$AveExpr))
    expect_equal(outm$tests[[1]]$logFC, outm2$tests[[1]]$logFC, ignore_attr = TRUE)
    expect_gt(cor(outm$tests[[1]]$t[idx], outm2$tests[[1]]$t[idx]), 0.9)
    expect_lt(cor(outm$tests[[1]]$t[idx], outm2$tests[[1]]$t[idx]), 0.99)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(outm2$tests[[1]]$logFC + qt(p = 0.975, df = outm2$tests[[1]]$df.total) *
                     outm2$tests[[1]]$se.logFC,
                 outm2$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(outm2$tests[[1]]$logFC - qt(p = 0.975, df = outm2$tests[[1]]$df.total) *
                     outm2$tests[[1]]$se.logFC,
                 outm2$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(outm2$tests[[1]]$t),
                               outm2$tests[[1]]$df.total, lower.tail = FALSE),
                 outm2$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(outm2$tests[[1]]$logFC / outm2$tests[[1]]$se.logFC,
                 outm2$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## Merged groups, with batch effect
    args <- args0
    args$groupComposition <- list(rbc_adnp = c("RBC_ctrl", "Adnp"))
    args$comparisons <- list(c("Adnp", "RBC_ctrl"), c("rbc_adnp", "Chd4BF"))
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
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$design, "list")
    expect_named(out$design, c("RBC_ctrl_vs_Adnp", "Chd4BF_vs_rbc_adnp"))
    expect_type(out$design$RBC_ctrl_vs_Adnp, "list")
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_named(out$design$RBC_ctrl_vs_Adnp$sampleData, c("fc", "bc"))
    expect_equal(out$design$RBC_ctrl_vs_Adnp$contrast, c(0, 0, 0, 1))
    expect_null(out$design$RBC_ctrl_vs_Adnp$sampleWeights)
    expect_type(out$design$Chd4BF_vs_rbc_adnp, "list")
    expect_named(out$design$Chd4BF_vs_rbc_adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_named(out$design$Chd4BF_vs_rbc_adnp$sampleData, c("fc", "bc"))
    expect_equal(out$design$Chd4BF_vs_rbc_adnp$contrast, c(0, 0, 0, 1))
    expect_null(out$design$Chd4BF_vs_rbc_adnp$sampleWeights)
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_type(out$curveparams[[2]], "list")
    expect_equal(nrow(out$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[2]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(out$tests[[2]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(out$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(out$plottitles[[2]], "Chd4BF vs rbc_adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[2]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Dhx9", "RBM8", "Zmym3"), "logFC"],
                 c(-13.171775, -7.840171, -8.960746, -8.212135, -12.540088),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Dhx9", "RBM8", "Zmym3"), "t"],
                 c(-21.414593, -20.341038, -8.374547, -7.661423, -7.468584),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Pogz", "Zfp462.B1AWL2"), "logFC"],
                 c(12.704305, 14.285579, 10.024301, 11.030955),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Pogz", "Zfp462.B1AWL2"), "t"],
                 c(8.693193, 8.146834, 7.774483, 7.091121),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## Merged groups, with batch effect, with sample weights
    args <- args0
    args$groupComposition <- list(rbc_adnp = c("RBC_ctrl", "Adnp"))
    args$comparisons <- list(c("Adnp", "RBC_ctrl"), c("rbc_adnp", "Chd4BF"))
    args$sce$batch <- c("B1", "B2", "B3", "B1", "B2", "B3", "B1", "B2", "B3")
    args$sce$sampleweight <-
        c(Adnp_IP04 = 1, Adnp_IP05 = 6, Adnp_IP06 = 2,
          Chd4BF_IP07 = 6, Chd4BF_IP08 = 1, Chd4BF_IP09 = 5,
          RBC_ctrl_IP01 = 7, RBC_ctrl_IP02 = 1, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 9)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "curveparams", "topsets", "messages",
                        "design", "featureCollections"))
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$design, "list")
    expect_named(out$design, c("RBC_ctrl_vs_Adnp", "Chd4BF_vs_rbc_adnp"))
    expect_type(out$design$RBC_ctrl_vs_Adnp, "list")
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_named(out$design$RBC_ctrl_vs_Adnp$sampleData, c("fc", "bc"))
    expect_equal(out$design$RBC_ctrl_vs_Adnp$contrast, c(0, 0, 0, 1))
    expect_equal(out$design$RBC_ctrl_vs_Adnp$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$RBC_ctrl_vs_Adnp$sampleData)])
    expect_type(out$design$Chd4BF_vs_rbc_adnp, "list")
    expect_named(out$design$Chd4BF_vs_rbc_adnp, c("design", "sampleData", "contrast",
                                                  "sampleWeights"))
    expect_named(out$design$Chd4BF_vs_rbc_adnp$sampleData, c("fc", "bc"))
    expect_equal(out$design$Chd4BF_vs_rbc_adnp$contrast, c(0, 0, 0, 1))
    expect_equal(out$design$Chd4BF_vs_rbc_adnp$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$Chd4BF_vs_rbc_adnp$sampleData)])
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_type(out$curveparams[[2]], "list")
    expect_equal(nrow(out$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[2]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(out$tests[[2]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(out$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(out$plottitles[[2]], "Chd4BF vs rbc_adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[2]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Pogz", "Zfp462.B1AWL2"), "logFC"],
                 c(13.237967, 15.502178, 9.519576, 10.742178),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Mbd3", "Mta1.F8WHY8", "Pogz", "Zfp462.B1AWL2"), "t"],
                 c(10.793285, 10.408011, 8.845311, 8.634606),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## Merged groups, with batch effect, single fit
    args <- args0
    args$groupComposition <- list(rbc_adnp = c("RBC_ctrl", "Adnp"))
    args$comparisons <- list(c("Adnp", "RBC_ctrl"), c("rbc_adnp", "Chd4BF"))
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
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_named(out$design$sampleData, c("fc", "bc"))
    expect_equal(out$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 0, 1))
    expect_equal(out$design$contrasts$Chd4BF_vs_rbc_adnp, c(0, 0, 1, -0.5))
    expect_null(out$design$sampleWeights)
    expect_type(out$curveparams[[1]], "list")
    expect_equal(nrow(out$tests[[1]]), 150)
    expect_type(out$curveparams[[2]], "list")
    expect_equal(nrow(out$tests[[2]]), 150)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[1]])))
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$tests[[2]])))
    expect_equal(out$tests[[1]]$pid, rownames(sce_mq_final))
    expect_equal(out$tests[[2]]$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnotes[[1]], 1, 8), "df.prior")
    expect_equal(substr(out$plotnotes[[2]], 1, 8), "df.prior")
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma")
    expect_equal(out$plottitles[[2]], "Chd4BF vs rbc_adnp, limma")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$tests[[2]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "logFC"],
                 c(-13.171775, -7.840171, -12.540088, -8.960746, -8.212135),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "t"],
                 c(-22.414726, -19.340982, -8.067907, -7.801567, -7.774172),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Chd4", "Mta1.F8WHY8", "Mbd3", "Adnp"), "logFC"],
                 c(8.171809, 13.463680, 12.375919, 3.418901),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][c("Chd4", "Mta1.F8WHY8", "Mbd3", "Adnp"), "t"],
                 c(15.068185, 11.455142, 9.552536, 9.026281),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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

    ## -------------------------------------------------------------------------
    ## proDA, singleFit = TRUE
    args <- args0
    args$assayForTests <- "log2_iBAQ_withNA"
    args$testType <- "proDA"
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

    ## -------------------------------------------------------------------------
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
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Dhx9", "RBM8", "Zmym3"), "logFC"],
                 c(-13.171775, -7.840171, -8.960746, -8.212135, -12.540088),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Dhx9", "RBM8", "Zmym3"), "t"],
                 c(-21.414593, -20.341038, -8.374547, -7.661423, -7.468584),
                 tolerance = 0.001)

    ## -------------------------------------------------------------------------
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    expect_type(out$design, "list")
    expect_named(out$design, "RBC_ctrl_vs_Adnp")
    expect_named(out$design$RBC_ctrl_vs_Adnp, c("design", "sampleData", "contrast",
                                                "sampleWeights"))
    expect_equal(out$design$RBC_ctrl_vs_Adnp$contrast, c(0, 1))
    expect_null(out$design$RBC_ctrl_vs_Adnp$sampleWeights)
    expect_named(out$design$RBC_ctrl_vs_Adnp$sampleData, c("fc", "bc"))
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_named(out$design$sampleData, c("fc", "bc"))
    expect_named(out$design$contrasts, "RBC_ctrl_vs_Adnp")
    expect_equal(out$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 0, 1))
    expect_null(out$design$sampleWeights)
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
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "logFC"],
                 c(-13.171775, -7.840171, -12.540088, -8.960746, -8.212135),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "t"],
                 c(-22.414726, -19.340982, -8.067907, -7.801567, -7.774172),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## With batch effect, single fit, with sampleweights
    args <- args0
    args$sce$batch <- c("B1", "B2", "B1", "B2", "B1", "B2", "B1", "B2", "B1")
    args$singleFit <- TRUE
    args$sce$sampleweight <-
        c(Adnp_IP04 = 1, Adnp_IP05 = 6, Adnp_IP06 = 2,
          Chd4BF_IP07 = 6, Chd4BF_IP08 = 1, Chd4BF_IP09 = 5,
          RBC_ctrl_IP01 = 7, RBC_ctrl_IP02 = 1, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
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
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_true(is.matrix(out$design$design))
    expect_equal(colnames(out$design$design), c("(Intercept)", "bcB2", "fcChd4BF", "fcRBC_ctrl"))
    expect_named(out$design$sampleData, c("fc", "bc"))
    expect_equal(out$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 0, 1))
    expect_equal(out$design$sampleWeights,
                 args$sce$sampleweight[rownames(out$design$sampleData)])
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "logFC"],
                 c(-13.368100, -8.183052, -11.991345, -9.386267, -8.982293),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Chd4", "Adnp", "Zmym3", "Dhx9", "RBM8"), "t"],
                 c(-21.738660, -21.106371, -7.198843, -7.161905, -7.915411),
                 tolerance = 0.001)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Dhx9", "Ssb", "Baz2a", "Adnp", "Chd4"), "logFC"],
                 c(-8.960746, -6.804194, -6.700705, -7.840171, -13.171775),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Dhx9", "Ssb", "Baz2a", "Adnp", "Chd4"), "t"],
                 c(-5.865848, -4.711314, -4.472601, -12.534539, -9.757926),
                 tolerance = 0.001)

    ## -------------------------------------------------------------------------
    ## With batch effect, single fit, subtract baseline, with sample weights
    args <- args0
    args$sce$batch <- c("B1", "B2", "B3", "B1", "B2", "B3", "B1", "B2", "B3")
    args$sce$sampleweight <-
        c(Adnp_IP04 = 1, Adnp_IP05 = 6, Adnp_IP06 = 2,
          Chd4BF_IP07 = 6, Chd4BF_IP08 = 1, Chd4BF_IP09 = 5,
          RBC_ctrl_IP01 = 7, RBC_ctrl_IP02 = 1, RBC_ctrl_IP03 = 2)[colnames(args$sce)]
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][c("Dhx9", "Ssb", "Baz2a", "Adnp", "Chd4"), "logFC"],
                 c(-8.193482, -6.340554, -5.400601, -7.452521, -15.161861),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][c("Dhx9", "Ssb", "Baz2a", "Adnp", "Chd4"), "t"],
                 c(-6.600978, -5.536790, -4.328110, -13.293583, -14.716825),
                 tolerance = 0.001)

    ## -------------------------------------------------------------------------
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
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts", "sampleWeights"))
    expect_named(out$design$sampleData, c("fc", "bc"))
    expect_equal(out$design$contrasts$RBC_ctrl_vs_Adnp, c(0, 0, 1))
    expect_equal(colnames(out$design$design), c("(Intercept)", "fcChd4BF", "fcRBC_ctrl"))
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## With batch effect, single fit, single batch, minlFC > 0
    args <- args0
    args$sce$batch <- "B1"
    args$singleFit <- TRUE
    args$minlFC <- 1
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
    expect_equal(out$plottitles[[1]], "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 1)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$tests[[1]]$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values (not working out of the box here since minlFC != 0)
    # expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
    #                            out$tests[[1]]$df.total, lower.tail = FALSE),
    #              out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    # expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
    #              out$tests[[1]]$t, ignore_attr = TRUE)

    ## -----------------------------------------------------------------------
    ## PD data
    args0_pd <- list(
        sce = sce_pd_final,
        comparisons = list(c("HIS4KO", "WT"), c("WT", "MET6KO")),
        groupComposition = NULL,
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
        singleFit = FALSE,
        extraColumns = NULL
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
    expect_type(out$design, "list")
    expect_named(out$design, c("WT_vs_HIS4KO", "MET6KO_vs_WT"))
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast",
                                            "sampleWeights"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, "fc")
    expect_equal(colnames(out$design$WT_vs_HIS4KO$design), c("(Intercept)", "fcWT"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 1))
    expect_null(out$design$WT_vs_HIS4KO$sampleWeights)
    expect_type(out$design$MET6KO_vs_WT, "list")
    expect_named(out$design$MET6KO_vs_WT, c("design", "sampleData", "contrast",
                                            "sampleWeights"))
    expect_named(out$design$MET6KO_vs_WT$sampleData, "fc")
    expect_equal(colnames(out$design$MET6KO_vs_WT$design), c("(Intercept)", "fcMET6KO"))
    expect_equal(out$design$MET6KO_vs_WT$contrast, c(0, 1))
    expect_null(out$design$MET6KO_vs_WT$sampleWeights)
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    # expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
    #                            out$tests[[1]]$df.total, lower.tail = FALSE),
    #              out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    # expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
    #              out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## With batch column, write to file
    tmp <- sce_pd_final
    tmp$batch <- rep(c("b1", "b2"), 8)
    args <- args0_pd
    args$sce <- tmp
    args$baseFileName <- tempfile()
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
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(outb$tests[[1]]$logFC + qt(p = 0.975, df = outb$tests[[1]]$df.total) *
                     outb$tests[[1]]$se.logFC,
                 outb$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(outb$tests[[1]]$logFC - qt(p = 0.975, df = outb$tests[[1]]$df.total) *
                     outb$tests[[1]]$se.logFC,
                 outb$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    # expect_equal(2 * stats::pt(abs(outb$tests[[1]]$t),
    #                            outb$tests[[1]]$df.total, lower.tail = FALSE),
    #              outb$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    # expect_equal(outb$tests[[1]]$logFC / outb$tests[[1]]$se.logFC,
    #              outb$tests[[1]]$t, ignore_attr = TRUE)

    ## Check the saved file
    expect_true(file.exists(paste0(args$baseFileName, "_testres_WT_vs_HIS4KO.txt")))
    fl <- read.delim(paste0(args$baseFileName, "_testres_WT_vs_HIS4KO.txt"))
    expect_equal(nrow(fl), sum(outb$tests$WT_vs_HIS4KO$showInVolcano, na.rm = TRUE))
    expect_true(all(fl$pid %in% rownames(args$sce)))
    tmpsce <- args$sce[fl$pid, ]
    tmpres <- outb$tests$WT_vs_HIS4KO[match(fl$pid, outb$tests$WT_vs_HIS4KO$pid), ]
    expect_equal(fl$pid, tmpres$pid)
    expect_equal(fl$logFC, tmpres$logFC)
    expect_equal(fl$P.Value, tmpres$P.Value)
    expect_equal(fl$Abundance.WT_S14, assay(tmpsce, "Abundance")[, "WT_S14"],
                 ignore_attr = TRUE)
    expect_equal(fl$Abundance.WT.avg,
                 rowMeans(assay(tmpsce, "Abundance")[, c("WT_S13", "WT_S14", "WT_S15", "WT_S16")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$Abundance.HIS4KO.avg,
                 rowMeans(fl[, c("Abundance.HIS4KO_S05", "Abundance.HIS4KO_S06",
                                 "Abundance.HIS4KO_S07", "Abundance.HIS4KO_S08")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$Abundance.WT.sd,
                 sqrt(matrixStats::rowVars(assay(tmpsce, "Abundance")[, c("WT_S13", "WT_S14", "WT_S15", "WT_S16")],
                                           na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_Abundance.WT.avg,
                 rowMeans(log2(assay(tmpsce, "Abundance"))[, c("WT_S13", "WT_S14", "WT_S15", "WT_S16")],
                          na.rm = TRUE), tolerance = 1e-5, ignore_attr = TRUE)
    expect_equal(fl$log2_Abundance.WT.sd,
                 sqrt(matrixStats::rowVars(log2(assay(tmpsce, "Abundance"))[, c("WT_S13", "WT_S14", "WT_S15", "WT_S16")],
                                           na.rm = TRUE)), tolerance = 1e-5, ignore_attr = TRUE)

})
