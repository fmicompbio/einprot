test_that("testing works", {

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    args0 <- list(
        sce = sce_mq_final,
        comparison = c("Adnp", "RBC_ctrl"),
        testType = "limma",
        assayForTests = "log2_iBAQ",
        assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2,
        minlFC = 0,
        featureCollections = fcoll_mq_final,
        complexFDRThr = 0.1,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        baseFileName = NULL,
        seed = 123,
        nperm = 25,
        volcanoS0 = 0.1,
        addAbundanceValues = TRUE,
        aName = "iBAQ"
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
    args$comparison <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'comparison' must be of class 'character'")
    args$comparison <- "Adnp"
    expect_error(do.call(runTest, args),
                 "'comparison' must have length 2")
    args$comparison <- c("missing", "Adnp")
    expect_error(do.call(runTest, args),
                 "All values in 'comparison' must be one of")
    args$comparison <- list(c("Adnp", "RBC_ctrl"))
    expect_error(do.call(runTest, args),
                 "'comparison' must be of class 'character'")

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

    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    out <- do.call(runTest, args0)
    expect_type(out, "list")
    expect_length(out, 7)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "topSets", "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(args0$sce))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 0)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$res$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$res$pid[1:5], out$res$IDsForSTRING[1:5])

    ## t-test, write results to file
    args <- args0
    args$testType <- "ttest"
    args$baseFileName <- tempfile()
    args$addAbundanceValues <- FALSE
    out1 <- do.call(runTest, args)
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt")))
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp_camera_complexes.txt")))
    expect_type(out1, "list")
    expect_length(out1, 7)
    expect_named(out1, c("res", "plotnote", "plottitle", "plotsubtitle",
                         "topSets", "featureCollections", "curveparam"))
    expect_s3_class(out1$res, "data.frame")
    expect_type(out1$plotnote, "character")
    expect_type(out1$plottitle, "character")
    expect_type(out1$featureCollections, "list")
    expect_type(out1$curveparam, "list")
    expect_named(out1$curveparam, c("x", "ta", "s0", "df"))
    expect_equal(nrow(out1$res), 70)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out1$res)))
    expect_false("iBAQ.Adnp_IP04" %in% colnames(out1$res))
    expect_equal(out1$res$pid, rownames(args$sce))
    expect_equal(out1$plotnote, "")
    expect_equal(out1$plottitle, "RBC_ctrl vs Adnp, t-test")
    expect_s4_class(out1$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out1$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out1$featureCollections$complexes)))
    expect_equal(out1$res$pid[1:5], out1$res$IDsForSTRING[1:5])
    expect_false(any(grepl("iBAQ", colnames(out1$res))))
    expect_equal(out1$res$logFC[1], mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 7:9]) -
                     mean(SummarizedExperiment::assay(args$sce, args$assayForTests)[1, 1:3]))
    ## Should correlate with limma results
    expect_equal(out$res$pid, out1$res$pid)
    idx <- which(!is.na(out$res$logFC))
    expect_equal(out$res$AveExpr, out1$res$AveExpr, ignore_attr = TRUE)
    expect_equal(out$res$logFC, out1$res$logFC, ignore_attr = TRUE)
    expect_gt(cor(out$res$t[idx], out1$res$sam[idx]), 0.9)
    expect_lt(cor(out$res$t[idx], out1$res$sam[idx]), 0.99)

    ## Different assay
    args <- args0
    args$assayForTests <- "Intensity"
    args$aName <- "LFQ.intensity"
    out <- do.call(runTest, args)
    expect_type(out, "list")
    expect_length(out, 7)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "topSets", "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "LFQ.intensity.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(args$sce))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 0)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$res$LFQ.intensity.Adnp_IP04,
                 SummarizedExperiment::assay(args0$sce, "LFQ.intensity")[, "Adnp_IP04"],
                 ignore_attr = TRUE)
    expect_equal(out$res$pid[1:5], out$res$IDsForSTRING[1:5])

    ## With batch effect
    args <- args0
    args$sce$batch <- c("B1", "B2", "B1", "B2", "B1", "B2", "B1", "B2", "B1")
    out <- do.call(runTest, args0)
    expect_type(out, "list")
    expect_length(out, 7)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "topSets", "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(sce_mq_final))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 0)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$res$iBAQ.Adnp_IP04,
                 SummarizedExperiment::assay(args$sce, "iBAQ")[, "Adnp_IP04"],
                 ignore_attr = TRUE)

    ## -----------------------------------------------------------------------
    ## PD data
    args0_pd <- list(
        sce = sce_pd_final,
        comparison = c("HIS4KO", "WT"),
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
        aName = "Abundance"
    )
    out <- do.call(runTest, args0_pd)
    expect_type(out, "list")
    expect_length(out, 7)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "topSets", "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "Abundance.HIS4KO_S07",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(sce_pd_final))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.5)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("WT_vs_HIS4KO_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))
    expect_equal(out$res$Abundance.HIS4KO_S05,
                 SummarizedExperiment::assay(args0_pd$sce, "Abundance")[, "HIS4KO_S05"],
                 ignore_attr = TRUE)

    ## With batch column
    tmp <- sce_pd_final
    tmp$batch <- rep(c("b1", "b2"), 8)
    args <- args0_pd
    args$sce <- tmp
    outb <- do.call(runTest, args)
    expect_type(outb, "list")
    expect_length(outb, 7)
    expect_named(outb, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "topSets", "featureCollections", "curveparam"))
    expect_s3_class(outb$res, "data.frame")
    expect_type(outb$plotnote, "character")
    expect_type(outb$plottitle, "character")
    expect_type(outb$featureCollections, "list")
    expect_type(outb$curveparam, "list")
    expect_equal(nrow(outb$res), 70)
    expect_true(all(c("adj.P.Val", "Abundance.HIS4KO_S07",
                      "showInVolcano", "IDsForSTRING") %in% colnames(outb$res)))
    expect_equal(outb$res$pid, rownames(sce_pd_final))
    expect_equal(substr(outb$plotnote, 1, 8), "df.prior")
    expect_equal(outb$plottitle, "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.5)")
    expect_s4_class(outb$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(outb$featureCollections$complexes), "DFrame")
    expect_true("WT_vs_HIS4KO_FDR" %in%
                    colnames(S4Vectors::mcols(outb$featureCollections$complexes)))
    expect_equal(outb$res$Abundance.HIS4KO_S05,
                 SummarizedExperiment::assay(args0_pd$sce, "Abundance")[, "HIS4KO_S05"],
                 ignore_attr = TRUE)
    ## Check that test results are different compared to without batch
    expect_equal(out$res$pid, outb$res$pid)
    expect_equal(which(!is.na(out$res$logFC)), which(!is.na(outb$res$logFC)))
    nonaidx <- which(!is.na(outb$res$logFC))
    expect_false(all(round(outb$res$P.Value[nonaidx], 13) ==
                         round(out$res$P.Value[nonaidx], 13)))
    expect_true(all(round(outb$res$logFC[nonaidx], 13) ==
                        round(out$res$logFC[nonaidx], 13)))
    expect_true(all(round(outb$res$AveExpr[nonaidx], 13) ==
                        round(out$res$AveExpr[nonaidx], 13)))
    expect_true(all(round(outb$res$Abundance.HIS4KO_S05[nonaidx], 13) ==
                        round(out$res$Abundance.HIS4KO_S05[nonaidx], 13)))
    expect_true(all(round(outb$res$log2_Abundance.WT.sd[nonaidx], 13) ==
                        round(out$res$log2_Abundance.WT.sd[nonaidx], 13)))
})
