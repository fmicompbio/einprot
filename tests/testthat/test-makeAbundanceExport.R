test_that("makeAbundanceExport works", {
    expect_error(makeAbundanceExport(testresList = 1, abundancePrefix = "pref"),
                 "'testresList' must be of class 'list'")
    expect_error(makeAbundanceExport(testresList = list(), abundancePrefix = 1),
                 "'abundancePrefix' must be of class 'character'")
    expect_error(makeAbundanceExport(testresList = list(),
                                     abundancePrefix = c("pref", "pref2")),
                 "'abundancePrefix' must have length 1")

    out_limma <- runTest(
        sce = sce_mq_final,
        comparisons = list(c("Adnp", "RBC_ctrl"), c("RBC_ctrl", "Chd4BF")),
        testType = "limma",
        assayForTests = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2, minlFC = 0, featureCollections = fcoll_mq_final,
        complexFDRThr = 0.8, volcanoAdjPvalThr = 0.05, volcanoLog2FCThr = 1,
        baseFileName = NULL, seed = 123, nperm = 25, volcanoS0 = 0.1,
        addAbundanceValues = TRUE, aName = "iBAQ", singleFit = FALSE
    )
    abExp <- makeAbundanceExport(testresList = out_limma$tests,
                                 abundancePrefix = "iBAQ")
    expect_s3_class(abExp, "data.frame")
    expect_equal(nrow(abExp), sum(out_limma$tests$RBC_ctrl_vs_Adnp$showInVolcano |
                                      out_limma$tests$Chd4BF_vs_RBC_ctrl$showInVolcano,
                                  na.rm = TRUE))
    expect_equal(ncol(abExp), 18L)
    expect_equal(abExp[, paste0("iBAQ.", colnames(sce_mq_final))],
                 as.data.frame(SummarizedExperiment::assay(sce_mq_final, "iBAQ")[abExp$pid, ]),
                 ignore_attr = TRUE)
    expect_equal(sort(abExp$pid[abExp$showInVolcano_Chd4BF_vs_RBC_ctrl]),
                 sort(rownames(out_limma$tests$Chd4BF_vs_RBC_ctrl)[out_limma$tests$Chd4BF_vs_RBC_ctrl$showInVolcano]))
    expect_equal(abExp$iBAQ.Adnp.avg,
                 rowMeans(abExp[, c(grep("iBAQ.Adnp_", colnames(abExp)))],
                          na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(abExp$iBAQ.Adnp.sd,
                 rowSds(abExp[, c(grep("iBAQ.Adnp_", colnames(abExp)))],
                        na.rm = TRUE),
                 ignore_attr = TRUE)

    ## Only one comparison
    ## -------------------------------------------------------------------------
    abExp <- makeAbundanceExport(testresList = out_limma$tests["RBC_ctrl_vs_Adnp"],
                                 abundancePrefix = "iBAQ")
    expect_s3_class(abExp, "data.frame")
    expect_equal(nrow(abExp), sum(out_limma$tests$RBC_ctrl_vs_Adnp$showInVolcano,
                                  na.rm = TRUE))
    expect_equal(ncol(abExp), 12L)
    expect_equal(abExp[, paste0("iBAQ.", c("Adnp_IP04", "Adnp_IP05", "RBC_ctrl_IP02"))],
                 as.data.frame(SummarizedExperiment::assay(sce_mq_final, "iBAQ")[abExp$pid, c("Adnp_IP04", "Adnp_IP05", "RBC_ctrl_IP02")]),
                 ignore_attr = TRUE)
    expect_equal(sort(abExp$pid[abExp$showInVolcano_RBC_ctrl_vs_Adnp]),
                 sort(rownames(out_limma$tests$RBC_ctrl_vs_Adnp)[out_limma$tests$RBC_ctrl_vs_Adnp$showInVolcano]))
    expect_equal(abExp$iBAQ.Adnp.avg,
                 rowMeans(abExp[, c(grep("iBAQ.Adnp_", colnames(abExp)))],
                          na.rm = TRUE),
                 ignore_attr = TRUE)
    expect_equal(abExp$iBAQ.Adnp.sd,
                 rowSds(abExp[, c(grep("iBAQ.Adnp_", colnames(abExp)))],
                        na.rm = TRUE),
                 ignore_attr = TRUE)
})
