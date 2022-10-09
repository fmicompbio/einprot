test_that("defineAssayNames works", {
    expect_error(defineAssayNames(aName = 1, normMethod = "none",
                                  doBatchCorr = TRUE),
                 "'aName' must be of class 'character'")
    expect_error(defineAssayNames(aName = c("a", "b"), normMethod = "none",
                                  doBatchCorr = TRUE),
                 "'aName' must have length 1")
    expect_error(defineAssayNames(aName = "iBAQ", normMethod = 1,
                                  doBatchCorr = TRUE),
                 "'normMethod' must be of class 'character'")
    expect_error(defineAssayNames(aName = "iBAQ", normMethod = c("m1", "m2"),
                                  doBatchCorr = TRUE),
                 "'normMethod' must have length 1")
    expect_error(defineAssayNames(aName = "iBAQ", normMethod = "none",
                                  doBatchCorr = 1),
                 "'doBatchCorr' must be of class 'logical'")
    expect_error(defineAssayNames(aName = "iBAQ", normMethod = "none",
                                  doBatchCorr = c(TRUE, FALSE)),
                 "'doBatchCorr' must have length 1")

    nms <- defineAssayNames(aName = "iBAQ", normMethod = "none",
                            doBatchCorr = FALSE)
    expect_type(nms, "list")
    expect_equal(nms$assayInput, "iBAQ")
    expect_equal(nms$assayLog2WithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputIndic, "imputed_iBAQ")
    expect_equal(nms$assayLog2NormWithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputed, "log2_iBAQ")
    expect_equal(nms$assayBatchCorr, "log2_iBAQ")

    nms <- defineAssayNames(aName = "iBAQ", normMethod = "center.median",
                            doBatchCorr = FALSE)
    expect_type(nms, "list")
    expect_equal(nms$assayInput, "iBAQ")
    expect_equal(nms$assayLog2WithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputIndic, "imputed_iBAQ")
    expect_equal(nms$assayLog2NormWithNA, "log2_iBAQ_withNA_norm")
    expect_equal(nms$assayImputed, "log2_iBAQ_norm")
    expect_equal(nms$assayBatchCorr, "log2_iBAQ_norm")

    nms <- defineAssayNames(aName = "iBAQ", normMethod = "none",
                            doBatchCorr = TRUE)
    expect_type(nms, "list")
    expect_equal(nms$assayInput, "iBAQ")
    expect_equal(nms$assayLog2WithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputIndic, "imputed_iBAQ")
    expect_equal(nms$assayLog2NormWithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputed, "log2_iBAQ")
    expect_equal(nms$assayBatchCorr, "log2_iBAQ_batchCorr")

    nms <- defineAssayNames(aName = "iBAQ", normMethod = "center.median",
                            doBatchCorr = TRUE)
    expect_type(nms, "list")
    expect_equal(nms$assayInput, "iBAQ")
    expect_equal(nms$assayLog2WithNA, "log2_iBAQ_withNA")
    expect_equal(nms$assayImputIndic, "imputed_iBAQ")
    expect_equal(nms$assayLog2NormWithNA, "log2_iBAQ_withNA_norm")
    expect_equal(nms$assayImputed, "log2_iBAQ_norm")
    expect_equal(nms$assayBatchCorr, "log2_iBAQ_norm_batchCorr")
})
