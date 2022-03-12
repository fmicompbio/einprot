test_that("imputation works", {
    expect_error(doImputation(sce = 1, method = "MinProb", assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(doImputation(sce = sce_mq_preimputation, method = 1,
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'method' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation,
                              method = c("impSeqRob", "MinProb"),
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'method' must have length 1")
    expect_error(doImputation(sce = sce_mq_preimputation,
                              method = c("missing"),
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "All values in 'method' must be one of")
    expect_error(doImputation(sce = sce_mq_preimputation, method = "MinProb",
                              assayName = 1,
                              imputedAssayName = "imputed"),
                 "'assayName' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation, method = "MinProb",
                              assayName = c("iBAQ", "iBAQ"),
                              imputedAssayName = "imputed"),
                 "'assayName' must have length 1")
    expect_error(doImputation(sce = sce_mq_preimputation, method = "MinProb",
                              assayName = "missing",
                              imputedAssayName = "imputed"),
                 "All values in 'assayName' must be one of")
    expect_error(doImputation(sce = sce_mq_preimputation, method = "MinProb",
                              assayName = "iBAQ",
                              imputedAssayName = 1),
                 "'imputedAssayName' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation, method = "MinProb",
                              assayName = "iBAQ",
                              imputedAssayName = c("imputed", "imputed")),
                 "'imputedAssayName' must have length 1")

    impout <- doImputation(sce = sce_mq_preimputation, method = "MinProb",
                           assayName = "iBAQ",
                           imputedAssayName = "imputedAssay")
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay" %in% SummarizedExperiment::assayNames(impout))

    impout <- doImputation(sce = sce_mq_preimputation, method = "impSeqRob",
                           assayName = "iBAQ",
                           imputedAssayName = "imputedAssay2")
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay2" %in% SummarizedExperiment::assayNames(impout))

})
