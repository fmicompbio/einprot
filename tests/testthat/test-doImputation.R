test_that("minProbGlobalFun works", {
    set.seed(123L)
    mat <- matrix(runif(70, min = 10, max = 15), nrow = 10)
    missingpos <- cbind(row = c(1, 1, 1, 4, 6, 9), col = c(2, 3, 7, 1, 6, 7))
    mat[missingpos] <- NA
    expect_equal(sum(is.na(mat)), nrow(missingpos))
    for (i in seq_len(nrow(missingpos))) {
        expect_equal(mat[missingpos[i, 1], missingpos[i, 2]], NA_real_)
    }

    expect_error(.minProbGlobalFun(mat = 1, lowQuantile = 0.01,
                                   multSigma = 1),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbGlobalFun(mat = data.frame(mat), lowQuantile = 0.01,
                                   multSigma = 1),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbGlobalFun(mat = mat, lowQuantile = 1,
                                   multSigma = 1),
                 "'lowQuantile' must be within")
    expect_error(.minProbGlobalFun(mat = mat, lowQuantile = "0.5",
                                   multSigma = 1),
                 "'lowQuantile' must be of class 'numeric'")
    expect_error(.minProbGlobalFun(mat = mat, lowQuantile = c(0.2, 0.3),
                                   multSigma = 1),
                 "'lowQuantile' must have length 1")
    expect_error(.minProbGlobalFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = "1"),
                 "'multSigma' must be of class 'numeric'")
    expect_error(.minProbGlobalFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = c(1, 2)),
                 "'multSigma' must have length 1")

    set.seed(1L)
    out <- .minProbGlobalFun(mat, lowQuantile = 0.01, multSigma = 1)
    trueq <- stats::quantile(mat[!is.na(mat)], 0.01)
    expect_equal(mat[!is.na(mat)], out[!is.na(mat)])
    imputed_values <- out[missingpos]
    expect_lte(abs(mean(imputed_values) - trueq), 0.06)

    set.seed(1L)
    out <- .minProbGlobalFun(mat, lowQuantile = 0.8, multSigma = 1)
    trueq <- stats::quantile(mat[!is.na(mat)], 0.8)
    expect_equal(mat[!is.na(mat)], out[!is.na(mat)])
    imputed_values <- out[missingpos]
    expect_lte(abs(mean(imputed_values) - trueq), 0.06)
})

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

    expect_true(sum(is.na(SummarizedExperiment::assay(sce_mq_preimputation,
                                                      "iBAQ"))) == 507)

    ## MinProb
    set.seed(42L)
    impout <- doImputation(sce = sce_mq_preimputation, method = "MinProb",
                           assayName = "log2_iBAQ_withNA",
                           imputedAssayName = "imputedAssay")
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay" %in% SummarizedExperiment::assayNames(impout))
    expect_true(sum(is.na(SummarizedExperiment::assay(impout,
                                                      "log2_iBAQ_withNA"))) == 507)
    expect_true(all(!is.na(SummarizedExperiment::assay(impout,
                                                       "imputedAssay"))))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"),
                 SummarizedExperiment::assay(sce_mq_preimputation, "log2_iBAQ_withNA"))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[!is.na(
        SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))],
                 SummarizedExperiment::assay(impout, "imputedAssay")[!is.na(
                     SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))])
    for (i in seq_len(9)) {
        imputed_values <- SummarizedExperiment::assay(impout, "imputedAssay")[is.na(
            SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[, i]), i]
        observed_values <- SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[!is.na(
            SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[, i]), i]
        ## Relative error should be less than 2.5%
        expect_lte(abs(mean(imputed_values) - stats::quantile(observed_values, 0.01)) /
                       stats::quantile(observed_values, 0.01), 0.025)
    }

    ## MinProbGlobal
    set.seed(42L)
    impout <- doImputation(sce = sce_mq_preimputation, method = "MinProbGlobal",
                           assayName = "log2_iBAQ_withNA",
                           imputedAssayName = "imputedAssay")
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay" %in% SummarizedExperiment::assayNames(impout))
    expect_true(sum(is.na(SummarizedExperiment::assay(impout,
                                                      "log2_iBAQ_withNA"))) == 507)
    expect_true(all(!is.na(SummarizedExperiment::assay(impout,
                                                       "imputedAssay"))))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"),
                 SummarizedExperiment::assay(sce_mq_preimputation, "log2_iBAQ_withNA"))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[!is.na(
        SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))],
        SummarizedExperiment::assay(impout, "imputedAssay")[!is.na(
            SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))])
    imputed_values <- SummarizedExperiment::assay(impout, "imputedAssay")[is.na(
        SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))]
    observed_values <- SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[!is.na(
        SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))]
    expect_lte(abs(mean(imputed_values) - stats::quantile(observed_values, 0.01)) /
                   stats::quantile(observed_values, 0.01), 0.01)

    ## impSeqRob
    impout <- doImputation(sce = sce_mq_preimputation, method = "impSeqRob",
                           assayName = "iBAQ",
                           imputedAssayName = "imputedAssay2")
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay2" %in% SummarizedExperiment::assayNames(impout))
    expect_true(sum(is.na(SummarizedExperiment::assay(impout,
                                                      "iBAQ"))) == 507)
    expect_true(all(!is.na(SummarizedExperiment::assay(impout,
                                                       "imputedAssay2"))))
    expect_equal(SummarizedExperiment::assay(impout, "iBAQ"),
                 SummarizedExperiment::assay(sce_mq_preimputation, "iBAQ"))
    expect_equal(SummarizedExperiment::assay(impout, "iBAQ")[!is.na(
        SummarizedExperiment::assay(impout, "iBAQ"))],
                 SummarizedExperiment::assay(impout, "imputedAssay2")[!is.na(
                     SummarizedExperiment::assay(impout, "iBAQ"))])
})
