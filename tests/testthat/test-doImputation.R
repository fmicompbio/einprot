test_that("minProbGlobalFun works", {
    set.seed(123L)
    mat <- matrix(runif(70, min = 10, max = 15), nrow = 10)
    mat0 <- mat
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
    out <- .minProbGlobalFun(mat0, lowQuantile = 0.01, multSigma = 1)
    expect_identical(out, mat0)

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

test_that("minProbGlobalOffsetFun works", {
    set.seed(123L)
    mat <- matrix(runif(70, min = 10, max = 15), nrow = 10)
    mat0 <- mat
    missingpos <- cbind(row = c(1, 1, 1, 4, 6, 9), col = c(2, 3, 7, 1, 6, 7))
    mat[missingpos] <- NA
    expect_equal(sum(is.na(mat)), nrow(missingpos))
    for (i in seq_len(nrow(missingpos))) {
        expect_equal(mat[missingpos[i, 1], missingpos[i, 2]], NA_real_)
    }

    expect_error(.minProbGlobalOffsetFun(mat = 1, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = NULL),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbGlobalOffsetFun(mat = data.frame(mat), lowQuantile = 0.01,
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = NULL),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 1,
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = NULL),
                 "'lowQuantile' must be within")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = "0.5",
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = NULL),
                 "'lowQuantile' must be of class 'numeric'")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = c(0.2, 0.3),
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = NULL),
                 "'lowQuantile' must have length 1")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = "1", constOffset = NULL,
                                         relOffset = NULL),
                 "'multSigma' must be of class 'numeric'")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = c(1, 2), constOffset = NULL,
                                         relOffset = NULL),
                 "'multSigma' must have length 1")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = "1",
                                         relOffset = NULL),
                 "'constOffset' must be of class 'numeric'")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = c(1, 2),
                                         relOffset = NULL),
                 "'constOffset' must have length 1")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = "1"),
                 "'relOffset' must be of class 'numeric'")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = NULL,
                                         relOffset = c(1, 2)),
                 "'relOffset' must have length 1")
    expect_error(.minProbGlobalOffsetFun(mat = mat, lowQuantile = 0.01,
                                         multSigma = 1, constOffset = -0.1,
                                         relOffset = -0.5),
                 "At most one of constOffset and relOffset can be specified")

    set.seed(1L)
    out <- .minProbGlobalOffsetFun(mat0, lowQuantile = 0.01, multSigma = 1,
                                   constOffset = -0.5, relOffset = NULL)
    expect_identical(out, mat0)

    set.seed(1L)
    out <- .minProbGlobalOffsetFun(mat, lowQuantile = 0.01, multSigma = 1,
                                   constOffset = NULL, relOffset = NULL)
    trueq <- stats::quantile(mat[!is.na(mat)], 0.01)
    expect_equal(mat[!is.na(mat)], out[!is.na(mat)])
    imputed_values <- out[missingpos]
    expect_lte(abs(mean(imputed_values) - trueq), 0.06)

    set.seed(1L)
    out <- .minProbGlobalOffsetFun(mat, lowQuantile = 0.01, multSigma = 1,
                                   constOffset = -0.4, relOffset = NULL)
    trueq <- stats::quantile(mat[!is.na(mat)], 0.01) - 0.4
    expect_equal(mat[!is.na(mat)], out[!is.na(mat)])
    imputed_values <- out[missingpos]
    expect_lte(abs(mean(imputed_values) - trueq), 0.06)

    set.seed(1L)
    out <- .minProbGlobalOffsetFun(mat, lowQuantile = 0.2, multSigma = 1,
                                   constOffset = NULL, relOffset = -0.4)
    trueq <- stats::quantile(mat[!is.na(mat)], 0.2) - 0.4 * 1.552274
    expect_equal(mat[!is.na(mat)], out[!is.na(mat)])
    imputed_values <- out[missingpos]
    expect_lte(abs(mean(imputed_values) - trueq), 0.06)
})

test_that("minProbOffsetFun works", {
    set.seed(123L)
    mat <- matrix(runif(70, min = 10, max = 15), nrow = 10)
    mat0 <- mat
    missingpos <- cbind(row = c(1, 1, 1, 4, 6, 9), col = c(2, 3, 7, 1, 6, 7))
    mat[missingpos] <- NA
    expect_equal(sum(is.na(mat)), nrow(missingpos))
    for (i in seq_len(nrow(missingpos))) {
        expect_equal(mat[missingpos[i, 1], missingpos[i, 2]], NA_real_)
    }

    expect_error(.minProbOffsetFun(mat = 1, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbOffsetFun(mat = data.frame(mat), lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'mat' must be of class 'matrix'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 1,
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'lowQuantile' must be within")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = "0.5",
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'lowQuantile' must be of class 'numeric'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = c(0.2, 0.3),
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'lowQuantile' must have length 1")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = "1", margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'multSigma' must be of class 'numeric'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = c(1, 2), margin = 2L, constOffset = NULL,
                                   relOffset = NULL),
                 "'multSigma' must have length 1")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = "2", constOffset = NULL,
                                   relOffset = NULL),
                 "'margin' must be of class 'numeric'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 3L, constOffset = NULL,
                                   relOffset = NULL),
                 "'margin' must be one of")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = "1",
                                   relOffset = NULL),
                 "'constOffset' must be of class 'numeric'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = c(1, 2),
                                   relOffset = NULL),
                 "'constOffset' must have length 1")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = "1"),
                 "'relOffset' must be of class 'numeric'")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = NULL,
                                   relOffset = c(1, 2)),
                 "'relOffset' must have length 1")
    expect_error(.minProbOffsetFun(mat = mat, lowQuantile = 0.01,
                                   multSigma = 1, margin = 2L, constOffset = -0.1,
                                   relOffset = -0.5),
                 "At most one of constOffset and relOffset can be specified")

    set.seed(1L)
    out <- .minProbOffsetFun(mat0, lowQuantile = 0.01, multSigma = 1,
                             margin = 2L, constOffset = NULL, relOffset = 0.2)
    expect_identical(out, mat0)

    set.seed(1L)
    out <- .minProbOffsetFun(mat, lowQuantile = 0.01, multSigma = 1,
                             margin = 2L, constOffset = NULL, relOffset = NULL)
    set.seed(1L)
    out2 <- MsCoreUtils::impute_MinProb(mat, q = 0.01, sigma = 1, MARGIN = 2L)
    expect_equal(out, out2)

    set.seed(1L)
    out <- .minProbOffsetFun(mat, lowQuantile = 0.01, multSigma = 1,
                             margin = 1L, constOffset = NULL, relOffset = NULL)
    set.seed(1L)
    out2 <- MsCoreUtils::impute_MinProb(mat, q = 0.01, sigma = 1, MARGIN = 1L)
    expect_equal(out, out2)
})

test_that("imputation works", {
    expect_error(doImputation(sce = 1, imputeMethod = "MinProb", assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = 1,
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'imputeMethod' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation,
                              imputeMethod = c("impSeqRob", "MinProb"),
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "'imputeMethod' must have length 1")
    expect_error(doImputation(sce = sce_mq_preimputation,
                              imputeMethod = c("missing"),
                              assayName = "iBAQ",
                              imputedAssayName = "imputed"),
                 "All values in 'imputeMethod' must be one of")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
                              assayName = 1,
                              imputedAssayName = "imputed"),
                 "'assayName' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
                              assayName = c("iBAQ", "iBAQ"),
                              imputedAssayName = "imputed"),
                 "'assayName' must have length 1")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
                              assayName = "missing",
                              imputedAssayName = "imputed"),
                 "All values in 'assayName' must be one of")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
                              assayName = "iBAQ",
                              imputedAssayName = 1),
                 "'imputedAssayName' must be of class 'character'")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
                              assayName = "iBAQ",
                              imputedAssayName = c("imputed", "imputed")),
                 "'imputedAssayName' must have length 1")
    expect_error(doImputation(sce = sce_mq_preimputation, imputeMethod = "custom",
                              assayName = "iBAQ",
                              imputedAssayName = "imputed", FUN = "1"),
                 "Please specify an imputation method")

    expect_true(sum(is.na(SummarizedExperiment::assay(sce_mq_preimputation,
                                                      "iBAQ"))) == 507)

    ## MinProb
    set.seed(42L)
    impout <- doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProb",
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
    set.seed(42L)
    impout2 <- doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProbOffset",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            lowQuantile = 0.01, multSigma = 1,
                            margin = 2L, constOffset = NULL, relOffset = NULL)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout2, "imputedAssay"))
    # with constant offset
    set.seed(42L)
    impout2 <- doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProbOffset",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            lowQuantile = 0.01, multSigma = 1,
                            margin = 2L, constOffset = -0.5, relOffset = NULL)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout2, "imputedAssay") +
                     0.5 * is.na(SummarizedExperiment::assay(impout2, "log2_iBAQ_withNA")))
    # with relative offset
    set.seed(42L)
    impout2 <- doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProbOffset",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            lowQuantile = 0.01, multSigma = 1,
                            margin = 2L, constOffset = NULL, relOffset = -1)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout2, "imputedAssay") +
                     1.159449 * is.na(SummarizedExperiment::assay(impout2, "log2_iBAQ_withNA")),
                 tolerance = 1e-5)

    ## MinProb, by row
    # ... first make sure that there are no features with all NA values (removes 5 rows)
    tmp <- sce_mq_preimputation[rowSums(!is.na(SummarizedExperiment::assay(sce_mq_preimputation))) > 0, ]
    set.seed(42L)
    impout <- doImputation(sce = tmp, imputeMethod = "MinProb",
                           assayName = "log2_iBAQ_withNA", MARGIN = 1L,
                           imputedAssayName = "imputedAssay", sigma = 0.25)
    expect_s4_class(impout, "SummarizedExperiment")
    expect_true("imputedAssay" %in% SummarizedExperiment::assayNames(impout))
    expect_true(sum(is.na(SummarizedExperiment::assay(impout,
                                                      "log2_iBAQ_withNA"))) == 507 - 5 * 9)
    expect_true(all(!is.na(SummarizedExperiment::assay(impout,
                                                       "imputedAssay"))))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"),
                 SummarizedExperiment::assay(tmp, "log2_iBAQ_withNA"))
    expect_equal(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[!is.na(
        SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))],
        SummarizedExperiment::assay(impout, "imputedAssay")[!is.na(
            SummarizedExperiment::assay(impout, "log2_iBAQ_withNA"))])
    for (i in seq_len(145)) {
        imputed_values <- SummarizedExperiment::assay(impout, "imputedAssay")[i, is.na(
            SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[i, ])]
        if (length(imputed_values) > 0) {
            observed_values <- SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[i, !is.na(
                SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")[i, ])]
            ## Relative error should be less than 10%
            expect_lte(abs(mean(imputed_values) - stats::quantile(observed_values, 0.01)) /
                           stats::quantile(observed_values, 0.01), 0.1)
        }
    }
    set.seed(42L)
    impout2 <- doImputation(sce = tmp, imputeMethod = "MinProbOffset",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            lowQuantile = 0.01, multSigma = 0.25,
                            margin = 1L, constOffset = NULL, relOffset = NULL)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout2, "imputedAssay"))
    set.seed(42L)
    impout3 <- doImputation(sce = tmp, imputeMethod = "custom",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            FUN = einprot:::.minProbOffsetFun,
                            lowQuantile = 0.01, multSigma = 0.25,
                            margin = 1L, constOffset = NULL, relOffset = NULL)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout3, "imputedAssay"))

    ## MinProbGlobal
    set.seed(42L)
    impout <- doImputation(sce = sce_mq_preimputation, imputeMethod = "MinProbGlobal",
                           assayName = "log2_iBAQ_withNA",
                           imputedAssayName = "imputedAssay", lowQuantile = 0.05)
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
    expect_lte(abs(mean(imputed_values) - stats::quantile(observed_values, 0.05)) /
                   stats::quantile(observed_values, 0.05), 0.01)
    set.seed(42L)
    impout2 <- doImputation(sce = sce_mq_preimputation, imputeMethod = "custom",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            FUN = einprot:::.minProbGlobalFun, lowQuantile = 0.05)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay"),
                 SummarizedExperiment::assay(impout2, "imputedAssay"))

    ## MinProbGlobalOffset
    set.seed(42L)
    impout <- doImputation(sce = sce_mq_preimputation,
                           imputeMethod = "MinProbGlobalOffset",
                           assayName = "log2_iBAQ_withNA",
                           imputedAssayName = "imputedAssay", constOffset = -0.5)
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
                   stats::quantile(observed_values, 0.01) - 0.5, 0.01)
    set.seed(42L)
    impout2 <- doImputation(sce = sce_mq_preimputation, imputeMethod = "custom",
                            assayName = "log2_iBAQ_withNA",
                            imputedAssayName = "imputedAssay",
                            FUN = einprot:::.minProbGlobalOffsetFun,
                            lowQuantile = 0.01, constOffset = -0.7)
    expect_equal(SummarizedExperiment::assay(impout, "imputedAssay") -
                     0.2 * is.na(SummarizedExperiment::assay(impout, "log2_iBAQ_withNA")),
                 SummarizedExperiment::assay(impout2, "imputedAssay"))

    ## impSeqRob
    impout <- doImputation(sce = sce_mq_preimputation, imputeMethod = "impSeqRob",
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
