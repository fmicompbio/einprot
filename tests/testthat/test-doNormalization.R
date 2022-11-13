test_that("normalization works", {
    expect_error(doNormalization(sce = 1, method = "center.median",
                                 assayName = "iBAQ",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(doNormalization(sce = sce_mq_preimputation, method = 1,
                                 assayName = "iBAQ",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "'method' must be of class 'character'")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = c("center.median", "div.mean"),
                                 assayName = "iBAQ",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "'method' must have length 1")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = c("missing"),
                                 assayName = "iBAQ",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "All values in 'method' must be one of")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = 1,
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "'assayName' must be of class 'character'")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = c("iBAQ", "iBAQ"),
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "'assayName' must have length 1")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = "missing",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = NULL),
                 "All values in 'assayName' must be one of")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = "iBAQ",
                                 normalizedAssayName = 1,
                                 spikeFeatures = NULL),
                 "'normalizedAssayName' must be of class 'character'")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = "iBAQ",
                                 normalizedAssayName = c("normalized", "normalized"),
                                 spikeFeatures = NULL),
                 "'normalizedAssayName' must have length 1")
    expect_error(doNormalization(sce = sce_mq_preimputation,
                                 method = "center.median",
                                 assayName = "iBAQ",
                                 normalizedAssayName = "normalized",
                                 spikeFeatures = 1),
                 "'spikeFeatures' must be of class 'character'")

    expect_true(sum(is.na(SummarizedExperiment::assay(sce_mq_preimputation,
                                                      "iBAQ"))) == 507)

    ## center.median
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "center.median",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = NULL)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    expect_equal(apply(assay(normout, "normalizedAssay"), 2, stats::median,
                       na.rm = TRUE), rep(0L, ncol(normout)),
                 ignore_attr = TRUE)

    ## center.mean
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "center.mean",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = NULL)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    expect_equal(apply(assay(normout, "normalizedAssay"), 2, mean,
                       na.rm = TRUE), rep(0L, ncol(normout)),
                 ignore_attr = TRUE)

    ## div.mean
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "div.mean",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = NULL)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    expect_equal(apply(assay(normout, "normalizedAssay"), 2, mean,
                       na.rm = TRUE), rep(1L, ncol(normout)),
                 ignore_attr = TRUE)

    ## spike-in features, center.mean
    spikes <- rownames(sce_mq_preimputation)[c(8, 90, 123)]
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "center.mean",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = spikes)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    ## mean across spikes should be the same in all samples
    spikeMeans <- apply(assay(normout, "normalizedAssay")[spikes, ], 2, mean,
                        na.rm = TRUE)
    expect_true(max(spikeMeans) - min(spikeMeans) < 1e-8)

    ## spike-in features, center.median
    spikes <- rownames(sce_mq_preimputation)[c(8, 90, 123)]
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "center.median",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = spikes)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    ## median across spikes should be the same in all samples
    spikeMedians <- apply(assay(normout, "normalizedAssay")[spikes, ], 2,
                          stats::median, na.rm = TRUE)
    expect_true(max(spikeMedians) - min(spikeMedians) < 1e-8)

    ## spike-in features, div.mean
    spikes <- rownames(sce_mq_preimputation)[c(8, 90, 123)]
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "div.mean",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = spikes)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    ## mean across spikes should be the same in all samples
    spikeMeans <- apply(assay(normout, "normalizedAssay")[spikes, ], 2, mean,
                        na.rm = TRUE)
    expect_true(max(spikeMeans) - min(spikeMeans) < 1e-8)

    ## spike-in features, div.median
    spikes <- rownames(sce_mq_preimputation)[c(8, 90, 123)]
    normout <- doNormalization(sce = sce_mq_preimputation,
                               method = "div.median",
                               assayName = "iBAQ",
                               normalizedAssayName = "normalizedAssay",
                               spikeFeatures = spikes)
    expect_s4_class(normout, "SummarizedExperiment")
    expect_true("normalizedAssay" %in% SummarizedExperiment::assayNames(normout))
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "iBAQ"))) == 507)
    expect_true(sum(is.na(SummarizedExperiment::assay(normout,
                                                      "normalizedAssay"))) == 507)
    ## mean across spikes should be the same in all samples
    spikeMedians <- apply(assay(normout, "normalizedAssay")[spikes, ], 2,
                          stats::median, na.rm = TRUE)
    expect_true(max(spikeMedians) - min(spikeMedians) < 1e-8)

    ## spike-in features, unsupported method
    expect_error(
        doNormalization(sce = sce_mq_preimputation,
                        method = "quantiles",
                        assayName = "iBAQ",
                        normalizedAssayName = "normalizedAssay",
                        spikeFeatures = spikes),
        "Unsupported normalization method with spike features"
    )
})
