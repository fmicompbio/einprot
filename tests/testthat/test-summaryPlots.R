test_that("generating summary plots works", {

    ## Fails with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## makeIntensityBoxplots
    expect_error(makeIntensityBoxplots(sce = 1, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = "lab"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(makeIntensityBoxplots(sce = sce_mq_initial,
                                       assayName = "iBAQ", doLog = TRUE,
                                       ylab = "ylab"),
                 'all(c("sample", "group") %in% colnames(', fixed = TRUE)
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = 1,
                                       doLog = TRUE, ylab = "lab"),
                 "'assayName' must be of class 'character'")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final,
                                       assayName = c("iBAQ", "log2_iBAQ"),
                                       doLog = TRUE, ylab = "lab"),
                 "'assayName' must have length 1")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = "missing",
                                       doLog = TRUE, ylab = "lab"),
                 "All values in 'assayName' must be one of")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = "log2_iBAQ",
                                       doLog = 1, ylab = "lab"),
                 "'doLog' must be of class 'logical'")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = "log2_iBAQ",
                                       doLog = c(TRUE, FALSE), ylab = "lab"),
                 "'doLog' must have length 1")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = 1),
                 "'ylab' must be of class 'character'")
    expect_error(makeIntensityBoxplots(sce = sce_mq_final, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = c("lab1", "lab2")),
                 "'ylab' must have length 1")

    ## makeMeanSDPlot
    expect_error(makeMeanSDPlot(sce = 1, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = "ylab"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = 1,
                                xlab = "xlab", ylab = "ylab"),
                 "'assayName' must be of class 'character'")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = c("iBAQ", "log2_iBAQ"),
                                xlab = "xlab", ylab = "ylab"),
                 "'assayName' must have length 1")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = "missing",
                                xlab = "xlab", ylab = "ylab"),
                 "All values in 'assayName' must be one of")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = 1),
                 "'ylab' must be of class 'character'")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = c("ylab1", "ylab2")),
                 "'ylab' must have length 1")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = "log2_iBAQ",
                                xlab = 1, ylab = "ylab"),
                 "'xlab' must be of class 'character'")
    expect_error(makeMeanSDPlot(sce = sce_mq_final, assayName = "log2_iBAQ",
                                xlab = c("xlab1", "xlab2"), ylab = "ylab"),
                 "'xlab' must have length 1")

    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    out <- makeIntensityBoxplots(sce_mq_final, assayName = "log2_iBAQ",
                                 doLog = FALSE, ylab = "lab")
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("col_id", "intensity", "sample", "group"))
    expect_equal(out$data$intensity[out$data$col_id == "Adnp_IP04"][4],
                 SummarizedExperiment::assay(sce_mq_final, "log2_iBAQ")[4, "Adnp_IP04"])
    expect_equal(out$data$group, gsub("_IP.*$", "", out$data$col_id))

    out <- makeIntensityBoxplots(sce_mq_final, assayName = "iBAQ",
                                 doLog = TRUE, ylab = "lab")
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("col_id", "intensity", "sample", "group"))
    expect_equal(out$data$intensity[out$data$col_id == "Adnp_IP04"][4],
                 SummarizedExperiment::assay(sce_mq_final, "iBAQ")[4, "Adnp_IP04"])
    expect_equal(out$data$group, gsub("_IP.*$", "", out$data$col_id))


    out <- makeMeanSDPlot(sce_mq_final, assayName = "log2_iBAQ",
                          xlab = "xlab", ylab = "ylab")
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 5)
    expect_named(out$data, c("pid", "col_id", "intensity", "mean_intensity",
                             "sd_intensity"))
    expect_equal(out$data$intensity[out$data$pid == "Zmym4" & out$data$col_id == "Adnp_IP04"],
                 SummarizedExperiment::assay(sce_mq_final, "log2_iBAQ")["Zmym4", "Adnp_IP04"])
    expect_equal(length(out$data$mean_intensity[out$data$pid == "Zmym4"]), 9)
    expect_equal(length(unique(out$data$mean_intensity[out$data$pid == "Zmym4"])), 1)
    expect_equal(unique(out$data$mean_intensity[out$data$pid == "Zmym4"]),
                 mean(SummarizedExperiment::assay(sce_mq_final, "log2_iBAQ")["Zmym4", ]))
    expect_equal(unique(out$data$sd_intensity[out$data$pid == "Zmym4"]),
                 stats::sd(SummarizedExperiment::assay(sce_mq_final, "log2_iBAQ")["Zmym4", ]))

    out <- makeMeanSDPlot(sce_mq_final, assayName = "iBAQ",
                          xlab = "xlab", ylab = "ylab")
    expect_s3_class(out, "ggplot")
    expect_equal(ncol(out$data), 5)
    expect_named(out$data, c("pid", "col_id", "intensity", "mean_intensity",
                             "sd_intensity"))
    expect_equal(out$data$intensity[out$data$pid == "Zmym4" & out$data$col_id == "Adnp_IP04"],
                 SummarizedExperiment::assay(sce_mq_final, "iBAQ")["Zmym4", "Adnp_IP04"])
    expect_equal(length(out$data$mean_intensity[out$data$pid == "Zmym4"]), 9)
    expect_equal(length(unique(out$data$mean_intensity[out$data$pid == "Zmym4"])), 1)
    expect_equal(unique(out$data$mean_intensity[out$data$pid == "Zmym4"]),
                 mean(SummarizedExperiment::assay(sce_mq_final, "iBAQ")["Zmym4", ]))
    expect_equal(unique(out$data$sd_intensity[out$data$pid == "Zmym4"]),
                 stats::sd(SummarizedExperiment::assay(sce_mq_final, "iBAQ")["Zmym4", ]))

})
