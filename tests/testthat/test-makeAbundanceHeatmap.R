test_that("makeAbundanceHeatmap works", {
    expect_error(makeAbundanceHeatmap(sce = 1, assayToPlot = "log2_iBAQ",
                                      doCenter = FALSE, settings = "report"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = 1,
                                      doCenter = FALSE, settings = "report"),
                 "'assayToPlot' must be of class 'character'")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = c("log2_iBAQ", "iBAQ"),
                                      doCenter = FALSE, settings = "report"),
                 "'assayToPlot' must have length 1")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "missing",
                                      doCenter = FALSE, settings = "report"),
                 "All values in 'assayToPlot' must be one of")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                                      doCenter = 1, settings = "report"),
                 "'doCenter' must be of class 'logical'")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                                      doCenter = c(TRUE, FALSE), settings = "report"),
                 "'doCenter' must have length 1")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                                      doCenter = FALSE, settings = TRUE),
                 "'settings' must be of class 'character'")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                                      doCenter = FALSE, settings = c("report", "export")),
                 "'settings' must have length 1")
    expect_error(makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                                      doCenter = FALSE, settings = "missing"),
                 "All values in 'settings' must be one of")

    ## settings = 'report'
    ht <- makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                               doCenter = FALSE, settings = "report")
    expect_s4_class(ht, "Heatmap")
    expect_equal(ht@name, "log2_iBAQ")
    expect_equal(dim(ht@matrix), dim(sce_mq_final))
    expect_equal(ht@row_names_param$labels, rownames(sce_mq_final))
    expect_equal(ht@row_names_param$show, FALSE)

    ## settings = 'export'
    ht <- makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                               doCenter = TRUE, settings = "export")
    expect_s4_class(ht, "Heatmap")
    expect_equal(ht@name, "log2_iBAQ\ncentered")
    expect_equal(dim(ht@matrix), dim(sce_mq_final))
    expect_equal(rowMeans(ht@matrix), rep(0, nrow(sce_mq_final)),
                 ignore_attr = TRUE)
    expect_equal(ht@row_names_param$labels, rownames(sce_mq_final))
    expect_equal(ht@row_names_param$show, TRUE)

    ## settings = NULL
    ht <- makeAbundanceHeatmap(sce = sce_mq_final, assayToPlot = "log2_iBAQ",
                               doCenter = FALSE, settings = NULL,
                               cluster_columns = TRUE, show_column_names = FALSE)
    expect_s4_class(ht, "Heatmap")
    expect_equal(ht@name, "log2_iBAQ")
    expect_equal(dim(ht@matrix), dim(sce_mq_final))
    expect_equal(ht@row_names_param$labels, rownames(sce_mq_final))
    expect_equal(ht@column_names_param$show, FALSE)
})
