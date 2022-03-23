test_that("plotting of imputation results works", {

    ## --------------------------------------------------------------------- ##
    ## plotImputationDistribution
    ## --------------------------------------------------------------------- ##
    expect_error(plotImputationDistribution(
        sce = 1, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        xlab = "log2 intensity"),
        "'sce' must be of class 'SummarizedExperiment'")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = 1, assayImputation = "imputed_iBAQ",
        xlab = "log2 intensity"),
        "'assayToPlot' must be of class 'character'")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = c("log2_iBAQ", "imputed_iBAQ"),
        assayImputation = "imputed_iBAQ",
        xlab = "log2 intensity"),
        "'assayToPlot' must have length 1")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "missing", assayImputation = "imputed_iBAQ",
        xlab = "log2 intensity"),
        "All values in 'assayToPlot' must be one of")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ", assayImputation = 1,
        xlab = "log2 intensity"),
        "'assayImputation' must be of class 'character'")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ",
        assayImputation = c("log2_iBAQ", "imputed_iBAQ"),
        xlab = "log2 intensity"),
        "'assayImputation' must have length 1")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "missing",
        xlab = "log2 intensity"),
        "All values in 'assayImputation' must be one of")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        xlab = 1),
        "'xlab' must be of class 'character'")
    expect_error(plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        xlab = c("lab1", "lab2")),
        "'xlab' must have length 1")

    out <- plotImputationDistribution(
        sce = sce_mq_final, assayToPlot = "log2_iBAQ",
        assayImputation = "imputed_iBAQ", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_false(any(grepl("^iBAQ\\.", out$data$sample)))
    expect_equal(nrow(out$data), 150 * 9)
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("pid", "sample", "log2intensity", "imputed"))

    ## PD input
    out <- plotImputationDistribution(
        sce = sce_pd_final, assayToPlot = "log2_Abundance",
        assayImputation = "imputed_Abundance", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_false(any(grepl("^Abundance\\.", out$data$sample)))
    expect_equal(nrow(out$data), 70 * 16)
    expect_equal(ncol(out$data), 4)
    expect_named(out$data, c("pid", "sample", "log2intensity", "imputed"))
})
