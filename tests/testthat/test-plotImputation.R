test_that("plotting of imputation results works", {

    ## --------------------------------------------------------------------- ##
    ## plotImputationDistribution
    ## --------------------------------------------------------------------- ##
    expect_error(plotImputationDistribution(
        qft = 1, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'qft' must be of class 'QFeatures'")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = 1, assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayToPlot' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = c("log2_iBAQ", "imputed_iBAQ"),
        assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayToPlot' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "missing", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "All values in 'assayToPlot' must be one of")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = 1,
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayImputation' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ",
        assayImputation = c("log2_iBAQ", "imputed_iBAQ"),
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayImputation' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "missing",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "All values in 'assayImputation' must be one of")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = 1, xlab = "log2 intensity"),
        "'iColPattern' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = c("^iBAQ\\.", "^Intensity\\."), xlab = "log2 intensity"),
        "'iColPattern' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = 1),
        "'xlab' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft_mq_final, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = c("lab1", "lab2")),
        "'xlab' must have length 1")

    out <- plotImputationDistribution(qft = qft_mq_final, assayToPlot = "log2_iBAQ",
                                      assayImputation = "imputed_iBAQ",
                                      iColPattern = "^iBAQ\\.", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_false(any(grepl("^iBAQ\\.", out$data$sample)))

    out <- plotImputationDistribution(qft = qft_mq_final, assayToPlot = "log2_iBAQ",
                                      assayImputation = "imputed_iBAQ",
                                      iColPattern = "", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_true(all(grepl("^iBAQ\\.", out$data$sample)))
})
