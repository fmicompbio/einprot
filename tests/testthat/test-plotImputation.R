test_that("plotting of imputation results works", {
    ## Preparation
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
                 "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
                 "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
    ecol <- paste0("iBAQ.", samples)
    qft <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "iBAQ",
                                    sep = "\t", nrows = 70)
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    qft <- addSampleAnnots(qft, iColPattern = "^iBAQ\\.",
                           sampleAnnot = sampleAnnot, mergeGroups = list())
    qft <- fixFeatureIds(qft)
    qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ")
    qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ_withNA")
    tmp <- qft[["log2_iBAQ"]]
    SummarizedExperiment::assay(tmp) <- !is.finite(SummarizedExperiment::assay(tmp))
    qft <- QFeatures::addAssay(qft, tmp, name = "imputed_iBAQ")
    qft <- QFeatures::zeroIsNA(qft, "iBAQ")
    qft <- QFeatures::infIsNA(qft, "log2_iBAQ")
    qft <- QFeatures::infIsNA(qft, "log2_iBAQ_withNA")
    nbr_na <- QFeatures::nNA(qft, i = seq_along(qft))
    set.seed(123)
    qft <- QFeatures::impute(qft, method = "MinProb", i = "log2_iBAQ")

    ## --------------------------------------------------------------------- ##
    ## plotImputationDistribution
    ## --------------------------------------------------------------------- ##
    expect_error(plotImputationDistribution(
        qft = 1, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'qft' must be of class 'QFeatures'")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = 1, assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayToPlot' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = c("log2_iBAQ", "imputed_iBAQ"),
        assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayToPlot' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "missing", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "All values in 'assayToPlot' must be one of")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = 1,
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayImputation' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ",
        assayImputation = c("log2_iBAQ", "imputed_iBAQ"),
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "'assayImputation' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = "missing",
        iColPattern = "^iBAQ\\.", xlab = "log2 intensity"),
        "All values in 'assayImputation' must be one of")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = 1, xlab = "log2 intensity"),
        "'iColPattern' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = c("^iBAQ\\.", "^Intensity\\."), xlab = "log2 intensity"),
        "'iColPattern' must have length 1")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = 1),
        "'xlab' must be of class 'character'")
    expect_error(plotImputationDistribution(
        qft = qft, assayToPlot = "log2_iBAQ", assayImputation = "imputed_iBAQ",
        iColPattern = "^iBAQ\\.", xlab = c("lab1", "lab2")),
        "'xlab' must have length 1")

    out <- plotImputationDistribution(qft = qft, assayToPlot = "log2_iBAQ",
                                      assayImputation = "imputed_iBAQ",
                                      iColPattern = "^iBAQ\\.", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_false(any(grepl("^iBAQ\\.", out$data$sample)))

    out <- plotImputationDistribution(qft = qft, assayToPlot = "log2_iBAQ",
                                      assayImputation = "imputed_iBAQ",
                                      iColPattern = "", xlab = "lab")
    expect_s3_class(out, "ggplot")
    expect_true(any(out$data$imputed))
    expect_true(any(!out$data$imputed))
    expect_true(all(grepl("^iBAQ\\.", out$data$sample)))
})
