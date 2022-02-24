test_that("generating summary plots works", {
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

    ## Fails with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## makeIntensityBoxplots
    expect_error(makeIntensityBoxplots(qft = 1, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = "lab"),
                 "'qft' must be of class 'QFeatures'")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = 1,
                                       doLog = TRUE, ylab = "lab"),
                 "'assayName' must be of class 'character'")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = c("iBAQ", "log2_iBAQ"),
                                       doLog = TRUE, ylab = "lab"),
                 "'assayName' must have length 1")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = "missing",
                                       doLog = TRUE, ylab = "lab"),
                 "All values in 'assayName' must be one of")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = "log2_iBAQ",
                                       doLog = 1, ylab = "lab"),
                 "'doLog' must be of class 'logical'")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = "log2_iBAQ",
                                       doLog = c(TRUE, FALSE), ylab = "lab"),
                 "'doLog' must have length 1")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = 1),
                 "'ylab' must be of class 'character'")
    expect_error(makeIntensityBoxplots(qft = qft, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = c("lab1", "lab2")),
                 "'ylab' must have length 1")

    ## makeMeanSDPlot
    expect_error(makeMeanSDPlot(qft = 1, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = "ylab"),
                 "'qft' must be of class 'QFeatures'")
    expect_error(makeMeanSDPlot(qft = qft, assayName = 1,
                                xlab = "xlab", ylab = "ylab"),
                 "'assayName' must be of class 'character'")
    expect_error(makeMeanSDPlot(qft = qft, assayName = c("iBAQ", "log2_iBAQ"),
                                xlab = "xlab", ylab = "ylab"),
                 "'assayName' must have length 1")
    expect_error(makeMeanSDPlot(qft = qft, assayName = "missing",
                                xlab = "xlab", ylab = "ylab"),
                 "All values in 'assayName' must be one of")
    expect_error(makeMeanSDPlot(qft = qft, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = 1),
                 "'ylab' must be of class 'character'")
    expect_error(makeMeanSDPlot(qft = qft, assayName = "log2_iBAQ",
                                xlab = "xlab", ylab = c("ylab1", "ylab2")),
                 "'ylab' must have length 1")
    expect_error(makeMeanSDPlot(qft = qft, assayName = "log2_iBAQ",
                                xlab = 1, ylab = "ylab"),
                 "'xlab' must be of class 'character'")
    expect_error(makeMeanSDPlot(qft = qft, assayName = "log2_iBAQ",
                                xlab = c("xlab1", "xlab2"), ylab = "ylab"),
                 "'xlab' must have length 1")

    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    expect_s3_class(makeIntensityBoxplots(qft, assayName = "log2_iBAQ",
                                          doLog = FALSE, ylab = "lab"),
                    "ggplot")
    expect_s3_class(makeIntensityBoxplots(qft, assayName = "iBAQ",
                                          doLog = TRUE, ylab = "lab"),
                    "ggplot")
    expect_s3_class(makeMeanSDPlot(qft, assayName = "log2_iBAQ",
                                   xlab = "xlab", ylab = "ylab"),
                    "ggplot")
    expect_s3_class(makeMeanSDPlot(qft, assayName = "iBAQ",
                                   xlab = "xlab", ylab = "ylab"),
                    "ggplot")
})
