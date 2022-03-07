test_that("generating summary plots works", {

    ## Fails with wrong arguments
    ## --------------------------------------------------------------------- ##
    ## makeIntensityBoxplots
    expect_error(makeIntensityBoxplots(sce = 1, assayName = "log2_iBAQ",
                                       doLog = TRUE, ylab = "lab"),
                 "'sce' must be of class 'SummarizedExperiment'")
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
    expect_s3_class(makeIntensityBoxplots(sce_mq_final, assayName = "log2_iBAQ",
                                          doLog = FALSE, ylab = "lab"),
                    "ggplot")
    expect_s3_class(makeIntensityBoxplots(sce_mq_final, assayName = "iBAQ",
                                          doLog = TRUE, ylab = "lab"),
                    "ggplot")
    expect_s3_class(makeMeanSDPlot(sce_mq_final, assayName = "log2_iBAQ",
                                   xlab = "xlab", ylab = "ylab"),
                    "ggplot")
    expect_s3_class(makeMeanSDPlot(sce_mq_final, assayName = "iBAQ",
                                   xlab = "xlab", ylab = "ylab"),
                    "ggplot")
})
