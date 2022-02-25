test_that("missing value plots work", {
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

    ## --------------------------------------------------------------------- ##
    ## plotMissingValuesHeatmap
    ## --------------------------------------------------------------------- ##
    expect_error(plotMissingValuesHeatmap(
        qft = 1, assayMissing = "imputed_iBAQ"),
        "'qft' must be of class 'QFeatures'")
    expect_error(plotMissingValuesHeatmap(
        qft = qft, assayMissing = 1),
        "'assayMissing' must be of class 'character'")
    expect_error(plotMissingValuesHeatmap(
        qft = qft, assayMissing = c("log2_iBAQ", "imputed_iBAQ")),
        "'assayMissing' must have length 1")
    expect_error(plotMissingValuesHeatmap(
        qft = qft, assayMissing = "missing"),
        "All values in 'assayMissing' must be one of")

    out <- plotMissingValuesHeatmap(qft = qft, assayMissing = "imputed_iBAQ")
    expect_s4_class(out, "Heatmap")

    ## --------------------------------------------------------------------- ##
    ## plotFractionDetectedPerSample
    ## --------------------------------------------------------------------- ##
    expect_error(plotFractionDetectedPerSample(
        dfNA = 1, aName = "iBAQ"),
        "'dfNA' must be of class 'DFrame'")
    dfNA <- nbr_na$nNAcols
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, aName = "iBAQ"),
        'all(c("name", "pNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na$nNAcols, aName = 1),
        "'aName' must be of class 'character'")
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na$nNAcols, aName = c("iBAQ", "imputed_iBAQ")),
        "'aName' must have length 1")
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na$nNAcols, aName = "missing"),
        "All values in 'aName' must be one of")

    out <- plotFractionDetectedPerSample(dfNA = nbr_na$nNAcols,
                                         aName = "iBAQ")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("assay", "name", "nNA", "pNA"))

    ## --------------------------------------------------------------------- ##
    ## plotDetectedInSamples
    ## --------------------------------------------------------------------- ##
    expect_error(plotDetectedInSamples(
        dfNA = 1, aName = "iBAQ"),
        "'dfNA' must be of class 'DFrame'")
    dfNA <- nbr_na$nNArows
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotDetectedInSamples(
        dfNA = dfNA, aName = "iBAQ"),
        'all(c("name", "pNA", "nNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na$nNArows, aName = 1),
        "'aName' must be of class 'character'")
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na$nNArows, aName = c("iBAQ", "imputed_iBAQ")),
        "'aName' must have length 1")
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na$nNArows, aName = "missing"),
        "All values in 'aName' must be one of")

    out <- plotDetectedInSamples(dfNA = nbr_na$nNArows,
                                 aName = "iBAQ")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n"))
})
