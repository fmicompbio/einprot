test_that("missing value plots work", {

    ## --------------------------------------------------------------------- ##
    ## plotMissingValuesHeatmap
    ## --------------------------------------------------------------------- ##
    expect_error(plotMissingValuesHeatmap(
        sce = 1, assayMissing = "imputed_iBAQ"),
        "'sce' must be of class 'SummarizedExperiment'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = 1),
        "'assayMissing' must be of class 'character'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = c("log2_iBAQ", "imputed_iBAQ")),
        "'assayMissing' must have length 1")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "missing"),
        "All values in 'assayMissing' must be one of")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "log2_iBAQ_withNA"),
        "Assay contains missing values")

    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation,
                                    assayMissing = "imputed_iBAQ")
    expect_s4_class(out, "Heatmap")

    ## --------------------------------------------------------------------- ##
    ## plotFractionDetectedPerSample
    ## --------------------------------------------------------------------- ##
    expect_error(plotFractionDetectedPerSample(
        dfNA = 1, aName = "iBAQ"),
        "'dfNA' must be of class 'DFrame'")
    dfNA <- nbr_na_mq$nNAcols
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, aName = "iBAQ"),
        'all(c("name", "pNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na_mq$nNAcols, aName = 1),
        "'aName' must be of class 'character'")
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na_mq$nNAcols, aName = c("iBAQ", "iBAQ")),
        "'aName' must have length 1")
    expect_error(plotFractionDetectedPerSample(
        dfNA = nbr_na_mq$nNAcols, aName = "missing"),
        "All values in 'aName' must be one of")

    out <- plotFractionDetectedPerSample(dfNA = nbr_na_mq$nNAcols,
                                         aName = "iBAQ")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("name", "nNA", "pNA", "assay"))

    ## --------------------------------------------------------------------- ##
    ## plotDetectedInSamples
    ## --------------------------------------------------------------------- ##
    expect_error(plotDetectedInSamples(
        dfNA = 1, aName = "iBAQ"),
        "'dfNA' must be of class 'DFrame'")
    dfNA <- nbr_na_mq$nNArows
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotDetectedInSamples(
        dfNA = dfNA, aName = "iBAQ"),
        'all(c("name", "pNA", "nNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na_mq$nNArows, aName = 1),
        "'aName' must be of class 'character'")
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na_mq$nNArows, aName = c("iBAQ", "iBAQ")),
        "'aName' must have length 1")
    expect_error(plotDetectedInSamples(
        dfNA = nbr_na_mq$nNArows, aName = "missing"),
        "All values in 'aName' must be one of")

    out <- plotDetectedInSamples(dfNA = nbr_na_mq$nNArows,
                                 aName = "iBAQ")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n", "nObs"))
    for (i in c(0, seq_len(9))) {
        expect_equal(sum(nbr_na_mq$nNArows$nNA == 9 - i),
                     out$data$n[out$data$nObs == i])
    }
    expect_equal(levels(out$data$nObs), as.character(c(0, seq_len(9))))

    ## PD data
    out <- plotDetectedInSamples(dfNA = nbr_na_pd$nNArows,
                                 aName = "Abundance")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n", "nObs"))
    for (i in c(0, seq_len(16))) {
        expect_equal(sum(nbr_na_pd$nNArows$nNA == 16 - i),
                     out$data$n[out$data$nObs == i])
    }
    expect_equal(levels(out$data$nObs), as.character(c(0, seq_len(16))))
})
