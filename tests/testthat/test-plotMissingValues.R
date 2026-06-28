test_that("missing value plots work", {

    ## -------------------------------------------------------------------------
    ## plotMissingValuesHeatmap
    ## -------------------------------------------------------------------------
    expect_error(plotMissingValuesHeatmap(
        sce = 1, assayMissing = "imputed_iBAQ", onlyRowsWithMissing = FALSE,
        settings = "clustered"),
        "'sce' must be of class 'SummarizedExperiment'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = 1,
        onlyRowsWithMissing = FALSE, settings = "clustered"),
        "'assayMissing' must be of class 'character'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = c("log2_iBAQ", "imputed_iBAQ"),
        onlyRowsWithMissing = FALSE, settings = "clustered"),
        "'assayMissing' must have length 1")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "missing",
        onlyRowsWithMissing = FALSE, settings = "clustered"),
        "All values in 'assayMissing' must be one of")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "log2_iBAQ_withNA",
        onlyRowsWithMissing = FALSE, settings = "clustered"),
        "Assay contains missing values")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "imputed_iBAQ",
        onlyRowsWithMissing = "FALSE", settings = "clustered"),
        "'onlyRowsWithMissing' must be of class 'logical'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "imputed_iBAQ",
        onlyRowsWithMissing = c(TRUE, FALSE), settings = "clustered"),
        "'onlyRowsWithMissing' must have length 1")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "imputed_iBAQ",
        onlyRowsWithMissing = FALSE, settings = 1),
        "'settings' must be of class 'character'")
    expect_error(plotMissingValuesHeatmap(
        sce = sce_mq_preimputation, assayMissing = "imputed_iBAQ",
        onlyRowsWithMissing = FALSE, settings = c("clustered", "clustered")),
        "'settings' must have length 1")

    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation,
                                    assayMissing = "imputed_iBAQ",
                                    onlyRowsWithMissing = FALSE,
                                    settings = "clustered")
    expect_s4_class(out, "Heatmap")
    expect_length(out@column_names_param$labels, 9L)
    expect_true(out@column_names_param$show)
    expect_length(out@row_names_param$labels, 150L)
    expect_equal(out@column_title, "Missing value pattern (white = missing)")
    expect_equal(out@name, "imputed")

    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation,
                                    assayMissing = "imputed_iBAQ",
                                    onlyRowsWithMissing = FALSE,
                                    settings = "sorted")
    expect_s4_class(out, "Heatmap")
    expect_length(out@column_names_param$labels, 9L)
    expect_true(out@column_names_param$show)
    expect_length(out@row_names_param$labels, 150L)
    expect_equal(out@column_title, "Missing value pattern (white = missing)")
    expect_equal(out@name, "imputed")

    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation,
                                    assayMissing = "imputed_iBAQ",
                                    onlyRowsWithMissing = FALSE,
                                    settings = NULL,
                                    cluster_rows = FALSE, show_column_names = FALSE)
    expect_s4_class(out, "Heatmap")
    expect_length(out@column_names_param$labels, 9L)
    expect_false(out@column_names_param$show)
    expect_length(out@row_names_param$labels, 150L)
    expect_equal(out@column_title, character(0))
    expect_equal(substr(out@name, 1, 7), "matrix_")

    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation,
                                    assayMissing = "imputed_iBAQ",
                                    onlyRowsWithMissing = TRUE,
                                    settings = "clustered")
    expect_s4_class(out, "Heatmap")
    expect_length(out@column_names_param$labels, 9L)
    expect_true(out@column_names_param$show)
    expect_length(out@row_names_param$labels, 103L)
    expect_equal(out@column_title, "Missing value pattern (white = missing)")
    expect_equal(out@name, "imputed")

    kp <- which(rowSums(assay(sce_mq_preimputation, "imputed_iBAQ")) == 0)
    out <- plotMissingValuesHeatmap(sce = sce_mq_preimputation[kp, ],
                                    assayMissing = "imputed_iBAQ",
                                    onlyRowsWithMissing = TRUE,
                                    settings = "clustered")
    expect_s4_class(out, "Heatmap")
    expect_length(out@column_names_param$labels, 9L)
    expect_true(out@column_names_param$show)
    expect_length(out@row_names_param$labels, 0L)
    expect_equal(out@column_title, "Missing value pattern (white = missing)")
    expect_equal(out@name, "imputed")

    ## -------------------------------------------------------------------------
    ## plotFractionDetectedPerSample
    ## -------------------------------------------------------------------------
    expect_error(plotFractionDetectedPerSample(
        dfNA = 1),
        "'dfNA' must be of class 'data.frame'")
    dfNA <- as.data.frame(nbr_na_mq$nNAcols) %>%
        dplyr::rename(sample = name)
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA),
        'all(c("sample", "pNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)
    dfNA <- as.data.frame(nbr_na_mq$nNAcols) %>%
        dplyr::rename(sample = name)
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, valueType = 1
    ), "'valueType' must be of class 'character'")
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, valueType = c("fraction", "percentage")
    ), "'valueType' must have length 1")
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, valueType = "error"
    ), "'valueType' must be one of")
    dfNA$pNA <- 1 + dfNA$pNA
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, valueType = "fraction"
    ), "'$dfNApNA' must be within", fixed = TRUE)
    dfNA$pNA <- 99 + dfNA$pNA
    expect_error(plotFractionDetectedPerSample(
        dfNA = dfNA, valueType = "percentage"
    ), "'$dfNApNA' must be within", fixed = TRUE)

    expect_message(
        out <- plotFractionDetectedPerSample(
            dfNA = as.data.frame(nbr_na_mq$nNAcols) %>%
                dplyr::rename(sample = name)),
        "Inferred valueType: fraction")
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("sample", "nNA", "pNA", "assay"))

    out2 <- plotFractionDetectedPerSample(
        dfNA = DataFrame(as.data.frame(nbr_na_mq$nNAcols) %>%
            dplyr::rename(sample = name)),
        valueType = "fraction")
    expect_s3_class(out2, "ggplot")
    expect_named(out2$data, c("sample", "nNA", "pNA", "assay"))
    expect_identical(out$data, out2$data)

    expect_message(
        out3 <- plotFractionDetectedPerSample(
            dfNA = DataFrame(as.data.frame(nbr_na_mq$nNAcols) %>%
                                 dplyr::rename(sample = name) %>%
                                 dplyr::mutate(pNA = 100 * pNA))),
        "Inferred valueType: percentage")
    expect_s3_class(out3, "ggplot")
    expect_named(out3$data, c("sample", "nNA", "pNA", "assay"))

    out4 <- plotFractionDetectedPerSample(
        dfNA = DataFrame(as.data.frame(nbr_na_mq$nNAcols) %>%
                             dplyr::rename(sample = name) %>%
                             dplyr::mutate(pNA = 100 * pNA)),
        valueType = "percentage")
    expect_s3_class(out4, "ggplot")
    expect_named(out4$data, c("sample", "nNA", "pNA", "assay"))

    ## -------------------------------------------------------------------------
    ## plotDetectedInSamples
    ## -------------------------------------------------------------------------
    expect_error(plotDetectedInSamples(
        dfNA = 1),
        "'dfNA' must be of class 'data.frame'")
    dfNA <- as.data.frame(nbr_na_mq$nNArows) %>%
        dplyr::rename(sample = name)
    colnames(dfNA) <- paste0("a", colnames(dfNA))
    expect_error(plotDetectedInSamples(
        dfNA = dfNA),
        'all(c("pNA", "nNA") %in% colnames(dfNA)) is not TRUE',
        fixed = TRUE)

    out <- plotDetectedInSamples(dfNA = as.data.frame(nbr_na_mq$nNArows))
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n", "nObs"))
    for (i in c(0, seq_len(9))) {
        expect_equal(sum(nbr_na_mq$nNArows$nNA == 9 - i),
                     out$data$n[out$data$nObs == i])
    }
    expect_equal(levels(out$data$nObs), as.character(c(0, seq_len(9))))

    out <- plotDetectedInSamples(dfNA = DataFrame(nbr_na_mq$nNArows))
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n", "nObs"))
    for (i in c(0, seq_len(9))) {
        expect_equal(sum(nbr_na_mq$nNArows$nNA == 9 - i),
                     out$data$n[out$data$nObs == i])
    }
    expect_equal(levels(out$data$nObs), as.character(c(0, seq_len(9))))

    ## PD data
    out <- plotDetectedInSamples(dfNA = as.data.frame(nbr_na_pd$nNArows))
    expect_s3_class(out, "ggplot")
    expect_named(out$data, c("nNA", "n", "nObs"))
    for (i in c(0, seq_len(16))) {
        expect_equal(sum(nbr_na_pd$nNArows$nNA == 16 - i),
                     out$data$n[out$data$nObs == i])
    }
    expect_equal(levels(out$data$nObs), as.character(c(0, seq_len(16))))
})
