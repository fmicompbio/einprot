test_that("seqLogoApp works", {
    testdf <- data.frame(einprotId = LETTERS[seq_len(5)],
                         seqWindow = c("AAAAA", "AEKAP", "AK-PA", "-EKAP", "AEPAK"))
    tmpf <- tempfile(fileext = ".csv")
    write.csv(testdf, file = tmpf, row.names = TRUE, quote = TRUE)
    tmpftxt <- tempfile(fileext = ".txt")
    write.csv(testdf, file = tmpftxt, row.names = TRUE, quote = TRUE)

    testdf2 <- data.frame(einprotId = LETTERS[seq_len(5)],
                         seqWindow2 = c("AAAAA", "AEKAP", "AK-PA", "-EKAP", "AEPAK"))
    tmpf2 <- tempfile(fileext = ".csv")
    write.csv(testdf2, file = tmpf2, row.names = TRUE, quote = TRUE)

    testdf3 <- data.frame(einprotId = LETTERS[seq_len(5)],
                          seqWindow = c("AAAAA", "AEKP", "AK-PA", "-EKAP", "AEPAK"))
    tmpf3 <- tempfile(fileext = ".csv")
    write.csv(testdf3, file = tmpf3, row.names = TRUE, quote = TRUE)

    expect_error(seqLogoApp(seqTableCsv = 1, exportName = "test"),
                 "'seqTableCsv' must be of class 'character'")
    expect_error(seqLogoApp(seqTableCsv = c(tmpf, tmpf), exportName = "test"),
                 "'seqTableCsv' must have length 1")
    expect_error(seqLogoApp(seqTableCsv = "missing", exportName = "test"),
                 "file.exists(seqTableCsv) is not TRUE", fixed = TRUE)
    expect_error(seqLogoApp(seqTableCsv = tmpf2, exportName = "test"),
                 '"seqWindow" %in% colnames(df) is not TRUE', fixed = TRUE)
    expect_error(seqLogoApp(seqTableCsv = tmpf3, exportName = "test"),
                 "All sequences in the seqWindow column")
    expect_error(seqLogoApp(seqTableCsv = tmpftxt, exportName = "test"),
                 'tools::file_ext(seqTableCsv) == "csv" is not TRUE', fixed = TRUE)
    expect_error(seqLogoApp(seqTableCsv = tmpf, exportName = 1),
                 "'exportName' must be of class 'character'")
    expect_error(seqLogoApp(seqTableCsv = tmpf, exportName = c("test", "test2")),
                 "'exportName' must have length 1")

    ## -------------------------------------------------------------------------
    app <- seqLogoApp(seqTableCsv = tmpf, exportName = "test")
    expect_s3_class(app, "shiny.appobj")
})
