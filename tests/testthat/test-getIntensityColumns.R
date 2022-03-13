test_that("getting column names of mq file works", {
    expect_error(getColumnNames(1),
                 "'inFile' must be of class 'character'")
    expect_error(getColumnNames(c("1", "2")),
                 "'inFile' must have length 1")
    expect_error(getColumnNames("missing"),
                 "file.exists(inFile) is not TRUE", fixed = TRUE)

    mqf <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                       package = "einprot")
    cnm <- getColumnNames(inFile = mqf)
    expect_type(cnm, "character")
    expect_length(cnm, 116)
    expect_equal(cnm[1], "Protein.IDs")
    expect_equal(cnm[2], "Majority.protein.IDs")
    expect_true(any(grepl("^LFQ\\.intensity\\.", cnm)))
    expect_true(any(grepl("^iBAQ\\.", cnm)))
})

test_that("getting column names of pd file works", {
    pdf <- system.file("extdata", "pdtmt_example",
                       "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                       package = "einprot")
    cnp <- getColumnNames(inFile = pdf)
    expect_type(cnp, "character")
    expect_length(cnp, 142)
    expect_equal(cnp[1], "Proteins.Unique.Sequence.ID")
    expect_equal(cnp[2], "Checked")
    expect_true(any(grepl("^Abundance\\.", cnp)))
})

test_that("getting intensity columns works (MQ)", {
    mqf <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                       package = "einprot")

    ## Mis-specified arguments
    expect_error(getIntensityColumns(1),
                 "'inFile' must be of class 'character'")
    expect_error(getIntensityColumns("missing"),
                 "file.exists(inFile) is not TRUE", fixed = TRUE)
    expect_error(getIntensityColumns(mqf, iColPattern = 1),
                 "'iColPattern' must be of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = c("a", "b")),
                 "'iColPattern' must have length 1")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     includeOnlySamples = 1),
                 "'includeOnlySamples' must be of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     excludeSamples = 1),
                 "'excludeSamples' must be of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     stopIfEmpty = 1),
                 "'stopIfEmpty' must be of class 'logical'")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     stopIfEmpty = c(TRUE, FALSE)),
                 "'stopIfEmpty' must have length 1")

    ## No matching samples - stop
    expect_error(getIntensityColumns(mqf, iColPattern = "missing",
                                     stopIfEmpty = TRUE),
                 "No samples were found matching")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     includeOnlySamples = "missing",
                                     stopIfEmpty = TRUE),
                 "No samples were retained")

    ## No matching samples - empty output
    expect_true(all(lengths(getIntensityColumns(mqf, iColPattern = "missing",
                                                stopIfEmpty = FALSE)) == 0))
    expect_equal(length(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                            includeOnlySamples = "missing",
                                            stopIfEmpty = FALSE)$iColsAll), 9)
    expect_equal(length(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                            includeOnlySamples = "missing",
                                            stopIfEmpty = FALSE)$iCols), 0)

    ## Keep all samples
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.")
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 9)
    ic <- getIntensityColumns(mqf, iColPattern = "^iBAQ\\.")
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 9)

    ## Remove samples
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              excludeSamples = c("LFQ.intensity.Adnp_IP04",
                                                 "LFQ.intensity.Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 7)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_false(any(c("LFQ.intensity.Adnp_IP04",
                       "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Specify only sample name
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              excludeSamples = c("Adnp_IP04", "Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 7)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_false(any(c("LFQ.intensity.Adnp_IP04",
                       "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Specify only group name
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              excludeSamples = c("Adnp", "Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 5)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_false(any(c("LFQ.intensity.Adnp_IP04",
                       "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Retain samples
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              includeOnlySamples = c("LFQ.intensity.Adnp_IP04",
                                                     "LFQ.intensity.Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 2)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Specify only sample name
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              includeOnlySamples = c("Adnp_IP04", "Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 2)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Specify only group name
    ic <- getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                              includeOnlySamples = c("Adnp", "Chd4BF_IP07"))
    expect_equal(length(ic$iColsAll), 9)
    expect_equal(length(ic$iCols), 4)
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iColsAll))
    expect_true(all(c("LFQ.intensity.Adnp_IP04",
                      "LFQ.intensity.Chd4BF_IP07") %in% ic$iCols))

    ## Specify both include and exclude -> error
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     includeOnlySamples = c("LFQ.intensity.Adnp_IP04",
                                                            "LFQ.intensity.Chd4BF_IP07"),
                                     excludeSamples = c("Adnp")),
                 "Please specify max one")
})

test_that("getting intensity columns works (PD)", {
    pdf <- system.file("extdata", "pdtmt_example",
                       "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                       package = "einprot")

    ## No matching samples - stop
    expect_error(getIntensityColumns(pdf, iColPattern = "missing",
                                     stopIfEmpty = TRUE),
                 "No samples were found matching")
    expect_error(getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = "missing",
        stopIfEmpty = TRUE),
        "No samples were retained")

    ## No matching samples - empty output
    expect_true(all(lengths(getIntensityColumns(
        pdf, iColPattern = "missing", stopIfEmpty = FALSE)) == 0))
    expect_equal(length(getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = "missing",
        stopIfEmpty = FALSE)$iColsAll), 16)
    expect_equal(length(getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = "missing",
        stopIfEmpty = FALSE)$iCols), 0)

    ## Keep all samples
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.")
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 16)
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundances\\.Count\\.F.+\\.Sample\\.")
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 16)

    ## Remove samples
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        excludeSamples = c("Abundance.F12.128C.Sample.HIS4KO_S05",
                           "Abundance.F12.131N.Sample.URA2KO_S10"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 14)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_false(any(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                       "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Specify only sample name
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        excludeSamples = c("HIS4KO_S05", "URA2KO_S10"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 14)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_false(any(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                       "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Specify only group name
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        excludeSamples = c("HIS4KO", "URA2KO"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 8)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_false(any(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                       "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Retain samples
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = c("Abundance.F12.128C.Sample.HIS4KO_S05",
                               "Abundance.F12.131N.Sample.URA2KO_S10"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 2)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Specify only sample name
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = c("HIS4KO_S05", "URA2KO_S10", "missing"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 2)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Specify only group name
    ic <- getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = c("HIS4KO", "URA2KO"))
    expect_equal(length(ic$iColsAll), 16)
    expect_equal(length(ic$iCols), 8)
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iColsAll))
    expect_true(all(c("Abundance.F12.128C.Sample.HIS4KO_S05",
                      "Abundance.F12.131N.Sample.URA2KO_S10") %in% ic$iCols))

    ## Specify both include and exclude -> error
    expect_error(getIntensityColumns(
        pdf, iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        includeOnlySamples = c("Abundance.F12.128C.Sample.HIS4KO_S05",
                               "Abundance.F12.131N.Sample.URA2KO_S10"),
        excludeSamples = "URA2KO"),
        "Please specify max one")
})
