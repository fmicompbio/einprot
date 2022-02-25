test_that("getting column names of mq file works", {
    expect_error(getColumnNames(1), "of class 'character'")
    expect_error(getColumnNames(c("1", "2")), "must have length 1")
    expect_error(getColumnNames("missing"), "is not TRUE")

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

test_that("getting intensity columns works", {
    mqf <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                       package = "einprot")

    ## Mis-specified arguments
    expect_error(getIntensityColumns(1), "of class 'character'")
    expect_error(getIntensityColumns("missing"), "is not TRUE")
    expect_error(getIntensityColumns(mqf, iColPattern = 1),
                 "of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = c("a", "b")),
                 "must have length 1")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     includeOnlySamples = 1),
                 "of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     excludeSamples = 1),
                 "of class 'character'")
    expect_error(getIntensityColumns(mqf, iColPattern = "^LFQ\\.intensity\\.",
                                     stopIfEmpty = 1),
                 "of class 'logical'")

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
