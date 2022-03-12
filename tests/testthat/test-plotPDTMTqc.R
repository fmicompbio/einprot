test_that("plotTDTMTqc works", {
    pdOutputFolder <- system.file("extdata", "pdtmt_example",
                                  package = "einprot")
    pdResultName <- "Fig2_m23139_RTS_QC_varMods"

    ## Fails with wrong arguments
    expect_error(plotPDTMTqc(pdOutputFolder = 1,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'pdOutputFolder' must be of class 'character'")
    expect_error(plotPDTMTqc(pdOutputFolder = c(pdOutputFolder, pdOutputFolder),
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'pdOutputFolder' must have length 1")
    expect_error(plotPDTMTqc(pdOutputFolder = "missing",
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "Missing files")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = 1,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'pdResultName' must be of class 'character'")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = c(pdResultName, pdResultName),
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'pdResultName' must have length 1")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = "missing",
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "Missing files")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = 1, poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'masterOnly' must be of class 'logical'")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = c(TRUE, FALSE), poiText = "",
                             doPlot = TRUE, textSize = 4),
                 "'masterOnly' must have length 1")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = 1,
                             doPlot = TRUE, textSize = 4),
                 "'poiText' must be of class 'character'")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = c("1", "2"),
                             doPlot = TRUE, textSize = 4),
                 "'poiText' must have length 1")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = 1, textSize = 4),
                 "'doPlot' must be of class 'logical'")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = c(TRUE, FALSE), textSize = 4),
                 "'doPlot' must have length 1")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = "4"),
                 "'textSize' must be of class 'numeric'")
    expect_error(plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                             pdResultName = pdResultName,
                             masterOnly = FALSE, poiText = "",
                             doPlot = TRUE, textSize = c(4, 5)),
                 "'textSize' must have length 1")

    ## Works with correct arguments
    out <- plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                       pdResultName = pdResultName,
                       masterOnly = FALSE, poiText = "",
                       doPlot = FALSE, textSize = 4)
    expect_type(out, "list")
    expect_length(out, 12)
    for (i in seq_len(12)) {
        expect_s3_class(out[[i]], "ggplot")
    }

    ## masterOnly = TRUE
    out <- plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                       pdResultName = pdResultName,
                       masterOnly = TRUE, poiText = "",
                       doPlot = FALSE, textSize = 4)
    expect_type(out, "list")
    expect_length(out, 12)
    for (i in seq_len(12)) {
        expect_s3_class(out[[i]], "ggplot")
    }

    ## doPlot = TRUE
    pdf(tempfile(fileext = ".pdf"))
    out <- plotPDTMTqc(pdOutputFolder = pdOutputFolder,
                       pdResultName = pdResultName,
                       masterOnly = TRUE, poiText = "",
                       doPlot = TRUE, textSize = 4)
    expect_type(out, "list")
    expect_length(out, 12)
    for (i in seq_len(12)) {
        expect_s3_class(out[[i]], "ggplot")
    }
    dev.off()

})
