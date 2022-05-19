test_that("subtracting baseline abundance works", {
    sce <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:60, nrow = 5)),
        colData = S4Vectors::DataFrame(batch = rep(c("B1", "B2"), each = 6),
                                       group = rep(c("G1", "G2", "G3"), 4))
    )
    rownames(sce) <- paste0("V", 1:5)
    colnames(sce) <- paste0("C", 1:12)
    ## Missing batch column
    scenb <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:60, nrow = 5)),
        colData = S4Vectors::DataFrame(group = rep(c("G1", "G2", "G3"), 4))
    )
    ## Missing group column
    sceng <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:60, nrow = 5)),
        colData = S4Vectors::DataFrame(batch = rep(c("B1", "B2"), each = 6))
    )
    ## Different batch levels
    scedb <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = matrix(1:60, nrow = 5)),
        colData = S4Vectors::DataFrame(batch = rep(c("B1", "B3"), each = 6),
                                       group = rep(c("G1", "G2", "G3"), 4))
    )
    ## Different assay name
    sceda <- SummarizedExperiment::SummarizedExperiment(
        assays = list(exprs = matrix(1:60, nrow = 5)),
        colData = S4Vectors::DataFrame(batch = rep(c("B1", "B2"), each = 6),
                                       group = rep(c("G1", "G2", "G3"), 4))
    )


    ## Fails with wrong input
    expect_error(getMatSubtractedBaseline(sce = 1, assayName = "counts",
                                          baselineGroup = "", sceFull = sce),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "", sceFull = 1),
                 "'sceFull' must be of class 'SummarizedExperiment'")
    expect_error(getMatSubtractedBaseline(sce = sceng, assayName = "counts",
                                          baselineGroup = "", sceFull = sce),
                 "all(!is.null(sce$group),", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = scenb, assayName = "counts",
                                          baselineGroup = "", sceFull = sce),
                 "all(!is.null(sce$group),", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "", sceFull = sceng),
                 "all(!is.null(sce$group),", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "", sceFull = scenb),
                 "all(!is.null(sce$group),", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = 1, sceFull = sce),
                 "'baselineGroup' must be of class 'character'")
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = c("G1", "G2"),
                                          sceFull = sce),
                 "'baselineGroup' must have length 1")
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "", sceFull = sce),
                 "All values in 'baselineGroup' must be one of")
    expect_error(getMatSubtractedBaseline(sce = scedb, assayName = "counts",
                                          baselineGroup = "G1", sceFull = sce),
                 "all(sce$batch %in% sceFull$batch", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "G1", sceFull = scedb),
                 "all(sce$batch %in% sceFull$batch", fixed = TRUE)
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = 1,
                                          baselineGroup = "G1", sceFull = sce),
                 "'assayName' must be of class 'character'")
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "missing",
                                          baselineGroup = "G1", sceFull = sce),
                 "All values in 'assayName' must be one of")
    expect_error(getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                          baselineGroup = "G1", sceFull = sceda),
                 "All values in 'assayName' must be one of")

    ## Works with correct arguments
    mat <- getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                    baselineGroup = "G1",
                                    sceFull = sce)
    expect_true(is(mat, "matrix"))
    expect_equal(nrow(mat), nrow(sce))
    expect_equal(ncol(mat), ncol(sce))
    expect_equal(rownames(mat), rownames(sce))
    expect_equal(colnames(mat), colnames(sce))
    for (i in seq_len(ncol(mat))) {
        expect_equal(mat[, i], assay(sce, "counts")[, i] -
                         rowMeans(assay(sce, "counts")[, sce$group == "G1" &
                                                           sce$batch == sce$batch[i], drop = FALSE]))
    }

    ## Subset as full
    scesub <- sce[, c(1:3, 7:9)]
    mat <- getMatSubtractedBaseline(sce = sce, assayName = "counts",
                                    baselineGroup = "G1",
                                    sceFull = scesub)
    expect_true(is(mat, "matrix"))
    expect_equal(nrow(mat), nrow(sce))
    expect_equal(ncol(mat), ncol(sce))
    expect_equal(rownames(mat), rownames(sce))
    expect_equal(colnames(mat), colnames(sce))
    for (i in seq_len(ncol(mat))) {
        expect_equal(mat[, i], assay(sce, "counts")[, i] -
                         rowMeans(assay(scesub, "counts")[, scesub$group == "G1" &
                                                              scesub$batch == sce$batch[i], drop = FALSE]))
    }

    ## Subset as target
    scesub <- sce[, c(1:3, 7:9)]
    mat <- getMatSubtractedBaseline(sce = scesub, assayName = "counts",
                                    baselineGroup = "G1",
                                    sceFull = sce)
    expect_true(is(mat, "matrix"))
    expect_equal(nrow(mat), nrow(scesub))
    expect_equal(ncol(mat), ncol(scesub))
    expect_equal(rownames(mat), rownames(scesub))
    expect_equal(colnames(mat), colnames(scesub))
    for (i in seq_len(ncol(mat))) {
        expect_equal(mat[, i], assay(scesub, "counts")[, i] -
                         rowMeans(assay(sce, "counts")[, sce$group == "G1" &
                                                           sce$batch == scesub$batch[i], drop = FALSE]))
    }

    ## Subset as both
    scesub <- sce[, c(1:3, 7:9)]
    mat <- getMatSubtractedBaseline(sce = scesub, assayName = "counts",
                                    baselineGroup = "G1",
                                    sceFull = scesub)
    expect_true(is(mat, "matrix"))
    expect_equal(nrow(mat), nrow(scesub))
    expect_equal(ncol(mat), ncol(scesub))
    expect_equal(rownames(mat), rownames(scesub))
    expect_equal(colnames(mat), colnames(scesub))
    for (i in seq_len(ncol(mat))) {
        expect_equal(mat[, i], assay(scesub, "counts")[, i] -
                         rowMeans(assay(scesub, "counts")[, scesub$group == "G1" &
                                                              scesub$batch == scesub$batch[i], drop = FALSE]))
    }
    for (i in which(scesub$group == "G1")) {
        expect_equal(mat[, i], rep(0, nrow(mat)), ignore_attr = TRUE)
    }

})
