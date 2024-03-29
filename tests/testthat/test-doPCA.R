test_that("doPCA works", {
    args0 <- list(
        sce = sce_mq_final,
        assayName = "log2_iBAQ",
        ncomponents = 4,
        ntop = Inf,
        plotpairs = list(c(1, 2)),
        maxNGroups = 10,
        maxTextWidthBarplot = NULL,
        colourBy = "group",
        subset_row = NULL,
        scale = FALSE
    )

    ## -------------------------------------------------------------------------
    ## Fails with wrong arguments
    ## -------------------------------------------------------------------------
    ## sce
    args <- args0
    args$sce <- 1
    expect_error(do.call(doPCA, args),
                 "'sce' must be of class 'SingleCellExperiment'")

    ## assayName
    args <- args0
    args$assayName <- 1
    expect_error(do.call(doPCA, args),
                 "'assayName' must be of class 'character'")
    args <- args0
    args$assayName <- c("log2_iBAQ", "log2_iBAQ")
    expect_error(do.call(doPCA, args),
                 "'assayName' must have length 1")
    args <- args0
    args$assayName <- "missing"
    expect_error(do.call(doPCA, args),
                 "All values in 'assayName' must be one of")

    ## ncomponents
    args <- args0
    args$ncomponents <- "1"
    expect_error(do.call(doPCA, args),
                 "'ncomponents' must be of class 'numeric'")
    args <- args0
    args$ncomponents <- c(1, 2)
    expect_error(do.call(doPCA, args),
                 "'ncomponents' must have length 1")

    ## ntop
    args <- args0
    args$ntop <- "1"
    expect_error(do.call(doPCA, args),
                 "'ntop' must be of class 'numeric'")
    args <- args0
    args$ntop <- c(1, 2)
    expect_error(do.call(doPCA, args),
                 "'ntop' must have length 1")
    args <- args0
    args$ntop <- 0
    expect_error(do.call(doPCA, args),
                 "'ntop' must be within [1,Inf] (inclusive)", fixed = TRUE)

    ## plotpairs
    args <- args0
    args$plotpairs <- c(1, 2)
    expect_error(do.call(doPCA, args),
                 "'plotpairs' must be of class 'list'")
    args <- args0
    args$plotpairs <- list(c(1, 2), c(1, 100))
    expect_error(do.call(doPCA, args),
                 "'plotpairs' requests components that will not be extracted")
    args <- args0
    args$plotpairs <- list(c(1, 2), c(1, 2, 3))
    expect_error(do.call(doPCA, args),
                 "'elm' must have length 2")

    ## maxNGroups
    args <- args0
    args$maxNGroups <- "1"
    expect_error(do.call(doPCA, args),
                 "'maxNGroups' must be of class 'numeric'")
    args <- args0
    args$maxNGroups <- c(1, 2)
    expect_error(do.call(doPCA, args),
                 "'maxNGroups' must have length 1")

    ## maxTextWidthBarplot
    args <- args0
    args$maxTextWidthBarplot <- "1"
    expect_error(do.call(doPCA, args),
                 "'maxTextWidthBarplot' must be of class 'numeric'")
    args <- args0
    args$maxTextWidthBarplot <- c(1, 2)
    expect_error(do.call(doPCA, args),
                 "'maxTextWidthBarplot' must have length 1")

    ## colourBy
    args <- args0
    args$colourBy <- 1
    expect_error(do.call(doPCA, args),
                 "'colourBy' must be of class 'character'")
    args <- args0
    args$colourBy <- c("group", "group")
    expect_error(do.call(doPCA, args),
                 "'colourBy' must have length 1")
    args <- args0
    args$colourBy <- "missing"
    expect_error(do.call(doPCA, args),
                 "All values in 'colourBy' must be one of")

    ## scale
    args <- args0
    args$scale <- 1
    expect_error(do.call(doPCA, args),
                 "'scale' must be of class 'logical'")
    args$scale <- c(TRUE, FALSE)
    expect_error(do.call(doPCA, args),
                 "'scale' must have length 1")

    ## -------------------------------------------------------------------------
    ## Manual calculation of PCs
    ## -------------------------------------------------------------------------
    pca_scaleF <- prcomp(t(assay(args0$sce, "log2_iBAQ")),
                         center = TRUE, scale. = FALSE)
    pca_scaleT <- prcomp(t(assay(args0$sce, "log2_iBAQ")),
                         center = TRUE, scale. = TRUE)
    set.seed(123)
    subrow <- sample(seq_len(150), 30)
    pca_scaleF_sub <- prcomp(t(assay(args0$sce, "log2_iBAQ")[subrow, ]),
                             center = TRUE, scale. = FALSE)
    pca_scaleT_sub <- prcomp(t(assay(args0$sce, "log2_iBAQ")[subrow, ]),
                             center = TRUE, scale. = TRUE)


    ## -------------------------------------------------------------------------
    ## Works with correct arguments
    ## -------------------------------------------------------------------------
    ## Defaults
    res <- do.call(doPCA, args0)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args0$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")
    expect_equal(rowData(res$sce)[, "PCA_log2_iBAQ_PC1"],
                 pca_scaleF$rotation[, 1], ignore_attr = TRUE)
    expect_equal(reducedDim(res$sce, "PCA_log2_iBAQ")[, 1],
                 pca_scaleF$x[, 1], ignore_attr = TRUE)

    ## Scale
    args <- args0
    args$scale <- TRUE
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")
    expect_equal(rowData(res$sce)[, "PCA_log2_iBAQ_PC1"],
                 pca_scaleT$rotation[, 1], ignore_attr = TRUE)
    expect_equal(reducedDim(res$sce, "PCA_log2_iBAQ")[, 1],
                 pca_scaleT$x[, 1], ignore_attr = TRUE)

    ## No scale, subset
    args <- args0
    args$scale <- FALSE
    args$subset_row <- subrow
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2, na.rm = TRUE),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")
    expect_equal(rowData(res$sce)[subrow, "PCA_log2_iBAQ_PC2"],
                 pca_scaleF_sub$rotation[, 2], ignore_attr = TRUE)
    expect_equal(rowData(res$sce)[setdiff(seq_len(nrow(res$sce)), subrow),
                                  "PCA_log2_iBAQ_PC1"],
                 rep(NA_real_, nrow(res$sce) - length(subrow)))
    expect_equal(reducedDim(res$sce, "PCA_log2_iBAQ")[, 2],
                 pca_scaleF_sub$x[, 2], ignore_attr = TRUE)

    ## Scale, subset
    args <- args0
    args$scale <- TRUE
    args$subset_row <- subrow
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2, na.rm = TRUE),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")
    expect_equal(rowData(res$sce)[subrow, "PCA_log2_iBAQ_PC2"],
                 pca_scaleT_sub$rotation[, 2], ignore_attr = TRUE)
    expect_equal(rowData(res$sce)[setdiff(seq_len(nrow(res$sce)), subrow),
                                  "PCA_log2_iBAQ_PC1"],
                 rep(NA_real_, nrow(res$sce) - length(subrow)))
    expect_equal(reducedDim(res$sce, "PCA_log2_iBAQ")[, 2],
                 pca_scaleT_sub$x[, 2], ignore_attr = TRUE)

    ## No colour
    args <- args0
    args["colourBy"] <- list(NULL)
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 3L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")

    ## Colour by different column
    args <- args0
    args$colourBy <- "sample"
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "sample"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")

    ## Two pairs, one rowData column already exists
    args <- args0
    args$maxTextWidthBarplot <- 2
    SummarizedExperiment::rowData(args$sce)[, "PCA_log2_iBAQ_PC1"] <-
        rep(2, nrow(args$sce))
    args$plotpairs <- list(c(1, 2), c(3, 1))
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC1"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, c("PC1_2", "PC3_1"))
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_s3_class(res$plotcoord$PC3_1, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, c("PC1_2", "PC3_1"))
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotcombined$PC3_1, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")

    ## Try to extract too many components -> will be reduced
    args <- args0
    args$ncomponents <- 15
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 8L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:8) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC3"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")

    ## Reduce maxNGroups - no legend
    args <- args0
    args$maxNGroups <- 1
    res <- do.call(doPCA, args)
    expect_type(res, "list")
    expect_named(res, c("sce", "plotcoord", "plotcombined", "plotpairs"))
    expect_s4_class(res$sce, "SingleCellExperiment")
    expect_true("PCA_log2_iBAQ" %in% SingleCellExperiment::reducedDimNames(res$sce))
    expect_equal(ncol(SingleCellExperiment::reducedDim(res$sce, "PCA_log2_iBAQ")),
                 4L)
    expect_true(all(paste0("PCA_log2_iBAQ_PC", 1:4) %in%
                        colnames(SummarizedExperiment::rowData(res$sce))))
    expect_equal(sum(SummarizedExperiment::rowData(res$sce)[, "PCA_log2_iBAQ_PC3"] ^ 2),
                 1)
    expect_type(res$plotcoord, "list")
    expect_named(res$plotcoord, "PC1_2")
    expect_s3_class(res$plotcoord$PC1_2, "ggplot")
    expect_equal(nrow(res$plotcoord$PC1_2$data), ncol(args$sce))
    expect_equal(ncol(res$plotcoord$PC1_2$data), 4L)
    expect_equal(colnames(res$plotcoord$PC1_2$data),
                 c("sampleLabel", "PC1", "PC2", "group"))
    expect_type(res$plotcombined, "list")
    expect_named(res$plotcombined, "PC1_2")
    expect_s3_class(res$plotcombined$PC1_2, "ggplot")
    expect_s3_class(res$plotpairs, "ggmatrix")

})
