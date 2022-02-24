test_that("testing works", {
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
    set.seed(123)
    qft <- QFeatures::impute(qft, method = "MinProb", i = "log2_iBAQ")
    fcoll <- prepareFeatureCollections(
        qft = qft, idCol = "Gene.names",
        includeFeatureCollections = "complexes",
        complexDbPath = system.file("extdata", "complexes",
                                    "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                    package = "einprot"),
        speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
        customComplexes = list(), minSizeToKeep = 2)

    ## Fail with wrong arguments
    ## --------------------------------------------------------------------- ##
    args0 <- list(
        qft = qft,
        comparison = c("Adnp", "RBC_ctrl"),
        testType = "limma",
        assayForTests = "log2_iBAQ",
        assayImputation = "imputed_iBAQ",
        minNbrValidValues = 2,
        minlFC = 0,
        featureCollections = fcoll,
        complexFDRThr = 0.1,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        baseFileName = NULL,
        seed = 123,
        nperm = 25,
        volcanoS0 = 0.1,
        addiBAQvalues = TRUE,
        iColPattern = "^iBAQ\\.",
        aName = "iBAQ"
    )

    ## qft
    args <- args0
    args$qft <- 1
    expect_error(do.call(runTest, args),
                 "'qft' must be of class 'QFeatures'")
    args <- args0
    args$qft$group <- as.numeric(as.factor(args$qft$group))
    expect_error(do.call(runTest, args),
                 "'$qftgroup' must be of class 'character'", fixed = TRUE)

    ## comparison
    args <- args0
    args$comparison <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'comparison' must be of class 'character'")
    args$comparison <- "Adnp"
    expect_error(do.call(runTest, args),
                 "'comparison' must have length 2")
    args$comparison <- c("missing", "Adnp")
    expect_error(do.call(runTest, args),
                 "All values in 'comparison' must be one of")
    args$comparison <- list(c("Adnp", "RBC_ctrl"))
    expect_error(do.call(runTest, args),
                 "'comparison' must be of class 'character'")

    ## testType
    args <- args0
    args$testType <- 1
    expect_error(do.call(runTest, args),
                 "'testType' must be of class 'character'")
    args$testType <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'testType' must be one of")
    args$testType <- c("ttest", "limma")
    expect_error(do.call(runTest, args),
                 "'testType' must have length 1")

    ## assayForTests
    args <- args0
    args$assayForTests <- 1
    expect_error(do.call(runTest, args),
                 "'assayForTests' must be of class 'character'")
    args$assayForTests <- c("log2_iBAQ", "imputed_iBAQ")
    expect_error(do.call(runTest, args),
                 "'assayForTests' must have length 1")
    args$assayForTests <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'assayForTests' must be one of")

    ## assayImputation
    args <- args0
    args$assayImputation <- 1
    expect_error(do.call(runTest, args),
                 "'assayImputation' must be of class 'character'")
    args$assayImputation <- c("log2_iBAQ", "imputed_iBAQ")
    expect_error(do.call(runTest, args),
                 "'assayImputation' must have length 1")
    args$assayImputation <- "missing"
    expect_error(do.call(runTest, args),
                 "All values in 'assayImputation' must be one of")

    ## minnbrValidValues
    args <- args0
    args$minNbrValidValues <- "2"
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runTest, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "2"
    expect_error(do.call(runTest, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runTest, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## featureCollections
    args <- args0
    args$featureCollections <- "2"
    expect_error(do.call(runTest, args),
                 "'featureCollections' must be of class 'list'")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "2"
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$complexFDRThr <- 2
    expect_error(do.call(runTest, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$volcanoAdjPvalThr <- 2
    expect_error(do.call(runTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(runTest, args),
                 "'baseFileName' must be of class 'character'")
    args$baseFileName <- c("pat1", "pat2")
    expect_error(do.call(runTest, args),
                 "'baseFileName' must have length 1")

    ## seed
    args <- args0
    args$seed <- "2"
    expect_error(do.call(runTest, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runTest, args),
                 "'seed' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## nperm
    args <- args0
    args$testType <- "ttest"
    args$nperm <- "2"
    expect_error(do.call(runTest, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runTest, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoS0
    args <- args0
    args$testType <- "ttest"
    args$volcanoS0 <- "2"
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(1, 2)
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runTest, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## addiBAQvalues
    args <- args0
    args$addiBAQvalues <- "2"
    expect_error(do.call(runTest, args),
                 "'addiBAQvalues' must be of class 'logical'")
    args$addiBAQvalues <- c(TRUE, FALSE)
    expect_error(do.call(runTest, args),
                 "'addiBAQvalues' must have length 1")

    ## iColPattern
    args <- args0
    args$addiBAQvalues <- TRUE
    args$iColPattern <- 1
    expect_error(do.call(runTest, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("pat1", "pat2")
    expect_error(do.call(runTest, args),
                 "'iColPattern' must have length 1")

    ## aName
    args <- args0
    args$addiBAQvalues <- TRUE
    args$aName <- 1
    expect_error(do.call(runTest, args),
                 "'aName' must be of class 'character'")
    args$aName <- c("pat1", "pat2")
    expect_error(do.call(runTest, args),
                 "'aName' must have length 1")

    ## Works with correct arguments
    ## --------------------------------------------------------------------- ##
    out <- do.call(runTest, args0)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(qft[[1]]))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 0)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))

    ## t-test, write results to file
    args <- args0
    args$testType <- "ttest"
    args$baseFileName <- tempfile()
    args$addiBAQvalues <- FALSE
    out <- do.call(runTest, args)
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp.txt")))
    expect_true(file.exists(paste0(args$baseFileName, "_testres_RBC_ctrl_vs_Adnp_camera_complexes.txt")))
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_named(out$curveparam, c("x", "ta", "s0", "df"))
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_false("iBAQ.Adnp_IP04" %in% colnames(out$res))
    expect_equal(out$res$pid, rownames(qft[[1]]))
    expect_equal(out$plotnote, "")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, t-test")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))

    ## With batch effect
    args <- args0
    args$qft$batch <- c("B1", "B2", "B1", "B2", "B1", "B2", "B1", "B2", "B1")
    out <- do.call(runTest, args0)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("res", "plotnote", "plottitle", "plotsubtitle",
                        "featureCollections", "curveparam"))
    expect_s3_class(out$res, "data.frame")
    expect_type(out$plotnote, "character")
    expect_type(out$plottitle, "character")
    expect_type(out$featureCollections, "list")
    expect_type(out$curveparam, "list")
    expect_equal(nrow(out$res), 70)
    expect_true(all(c("adj.P.Val", "iBAQ.Adnp_IP04",
                      "showInVolcano", "IDsForSTRING") %in% colnames(out$res)))
    expect_equal(out$res$pid, rownames(qft[[1]]))
    expect_equal(substr(out$plotnote, 1, 8), "df.prior")
    expect_equal(out$plottitle, "RBC_ctrl vs Adnp, limma treat (H0: |log2FC| <= 0)")
    expect_s4_class(out$featureCollections$complexes, "CharacterList")
    expect_s4_class(S4Vectors::mcols(out$featureCollections$complexes), "DFrame")
    expect_true("RBC_ctrl_vs_Adnp_FDR" %in%
                    colnames(S4Vectors::mcols(out$featureCollections$complexes)))

})
