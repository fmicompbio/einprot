test_that("runMaxQuantAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "MaxQuantAnalysis"
    args0 <- list(
        templateRmd = system.file("extdata/process_MaxQuant_template.Rmd",
                                  package = "einprot"),
        outputDir = outDir,
        outputBaseName = outBaseName,
        reportTitle = "MaxQuant LFQ data processing",
        reportAuthor = "",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "type1"),
        species = "mouse",
        mqFile = system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                             package = "einprot"),
        mqParameterFile = system.file("extdata", "mq_example", "1356_mqpar.xml",
                                      package = "einprot"),
        aName = "iBAQ",
        iColPattern = "^iBAQ\\\\.",
        sampleAnnot = data.frame(
            sample = c("Adnp_IP04", "Adnp_IP05",
                       "Adnp_IP06", "Chd4BF_IP07",
                       "Chd4BF_IP08", "Chd4BF_IP09",
                       "RBC_ctrl_IP01", "RBC_ctrl_IP02",
                       "RBC_ctrl_IP03"),
            group = c("Adnp", "Adnp", "Adnp", "Chd4BF", "Chd4BF",
                      "Chd4BF", "RBC_ctrl", "RBC_ctrl", "RBC_ctrl")),
        includeOnlySamples = "",
        excludeSamples = "",
        minScore = 10,
        minPeptides = 2,
        imputeMethod = "MinProb",
        mergeGroups = list(),
        comparisons = list(),
        ctrlGroup = "",
        allPairwiseComparisons = TRUE,
        normMethod = "none",
        stattest = "limma",
        minNbrValidValues = 2,
        minlFC = 0,
        nperm = 250,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        volcanoMaxFeatures = 25,
        volcanoS0 = 0.1,
        volcanoFeaturesToLabel = "",
        addInteractiveVolcanos = FALSE,
        complexFDRThr = 0.1,
        seed = 42,
        includeFeatureCollections = "complexes",
        customComplexes = list(),
        complexSpecies = "all",
        complexDbPath = system.file(
            "extdata", "complexes",
            "complexdb_einprot0.5.0_20220211_orthologs.rds",
            package = "einprot"),
        customYml = NULL,
        doRender = FALSE
    )

    ## Fail with wrong parameter values (essentially the same tests as for
    ## .checkArgumentsMaxQuant())
    ## --------------------------------------------------------------------- ##
    ## templateRmd
    args <- args0
    args$templateRmd <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'forceOverwrite' must have length 1")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'doRender' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'namesexperimentInfo' must not be NULL")

    ## species
    args <- args0
    args$species <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    ## mqFile
    args <- args0
    args$mqFile <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqFile' must be of class 'character'")
    args$mqFile <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqFile' must have length 1")
    args$mqFile <- "missing"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqFile' must point to an existing file")

    ## mqParameterFile
    args <- args0
    args$mqParameterFile <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqParameterFile' must be of class 'character'")
    args$mqParameterFile <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqParameterFile' must have length 1")
    args$mqParameterFile <- "missing"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mqParameterFile' must point to an existing file")

    ## aName
    args <- args0
    args$aName <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'aName' must be of class 'character'")
    args$aName <- c("name1", "name2")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'aName' must have length 1")

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("^LFQ\\\\.intensity\\\\.", "^iBAQ\\\\.")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'iColPattern' must have length 1")
    args$iColPattern <- c("^LFQ\\.intensity\\.")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'iColPattern' must be one of")

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$sample[1] <- paste0("iBAQ.",
                                         args$sampleAnnot$sample[1])
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "Please specify max one of includeOnlySamples")

    ## minScore
    args <- args0
    args$minScore <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minScore' must be of class 'numeric'")
    args$minScore <- c(1, 2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minScore' must have length 1")

    ## minPeptides
    args <- args0
    args$minPeptides <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minPeptides' must be of class 'numeric'")
    args$minPeptides <- c(1, 2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minPeptides' must have length 1")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'imputeMethod' must be one of")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4"))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "A given name can just be part of one merged group")

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "Each entry in 'comparisons' must have exactly two elements")

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'ctrlGroup' must be of class 'character'")
    # args$ctrlGroup <- c("name1", "name2")
    # expect_error(do.call(runMaxQuantAnalysis, args),
    #              "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'allPairwiseComparisons' must have length 1")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'normMethod' must be one of")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'stattest' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'namescustomComplexes' must not be NULL")

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'complexDbPath' must point to an existing file")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "'customYml' must point to an existing file")

    ## Test fix for backwards compatibility with mergeGroups/ctrlGroup
    args <- args0
    args$mergeGroups <- list(G1 = c("Adnp", "RBC_ctrl"))
    args$ctrlGroup <- c("Adnp", "RBC_ctrl")
    expect_error(do.call(runMaxQuantAnalysis, args),
                 "If 'mergeGroups' is specified, 'ctrlGroup' should not")

    args <- args0
    args$doRender <- FALSE
    args$ctrlGroup <- c("Adnp", "RBC_ctrl")
    res <- do.call(runMaxQuantAnalysis, args)
    expect_type(res, "character")
    tmp <- readLines(res)
    expect_true(length(
        grep('ctrlGroup <- \"Adnp.RBC_ctrl\"', tmp, fixed = TRUE)) > 0)
    expect_true(length(
        grep('mergeGroups <- list(Adnp.RBC_ctrl = c(\"Adnp\", \"RBC_ctrl\"))',
             tmp, fixed = TRUE)) > 0)


    ## Generate report
    ## --------------------------------------------------------------------- ##
    ## Without rendering
    args <- args0
    res <- do.call(runMaxQuantAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))

    ## Stop if forceOverwrite = FALSE
    args <- args0
    args$forceOverwrite <- FALSE
    expect_error(res <- do.call(runMaxQuantAnalysis, args))

    ## Message if forceOverwrite = TRUE
    args <- args0
    args$forceOverwrite <- TRUE
    expect_message(res <- do.call(runMaxQuantAnalysis, args),
                   "already exists but forceOverwrite = TRUE")


    ## With rendering
    skip_if(!capabilities()["X11"])
    args <- args0
    args$doRender <- TRUE
    res <- do.call(runMaxQuantAnalysis, args)

    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".html"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".html"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_iSEE.R"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_sce.rds"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_heatmap_centered.pdf"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_feature_info.txt"))))
})
