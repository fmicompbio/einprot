test_that("argument checking for PD-TMT works", {
    tmpd <- tempdir()
    write.table(data.frame(Abundance.F352.Sample.S1 = 1:3, Abundance.F543.Sample.S2 = 4:6),
                file = file.path(tmpd, "tmp_Proteins.txt"), sep = "\t")
    write.table(1:3, file = file.path(tmpd, "tmp_InputFiles.txt"))
    write.table(1:3, file = file.path(tmpd, "tmp_StudyInformation.txt"))
    write.table(1:3, file = file.path(tmpd, "tmp.pdAnalysis"))

    ## Working arguments
    args0 <- list(
        templateRmd = system.file("extdata", "process_PD_TMT_template.Rmd",
                                  package = "einprot"),
        outputDir = tempdir(),
        outputBaseName = "baseName",
        reportTitle = "reportTitle",
        reportAuthor = "reportAuthor",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "experiment",
                              "Sample name" = "example"),
        species = "roundworm",
        pdOutputFolder = tmpd,
        pdResultName = "tmp",
        pdAnalysisFile = file.path(tmpd, "tmp.pdAnalysis"),
        iColPattern = "^Abundance\\\\.F.+\\\\.Sample\\\\.",
        sampleAnnot = data.frame(sample = c("S1", "S2"),
                                 group = c("G1", "G2")),
        includeOnlySamples = "",
        excludeSamples = "",
        minScore = 2,
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
        nperm = 100,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        volcanoMaxFeatures = 10,
        volcanoS0 = 0.1,
        volcanoFeaturesToLabel = c(""),
        addInteractiveVolcanos = FALSE,
        complexFDRThr = 0.1,
        maxNbrComplexesToPlot = 10,
        seed = 123,
        includeFeatureCollections = "complexes",
        minSizeToKeepSet = 2,
        customComplexes = list(),
        complexSpecies = "all",
        complexDbPath = system.file("extdata", "complexes",
                                    "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                    package = "einprot"),
        customYml = NULL,
        doRender = TRUE,
        generateQCPlot = TRUE
    )

    ## Check that it works with these parameter values
    expect_null(do.call(.checkArgumentsPDTMT, args0))

    ## templateRmd
    args <- args0
    args$templateRmd <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'forceOverwrite' must have length 1")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'doRender' must have length 1")

    ## generateQCPlot
    args <- args0
    args$generateQCPlot <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'generateQCPlot' must be of class 'logical'")
    args$generateQCPlot <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'generateQCPlot' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'namesexperimentInfo' must not be NULL")
    args$experimentInfo <- list()
    expect_null(do.call(.checkArgumentsPDTMT, args))

    ## species
    args <- args0
    args$species <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)
    args$species <- "Mouse"
    expect_null(do.call(.checkArgumentsPDTMT, args))

    ## pdOutputFolder
    args <- args0
    args$pdOutputFolder <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdOutputFolder' must be of class 'character'")
    args$pdOutputFolder <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdOutputFolder' must have length 1")

    ## pdResultName
    args <- args0
    args$pdResultName <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdResultName' must be of class 'character'")
    args$pdResultName <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdResultName' must have length 1")

    args <- args0
    args$pdResultName <- "missing"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "missing_Proteins.txt doesn't exist")

    ## pdAnalysisFile
    args <- args0
    args$pdAnalysisFile <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdAnalysisFile' must be of class 'character'")
    args$pdAnalysisFile <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdAnalysisFile' must have length 1")
    args$pdAnalysisFile <- "missing"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'pdAnalysisFile' must point to an existing file")

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("^Abundance\\\\.F.+\\\\.Sample\\\\.", "^Abundance\\\\.F.+\\\\.Sample\\\\.")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'iColPattern' must have length 1")
    args$iColPattern <- c("^LFQ\\.intensity\\.")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'iColPattern' must be one of")

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$sample[1] <- paste0("iBAQ.",
                                         args$sampleAnnot$sample[1])
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "Please specify max one of includeOnlySamples")

    ## minScore
    args <- args0
    args$minScore <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minScore' must be of class 'numeric'")
    args$minScore <- c(1, 2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minScore' must have length 1")

    ## minPeptides
    args <- args0
    args$minPeptides <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minPeptides' must be of class 'numeric'")
    args$minPeptides <- c(1, 2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minPeptides' must have length 1")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'imputeMethod' must be one of")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4"))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "A given name can just be part of one merged group")
    args$mergeGroups <- list()
    expect_null(do.call(.checkArgumentsPDTMT, args))

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "Each entry in 'comparisons' must have exactly two elements")
    args$comparisons <- list()
    expect_null(do.call(.checkArgumentsPDTMT, args))

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'ctrlGroup' must be of class 'character'")
    args$ctrlGroup <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'allPairwiseComparisons' must have length 1")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'normMethod' must be one of")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'stattest' must be one of")
    args$stattest <- "ttest"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'ttest' is currently not supported")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'namescustomComplexes' must not be NULL")
    args$customComplexes <- list()
    expect_null(do.call(.checkArgumentsPDTMT, args))

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'complexDbPath' must point to an existing file")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(.checkArgumentsPDTMT, args),
                 "'customYml' must point to an existing file")
})