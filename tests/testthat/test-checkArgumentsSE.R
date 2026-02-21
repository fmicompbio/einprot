test_that("argument checking for SE works", {
    ## Working arguments
    args0 <- list(
        templateRmd = system.file("extdata", "process_basic_template.Rmd",
                                  package = "einprot"),
        outputDir = tempdir(),
        outputBaseName = "baseName",
        reportTitle = "reportTitle",
        reportAuthor = "reportAuthor",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "experiment",
                              "Sample name" = "example"),
        species = "roundworm",
        inFile = system.file("extdata", "mq_example", "1356_sce_init.rds",
                             package = "einprot"),
        summaryInfo = list(),
        aName = "LFQ.intensity",
        idCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                        combineWhen = "nonunique",
                                        splitSeparator = ";", joinSeparator = "."),
        labelCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                           combineWhen = "nonunique",
                                           splitSeparator = ";", joinSeparator = "."),
        geneIdCol = function(df) getFirstId(df, colName = "Gene.names",
                                            separator = ";"),
        proteinIdCol = "Majority.protein.IDs",
        stringIdCol = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                              combineWhen = "missing", splitSeparator = ";",
                                              joinSeparator = ".", makeUnique = FALSE),
        extraFeatureCols = NULL,
        filtersSE = einprotMQFilters,
        imputeMethod = "MinProb",
        assaysForExport = c("iBAQ", "Top3"),
        addAbundanceValues = TRUE,
        addHeatmaps = TRUE,
        mergeGroups = list(),
        comparisons = list(),
        ctrlGroup = "",
        allPairwiseComparisons = TRUE,
        singleFit = FALSE,
        subtractBaseline = FALSE,
        baselineGroup = "",
        normMethod = "none",
        spikeFeatures = NULL,
        stattest = "limma",
        minNbrValidValues = 2,
        minlFC = 0,
        samSignificance = TRUE,
        nperm = 100,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        volcanoMaxFeatures = 10,
        volcanoLabelSign = "both",
        volcanoS0 = 0.1,
        volcanoFeaturesToLabel = c(""),
        addInteractiveVolcanos = FALSE,
        interactiveDisplayColumns = NULL,
        interactiveGroupColumn = NULL,
        complexFDRThr = 0.1,
        maxNbrComplexesToPlot = Inf,
        maxComplexSimilarity = 1,
        seed = 123,
        includeFeatureCollections = "complexes",
        minSizeToKeepSet = 2,
        customComplexes = list(),
        complexSpecies = "all",
        complexDbPath = system.file(
            "extdata", "complexes",
            "complexdb_einprot0.5.0_20220323_orthologs.rds",
            package = "einprot"),
        stringVersion = "11.5",
        stringDir = "",
        linkTableColumns = c(),
        customYml = NULL,
        doRender = TRUE
    )

    ## Check that it works with these parameter values
    expect_null(do.call(.checkArgumentsSE, args0))

    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'forceOverwrite' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'namesexperimentInfo' must not be NULL")
    args$experimentInfo <- list()
    expect_null(do.call(.checkArgumentsSE, args))

    ## species
    args <- args0
    args$species <- 1
    expect_warning(do.call(.checkArgumentsSE, args),
                   "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(.checkArgumentsSE, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)
    args$species <- "Mouse"
    expect_null(do.call(.checkArgumentsSE, args))

    ## inFile
    args <- args0
    args$inFile <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'inFile' must be of class 'character'")
    args$inFile <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'inFile' must have length 1")
    args$inFile <- "missing"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'inFile' must point to an existing file")

    ## summaryInfo
    args <- args0
    args$summaryInfo <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'summaryInfo' must be of class 'list'")
    args$summaryInfo <- list("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'namessummaryInfo' must not be NULL")

    ## aName
    args <- args0
    args$aName <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'aName' must be of class 'character'")
    args$aName <- c("iBAQ", "Unique.peptides")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'aName' must have length 1")
    args$aName <- "missing"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'aName' must be one of")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'idCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'labelCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$geneIdCol <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'geneIdCol' must be of class 'character'")

    ## proteinIdCol
    args <- args0
    args$proteinIdCol <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'proteinIdCol' must be of class 'character'")

    ## stringIdCol
    args <- args0
    args$stringIdCol <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stringIdCol' must be of class 'character'")

    ## extraFeatureCols
    args <- args0
    args$extraFeatureCols <- function(df) df$einprotId
    expect_error(do.call(.checkArgumentsSE, args),
                 "'extraFeatureCols' must be of class 'list'")
    args$extraFeatureCols <- list(function(df) df$einprotId)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'namesextraFeatureCols' must not be NULL")
    args$extraFeatureCols <- list(newCol = 1)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'i' must be of class 'character'")
    args$extraFeatureCols <- list(newCol = function(x, y) "a")
    expect_error(do.call(.checkArgumentsSE, args),
                 "length(formals(i)) == 1 is not TRUE", fixed = TRUE)

    ## filtersSE
    args <- args0
    args$filtersSE <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'filtersSE' must be of class 'list'")
    args$filtersSE <- list(function(x) x)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'namesfiltersSE' must not be NULL")
    args$filtersSE <- list(f1 = function(x, y) x + y)
    expect_error(do.call(.checkArgumentsSE, args),
                 "is not TRUE")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'imputeMethod' must be one of")

    ## assaysForExport
    args <- args0
    args$assaysForExport <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'assaysForExport' must be of class 'character'")

    ## addAbundanceValues
    args <- args0
    args$addAbundanceValues <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addAbundanceValues' must be of class 'logical'")
    args$addAbundanceValues <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addAbundanceValues' must have length 1")

    ## addHeatmaps
    args <- args0
    args$addHeatmaps <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addHeatmaps' must be of class 'logical'")
    args$addHeatmaps <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addHeatmaps' must have length 1")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(.checkArgumentsSE, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(.checkArgumentsSE, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list()
    expect_null(do.call(.checkArgumentsSE, args))

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsSE, args),
                 "Each entry in 'comparisons' must have exactly two elements")
    args$comparisons <- list()
    expect_null(do.call(.checkArgumentsSE, args))

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'ctrlGroup' must be of class 'character'")
    args$ctrlGroup <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'baselineGroup' must be of class 'character'")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'normMethod' must be one of")

    ## spikeFeatures
    args <- args0
    args$spikeFeatures <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'spikeFeatures' must be of class 'character'")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'stattest' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## samSignificance
    args <- args0
    args$samSignificance <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'samSignificance' must have length 1")

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'interactiveGroupColumn' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## maxComplexSimilarity
    args <- args0
    args$maxComplexSimilarity <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'maxComplexSimilarity' must be of class 'numeric'")
    args$maxComplexSimilarity <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'maxComplexSimilarity' must have length 1")

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsSE, args),
                 "'namescustomComplexes' must not be NULL")
    args$customComplexes <- list()
    expect_null(do.call(.checkArgumentsSE, args))

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(.checkArgumentsSE, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'complexDbPath' must point to an existing file")

    ## stringVersion
    args <- args0
    args$stringVersion <- 11
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stringVersion' must be of class 'character'")
    args$stringVersion <- c("11.0", "11.5")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stringVersion' must have length 1")

    ## stringDir
    args <- args0
    args$stringDir <- 11
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stringDir' must be of class 'character'")
    args$stringDir <- c("", "")
    expect_error(do.call(.checkArgumentsSE, args),
                 "'stringDir' must have length 1")

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(.checkArgumentsSE, args),
                 "'customYml' must point to an existing file")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(.checkArgumentsSE, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsSE, args),
                 "'doRender' must have length 1")
})
