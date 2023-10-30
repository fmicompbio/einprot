test_that("runDIANNAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "DIANNAnalysis"
    args0 <- list(
        templateRmd = system.file("extdata", "process_basic_template.Rmd",
                                  package = "einprot"),
        outputDir = outDir,
        outputBaseName = outBaseName,
        reportTitle = "reportTitle",
        reportAuthor = "reportAuthor",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "experiment",
                              "Sample name" = "example"),
        species = "roundworm",
        diannFile = system.file("extdata", "diann_example",
                                "PXD028735.pg_matrix.tsv", package = "einprot"),
        diannFileType = "pg_matrix",
        outLevel = "pg",
        diannLogFile = system.file("extdata", "diann_example",
                                   "diann-output.log.txt", package = "einprot"),
        aName = "MaxLFQ",
        idCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                            "Majority.protein.IDs")),
        labelCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                               "Majority.protein.IDs")),
        geneIdCol = function(df) getFirstId(df, colName = "Gene.names"),
        proteinIdCol = "Majority.protein.IDs",
        stringIdCol = function(df) combineIds(df, combineCols = c("Gene.names",
                                                                  "Majority.protein.IDs"),
                                              combineWhen = "missing",
                                              makeUnique = FALSE),
        sampleAnnot = data.frame(
            sample = c("LFQ_Orbitrap_AIF_Condition_A_Sample_Beta_01",
                       "LFQ_Orbitrap_AIF_Condition_A_Sample_Gamma_01",
                       "LFQ_Orbitrap_AIF_Condition_A_Sample_Alpha_01",
                       "LFQ_Orbitrap_AIF_Condition_B_Sample_Beta_01",
                       "LFQ_Orbitrap_AIF_Condition_B_Sample_Gamma_01",
                       "LFQ_Orbitrap_AIF_Condition_B_Sample_Alpha_01"),
            group = c("A", "A", "A", "B", "B", "B")),
        includeOnlySamples = "",
        excludeSamples = "",
        minScore = 10,
        minPeptides = 2,
        imputeMethod = "MinProb",
        assaysForExport = NULL,
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
        doRender = FALSE
    )

    ## Fail with wrong parameter values (essentially the same tests as for
    ## .checkArgumentsDIANN())
    ## -------------------------------------------------------------------------
    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'forceOverwrite' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'namesexperimentInfo' must not be NULL")

    ## species
    args <- args0
    args$species <- 1
    expect_warning(do.call(runDIANNAnalysis, args),
                   "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(runDIANNAnalysis, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    ## diannFile
    args <- args0
    args$diannFile <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannFile' must be of class 'character'")
    args$diannFile <- "missing"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannFile' must point to an existing file")

    ## diannFileType
    args <- args0
    args$diannFileType <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannFileType' must be of class 'character'")
    args$diannFileType <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "values in 'diannFileType' must be one of")
    args$diannFileType = c("pr_matrix", "pg_matrix")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannFileType' must have length 1")

    ## outLevel
    args <- args0
    args$outLevel <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outLevel' must be of class 'character'")
    args$outLevel <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "values in 'outLevel' must be one of")
    args$outLevel = c("pr", "pg")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'outLevel' must have length 1")
    args$diannFileType <- "pg_matrix"
    args$outLevel <- "pr"
    expect_error(do.call(runDIANNAnalysis, args),
                 "Can't get pr information from pg_matrix")

    ## diannLogFile
    args <- args0
    args$diannLogFile <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannLogFile' must be of class 'character'")
    args$diannLogFile <- "missing"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'diannLogFile' must point to an existing file")

    ## aName
    args <- args0
    args$aName <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'aName' must be of class 'character'")
    args$aName <- NULL
    args <- c(args, list(aName = NULL))
    expect_error(do.call(runDIANNAnalysis, args),
                 "'aName' must not be NULL")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'idCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'labelCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$geneIdCol <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'geneIdCol' must be of class 'character'")

    ## proteinIdCol
    args <- args0
    args$proteinIdCol <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'proteinIdCol' must be of class 'character'")

    ## stringIdCol
    args <- args0
    args$stringIdCol <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stringIdCol' must be of class 'character'")

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(runDIANNAnalysis, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(runDIANNAnalysis, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(runDIANNAnalysis, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    # args$sampleAnnot <- args0$sampleAnnot
    # args$sampleAnnot$sample[1] <- paste0("iBAQ.",
    #                                      args$sampleAnnot$sample[1])
    # expect_error(do.call(runDIANNAnalysis, args),
    #              "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(runDIANNAnalysis, args),
                 "Please specify max one of includeOnlySamples")

    ## minScore
    args <- args0
    args$minScore <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minScore' must be of class 'numeric'")
    args$minScore <- c(1, 2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minScore' must have length 1")

    ## minPeptides
    args <- args0
    args$minPeptides <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minPeptides' must be of class 'numeric'")
    args$minPeptides <- c(1, 2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minPeptides' must have length 1")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'imputeMethod' must be one of")

    ## assaysForExport
    args <- args0
    args$assaysForExport <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'assaysForExport' must be of class 'character'")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(runDIANNAnalysis, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(runDIANNAnalysis, args),
                 "'mergeGroups' must be a named list")

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(runDIANNAnalysis, args),
                 "Each entry in 'comparisons' must have exactly two elements")

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'ctrlGroup' must be of class 'character'")
    # args$ctrlGroup <- c("name1", "name2")
    # expect_error(do.call(runDIANNAnalysis, args),
    #              "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'baselineGroup' must be of class 'character'")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'normMethod' must be one of")

    ## spikeFeatures
    args <- args0
    args$spikeFeatures <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'spikeFeatures' must be of class 'character'")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'stattest' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## samSignificance
    args <- args0
    args$samSignificance <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'samSignificance' must have length 1")

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'interactiveGroupColumn' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(runDIANNAnalysis, args),
                 "'namescustomComplexes' must not be NULL")

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(runDIANNAnalysis, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'complexDbPath' must point to an existing file")

    ## stringVersion
    args <- args0
    args$stringVersion <- 11
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stringVersion' must be of class 'character'")
    args$stringVersion <- c("11.0", "11.5")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stringVersion' must have length 1")

    ## stringDir
    args <- args0
    args$stringDir <- 11
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stringDir' must be of class 'character'")
    args$stringDir <- c("", "")
    expect_error(do.call(runDIANNAnalysis, args),
                 "'stringDir' must have length 1")

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(runDIANNAnalysis, args),
                 "'customYml' must point to an existing file")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(runDIANNAnalysis, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(runDIANNAnalysis, args),
                 "'doRender' must have length 1")

    ## Test fix for backwards compatibility with mergeGroups/ctrlGroup
    args <- args0
    args$mergeGroups <- list(G1 = c("A", "B"))
    args$ctrlGroup <- c("A", "B")
    expect_error(do.call(runDIANNAnalysis, args),
                 "If 'mergeGroups' is specified, 'ctrlGroup' should not")

    ## Generate report
    ## -------------------------------------------------------------------------
    ## Without rendering
    args <- args0
    res <- do.call(runDIANNAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))

    ## Stop if forceOverwrite = FALSE
    args <- args0
    args$forceOverwrite <- FALSE
    expect_error(res <- do.call(runDIANNAnalysis, args))

    ## Message if forceOverwrite = TRUE
    args <- args0
    args$forceOverwrite <- TRUE
    expect_message(res <- do.call(runDIANNAnalysis, args),
                   "already exists but forceOverwrite = TRUE")

    ## In new, non-existing directory and with custom yml
    args <- args0
    args$outputDir <- file.path(outDir, "new_diann_dir")
    args$customYml <- system.file("extdata", "custom.yml",
                                  package = "einprot")
    res <- do.call(runDIANNAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(args$outputDir, paste0(outBaseName, ".Rmd"))))
    tmp <- readLines(file.path(args$outputDir, paste0(outBaseName, ".Rmd")))
    # expect_true(grepl("theme: journal", tmp[4]))

    ## With rendering
    skip_if(!capabilities()["X11"])
    args <- args0
    args$doRender <- TRUE
    res <- do.call(runDIANNAnalysis, args)

    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".html"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".html"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_iSEE.R"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_sce.rds"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_heatmap_centered.pdf"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_feature_info.txt"))))


})
