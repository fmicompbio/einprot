test_that("argument checking for FP works", {
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
        fragpipeDir = system.file("extdata", "fp_example",
                                  package = "einprot"),
        idCol = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID")),
        labelCol = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID")),
        geneIdCol = function(df) getFirstId(df, colName = "Gene"),
        proteinIdCol = function(df) getFirstId(df, colName = "Protein.ID"),
        stringIdCol = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID"),
                                              combineWhen = "missing", makeUnique = FALSE),
        iColPattern = "\\.MaxLFQ\\.Intensity$",
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
        assaysForExport = NULL,
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
    expect_null(do.call(.checkArgumentsFragPipe, args0))

    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'forceOverwrite' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'namesexperimentInfo' must not be NULL")
    args$experimentInfo <- list()
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## species
    args <- args0
    args$species <- 1
    expect_warning(do.call(.checkArgumentsFragPipe, args),
                   "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)
    args$species <- "Mouse"
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## fragpipeDir
    args <- args0
    args$fragpipeDir <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'fragpipeDir' must be of class 'character'")
    args$fragpipeDir <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'fragpipeDir' must have length 1")
    args$fragpipeDir <- "missing"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "The file missing/combined_protein.tsv doesn't exist")

    ## Multiple log files
    dir.create(file.path(tempdir(), "fragpipe_dir_temp"))
    file.copy(system.file("extdata", "fp_example", package = "einprot"),
              file.path(tempdir(), "fragpipe_dir_temp"), recursive = TRUE)
    file.copy(file.path(tempdir(), "fragpipe_dir_temp", "fp_example",
                        "log_2023-04-12_20-12-46.txt"),
              file.path(tempdir(), "fragpipe_dir_temp", "fp_example",
                        "log_2024-04-12_20-12-46.txt"))
    args <- args0
    args$fragpipeDir <- file.path(tempdir(), "fragpipe_dir_temp", "fp_example")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "There are more than one log file")

    file.copy(file.path(tempdir(), "fragpipe_dir_temp", "fp_example",
                        "fragpipe.workflow"),
              file.path(tempdir(), "fragpipe_dir_temp", "fp_example",
                        "fragpipe2.workflow"))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "There are more than one workflow file")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'idCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'labelCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$geneIdCol <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'geneIdCol' must be of class 'character'")

    ## proteinIdCol
    args <- args0
    args$proteinIdCol <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'proteinIdCol' must be of class 'character'")
    args <- args0
    args$proteinIdCol <- "Protein.ID"
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## stringIdCol
    args <- args0
    args$stringIdCol <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stringIdCol' must be of class 'character'")

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("\\.MaxLFQ\\.Intensity$", "\\.MaxLFQ\\.Intensity$")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'iColPattern' must have length 1")
    args$iColPattern <- c("\\\\.MaxLFQ\\\\.Intensity$")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'iColPattern' must be one of")
    args$iColPattern <- c("^iBAQ\\.")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'iColPattern' must be one of")
    ## Works without escaped periods
    args <- args0
    args$iColPattern <- ".MaxLFQ.Intensity$"
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$sample[1] <- paste0("iBAQ.",
                                         args$sampleAnnot$sample[1])
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "Please specify max one of includeOnlySamples")

    ## minScore
    args <- args0
    args$minScore <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minScore' must be of class 'numeric'")
    args$minScore <- c(1, 2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minScore' must have length 1")

    ## minPeptides
    args <- args0
    args$minPeptides <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minPeptides' must be of class 'numeric'")
    args$minPeptides <- c(1, 2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minPeptides' must have length 1")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'imputeMethod' must be one of")

    ## assaysForExport
    args <- args0
    args$assaysForExport <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'assaysForExport' must be of class 'character'")

    ## addHeatmaps
    args <- args0
    args$addHeatmaps <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'addHeatmaps' must be of class 'logical'")
    args$addHeatmaps <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'addHeatmaps' must have length 1")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list()
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "Each entry in 'comparisons' must have exactly two elements")
    args$comparisons <- list()
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'ctrlGroup' must be of class 'character'")
    args$ctrlGroup <- c("name1", "name2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'baselineGroup' must be of class 'character'")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'normMethod' must be one of")

    ## spikeFeatures
    args <- args0
    args$spikeFeatures <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'spikeFeatures' must be of class 'character'")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'stattest' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## samSignificance
    args <- args0
    args$samSignificance <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'samSignificance' must have length 1")

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'interactiveGroupColumn' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'namescustomComplexes' must not be NULL")
    args$customComplexes <- list()
    expect_null(do.call(.checkArgumentsFragPipe, args))

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'complexDbPath' must point to an existing file")

    ## stringVersion
    args <- args0
    args$stringVersion <- 11
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stringVersion' must be of class 'character'")
    args$stringVersion <- c("11.0", "11.5")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stringVersion' must have length 1")

    ## stringDir
    args <- args0
    args$stringDir <- 11
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stringDir' must be of class 'character'")
    args$stringDir <- c("", "")
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'stringDir' must have length 1")

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'customYml' must point to an existing file")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(.checkArgumentsFragPipe, args),
                 "'doRender' must have length 1")
})
