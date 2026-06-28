test_that("runSpectronautAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "SpectronautAnalysis"
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
        spectronautFile = system.file("extdata", "spectronaut_example",
                                      "SN_unit_test_subset.tsv", package = "einprot"),
        spectronautFileType = "pg_pivot",
        outLevel = "pg",
        spectronautLogFile = system.file(
            "extdata", "spectronaut_example",
            "PG.QValRUN_005_PGPEP_075_Minlog2PrecInt_0_reduced_pivot_filtered.setup.txt",
            package = "einprot"),
        aName = "PG.Quantity",
        idCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                            "PG.ProteinGroups")),
        labelCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                               "PG.ProteinGroups")),
        geneIdCol = function(df) getFirstId(df, colName = "PG.Genes"),
        proteinIdCol = "PG.ProteinGroups",
        stringIdCol = function(df) combineIds(df, combineCols = c("PG.Genes",
                                                                  "PG.ProteinGroups"),
                                              combineWhen = "missing",
                                              makeUnique = FALSE),
        extraFeatureCols = NULL,
        iColPattern = ".PG.Quantity$",
        sampleAnnot = data.frame(
            sample = c("LFQ_24min_500ng_2Th_3ms_A_01",
                       "LFQ_24min_500ng_2Th_3ms_A_02",
                       "LFQ_24min_500ng_2Th_3ms_A_03",
                       "LFQ_24min_500ng_2Th_3ms_B_01",
                       "LFQ_24min_500ng_2Th_3ms_B_02",
                       "LFQ_24min_500ng_2Th_3ms_B_03"),
            group = c("A", "A", "A", "B", "B", "B")),
        includeOnlySamples = "",
        excludeSamples = "",
        filtersDF = list(),
        filtersSE = list(),
        imputeMethod = "MinProb",
        imputeArgs = list(),
        assaysForExport = NULL,
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
        doRender = FALSE
    )

    ## Fail with wrong parameter values (essentially the same tests as for
    ## .checkArgumentsSpectronaut())
    ## -------------------------------------------------------------------------
    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'forceOverwrite' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'namesexperimentInfo' must not be NULL")

    ## species
    args <- args0
    args$species <- c("Mouse", "Human")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    ## spectronautFile
    args <- args0
    args$spectronautFile <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautFile' must be of class 'character'")
    args$spectronautFile <- "missing"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautFile' must point to an existing file")

    ## spectronautFileType
    args <- args0
    args$spectronautFileType <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautFileType' must be of class 'character'")
    args$spectronautFileType <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "values in 'spectronautFileType' must be one of")
    args$spectronautFileType = c("pg_pivot", "pg_pivot")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautFileType' must have length 1")

    ## outLevel
    args <- args0
    args$outLevel <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outLevel' must be of class 'character'")
    args$outLevel <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "values in 'outLevel' must be one of")
    args$outLevel = c("pg", "pg")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'outLevel' must have length 1")

    ## spectronautLogFile
    args <- args0
    args$spectronautLogFile <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautLogFile' must be of class 'character'")
    args$spectronautLogFile <- "missing"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spectronautLogFile' must point to an existing file")

    ## aName
    # args <- args0
    # args$aName <- 1
    # expect_error(do.call(runSpectronautAnalysis, args),
    #              "'aName' must be of class 'character'")
    # args$aName <- NULL
    # args <- c(args, list(aName = NULL))
    # expect_error(do.call(runSpectronautAnalysis, args),
    #              "'aName' must not be NULL")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'idCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'labelCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$geneIdCol <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'geneIdCol' must be of class 'character'")

    ## proteinIdCol
    args <- args0
    args$proteinIdCol <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'proteinIdCol' must be of class 'character'")

    ## stringIdCol
    args <- args0
    args$stringIdCol <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stringIdCol' must be of class 'character'")

    ## extraFeatureCols
    args <- args0
    args$extraFeatureCols <- function(df) df$einprotId
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'extraFeatureCols' must be of class 'list'")
    args$extraFeatureCols <- list(function(df) df$einprotId)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'namesextraFeatureCols' must not be NULL")
    args$extraFeatureCols <- list(newCol = 1)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'i' must be of class 'character'")
    args$extraFeatureCols <- list(newCol = function(x, y) "a")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "length(formals(i)) == 1 is not TRUE", fixed = TRUE)

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("\\.PG\\.Quantity$", "\\.PG\\.Quantity$")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'iColPattern' must have length 1")
    # args$iColPattern <- c("\\\\.PG\\\\.Quantity$")
    # expect_error(do.call(runSpectronautAnalysis, args),
    #              "All values in 'iColPattern' must be one of")

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    # args$sampleAnnot <- args0$sampleAnnot
    # args$sampleAnnot$sample[1] <- paste0("iBAQ.",
    #                                      args$sampleAnnot$sample[1])
    # expect_error(do.call(runSpectronautAnalysis, args),
    #              "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "Please specify max one of includeOnlySamples")

    ## filtersDF
    args <- args0
    args$filtersDF <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'filtersDF' must be of class 'list'")
    args$filtersDF <- list(function(x) x)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'namesfiltersDF' must not be NULL")
    args$filtersDF <- list(f1 = function(x, y) x + y)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "is not TRUE")

    ## filtersSE
    args <- args0
    args$filtersSE <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'filtersSE' must be of class 'list'")
    args$filtersSE <- list(function(x) x)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'namesfiltersSE' must not be NULL")
    args$filtersSE <- list(f1 = function(x, y) x + y)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "is not TRUE")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'imputeMethod' must be one of")

    ## imputeArgs
    args <- args0
    args$imputeArgs <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'imputeArgs' must be of class 'list'")

    ## assaysForExport
    args <- args0
    args$assaysForExport <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'assaysForExport' must be of class 'character'")

    ## addAbundanceValues
    args <- args0
    args$addAbundanceValues <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addAbundanceValues' must be of class 'logical'")
    args$addAbundanceValues <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addAbundanceValues' must have length 1")

    ## addHeatmaps
    args <- args0
    args$addHeatmaps <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addHeatmaps' must be of class 'logical'")
    args$addHeatmaps <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addHeatmaps' must have length 1")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'mergeGroups' must be a named list")

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(runSpectronautAnalysis, args),
                 "Each entry in 'comparisons' must have exactly two elements")

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'ctrlGroup' must be of class 'character'")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'baselineGroup' must be of class 'character'")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'normMethod' must be one of")

    ## spikeFeatures
    args <- args0
    args$spikeFeatures <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'spikeFeatures' must be of class 'character'")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'stattest' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## samSignificance
    args <- args0
    args$samSignificance <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'samSignificance' must have length 1")

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'addInteractiveVolcanos' must have length 1")

    ## interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'interactiveGroupColumn' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## maxComplexSimilarity
    args <- args0
    args$maxComplexSimilarity <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'maxComplexSimilarity' must be of class 'numeric'")
    args$maxComplexSimilarity <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'maxComplexSimilarity' must have length 1")

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'namescustomComplexes' must not be NULL")

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'complexDbPath' must point to an existing file")

    ## stringVersion
    args <- args0
    args$stringVersion <- 11
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stringVersion' must be of class 'character'")
    args$stringVersion <- c("11.0", "11.5")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stringVersion' must have length 1")

    ## stringDir
    args <- args0
    args$stringDir <- 11
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stringDir' must be of class 'character'")
    args$stringDir <- c("", "")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'stringDir' must have length 1")

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'customYml' must point to an existing file")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(runSpectronautAnalysis, args),
                 "'doRender' must have length 1")


    ## Test fix for backwards compatibility with mergeGroups/ctrlGroup
    args <- args0
    args$mergeGroups <- list(G1 = c("A", "B"))
    args$ctrlGroup <- c("A", "B")
    expect_error(do.call(runSpectronautAnalysis, args),
                 "If 'mergeGroups' is specified, 'ctrlGroup' should not")

    ## Generate report
    ## -------------------------------------------------------------------------
    ## Without rendering
    args <- args0
    res <- do.call(runSpectronautAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))

    ## Stop if forceOverwrite = FALSE
    args <- args0
    args$forceOverwrite <- FALSE
    expect_error(res <- do.call(runSpectronautAnalysis, args))

    ## Message if forceOverwrite = TRUE
    args <- args0
    args$forceOverwrite <- TRUE
    expect_message(res <- do.call(runSpectronautAnalysis, args),
                   "already exists but forceOverwrite = TRUE")

    ## In new, non-existing directory and with custom yml
    args <- args0
    args$outputDir <- file.path(outDir, "new_Spectronaut_dir")
    args$customYml <- system.file("extdata", "custom.yml",
                                  package = "einprot")
    res <- do.call(runSpectronautAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(args$outputDir, paste0(outBaseName, ".Rmd"))))
    tmp <- readLines(file.path(args$outputDir, paste0(outBaseName, ".Rmd")))
    # expect_true(grepl("theme: journal", tmp[4]))

    ## With rendering
    skip_if(!capabilities()["X11"])
    args <- args0
    args$doRender <- TRUE
    res <- do.call(runSpectronautAnalysis, args)

    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".html"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".html"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_iSEE.R"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_sce.rds"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_heatmap_centered.pdf"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_feature_info.txt"))))
})
