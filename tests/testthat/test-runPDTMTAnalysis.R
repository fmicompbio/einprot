test_that("runPDTMTAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "PDTMTAnalysis"
    args0 <- list(
        templateRmd = system.file("extdata/process_basic_template.Rmd",
                                  package = "einprot"),
        outputDir = outDir,
        outputBaseName = outBaseName,
        reportTitle = "PD TMT data processing",
        reportAuthor = "",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "type1"),
        species = "baker's yeast",
        pdOutputFolder = system.file("extdata", "pdtmt_example",
                                     package = "einprot"),
        pdResultName = "Fig2_m23139_RTS_QC_varMods",
        inputLevel = "Proteins",
        pdAnalysisFile = system.file("extdata", "pdtmt_example",
                                     "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
                                     package = "einprot"),
        idCol = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession"),
                                        combineWhen = "nonunique",
                                        splitSeparator = ";", joinSeparator = "."),
        labelCol = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession"),
                                           combineWhen = "nonunique",
                                           splitSeparator = ";", joinSeparator = "."),
        geneIdCol = function(df) getFirstId(df, colName = "Gene.Symbol",
                                            separator = ";"),
        proteinIdCol = function(df) getFirstId(df, colName = "Accession",
                                               separator = ";"),
        stringIdCol = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession"),
                                              combineWhen = "missing", makeUnique = FALSE),
        extraFeatureCols = NULL,
        modificationsCol = "Modifications.in.Master.Proteins",
        excludeUnmodifiedPeptides = FALSE,
        keepModifications = NULL,
        iColPattern = "^Abundance\\.F.+\\.Sample\\.",
        sampleAnnot = data.frame(
            sample = c("HIS4KO_S05", "HIS4KO_S06", "HIS4KO_S07", "HIS4KO_S08",
                       "MET6KO_S01", "MET6KO_S02", "MET6KO_S03", "MET6KO_S04",
                       "URA2KO_S09", "URA2KO_S10", "URA2KO_S11", "URA2KO_S12",
                       "WT_S13", "WT_S14", "WT_S15", "WT_S16"),
            group = c(rep("HIS4KO", 4), rep("MET6KO", 4), rep("URA2KO", 4),
                      rep("WT", 4))),
        includeOnlySamples = "",
        excludeSamples = "",
        minScore = 10,
        minDeltaScore = 0.2,
        minPeptides = 2,
        minPSMs = 2,
        masterProteinsOnly = FALSE,
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
        normMethod = "center.median",
        spikeFeatures = NULL,
        stattest = "limma",
        minNbrValidValues = 2,
        minlFC = 0,
        samSignificance = FALSE,
        nperm = 250,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        volcanoMaxFeatures = 25,
        volcanoLabelSign = "both",
        volcanoS0 = 0.1,
        volcanoFeaturesToLabel = "",
        addInteractiveVolcanos = FALSE,
        interactiveDisplayColumns = NULL,
        interactiveGroupColumn = NULL,
        complexFDRThr = 0.1,
        maxNbrComplexesToPlot = Inf,
        maxComplexSimilarity = 1,
        seed = 42,
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
    ## runPDTMTAnalysis())
    ## -------------------------------------------------------------------------
    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'forceOverwrite' must have length 1")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'doRender' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'namesexperimentInfo' must not be NULL")

    ## species
    args <- args0
    args$species <- c("Mouse", "Human")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    ## pdOutputFolder
    args <- args0
    args$pdOutputFolder <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdOutputFolder' must be of class 'character'")
    args$pdOutputFolder <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdOutputFolder' must have length 1")

    ## pdResultName
    args <- args0
    args$pdResultName <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdResultName' must be of class 'character'")
    args$pdResultName <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdResultName' must have length 1")

    args <- args0
    args$pdResultName <- "missing"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "missing_Proteins.txt doesn't exist")

    ## inputLevel
    args <- args0
    args$inputLevel <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'inputLevel' must be of class 'character'")
    args$inputLevel <- c("Proteins", "PeptideGroups")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'inputLevel' must have length 1")

    args <- args0
    args$inputLevel <- "missing"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'inputLevel' must be one of")

    ## pdAnalysisFile
    args <- args0
    args$pdAnalysisFile <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdAnalysisFile' must be of class 'character'")
    args$pdAnalysisFile <- c("name1", "name2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdAnalysisFile' must have length 1")
    args$pdAnalysisFile <- "missing"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'pdAnalysisFile' must point to an existing file")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'idCol' must be of class 'character'")

    ## labelCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'labelCol' must be of class 'character'")

    ## geneIdCol
    args <- args0
    args$geneIdCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'geneIdCol' must be of class 'character'")

    ## proteinIdCol
    args <- args0
    args$proteinIdCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'proteinIdCol' must be of class 'character'")

    ## stringIdCol
    args <- args0
    args$stringIdCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stringIdCol' must be of class 'character'")

    ## extraFeatureCols
    args <- args0
    args$extraFeatureCols <- function(df) df$einprotId
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'extraFeatureCols' must be of class 'list'")
    args$extraFeatureCols <- list(function(df) df$einprotId)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'namesextraFeatureCols' must not be NULL")
    args$extraFeatureCols <- list(newCol = 1)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'i' must be of class 'character'")
    args$extraFeatureCols <- list(newCol = function(x, y) "a")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "length(formals(i)) == 1 is not TRUE", fixed = TRUE)

    ## modificationsCol
    args <- args0
    args$inputLevel <- "PeptideGroups"
    args$modificationsCol <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'modificationsCol' must be of class 'character'")
    args$modificationsCol <- c("Modifications", "Modifications.in.Master.Proteins")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'modificationsCol' must have length 1")

    ## excludeUnmodifiedPeptides
    args <- args0
    args$inputLevel <- "PeptideGroups"
    args$excludeUnmodifiedPeptides <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'excludeUnmodifiedPeptides' must be of class 'logical'")
    args$excludeUnmodifiedPeptides <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'excludeUnmodifiedPeptides' must have length 1")

    ## keepModifications
    args <- args0
    args$inputLevel <- "PeptideGroups"
    args$keepModifications <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'keepModifications' must be of class 'character'")

    ## iColPattern
    args <- args0
    args$iColPattern <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'iColPattern' must be of class 'character'")
    args$iColPattern <- c("^Abundance\\.F.+\\.Sample\\.",
                          "^Abundance\\.F.+\\.Sample\\.")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'iColPattern' must have length 1")
    args$iColPattern <- c("^LFQ\\\\.intensity\\\\.")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'iColPattern' must be one of")
    args$iColPattern <- c("^Abundance\\\\.F.+\\\\.Sample\\\\.")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'iColPattern' must be one of")

    ## sampleAnnot
    args <- args0
    colnames(args$sampleAnnot) <- c("sample", "condition")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    colnames(args$sampleAnnot) <- NULL
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'colnamessampleAnnot' must not be NULL")
    args$sampleAnnot <- as.matrix(args0$sampleAnnot)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'sampleAnnot' must be of class 'data.frame'")
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$group <- as.numeric(as.factor(args$sampleAnnot$group))
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'$sampleAnnotgroup' must be of class 'character'",
                 fixed = TRUE)
    args$sampleAnnot <- args0$sampleAnnot
    args$sampleAnnot$sample[1] <- paste0("iBAQ.",
                                         args$sampleAnnot$sample[1])
    expect_error(do.call(runPDTMTAnalysis, args),
                 "Not all sample names are available in the sample")

    ## includeOnlySamples
    args <- args0
    args$includeOnlySamples <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'includeOnlySamples' must be of class 'character'")

    ## excludeSamples
    args <- args0
    args$excludeSamples <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'excludeSamples' must be of class 'character'")
    args$includeOnlySamples <- "Adnp"
    args$excludeSamples <- "RBC_ctrl"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "Please specify max one of includeOnlySamples")

    ## minScore/minDeltaScore
    args <- args0
    args$minScore <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minScore' must be of class 'numeric'")
    args$minScore <- c(1, 2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minScore' must have length 1")

    args <- args0
    args$inputLevel <- "PeptideGroups"
    args$minDeltaScore <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minDeltaScore' must be of class 'numeric'")
    args$minDeltaScore <- c(1, 2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minDeltaScore' must have length 1")

    ## minPeptides/minPSMs
    args <- args0
    args$minPeptides <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minPeptides' must be of class 'numeric'")
    args$minPeptides <- c(1, 2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minPeptides' must have length 1")

    args <- args0
    args$inputLevel <- "PeptideGroups"
    args$minPSMs <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minPSMs' must be of class 'numeric'")
    args$minPSMs <- c(1, 2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minPSMs' must have length 1")

    ## masterProteinsOnly
    args <- args0
    args$masterProteinsOnly <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'masterProteinsOnly' must be of class 'logical'")
    args$masterProteinsOnly <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'masterProteinsOnly' must have length 1")

    ## imputeMethod
    args <- args0
    args$imputeMethod <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'imputeMethod' must be of class 'character'")
    args$imputeMethod <- c("MinProb", "impSeqRob")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'imputeMethod' must have length 1")
    args$imputeMethod <- "wrong"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'imputeMethod' must be one of")

    ## assaysForExport
    args <- args0
    args$assaysForExport <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'assaysForExport' must be of class 'character'")

    ## addHeatmaps
    args <- args0
    args$addHeatmaps <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'addHeatmaps' must be of class 'logical'")
    args$addHeatmaps <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'addHeatmaps' must have length 1")

    ## mergeGroups
    args <- args0
    args$mergeGroups <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'mergeGroups' must be of class 'list'")
    args$mergeGroups <- list(c("c1", "c2"))
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'mergeGroups' must be a named list")
    args$mergeGroups <- list(g1 = c("c1", "c2"),
                             g1 = c("c3", "c4"))
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'mergeGroups' must be a named list")

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(runPDTMTAnalysis, args),
                 "Each entry in 'comparisons' must have exactly two elements")

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'ctrlGroup' must be of class 'character'")
    # args$ctrlGroup <- c("name1", "name2")
    # expect_error(do.call(runPDTMTAnalysis, args),
    #              "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'baselineGroup' must be of class 'character'")

    ## normMethod
    args <- args0
    args$normMethod <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'normMethod' must be of class 'character'")
    args$normMethod <- c("none", "center.median")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'normMethod' must have length 1")
    args$normMethod <- "wrong"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'normMethod' must be one of")

    ## spikeFeatures
    args <- args0
    args$spikeFeatures <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'spikeFeatures' must be of class 'character'")

    ## stattest
    args <- args0
    args$stattest <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stattest' must be of class 'character'")
    args$stattest <- c("limma", "ttest")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stattest' must have length 1")
    args$stattest <- "wrong"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'stattest' must be one of")
    # args$stattest <- "ttest"
    # expect_error(do.call(runPDTMTAnalysis, args),
    #              "'ttest' is currently not supported")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## samSignificance
    args <- args0
    args$samSignificance <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'samSignificance' must be of class 'logical'")
    args$samSignificance <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'samSignificance' must have length 1")

    ## nperm
    args <- args0
    args$nperm <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'nperm' must be of class 'numeric'")
    args$nperm <- c(10, 20)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'nperm' must have length 1")
    args$nperm <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'nperm' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoS0
    args <- args0
    args$volcanoS0 <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoS0' must be of class 'numeric'")
    args$volcanoS0 <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoS0' must have length 1")
    args$volcanoS0 <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoS0' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'addInteractiveVolcanos' must have length 1")

    # interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'interactiveGroupColumn' must have length 1")

    ## complexFDRThr
    args <- args0
    args$complexFDRThr <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexFDRThr' must be of class 'numeric'")
    args$complexFDRThr <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexFDRThr' must have length 1")
    args$complexFDRThr <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexFDRThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## maxNbrComplexesToPlot
    args <- args0
    args$maxNbrComplexesToPlot <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'maxNbrComplexesToPlot' must be of class 'numeric'")
    args$maxNbrComplexesToPlot <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'maxNbrComplexesToPlot' must have length 1")
    args$maxNbrComplexesToPlot <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'maxNbrComplexesToPlot' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## maxComplexSimilarity
    args <- args0
    args$maxComplexSimilarity <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'maxComplexSimilarity' must be of class 'numeric'")
    args$maxComplexSimilarity <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'maxComplexSimilarity' must have length 1")

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## includeFeatureCollections
    args <- args0
    args$includeFeatureCollections <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'includeFeatureCollections' must be of class 'character'")
    args$includeFeatureCollections <- "wrong"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'includeFeatureCollections' must be one of")

    ## minSizeToKeepSet
    args <- args0
    args$minSizeToKeepSet <- "1"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minSizeToKeepSet' must be of class 'numeric'")
    args$minSizeToKeepSet <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minSizeToKeepSet' must have length 1")
    args$minSizeToKeepSet <- -1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'minSizeToKeepSet' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## customComplexes
    args <- args0
    args$customComplexes <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'customComplexes' must be of class 'list'")
    args$customComplexes <- list(c("c1", "c2"),
                                 c("c2", "c4", "c3"))
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'namescustomComplexes' must not be NULL")

    ## complexSpecies
    args <- args0
    args$complexSpecies <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexSpecies' must be of class 'character'")
    args$complexSpecies <- c("all", "current")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexSpecies' must have length 1")
    args$complexSpecies <- "wrong"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "All values in 'complexSpecies' must be one of")

    ## complexDbPath
    args <- args0
    args$complexDbPath <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexDbPath' must be of class 'character'")
    args$complexDbPath <- "missing_file"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'complexDbPath' must point to an existing file")

    ## stringVersion
    args <- args0
    args$stringVersion <- 11
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stringVersion' must be of class 'character'")
    args$stringVersion <- c("11.0", "11.5")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stringVersion' must have length 1")

    ## stringDir
    args <- args0
    args$stringDir <- 11
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stringDir' must be of class 'character'")
    args$stringDir <- c("", "")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'stringDir' must have length 1")

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(runPDTMTAnalysis, args),
                 "'customYml' must point to an existing file")

    ## Test fix for backwards compatibility with mergeGroups/ctrlGroup
    args <- args0
    args$mergeGroups <- list(G1 = c("URA2KO", "WT"))
    args$ctrlGroup <- c("URA2KO", "WT")
    expect_error(do.call(runPDTMTAnalysis, args),
                 "If 'mergeGroups' is specified, 'ctrlGroup' should not")

    ## Note that the tests below depend on the value of 'fcn' in
    ## generateConfigChunk ("deparse" or "dump")
    args <- args0
    args$doRender <- FALSE
    args$ctrlGroup <- c("URA2KO", "WT")
    res <- do.call(runPDTMTAnalysis, args)
    expect_type(res, "character")
    tmp <- readLines(res)
    ## if fcn = "dump"
    idx <- grep("ctrlGroup <-", tmp, fixed = TRUE)
    expect_length(length(idx), 1L)
    # expect_equal(tmp[idx + 1], '\"URA2KO.WT\"')
    idx <- grep("mergeGroups <-", tmp, fixed = TRUE)
    expect_length(length(idx), 1L)
    # expect_equal(tmp[idx + 1], 'list(URA2KO.WT = c(\"URA2KO\", \"WT\"))')
    ## if fcn = "deparse"
    # expect_true(length(
    #     grep('ctrlGroup <- \"URA2KO.WT\"', tmp, fixed = TRUE)) > 0)
    # expect_true(length(
    #     grep('mergeGroups <- list(URA2KO.WT = c(\"URA2KO\", \"WT\"))',
    #          tmp, fixed = TRUE)) > 0)


    ## Generate report
    ## -------------------------------------------------------------------------
    ## Without rendering
    args <- args0
    res <- do.call(runPDTMTAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))

    ## Stop if forceOverwrite = FALSE
    args <- args0
    args$forceOverwrite <- FALSE
    expect_error(res <- do.call(runPDTMTAnalysis, args))

    ## Message if forceOverwrite = TRUE
    args <- args0
    args$forceOverwrite <- TRUE
    expect_message(res <- do.call(runPDTMTAnalysis, args),
                   "already exists but forceOverwrite = TRUE")

    ## iColPattern without escaped period
    args <- args0
    args$forceOverwrite <- TRUE
    args$iColPattern <- "^Abundance.F.+.Sample."
    expect_message(res <- do.call(runPDTMTAnalysis, args),
                   "already exists but forceOverwrite = TRUE")

    ## In new, non-existing directory and with custom yml
    args <- args0
    args$outputDir <- file.path(outDir, "new_pd_dir")
    args$customYml <- system.file("extdata", "custom.yml",
                                  package = "einprot")
    res <- do.call(runPDTMTAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(args$outputDir, paste0(outBaseName, ".Rmd"))))
    tmp <- readLines(file.path(args$outputDir, paste0(outBaseName, ".Rmd")))
    # expect_true(grepl("theme: journal", tmp[4]))

    ## With rendering
    skip_if(!capabilities()["X11"])
    args <- args0
    args$doRender <- TRUE
    res <- do.call(runPDTMTAnalysis, args)

    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".html"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".html"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_iSEE.R"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_sce.rds"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_feature_info.txt"))))

})
