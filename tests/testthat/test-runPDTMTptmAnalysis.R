test_that("runPDTMTptmAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "PDTMTAnalysis"
    args0 <- list(
        templateRmd = system.file("extdata", "process_PD_TMT_PTM_template.Rmd",
                                  package = "einprot"),
        outputDir = outDir,
        outputBaseName = outBaseName,
        reportTitle = "reportTitle",
        reportAuthor = "reportAuthor",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "experiment",
                              "Sample name" = "example"),
        species = "roundworm",
        sceProteins = system.file("extdata", "pdtmt_example",
                                  "Fig2_m23139_RTS_QC_varMods_Proteins_sce.rds",
                                  package = "einprot"),
        scePeptides = system.file("extdata", "pdtmt_example",
                                  "Fig2_m23139_RTS_QC_varMods_PeptideGroups_sce.rds",
                                  package = "einprot"),
        assayForTests = "log2_Abundance_norm",
        assayImputation = "imputed_Abundance",
        idCol = function(df) combineIds(df, combineCols = c("Annotated.Sequence", "Positions.in.Proteins"),
                                        combineWhen = "nonunique", makeUnique = TRUE,
                                        splitSeparator = ";", joinSeparator = "."),
        labelCol = function(df) combineIds(df, combineCols = c("Annotated.Sequence", "Positions.in.Proteins"),
                                           combineWhen = "nonunique", makeUnique = FALSE,
                                           splitSeparator = ";", joinSeparator = "."),
        proteinIdColProteins = "einprotProtein",
        proteinIdColPeptides = "einprotProtein",
        comparisons = list(),
        ctrlGroup = "",
        allPairwiseComparisons = TRUE,
        singleFit = FALSE,
        subtractBaseline = FALSE,
        baselineGroup = "",
        testType = "interaction",
        minNbrValidValues = 2,
        minlFC = 0,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        volcanoMaxFeatures = 10,
        volcanoLabelSign = "both",
        volcanoFeaturesToLabel = c(""),
        addInteractiveVolcanos = FALSE,
        interactiveDisplayColumns = NULL,
        interactiveGroupColumn = NULL,
        seed = 123,
        linkTableColumns = c(),
        customYml = NULL,
        doRender = FALSE
    )

    ## Fail with wrong parameter values
    ## (essentially the same tests as for checkArgumentsPDTMTptm)
    ## -------------------------------------------------------------------------
    ## templateRmd
    args <- args0
    args$templateRmd <- c(args$templateRmd, args$templateRmd)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'templateRmd' must have length 1")
    args$templateRmd <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'templateRmd' must be of class 'character'")
    args$templateRmd <- "missing_file"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'templateRmd' must point to an existing file")

    ## outputDir
    args <- args0
    args$outputDir <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'outputDir' must be of class 'character'")
    args$outputDir <- c("dir1", "dir2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'outputDir' must have length 1")

    ## outputBaseName
    args <- args0
    args$outputBaseName <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'outputBaseName' must be of class 'character'")
    args$outputBaseName <- c("name1", "name2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'outputBaseName' must have length 1")

    ## reportTitle
    args <- args0
    args$reportTitle <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'reportTitle' must be of class 'character'")
    args$reportTitle <- c("name1", "name2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'reportTitle' must have length 1")

    ## reportAuthor
    args <- args0
    args$reportAuthor <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'reportAuthor' must be of class 'character'")
    args$reportAuthor <- c("name1", "name2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'reportAuthor' must have length 1")

    ## forceOverwrite
    args <- args0
    args$forceOverwrite <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'forceOverwrite' must be of class 'logical'")
    args$forceOverwrite <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'forceOverwrite' must have length 1")

    ## experimentInfo
    args <- args0
    args$experimentInfo <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'experimentInfo' must be of class 'list'")
    args$experimentInfo <- list(1, 2)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'namesexperimentInfo' must not be NULL")

    ## species
    args <- args0
    args$species <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "Unknown species 1")
    args$species <- c("Mouse", "Human")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "length(species) == 1 is not TRUE", fixed = TRUE)

    ## sceProteins
    args <- args0
    args$sceProteins <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'sceProteins' must be of class 'character'")
    args$sceProteins <- "missing"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'sceProteins' must point to an existing file")

    ## scePeptides
    args <- args0
    args$scePeptides <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'scePeptides' must be of class 'character'")
    args$scePeptides <- "missing"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'scePeptides' must point to an existing file")

    ## assayForTests
    args <- args0
    args$assayForTests <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'assayForTests' must be of class 'character'")
    args$assayForTests <- c("log2_Abundance_norm", "log2_Abundance_norm")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'assayForTests' must have length 1")

    ## assayImputation
    args <- args0
    args$assayImputation <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'assayImputation' must be of class 'character'")
    args$assayImputation <- c("imputed_Abundance", "imputed_Abundance")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'assayImputation' must have length 1")

    ## idCol
    args <- args0
    args$idCol <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'idCol' must be of class 'character'")

    ## labelCol
    args <- args0
    args$labelCol <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'labelCol' must be of class 'character'")

    ## proteinIdColProteins
    args <- args0
    args$proteinIdColProteins <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'proteinIdColProteins' must be of class 'character'")
    args$proteinIdColProteins <- function(x, y) x + y
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "length(formals(proteinIdColProteins)) == 1 is not TRUE", fixed = TRUE)

    ## proteinIdColPeptides
    args <- args0
    args$proteinIdColPeptides <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'proteinIdColPeptides' must be of class 'character'")
    args$proteinIdColPeptides <- function(x, y) x + y
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "length(formals(proteinIdColPeptides)) == 1 is not TRUE", fixed = TRUE)

    ## comparisons
    args <- args0
    args$comparisons <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(g1 = c("c1", "c2"),
                             g2 = c("c2", "c4", "c3"))
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "Each entry in 'comparisons' must have exactly two elements")

    ## ctrlGroup
    args <- args0
    args$ctrlGroup <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'ctrlGroup' must be of class 'character'")
    args$ctrlGroup <- c("name1", "name2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'ctrlGroup' must have length 1")

    ## allPairwiseComparisons
    args <- args0
    args$allPairwiseComparisons <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'allPairwiseComparisons' must be of class 'logical'")
    args$allPairwiseComparisons <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'allPairwiseComparisons' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$baselineGroup <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'baselineGroup' must be of class 'character'")

    ## testType
    args <- args0
    args$testType <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'testType' must be of class 'character'")
    args$testType <- c("interaction", "welch")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'testType' must have length 1")
    args$testType <- "wrong"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "All values in 'testType' must be one of")

    ## minNbrValidValues
    args <- args0
    args$minNbrValidValues <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(10, 20)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(10, 20)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoMaxFeatures
    args <- args0
    args$volcanoMaxFeatures <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoMaxFeatures' must be of class 'numeric'")
    args$volcanoMaxFeatures <- c(0.1, 0.2)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoMaxFeatures' must have length 1")
    args$volcanoMaxFeatures <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoMaxFeatures' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoLabelSign
    args <- args0
    args$volcanoLabelSign <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoLabelSign' must be of class 'character'")
    args$volcanoLabelSign <- c("both", "pos")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoLabelSign' must have length 1")
    args$volcanoLabelSign <- "missing"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "All values in 'volcanoLabelSign' must be one of")

    ## volcanoFeaturesToLabel
    args <- args0
    args$volcanoFeaturesToLabel <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'volcanoFeaturesToLabel' must be of class 'character'")

    ## addInteractiveVolcanos
    args <- args0
    args$addInteractiveVolcanos <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'addInteractiveVolcanos' must be of class 'logical'")
    args$addInteractiveVolcanos <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'addInteractiveVolcanos' must have length 1")

    # interactiveDisplayColumns
    args <- args0
    args$interactiveDisplayColumns <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'interactiveDisplayColumns' must be of class 'character'")

    # interactiveGroupColumn
    args <- args0
    args$interactiveGroupColumn <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'interactiveGroupColumn' must be of class 'character'")
    args$interactiveGroupColumn <- c("col1", "col2")
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'interactiveGroupColumn' must have length 1")

    ## seed
    args <- args0
    args$seed <- "1"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'seed' must be of class 'numeric'")
    args$seed <- c(4, 5)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'seed' must have length 1")
    args$seed <- -1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'seed' must be within [1,Inf] (inclusive)",
                 fixed = TRUE)

    ## linkTableColumns
    args <- args0
    args$linkTableColumns <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'linkTableColumns' must be of class 'character'")

    ## customYml
    args <- args0
    args$customYml <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'customYml' must be of class 'character'")
    args$customYml <- "missing_file"
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'customYml' must point to an existing file")

    ## doRender
    args <- args0
    args$doRender <- 1
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'doRender' must be of class 'logical'")
    args$doRender <- c(TRUE, FALSE)
    expect_error(do.call(runPDTMTptmAnalysis, args),
                 "'doRender' must have length 1")

    ## Generate report
    ## -------------------------------------------------------------------------
    ## Without rendering
    args <- args0
    res <- do.call(runPDTMTptmAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))

    ## Non-existing output directory
    args <- args0
    args$outputDir <- file.path(outDir, "new_directory_pdtmtptm")
    res <- do.call(runPDTMTptmAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(outDir, "new_directory_pdtmtptm",
                                      paste0(outBaseName, ".Rmd"))))

    ## Stop if forceOverwrite = FALSE
    args <- args0
    args$forceOverwrite <- FALSE
    expect_error(res <- do.call(runPDTMTptmAnalysis, args))

    ## Message if forceOverwrite = TRUE
    args <- args0
    args$forceOverwrite <- TRUE
    expect_message(res <- do.call(runPDTMTptmAnalysis, args),
                   "already exists but forceOverwrite = TRUE")

    ## In new, non-existing directory and with custom yml
    args <- args0
    args$outputDir <- file.path(outDir, "new_pd_dir")
    args$customYml <- system.file("extdata", "custom.yml",
                                  package = "einprot")
    res <- do.call(runPDTMTptmAnalysis, args)
    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".Rmd"))
    expect_true(file.exists(file.path(args$outputDir, paste0(outBaseName, ".Rmd"))))
    tmp <- readLines(file.path(args$outputDir, paste0(outBaseName, ".Rmd")))
    # expect_true(grepl("theme: journal", tmp[4]))

    ## With rendering
    skip_if(!capabilities()["X11"])
    args <- args0
    args$doRender <- TRUE
    res <- do.call(runPDTMTptmAnalysis, args)

    expect_type(res, "character")
    expect_equal(basename(res), paste0(outBaseName, ".html"))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".Rmd"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, ".html"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_iSEE.R"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_sce.rds"))))
    expect_true(file.exists(file.path(outDir, paste0(outBaseName, "_feature_info.txt"))))

})
