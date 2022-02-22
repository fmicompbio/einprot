test_that("runMaxQuantAnalysis works", {
    outDir <- tempdir()
    outBaseName <- "MaxQuantAnalysis"
    args <- list(
        templateRmd = system.file("extdata/process_MaxQuant_template.Rmd",
                                  package = "einprot"),
        outputDir = outDir, outputBaseName = outBaseName,
        reportTitle = "MaxQuant LFQ data processing", reportAuthor = "",
        forceOverwrite = TRUE,
        experimentInfo = list("Experiment type" = "type1"),
        species = "mouse",
        mqFile = system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                             package = "einprot"),
        mqParameterFile = system.file("extdata", "mq_example", "1356_mqpar.xml",
                                      package = "einprot"),
        aName = "iBAQ", iColPattern = "^iBAQ\\\\.",
        sampleAnnot = data.frame(sample = c("Adnp_IP04", "Adnp_IP05",
                                            "Adnp_IP06", "Chd4BF_IP07",
                                            "Chd4BF_IP08", "Chd4BF_IP09",
                                            "RBC_ctrl_IP01", "RBC_ctrl_IP02",
                                            "RBC_ctrl_IP03"),
                                 group = c("Adnp", "Adnp", "Adnp", "Chd4BF", "Chd4BF",
                                           "Chd4BF", "RBC_ctrl", "RBC_ctrl", "RBC_ctrl")),
        includeOnlySamples = "", excludeSamples = "",
        minScore = 10, minPeptides = 2, imputeMethod = "MinProb",
        mergeGroups = list(), comparisons = list(),
        ctrlGroup = "", allPairwiseComparisons = TRUE,
        normMethod = "none", stattest = "limma", minNbrValidValues = 2,
        minlFC = 0, nperm = 250, volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1, volcanoMaxFeatures = 25,
        volcanoS0 = 0.1, volcanoFeaturesToLabel = "",
        addInteractiveVolcanos = FALSE, complexFDRThr = 0.1, seed = 42,
        includeFeatureCollections = "complexes", customComplexes = list(),
        complexSpecies = "all",
        complexDbPath = system.file("extdata", "complexes",
                                    "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                    package = "einprot")
    )
    skip_if(!capabilities()["X11"])
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
