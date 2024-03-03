test_that("text snippet generation works", {
    ## testText
    expect_error(testText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(testText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(testText(testType = "missing"),
                 "All values in 'testType' must be one of")
    expect_error(testText(testType = "limma", minlFC = "1"),
                 "'minlFC' must be of class 'numeric'")
    expect_error(testText(testType = "limma", minlFC = c(1, 2)),
                 "'minlFC' must have length 1")
    expect_error(testText(testType = "limma", samSignificance = "1"),
                 "'samSignificance' must be of class 'logical'")
    expect_error(testText(testType = "limma", samSignificance = c(TRUE, FALSE)),
                 "'samSignificance' must have length 1")

    expect_type(testText(testType = "limma", minlFC = 1), "character")
    expect_equal(length(testText(testType = "limma")), 1)
    expect_true(grepl("the treat function", testText(testType = "limma", minlFC = 1)))

    expect_type(testText(testType = "limma", minlFC = 0), "character")
    expect_equal(length(testText(testType = "limma")), 1)
    expect_false(grepl("the treat function", testText(testType = "limma", minlFC = 0)))
    expect_true(grepl("limma", testText(testType = "limma", minlFC = 0)))

    expect_type(testText(testType = "ttest", samSignificance = TRUE), "character")
    expect_equal(length(testText(testType = "ttest", samSignificance = TRUE)), 1)
    expect_true(grepl("a Student's t-test", testText(testType = "ttest",
                                                     samSignificance = TRUE)))
    expect_true(grepl("Tusher", testText(testType = "ttest",
                                         samSignificance = TRUE)))
    expect_type(testText(testType = "ttest", samSignificance = FALSE), "character")
    expect_equal(length(testText(testType = "ttest", samSignificance = FALSE)), 1)
    expect_true(grepl("a Student's t-test", testText(testType = "ttest",
                                                     samSignificance = FALSE)))
    expect_false(grepl("Tusher", testText(testType = "ttest",
                                          samSignificance = FALSE)))

    expect_type(testText(testType = "proDA"), "character")
    expect_equal(length(testText(testType = "proDA")), 1)
    expect_true(grepl("proDA", testText(testType = "proDA")))

    expect_type(testText(testType = "none"), "character")
    expect_equal(length(testText(testType = "none")), 1)
    expect_true(grepl("no statistical testing", testText(testType = "none")))

    ## normText
    expect_error(normText(normMethod = 1),
                 "'normMethod' must be of class 'character'")
    expect_error(normText(normMethod = c("limma", "ttest")),
                 "'normMethod' must have length 1")

    expect_type(normText(normMethod = "none"), "character")
    expect_equal(length(normText(normMethod = "none")), 1)
    expect_true(grepl("are not normalized", normText(normMethod = "none")))
    expect_type(normText(normMethod = "center.median"), "character")
    expect_equal(length(normText(normMethod = "center.median")), 1)
    expect_true(grepl("using the center.median method",
                      normText(normMethod = "center.median")))

    ## saText
    expect_error(saText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(saText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(saText(testType = "missing"),
                 "All values in 'testType' must be one of")

    expect_type(saText(testType = "limma"), "character")
    expect_equal(length(saText(testType = "limma")), 1)
    expect_true(grepl("square root of the residual standard",
                      saText(testType = "limma")))

    expect_type(saText(testType = "ttest"), "character")
    expect_equal(length(saText(testType = "ttest")), 1)
    expect_equal(saText(testType = "ttest"), "")

    expect_type(saText(testType = "proDA"), "character")
    expect_equal(length(saText(testType = "proDA")), 1)
    expect_equal(saText(testType = "proDA"), "")

    expect_type(saText(testType = "none"), "character")
    expect_equal(length(saText(testType = "none")), 1)
    expect_equal(saText(testType = "none"), "")

    ## expDesignText
    expect_error(expDesignText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(expDesignText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(expDesignText(testType = "missing"),
                 "All values in 'testType' must be one of")

    expect_type(expDesignText(testType = "limma"), "character")
    expect_equal(length(expDesignText(testType = "limma")), 1)
    expect_true(grepl("experimental design", expDesignText(testType = "limma")))

    expect_type(expDesignText(testType = "ttest"), "character")
    expect_equal(length(expDesignText(testType = "ttest")), 1)
    expect_equal(expDesignText(testType = "ttest"), "")

    expect_type(expDesignText(testType = "proDA"), "character")
    expect_equal(length(expDesignText(testType = "proDA")), 1)
    expect_equal(expDesignText(testType = "proDA"), "")

    expect_type(expDesignText(testType = "none"), "character")
    expect_equal(length(expDesignText(testType = "none")), 1)
    expect_equal(expDesignText(testType = "none"), "")

    ## introText
    expect_error(introText(expType = 1),
                 "'expType' must be of class 'character'")
    expect_error(introText(expType = c("MaxQuant", "FragPipe")),
                 "'expType' must have length 1")
    expect_error(introText(expType = "missing"),
                 "All values in 'expType' must be one of")

    expect_type(introText(expType = "MaxQuant"), "character")
    expect_equal(length(introText(expType = "MaxQuant")), 1)
    expect_true(grepl("[MaxQuant](https://www.maxquant.org/)",
                      introText(expType = "MaxQuant"), fixed = TRUE))

    expect_type(introText(expType = "FragPipe"), "character")
    expect_equal(length(introText(expType = "FragPipe")), 1)
    expect_true(grepl("[FragPipe](https://fragpipe.nesvilab.org/)",
                      introText(expType = "FragPipe"), fixed = TRUE))

    expect_type(introText(expType = "ProteomeDiscoverer"), "character")
    expect_equal(length(introText(expType = "ProteomeDiscoverer")), 1)
    expect_true(grepl("[Proteome Discoverer](https://www.thermofisher.com",
                      introText(expType = "ProteomeDiscoverer"), fixed = TRUE))

    expect_type(introText(expType = "DIANN"), "character")
    expect_equal(length(introText(expType = "DIANN")), 1)
    expect_true(grepl("[DIA-NN](https://github.com/vdemichev/DiaNN)",
                      introText(expType = "DIANN"), fixed = TRUE))

    expect_type(introText(expType = "Spectronaut"), "character")
    expect_equal(length(introText(expType = "Spectronaut")), 1)
    expect_true(grepl("[Spectronaut](https://biognosys.com/",
                      introText(expType = "Spectronaut"), fixed = TRUE))

    ## inputText
    expect_error(inputText(expTypeLevel = 1),
                 "'expTypeLevel' must be of class 'character'")
    expect_error(inputText(expTypeLevel = c("MaxQuant", "FragPipe")),
                 "'expTypeLevel' must have length 1")
    expect_error(inputText(expTypeLevel = "missing"),
                 "All values in 'expTypeLevel' must be one of")

    expect_type(inputText(expTypeLevel = "MaxQuant"), "character")
    expect_equal(length(inputText(expTypeLevel = "MaxQuant")), 1)
    expect_true(grepl("The input to this workflow is a `proteinGroups.txt` file",
                      inputText(expTypeLevel = "MaxQuant"), fixed = TRUE))

    expect_type(inputText(expTypeLevel = "FragPipe"), "character")
    expect_equal(length(inputText(expTypeLevel = "FragPipe")), 1)
    expect_true(grepl("The input to this workflow is a `combined_protein` file",
                      inputText(expTypeLevel = "FragPipe"), fixed = TRUE))

    expect_type(inputText(expTypeLevel = "ProteomeDiscovererProteins"), "character")
    expect_equal(length(inputText(expTypeLevel = "ProteomeDiscovererProteins")), 1)
    expect_true(grepl("The input to this workflow is a `Proteins.txt` file from",
                      inputText(expTypeLevel = "ProteomeDiscovererProteins"), fixed = TRUE))

    expect_type(inputText(expTypeLevel = "ProteomeDiscovererPeptideGroups"), "character")
    expect_equal(length(inputText(expTypeLevel = "ProteomeDiscovererPeptideGroups")), 1)
    expect_true(grepl("The input to this workflow is a `PeptideGroups.txt` file from",
                      inputText(expTypeLevel = "ProteomeDiscovererPeptideGroups"), fixed = TRUE))

    expect_type(inputText(expTypeLevel = "DIANN"), "character")
    expect_equal(length(inputText(expTypeLevel = "DIANN")), 1)
    expect_true(grepl("The input to this workflow is a `pg_matrix.tsv`,",
                      inputText(expTypeLevel = "DIANN"), fixed = TRUE))

    expect_type(inputText(expTypeLevel = "Spectronaut"), "character")
    expect_equal(length(inputText(expTypeLevel = "Spectronaut")), 1)
    expect_true(grepl("The input to this workflow is a `Report.tsv`",
                      inputText(expTypeLevel = "Spectronaut"), fixed = TRUE))

    ## emptySampleText
    sce0 <- sce_mq_final
    expect_error(emptySampleText(sce = 1, assayName = "log2_iBAQ_withNA"),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(emptySampleText(sce = sce0, assayName = 1),
                 "'assayName' must be of class 'character'")
    expect_error(emptySampleText(sce = sce0,
                                 assayName = c("log2_iBAQ_withNA",
                                               "log2_iBAQ_withNA")),
                 "'assayName' must have length 1")
    expect_error(emptySampleText(sce = sce0, assayName = "missing"),
                 "assayName %in% SummarizedExperiment::assayNames(sce)",
                 fixed = TRUE)

    expect_type(emptySampleText(sce = sce0, assayName = "log2_iBAQ_withNA"),
                "character")
    expect_equal(emptySampleText(sce = sce0, assayName = "log2_iBAQ_withNA"),
                 "")
    SummarizedExperiment::assay(sce0, "log2_iBAQ_withNA")[, c(2, 5)] <- NA
    expect_type(emptySampleText(sce = sce0, assayName = "log2_iBAQ_withNA"),
                "character")
    expect_equal(emptySampleText(sce = sce0, assayName = "log2_iBAQ_withNA"),
                 paste0("The following sample(s) do not have any detected ",
                        "features and will be removed from further analysis: ",
                        "Adnp_IP05, Chd4BF_IP08."))

})
