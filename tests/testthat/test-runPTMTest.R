test_that("runPTMTest works", {
    sce_proteins <- readRDS(system.file("extdata", "pdtmt_example",
                                        "Fig2_m23139_RTS_QC_varMods_Proteins_sce.rds",
                                        package = "einprot"))
    rownames(sce_proteins) <- make.unique(rownames(sce_proteins))
    sce_peptides <- readRDS(system.file("extdata", "pdtmt_example",
                                        "Fig2_m23139_RTS_QC_varMods_PeptideGroups_sce.rds",
                                        package = "einprot"))

    ## =========================================================================
    ## Fails with wrong arguments
    ## =========================================================================
    args0 <- list(
        sceProteins = sce_proteins,
        scePeptides = sce_peptides,
        matchColProteins = "einprotMatchProt",
        matchColPeptides = "einprotMatchProt",
        testType = "interaction",
        comparisons = list(c("HIS4KO", "WT")),
        groupComposition = NULL,
        assayForTests = "log2_Abundance_norm",
        assayImputation = "imputed_Abundance",
        minNbrValidValues = 2,
        minlFC = 0,
        volcanoAdjPvalThr = 0.05,
        volcanoLog2FCThr = 1,
        baseFileName = NULL,
        singleFit = FALSE,
        subtractBaseline = FALSE,
        baselineGroup = ""
    )

    ## sceProteins
    args <- args0
    args$sceProteins <- 1
    expect_error(do.call(runPTMTest, args),
                 "'sceProteins' must be of class 'SummarizedExperiment'")
    args <- args0
    args$sceProteins$group <- as.numeric(as.factor(args$sceProteins$group))
    expect_error(do.call(runPTMTest, args),
                 "'$sceProteinsgroup' must be of class 'character'", fixed = TRUE)
    args <- args0
    rownames(args$sceProteins)[2] <- rownames(args$sceProteins)[1]
    expect_error(do.call(runPTMTest, args),
                 "The row names of sceProteins cannot contain duplicated entries.")

    ## scePeptides
    args <- args0
    args$scePeptides <- 1
    expect_error(do.call(runPTMTest, args),
                 "'scePeptides' must be of class 'SummarizedExperiment'")
    args <- args0
    args$scePeptides$group <- as.numeric(as.factor(args$scePeptides$group))
    expect_error(do.call(runPTMTest, args),
                 "'$scePeptidesgroup' must be of class 'character'", fixed = TRUE)
    args <- args0
    rownames(args$scePeptides)[2] <- rownames(args$scePeptides)[1]
    expect_error(do.call(runPTMTest, args),
                 "The row names of scePeptides cannot contain duplicated entries.")

    ## matchColProteins
    args <- args0
    args$matchColProteins <- 1
    expect_error(do.call(runPTMTest, args),
                 "'matchColProteins' must be of class 'character'")
    args <- args0
    args$matchColProteins <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "All values in 'matchColProteins' must be one of")
    args <- args0
    args$matchColProteins <- c("einprotMatchProt", "einprotMatchProt")
    expect_error(do.call(runPTMTest, args),
                 "'matchColProteins' must have length 1")

    ## matchColPeptides
    args <- args0
    args$matchColPeptides <- 1
    expect_error(do.call(runPTMTest, args),
                 "'matchColPeptides' must be of class 'character'")
    args <- args0
    args$matchColPeptides <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "All values in 'matchColPeptides' must be one of")
    args <- args0
    args$matchColPeptides <- c("einprotMatchProt", "einprotMatchProt")
    expect_error(do.call(runPTMTest, args),
                 "'matchColPeptides' must have length 1")

    ## testType
    args <- args0
    args$testType <- 1
    expect_error(do.call(runPTMTest, args),
                 "'testType' must be of class 'character'")
    args$testType <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "All values in 'testType' must be one of")
    args$testType <- c("interaction", "welch")
    expect_error(do.call(runPTMTest, args),
                 "'testType' must have length 1")

    ## comparisons
    args <- args0
    args$comparisons <- c(1, 2)
    expect_error(do.call(runPTMTest, args),
                 "'comparisons' must be of class 'list'")
    args$comparisons <- list(c("HIS4KO"))
    expect_error(do.call(runPTMTest, args),
                 "'comparison' must have length 2")
    args$comparisons <- list(c("missing", "HIS4KO"))
    expect_error(do.call(runPTMTest, args),
                 "All values in 'comparison' must be one of")

    ## groupComposition
    args <- args0
    args$groupComposition <- 1
    expect_error(do.call(runPTMTest, args),
                 "'groupComposition' must be of class 'list'")

    ## assayForTests
    args <- args0
    args$assayForTests <- 1
    expect_error(do.call(runPTMTest, args),
                 "'assayForTests' must be of class 'character'")
    args$assayForTests <- c("log2_Abundance_norm", "imputed_Abundance")
    expect_error(do.call(runPTMTest, args),
                 "'assayForTests' must have length 1")
    args$assayForTests <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "All values in 'assayForTests' must be one of")

    ## assayImputation
    args <- args0
    args$assayImputation <- 1
    expect_error(do.call(runPTMTest, args),
                 "'assayImputation' must be of class 'character'")
    args$assayImputation <- c("log2_Abundance_norm", "imputed_Abundance")
    expect_error(do.call(runPTMTest, args),
                 "'assayImputation' must have length 1")
    args$assayImputation <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "All values in 'assayImputation' must be one of")

    ## minnbrValidValues
    args <- args0
    args$minNbrValidValues <- "2"
    expect_error(do.call(runPTMTest, args),
                 "'minNbrValidValues' must be of class 'numeric'")
    args$minNbrValidValues <- c(1, 2)
    expect_error(do.call(runPTMTest, args),
                 "'minNbrValidValues' must have length 1")
    args$minNbrValidValues <- -1
    expect_error(do.call(runPTMTest, args),
                 "'minNbrValidValues' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## minlFC
    args <- args0
    args$minlFC <- "2"
    expect_error(do.call(runPTMTest, args),
                 "'minlFC' must be of class 'numeric'")
    args$minlFC <- c(1, 2)
    expect_error(do.call(runPTMTest, args),
                 "'minlFC' must have length 1")
    args$minlFC <- -1
    expect_error(do.call(runPTMTest, args),
                 "'minlFC' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## volcanoAdjPvalThr
    args <- args0
    args$volcanoAdjPvalThr <- "2"
    expect_error(do.call(runPTMTest, args),
                 "'volcanoAdjPvalThr' must be of class 'numeric'")
    args$volcanoAdjPvalThr <- c(0.1, 0.2)
    expect_error(do.call(runPTMTest, args),
                 "'volcanoAdjPvalThr' must have length 1")
    args$volcanoAdjPvalThr <- -1
    expect_error(do.call(runPTMTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)
    args$volcanoAdjPvalThr <- 2
    expect_error(do.call(runPTMTest, args),
                 "'volcanoAdjPvalThr' must be within [0,1] (inclusive)",
                 fixed = TRUE)

    ## volcanoLog2FCThr
    args <- args0
    args$volcanoLog2FCThr <- "2"
    expect_error(do.call(runPTMTest, args),
                 "'volcanoLog2FCThr' must be of class 'numeric'")
    args$volcanoLog2FCThr <- c(0.1, 0.2)
    expect_error(do.call(runPTMTest, args),
                 "'volcanoLog2FCThr' must have length 1")
    args$volcanoLog2FCThr <- -1
    expect_error(do.call(runPTMTest, args),
                 "'volcanoLog2FCThr' must be within [0,Inf] (inclusive)",
                 fixed = TRUE)

    ## baseFileName
    args <- args0
    args$baseFileName <- 1
    expect_error(do.call(runPTMTest, args),
                 "'baseFileName' must be of class 'character'")
    args$baseFileName <- c("pat1", "pat2")
    expect_error(do.call(runPTMTest, args),
                 "'baseFileName' must have length 1")

    ## singleFit
    args <- args0
    args$singleFit <- "2"
    expect_error(do.call(runPTMTest, args),
                 "'singleFit' must be of class 'logical'")
    args$singleFit <- c(TRUE, FALSE)
    expect_error(do.call(runPTMTest, args),
                 "'singleFit' must have length 1")

    ## subtractBaseline
    args <- args0
    args$subtractBaseline <- "1"
    expect_error(do.call(runPTMTest, args),
                 "'subtractBaseline' must be of class 'logical'")
    args$subtractBaseline <- c(TRUE, FALSE)
    expect_error(do.call(runPTMTest, args),
                 "'subtractBaseline' must have length 1")

    ## baselineGroup
    args <- args0
    args$subtractBaseline <- TRUE
    args$baselineGroup <- 1
    expect_error(do.call(runPTMTest, args),
                 "'baselineGroup' must be of class 'character'")
    args$baselineGroup <- "missing"
    expect_error(do.call(runPTMTest, args),
                 "baselineGroup %in% intersect(sceProteins$group, scePeptides$group) is not TRUE", fixed = TRUE)
    args$baselineGroup <- "WT"
    expect_error(do.call(runPTMTest, args),
                 '"batch" %in% colnames(SummarizedExperiment::colData(sceProteins)) is not TRUE', fixed = TRUE)

    ## Fails if duplicated comparison names
    args <- args0
    args$comparisons <- list(c("HIS4KO", "WT"),
                             WT_vs_HIS4KO = c("MET6KO", "WT"))
    expect_error(do.call(runPTMTest, args),
                 "Duplicated comparison names not allowed")

    ## Fails if groupComposition contains nonexistent group
    args <- args0
    args$groupComposition <- list(tmpgroup = c("HIS4KO", "missing"))
    args$comparisons <- list(c("MET6KO", "tmpgroup"))
    expect_error(do.call(runPTMTest, args), "Missing group(s) in sceProteins/scePeptides$groups", fixed = TRUE)

    ## Fails if the groups overlap
    args <- args0
    args$groupComposition <- list(tmpgroup = c("HIS4KO", "WT"),
                                  HIS4KO = "HIS4KO")
    args$comparisons <- list(c("tmpgroup", "HIS4KO"))
    expect_error(do.call(runPTMTest, args),
                 "The same original group is part of both groups")

    ## Fails for specific combinations of input arguments
    args <- args0
    args$testType <- "welch"
    args$minlFC <- 1
    expect_error(do.call(runPTMTest, args),
                 "The 'welch' test is currently not implemented for minlFC > 0")

    ## =========================================================================
    ## Works with correct arguments
    ## =========================================================================
    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0, no batch, singleFit = FALSE, no merged groups
    ## -------------------------------------------------------------------------
    out <- do.call(runPTMTest, args0)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, "WT_vs_HIS4KO")
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.037949704),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.8816754, 0.1258899, 1.0522788),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.09583327, 0.90285744, 0.32275149),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test
    ## same as above, but reorder columns in protein SCE
    ## minlFC = 0, no batch, singleFit = FALSE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    set.seed(123L)
    args$sceProteins <- args$sceProteins[, sample(seq_len(ncol(args$sceProteins)),
                                                  ncol(args$sceProteins),
                                                  replace = FALSE)]
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, "WT_vs_HIS4KO")
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.037949704),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.8816754, 0.1258899, 1.0522788),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.09583327, 0.90285744, 0.32275149),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checking
    ## partly different samples in protein and peptide SCE
    ## minlFC = 0, no batch, singleFit = FALSE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$sceProteins <- args$sceProteins[, 2:16] ## remove first HIS4KO
    args$scePeptides <- args$scePeptides[, 1:15] ## remove last WT
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, "WT_vs_HIS4KO")
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.03533889, 0.01606748, 0.04301649),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(0.7420381, 0.3548616, 1.1122455),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.4799475, 0.7321567, 0.2993434),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0, no batch, singleFit = TRUE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$singleFit <- TRUE
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.037949704),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.5535080, 0.1271131, 1.3242122),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.1423730, 0.9006392, 0.2064207),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checking
    ## minlFC = 0, no batch, singleFit = TRUE, no merged groups
    ## different samples in prot and pept
    ## -------------------------------------------------------------------------
    args <- args0
    args$sceProteins <- args$sceProteins[, 2:16] ## remove first HIS4KO
    args$scePeptides <- args$scePeptides[, 1:15] ## remove last WT
    args$singleFit <- TRUE
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.03533889, 0.01606748, 0.04301649),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(0.6374745, 0.3643189, 1.4379420),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.5353549, 0.7217279, 0.1750898),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - testing
    ## minlFC = 0, no batch, singleFit = TRUE, no merged groups
    ## don't filter by imputation status (assayImputation = NULL)
    ## -------------------------------------------------------------------------
    args <- args0
    args$singleFit <- TRUE
    args$assayImputation <- NULL
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.105294187),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.4841831, 0.1212044, 1.7184080),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.15927558, 0.90520093, 0.10709626),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 80)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 0)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - testing
    ## minlFC = 0, no batch, singleFit = TRUE, no merged groups
    ## don't filter by imputation status (minNbrValidValues = 0)
    ## -------------------------------------------------------------------------
    args <- args0
    args$singleFit <- TRUE
    args$minNbrValidValues <- 0
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.105294187),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.4841831, 0.1212044, 1.7184080),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.15927558, 0.90520093, 0.10709626),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 80)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 0)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0, batch (but will not be used), singleFit = TRUE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$singleFit <- TRUE
    tmpprot <- sce_proteins
    tmpprot$batch <- sample(LETTERS[1:2], ncol(tmpprot), replace = TRUE)
    tmppept <- sce_peptides
    tmppept$batch <- sample(LETTERS[1:2], ncol(tmppept), replace = TRUE)
    args$scePeptides <- tmppept
    args$sceProteins <- tmpprot
    expect_warning({
        out <- do.call(runPTMTest, args)
    }, "The 'interaction' test will currently not use the batch information")
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.037949704),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.5535080, 0.1271131, 1.3242122),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.1423730, 0.9006392, 0.2064207),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0.04, no batch, singleFit = TRUE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$singleFit <- TRUE
    args$minlFC <- 0.04
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_type(out$design$contrasts, "list")
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$WT_vs_HIS4KO, c(rep(0, 16), -1, 0, 0, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.04)")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.105294187),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "t"],
                 c(0.7220503, 0.0000000, 1.0905134),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.2568110, 0.9472806, 0.1614355),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    # expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
    #                            out$tests[[1]]$df.total, lower.tail = FALSE),
    #              out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    # expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
    #              out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0, no batch, singleFit = FALSE, merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$groupComposition <- list(his_met = c("HIS4KO", "MET6KO"))
    args$comparisons <- list(c("HIS4KO", "WT"), c("WT", "his_met"))
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 2)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$plotsubtitles[[2]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("WT_vs_HIS4KO", "his_met_vs_WT"))
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_type(out$design$his_met_vs_WT, "list")
    expect_named(out$design$his_met_vs_WT, c("design", "sampleData", "contrast"))
    expect_named(out$design$his_met_vs_WT$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$his_met_vs_WT$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_equal(nrow(out$tests[[2]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[2]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$tests[[2]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plotnotes[[2]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    expect_equal(out$plottitles[[2]], "his_met vs WT, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.037949704),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.8816754, 0.1258899, 1.0522788),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.09583327, 0.90285744, 0.32275149),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # second contrast
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[2]]$pid), "logFC"],
                 c(-0.08803610, -0.01335215, -0.04029302),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[2]]$pid), "t"],
                 c(-2.3268713, -0.4193856, -1.4923621),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[2]]$pid), "P.Value"],
                 c(0.03800047, 0.68223814, 0.16105265),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[2]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[2]]$t)), 28)
    expect_equal(sum(out$tests[[2]]$showInVolcano, na.rm = TRUE), 1)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0, no batch, singleFit = TRUE, merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$groupComposition <- list(his_met = c("HIS4KO", "MET6KO"))
    args$comparisons <- list(c("HIS4KO", "MET6KO"), c("WT", "his_met"))
    args$singleFit <- TRUE
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 2)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_s3_class(out$tests[[2]], "data.frame")
    expect_type(out$plotnotes[[2]], "character")
    expect_type(out$plottitles[[2]], "character")
    expect_type(out$plotsubtitles[[2]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, c("design", "sampleData", "contrasts"))
    expect_named(out$design$sampleData, c("sample", "group", "dataLevel"))
    expect_equal(out$design$contrasts$MET6KO_vs_HIS4KO, c(rep(0, 16), -1, 1, 0, 0))
    expect_equal(out$design$contrasts$his_met_vs_WT, c(rep(0, 16), 0.5, 0.5, 0, -1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_equal(nrow(out$tests[[2]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[2]])))
    expect_equal(out$tests[[2]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "MET6KO vs HIS4KO, limma")
    expect_equal(out$plotnotes[[2]], "")
    expect_equal(out$plottitles[[2]], "his_met vs WT, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(-0.026599011, -0.017582252, -0.004686624),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(-0.5528988, -0.4900067, -0.1635345),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.5889597, 0.6316376, 0.8724099),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 0)
    ## second contrast
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(-0.08803610, -0.01335215, -0.04029302),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "t"],
                 c(-2.1130528, -0.4296831, -1.6234853),
                 tolerance = 0.001)
    expect_equal(out$tests[[2]][match(c("P0CX40.42", "P26785.126", "P05750.411"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.05282949, 0.67389594, 0.12654154),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[2]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[2]]$t)), 28)
    expect_equal(sum(out$tests[[2]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df.total, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## interaction test - checked
    ## minlFC = 0.04, no batch, singleFit = FALSE, no merged groups
    ## -------------------------------------------------------------------------
    args <- args0
    args$minlFC <- 0.04
    out <- do.call(runPTMTest, args)
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_named(out$design, "WT_vs_HIS4KO")
    expect_type(out$design$WT_vs_HIS4KO, "list")
    expect_named(out$design$WT_vs_HIS4KO, c("design", "sampleData", "contrast"))
    expect_named(out$design$WT_vs_HIS4KO$sampleData, c("sample", "dataLevel", "group"))
    expect_equal(out$design$WT_vs_HIS4KO$contrast, c(0, 0, 0, 0, 0, 0, 0, 0, -1, 1))
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma treat (H0: |log2FC| <= 0.04)")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.105294187),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "t"],
                 c(0.8745782, 0.0000000, 1.0825026),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.2132518, 0.9484406, 0.1759112),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    ## Check consistency of values
    ## logFC +/- t * se = CI.R/CI.L
    expect_equal(out$tests[[1]]$logFC + qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.R, ignore_attr = TRUE)
    expect_equal(out$tests[[1]]$logFC - qt(p = 0.975, df = out$tests[[1]]$df.total) *
                     out$tests[[1]]$se.logFC,
                 out$tests[[1]]$CI.L, ignore_attr = TRUE)
    ## p-values
    # expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
    #                            out$tests[[1]]$df.total, lower.tail = FALSE),
    #              out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    # expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
    #              out$tests[[1]]$t, ignore_attr = TRUE)

    ## -------------------------------------------------------------------------
    ## welch test - checked
    ## minlFC = 0, no batch, singleFit = FALSE, no merged groups
    ## write results to file
    ## -------------------------------------------------------------------------
    args <- args0
    args$testType <- "welch"
    args$baseFileName <- tempfile()
    out <- do.call(runPTMTest, args)
    expect_true(file.exists(paste0(args$baseFileName, "_ptmtestres_WT_vs_HIS4KO.txt")))
    expect_type(out, "list")
    expect_length(out, 6)
    expect_named(out, c("plottitles", "plotsubtitles", "plotnotes",
                        "tests", "messages", "design"))
    expect_length(out$tests, 1)
    expect_s3_class(out$tests[[1]], "data.frame")
    expect_type(out$plotnotes[[1]], "character")
    expect_type(out$plottitles[[1]], "character")
    expect_type(out$plotsubtitles[[1]], "character")
    expect_type(out$design, "list")
    expect_equal(out$design, list())
    expect_equal(nrow(out$tests[[1]]), 80)
    expect_true(all(c("adj.P.Val", "showInVolcano") %in%
                        colnames(out$tests[[1]])))
    expect_equal(out$tests[[1]]$pid, rownames(args0$scePeptides))
    expect_equal(out$plotnotes[[1]], "")
    expect_equal(out$plottitles[[1]], "WT vs HIS4KO, limma")
    ## Compare to values calculated manually
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "logFC"],
                 c(0.074736597, 0.004561029, 0.105294187),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "t"],
                 c(1.13347625, 0.09635668, 2.05349057),
                 tolerance = 0.001)
    expect_equal(out$tests[[1]][match(c("P0CX40.42", "P26785.126", "Q02326.371"),
                                      out$tests[[1]]$pid), "P.Value"],
                 c(0.26273592, 0.92365081, 0.04561944),
                 tolerance = 0.001)
    expect_equal(sum(!is.na(out$tests[[1]]$t)), 52)
    expect_equal(sum(is.na(out$tests[[1]]$t)), 28)
    expect_equal(sum(out$tests[[1]]$showInVolcano, na.rm = TRUE), 1)
    # ## Check consistency of values
    ## p-values
    expect_equal(2 * stats::pt(abs(out$tests[[1]]$t),
                               out$tests[[1]]$df, lower.tail = FALSE),
                 out$tests[[1]]$P.Value, ignore_attr = TRUE)
    ## t-statistics
    expect_equal(out$tests[[1]]$logFC / out$tests[[1]]$se.logFC,
                 out$tests[[1]]$t, ignore_attr = TRUE)

})
