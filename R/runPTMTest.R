#' @keywords internal
#' @noRd
#' @importFrom SummarizedExperiment rowData colData assayNames
.checkArgumentsPTMTest <- function(sceProteins, scePeptides, matchColProteins,
                                   matchColPeptides, testType, comparisons,
                                   groupComposition, assayForTests,
                                   assayImputation, minNbrValidValues,
                                   minlFC, volcanoAdjPvalThr,
                                   volcanoLog2FCThr, baseFileName,
                                   singleFit, subtractBaseline,
                                   baselineGroup) {
    .assertVector(x = sceProteins, type = "SummarizedExperiment")
    .assertVector(x = scePeptides, type = "SummarizedExperiment")
    .assertScalar(x = matchColProteins, type = "character",
                  validValues = colnames(SummarizedExperiment::rowData(sceProteins)))
    .assertScalar(x = matchColPeptides, type = "character",
                  validValues = colnames(SummarizedExperiment::rowData(scePeptides)))
    stopifnot("group" %in% colnames(SummarizedExperiment::colData(sceProteins)),
              "group" %in% colnames(SummarizedExperiment::colData(scePeptides)))
    .assertVector(x = sceProteins$group, type = "character")
    .assertVector(x = scePeptides$group, type = "character")
    .assertVector(x = comparisons, type = "list")
    .assertVector(x = groupComposition, type = "list", allowNULL = TRUE)
    for (comparison in comparisons) {
        .assertVector(x = comparison, type = "character", len = 2,
                      validValues = unique(c(intersect(sceProteins$group,
                                                       scePeptides$group),
                                             names(groupComposition))))
    }
    .assertScalar(x = testType, type = "character",
                  validValues = c("interaction", "welch"))
    .assertScalar(x = assayForTests, type = "character",
                  validValues = intersect(
                      SummarizedExperiment::assayNames(sceProteins),
                      SummarizedExperiment::assayNames(scePeptides)))
    .assertScalar(x = assayImputation, type = "character",
                  validValues = intersect(
                      SummarizedExperiment::assayNames(sceProteins),
                      SummarizedExperiment::assayNames(scePeptides)),
                  allowNULL = TRUE)
    .assertScalar(x = minNbrValidValues, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = minlFC, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = volcanoAdjPvalThr, type = "numeric", rngIncl = c(0, 1))
    .assertScalar(x = volcanoLog2FCThr, type = "numeric", rngIncl = c(0, Inf))
    .assertScalar(x = baseFileName, type = "character", allowNULL = TRUE)
    .assertScalar(x = singleFit, type = "logical")
    .assertScalar(x = subtractBaseline, type = "logical")
    .assertScalar(x = baselineGroup, type = "character")
    if (subtractBaseline) {
        stopifnot(baselineGroup %in% intersect(sceProteins$group,
                                               scePeptides$group))
        stopifnot("batch" %in% colnames(SummarizedExperiment::colData(sceProteins)),
                  "batch" %in% colnames(SummarizedExperiment::colData(scePeptides)))
    }

    ## rownames(sceProteins) and rownames(scePeptides) must be unique
    ## (otherwise limma will remove the row names from the result table)
    if (any(duplicated(rownames(sceProteins)))) {
        stop("The row names of sceProteins cannot contain duplicated entries.")
    }
    if (any(duplicated(rownames(scePeptides)))) {
        stop("The row names of scePeptides cannot contain duplicated entries.")
    }

    ## If comparisons is not named, add names
    comparisons <- .assignNamesToComparisons(comparisons)

    ## If there are entries in unlist(comparisons) that are not defined in
    ## groupComposition, add them to the latter
    stdf <- setdiff(unlist(comparisons), names(groupComposition))
    groupComposition <- c(groupComposition,
                          setNames(as.list(stdf), stdf))

    ## Check that all entries in groupComposition[comparisons] are in the
    ## group column
    stdf <- setdiff(unlist(groupComposition[unlist(comparisons)]),
                    intersect(sceProteins$group, scePeptides$group))
    if (length(stdf) > 0) {
        stop("Missing group(s) in sceProteins/scePeptides$groups: ",
             paste(stdf, collapse = ", "))
    }

    if (minlFC > 0 && testType == "welch") {
        stop("The 'welch' test is currently not implemented for minlFC > 0")
    }
    if (("batch" %in% colnames(SummarizedExperiment::colData(sceProteins)) ||
         "batch" %in% colnames(SummarizedExperiment::colData(scePeptides))) &&
        testType == "interaction") {
        warning("The 'interaction' test will currently not use the batch information")
    }

    ## Return (possibly modified) comparisons and groupComposition
    return(list(comparisons = comparisons, groupComposition = groupComposition))
}

#' Perform PTM test
#'
#' Perform a test for post-translational modifications. If \code{testType} is
#' "welch", the functions applies a Welch t-test to the log-fold changes
#' obtained from independent tests on the peptide and protein levels.
#' This effectively 'adjusts' the changes on the peptide level for the changes
#' seen in the corresponding protein. This approach is similar to the method
#' implemented in the \code{MSstatsPTM} package (Kohler et al 2022).
#' If \code{testType} is "interaction", data from the peptide and protein
#' level are concatenated, and a model is fit to test for the significance
#' of the interaction between the "value type" and the condition, i.e.,
#' whether the difference between the groups depend on whether we are
#' considering the peptide or the protein-level abundances.
#'
#' @param sceProteins A \code{SummarizedExperiment} object (or a derivative)
#'     with protein-level abundances.
#' @param scePeptides A \code{SummarizedExperiment} object (or a derivative)
#'     with peptide-level abundances.
#' @param matchColProteins,matchColPeptides Character scalars indicating
#'     columns of \code{rowData(sceProteins)} and \code{rowData(scePeptides)},
#'     respectively, that will be used to extract matching record pairs.
#'     Typically, this will be a column with the protein identifier.
#' @param testType Either "welch" or "interaction", the type of test to
#'     perform. See Details for a description.
#' @param comparisons A list of character vectors of length 2, each giving the
#'     two groups to be compared.
#' @param groupComposition A list providing the composition of each group
#'     used in any of the comparisons. If \code{NULL}, assumes that each
#'     group used in \code{comparisons} consists of a single group in the
#'     \code{group} column of \code{colData(sceProteins)} and
#'     \code{colData(scePeptides)}.
#' @param assayForTests Character scalar, the name of an assay of the
#'     \code{SummarizedExperiment} object with values that will be used to
#'     perform the test.
#' @param assayImputation Character scalar, the name of an assay of
#'     \code{sce} with logical values indicating whether an entry was imputed
#'     or not.
#' @param minNbrValidValues Numeric scalar, the minimum number of valid
#'     (non-imputed) values that must be present for a features to include it
#'     in the result table.
#' @param minlFC Non-negative numeric scalar, the logFC threshold to use for
#'     limma-treat. If \code{minlFC} = 0, \code{limma::eBayes} is used instead.
#' @param volcanoAdjPvalThr Numeric scalar giving the FDR threshold for
#'     significance (for later use in volcano plots).
#' @param volcanoLog2FCThr Numeric scalar giving the logFC threshold for
#'     significance (for later use in volcano plots).
#' @param baseFileName Character scalar or \code{NULL}, the base file name of
#'     the output text files. If \code{NULL}, no result files are generated.
#' @param singleFit Logical scalar, whether to fit a single model to the full
#'     data set and extract relevant results using contrasts. If \code{FALSE},
#'     the data set will be subset for each comparison to only the relevant
#'     samples.
#' @param subtractBaseline Logical scalar, whether to subtract the background/
#'     reference value for each feature in each batch before fitting the
#'     model. If \code{TRUE}, requires that a 'batch' column is available.
#' @param baselineGroup Character scalar representing the reference group.
#'     Only used if \code{subtractBaseline} is \code{TRUE}, in which case the
#'     abundance values for a given sample will be adjusted by subtracting the
#'     average value across all samples in the \code{baselineGroup} from the
#'     same batch as the original sample.
#'
#' @references
#' Kohler D, Tsai T-H, Vershueren E, Huang T, Hinkle T, Phu L, Choi M,
#' Vitek O: MSstatsPTM: Statistical relative quantification of post-translational
#' modifiations in bottom-up mass spectrometry-based proteomics. Molecular and
#' Cellular Proteomics (2022).
#'
#' @export
#' @author Charlotte Soneson
#'
#' @importFrom SummarizedExperiment assay colData rowData
#' @importFrom dplyr select mutate left_join
#' @importFrom tibble rownames_to_column
#' @importFrom limma lmFit contrasts.fit eBayes topTable topTreat treat
#' @importFrom S4Vectors combineCols
#' @importFrom tidyselect any_of
#'
runPTMTest <- function(sceProteins, scePeptides, matchColProteins,
                       matchColPeptides, testType, comparisons,
                       groupComposition = NULL, assayForTests,
                       assayImputation = NULL, minNbrValidValues = 0,
                       minlFC = 0, volcanoAdjPvalThr = 0.05,
                       volcanoLog2FCThr = 1, baseFileName = NULL,
                       singleFit = FALSE, subtractBaseline = FALSE,
                       baselineGroup = "") {
    ## --------------------------------------------------------------------- ##
    ## Pre-flight checks
    ## --------------------------------------------------------------------- ##
    chk <- .checkArgumentsPTMTest(
        sceProteins = sceProteins, scePeptides = scePeptides,
        matchColProteins = matchColProteins, matchColPeptides = matchColPeptides,
        testType = testType, comparisons = comparisons,
        groupComposition = groupComposition, assayForTests = assayForTests,
        assayImputation = assayImputation, minNbrValidValues = minNbrValidValues,
        minlFC = minlFC, volcanoAdjPvalThr = volcanoAdjPvalThr,
        volcanoLog2FCThr = volcanoLog2FCThr, baseFileName = baseFileName,
        singleFit = singleFit, subtractBaseline = subtractBaseline,
        baselineGroup = baselineGroup)
    comparisons <- chk$comparisons
    groupComposition <- chk$groupComposition

    ## --------------------------------------------------------------------- ##
    ## Create row-matched objects
    ## --------------------------------------------------------------------- ##
    ## Store feature names from peptide SE - will be used to initiate the
    ## final result table
    all_peptide_ids <- rownames(scePeptides)

    shared <- intersect(
        SummarizedExperiment::rowData(scePeptides)[[matchColPeptides]],
        SummarizedExperiment::rowData(sceProteins)[[matchColProteins]])
    scePeptides <- scePeptides[which(
        SummarizedExperiment::rowData(scePeptides)[[matchColPeptides]] %in%
            shared), ]
    sceProteins <- sceProteins[match(
        SummarizedExperiment::rowData(scePeptides)[[matchColPeptides]],
        SummarizedExperiment::rowData(sceProteins)[[matchColPeptides]]), ]

    ## --------------------------------------------------------------------- ##
    ## Initialize result lists
    ## --------------------------------------------------------------------- ##
    returndesign <- list()
    tests <- list()
    plottitles <- list()
    plotsubtitles <- list()
    plotnotes <- list()
    messages <- list()

    ## --------------------------------------------------------------------- ##
    ## Run test
    ## --------------------------------------------------------------------- ##
    if (testType == "interaction") {
        ## Interaction test - create a merged object
        colnames(scePeptides) <- paste0(colnames(scePeptides), "_peptide")
        scePeptides$dataLevel <- "peptide"
        colnames(sceProteins) <- paste0(colnames(sceProteins), "_protein")
        sceProteins$dataLevel <- "protein"
        suppressWarnings({
            sceMerged <- S4Vectors::combineCols(scePeptides, sceProteins,
                                                use.names = FALSE, delayed = FALSE,
                                                fill = NA)
        })
        exprvals <- SummarizedExperiment::assay(sceMerged, assayForTests)
        df <- SummarizedExperiment::colData(sceMerged)

        if (singleFit) {
            ## If singleFit, then fit a model to the entire data set, with a
            ## column for each sample, and a separate 'peptide effect' for
            ## each group
            design <- model.matrix(~ sample, data = df)
            for (gr in unique(df$group)) {
                design <- cbind(design, as.numeric(df$group == gr &
                                                       df$dataLevel == "peptide"))
                colnames(design)[ncol(design)] <- paste0(gr, "_pept")
            }
            returndesign <- list(design = design,
                                 sampleData = df[, c("sample", "group", "dataLevel")],
                                 contrasts = list())
            fit0 <- limma::lmFit(exprvals, design)
        }
    } else if (testType == "welch") {
        rownames(sceProteins) <- make.unique(rownames(sceProteins))
        rownames(scePeptides) <- make.unique(rownames(scePeptides))
        res_proteins <- runTest(
            sce = sceProteins, comparisons = comparisons,
            groupComposition = groupComposition, testType = "limma",
            assayForTests = assayForTests, assayImputation = assayImputation,
            minNbrValidValues = minNbrValidValues, minlFC = 0,
            featureCollections = list(), complexFDRThr = 0.1,
            volcanoAdjPvalThr = volcanoAdjPvalThr,
            volcanoLog2FCThr = volcanoLog2FCThr, baseFileName = NULL, seed = 1,
            samSignificance = FALSE, nperm = 0, volcanoS0 = 1,
            addAbundanceValues = FALSE, aName = assayForTests,
            singleFit = singleFit, subtractBaseline = subtractBaseline,
            baselineGroup = baselineGroup)
        res_peptides <- runTest(
            sce = scePeptides, comparisons = comparisons,
            groupComposition = groupComposition, testType = "limma",
            assayForTests = assayForTests, assayImputation = assayImputation,
            minNbrValidValues = minNbrValidValues, minlFC = 0,
            featureCollections = list(), complexFDRThr = 0.1,
            volcanoAdjPvalThr = volcanoAdjPvalThr,
            volcanoLog2FCThr = volcanoLog2FCThr, baseFileName = NULL, seed = 1,
            samSignificance = FALSE, nperm = 0, volcanoS0 = 1,
            addAbundanceValues = FALSE, aName = assayForTests,
            singleFit = singleFit, subtractBaseline = subtractBaseline,
            baselineGroup = baselineGroup)
    } else {
        ## Should never end up here
        #nocov start
        stop("Invalid testType: ", testType)
        #nocov end
    }

    ## Perform each comparison.
    ## If testType == "interaction" and singleFit == TRUE, the model is fit above
    ## If testType == "interaction" and singleFit == FALSE, fit the models below
    ## If testType == "welch", fit the separate models above
    for (comparisonName in names(comparisons)) {
        if (testType == "interaction") {
            comparison <- comparisons[[comparisonName]]
            sceProteinsSub <- sceProteins[, sceProteins$group %in%
                                              unlist(groupComposition[comparison])]
            scePeptidesSub <- scePeptides[, scePeptides$group %in%
                                              unlist(groupComposition[comparison])]
            sceMergedSub <- sceMerged[, sceMerged$group %in%
                                          unlist(groupComposition[comparison])]
            if (nrow(sceProteinsSub) != nrow(sceMergedSub) ||
                nrow(scePeptidesSub) != nrow(sceMergedSub)) {
                ## Should never end up in here
                #nocov start
                stop("Something went wrong - subsetted proteins and peptides ",
                     "SCEs are not the same size as subsetted merged SCE")
                #nocov end
            }
            ## Only consider features with at least a certain number of valid values
            if (!is.null(assayImputation)) {
                imputedvalsPeptides <- SummarizedExperiment::assay(
                    scePeptidesSub, assayImputation, withDimnames = TRUE)
                imputedvalsProteins <- SummarizedExperiment::assay(
                    sceProteinsSub, assayImputation, withDimnames = TRUE)
                keep <- (rowSums(!imputedvalsPeptides) >= minNbrValidValues) &
                    (rowSums(!imputedvalsProteins) >= minNbrValidValues)
            } else {
                keep <- rep(TRUE, nrow(sceMergedSub))
            }

            if (singleFit) {
                fit <- fit0[keep, ]
                contrast <-
                    (1/length(groupComposition[[comparison[2]]])) *
                    (colnames(design) %in%
                         paste0(groupComposition[[comparison[2]]], "_pept")) -
                    (1/length(groupComposition[[comparison[1]]])) *
                    (colnames(design) %in%
                         paste0(groupComposition[[comparison[1]]], "_pept"))
                returndesign$contrasts[[comparisonName]] <- contrast
                fit <- limma::contrasts.fit(fit, contrasts = contrast)
                if (minlFC == 0) {
                    fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                    res <- limma::topTable(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                } else {
                    fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                    res <- limma::topTreat(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                }
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::mutate(s2.prior = fit$s2.prior,
                                  weights = fit$weights,
                                  sigma = fit$sigma,
                                  se.logFC = sqrt(fit$s2.post) *
                                      fit$stdev.unscaled[, 1],
                                  df.total = fit$df.total)
            } else {
                ## Create a new vector with the "merged" group names
                ## Check that it has length 2
                ## Also check that no column is in both groups
                c1 <- which(sceMergedSub$group %in% groupComposition[[comparison[1]]])
                c2 <- which(sceMergedSub$group %in% groupComposition[[comparison[2]]])
                if (length(intersect(c1, c2)) > 0) {
                    stop("The same original group is part of both groups ",
                         "to be compared")
                }
                if (!all(sort(union(c1, c2)) == seq_len(ncol(sceMergedSub)))) {
                    #nocov start
                    stop("Subsetting error - not all samples seem to have ",
                         "an assigned group")
                    #nocov end
                }
                fc <- rep(NA_character_, ncol(sceMergedSub))
                fc[c1] <- comparison[1]
                fc[c2] <- comparison[2]
                exprvals <- SummarizedExperiment::assay(sceMergedSub, assayForTests)
                exprvals <- exprvals[keep, , drop = FALSE]
                df <- SummarizedExperiment::colData(sceMergedSub)
                design <- model.matrix(~ sample, data = df)
                for (gr in comparison) {
                    design <- cbind(design, as.numeric(df$group %in% groupComposition[[gr]] &
                                                           df$dataLevel == "peptide"))
                    colnames(design)[ncol(design)] <- paste0(gr, "_pept")
                }

                contrast <- (colnames(design) == paste0(comparison[2], "_pept")) -
                    (colnames(design) == paste0(comparison[1], "_pept"))
                returndesign[[comparisonName]] <-
                    list(design = design,
                         sampleData = data.frame(df[, c("sample", "dataLevel")],
                                                 group = fc),
                         contrast = contrast)

                ## ------------------------------------------------------------- ##
                ## Run test
                ## ------------------------------------------------------------- ##
                fit <- limma::lmFit(exprvals, design)
                fit <- limma::contrasts.fit(fit, contrasts = contrast)
                if (minlFC == 0) {
                    fit <- limma::eBayes(fit, trend = TRUE, robust = FALSE)
                    res <- limma::topTable(fit, coef = 1, number = Inf,
                                           confint = TRUE, sort.by = "none")
                } else {
                    fit <- limma::treat(fit, fc = 2^minlFC, trend = TRUE, robust = FALSE)
                    res <- limma::topTreat(fit, coef = 1, confint = TRUE,
                                           number = Inf, sort.by = "none")
                }
                res <- res %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::mutate(s2.prior = fit$s2.prior,
                                  weights = fit$weights,
                                  sigma = fit$sigma,
                                  se.logFC = sqrt(fit$s2.post) *
                                      fit$stdev.unscaled[, 1],
                                  df.total = fit$df.total)
            }
        } else if (testType == "welch") {
            resPeptide <- res_peptides$tests[[comparisonName]]
            resProtein <- res_proteins$tests[[comparisonName]]
            res <- data.frame(pid = resPeptide$pid)
            res$logFC <- resPeptide$logFC - resProtein$logFC
            s2_1 <- resPeptide$se.logFC ^ 2
            s2_2 <- resProtein$se.logFC ^ 2
            res$se.logFC <- sqrt(s2_1 + s2_2)
            numer <- (s2_1 + s2_2) ^ 2
            denom <- (s2_1 ^ 2 / resPeptide$df.total +
                          s2_2 ^ 2 / resProtein$df.total)
            res$df <- numer / denom
            res$t <- res$logFC / res$se.logFC
            res$P.Value <- 2 * stats::pt(abs(res$t), res$df, lower.tail = FALSE)
            res$adj.P.Val <- p.adjust(res$P.Value, method = "BH")
        } else {
            ## Should never end up here
            #nocov start
            stop("Invalid testType: ", testType)
            #nocov end
        }

        ## Calculate -log10(p). For features with p-value = 0, use half the
        ## smallest non-zero p-value as a proxy to be able to make volcano plots.
        res <- res %>%
            dplyr::mutate(mlog10p = -log10(.data$P.Value)) %>%
            dplyr::mutate(mlog10p = replace(
                .data$mlog10p, .data$P.Value == 0,
                -log10(min(.data$P.Value[which(.data$P.Value > 0)])/2))) %>%
            dplyr::left_join(as.data.frame(
                SummarizedExperiment::rowData(scePeptides)) %>%
                    tibble::rownames_to_column("pid") %>%
                    dplyr::select(tidyselect::any_of(
                        c("pid", "einprotGene", "einprotProtein",
                          "einprotLabel"))),
                by = "pid")

        ## ----------------------------------------------------------------- ##
        ## Determine significant features
        ## ----------------------------------------------------------------- ##
        res <- data.frame(pid = all_peptide_ids) %>%
            dplyr::left_join(res, by = "pid")
        res$showInVolcano <- res$adj.P.Val <= volcanoAdjPvalThr &
            abs(res$logFC) >= volcanoLog2FCThr

        ## ----------------------------------------------------------------- ##
        ## Write results to file
        ## ----------------------------------------------------------------- ##
        if (!is.null(baseFileName)) {
            write.table(res %>%
                            dplyr::filter(.data$showInVolcano) %>%
                            dplyr::arrange(desc(.data$logFC)),
                        file = paste0(baseFileName, "_ptmtestres_", comparisonName,
                                      ".txt"),
                        row.names = FALSE, col.names = TRUE,
                        quote = FALSE, sep = "\t")
        }

        ## ----------------------------------------------------------------- ##
        ## Generate return values
        ## ----------------------------------------------------------------- ##
        if (minlFC == 0) {
            plottitle <- paste0(sub("_vs_", " vs ", comparisonName), ", limma")
        } else {
            plottitle <- paste0(sub("_vs_", " vs ", comparisonName),
                                ", limma treat (H0: |log2FC| <= ", minlFC, ")")
        }
        # plotnote <- paste0("df.prior = ", round(fit$df.prior, digits = 2))
        plotnote <- ""
        plotsubtitle <- paste0("Adj.p threshold = ", volcanoAdjPvalThr,
                               ", |log2FC| threshold = ", volcanoLog2FCThr)

        ## ----------------------------------------------------------------- ##
        ## Populate result lists
        ## ----------------------------------------------------------------- ##
        plottitles[[comparisonName]] <- plottitle
        plotsubtitles[[comparisonName]] <- plotsubtitle
        plotnotes[[comparisonName]] <- plotnote
        tests[[comparisonName]] <- res
    }

    return(list(plottitles = plottitles, plotsubtitles = plotsubtitles,
                plotnotes = plotnotes, tests = tests,
                messages = messages, design = returndesign))
}
