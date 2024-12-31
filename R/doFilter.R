#' @keywords internal
#' @noRd
#'
#' @importFrom ggplot2 ggplot aes geom_col labs
#' @importFrom cowplot theme_cowplot
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr summarize everything
#' @importFrom ComplexUpset upset
#'
.makeFilterPlot <- function(filtdf, plotUpset) {
    if (plotUpset && any(rowSums(filtdf) > 0)) {
        if (ncol(filtdf) > 1) {
            print(ComplexUpset::upset(filtdf[rowSums(filtdf) > 0, , drop = FALSE],
                                      intersect = colnames(filtdf)))
        } else {
            print(ggplot2::ggplot(
                data = filtdf %>%
                    dplyr::summarize(across(dplyr::everything(),
                                            function(x) length(which(x > 0)))) %>%
                    tidyr::pivot_longer(cols = dplyr::everything(),
                                        names_to = "criterion", values_to = "number"),
                ggplot2::aes(x = .data$criterion, y = .data$number)) +
                    ggplot2::geom_col() +
                    cowplot::theme_cowplot() +
                    ggplot2::labs(x = "", y = "Number of excluded features")
            )
        }
    }
}

#' Filter out features in MaxQuant data
#'
#' Exclude features with 'Score' below \code{minScore}, 'Peptides' below
#' \code{minPeptides}, or identified as either 'Reverse',
#' 'Potential.contaminant' or 'Only.identified.by.site' by \code{MaxQuant}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param minScore Numeric scalar, the minimum allowed value in the 'Score'
#'     column in order to retain the feature.
#' @param minPeptides Numeric scalar, the minimum allowed value in the
#'     'Peptides' column in order to retain the feature.
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @returns A filtered object of the same type as \code{sce}.
#'
#' @examples
#' sce <- importExperiment(inFile = system.file("extdata", "mq_example",
#'                                              "1356_proteinGroups.txt",
#'                                              package = "einprot"),
#'                         iColPattern = "^LFQ.intensity.")$sce
#'
#' dim(sce)
#' sce <- filterMaxQuant(sce = sce, minScore = 2, minPeptides = 2,
#'                       plotUpset = TRUE)
#' dim(sce)
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate across
#' @importFrom ComplexUpset upset
#' @importFrom rlang .data
#'
filterMaxQuant <- function(sce, minScore, minPeptides, plotUpset = TRUE,
                           exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = minScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    ## Make sure that the columns used for filtering are character vectors
    if ("Reverse" %in% colnames(rowData(sce))) {
        rowData(sce)$Reverse[is.na(rowData(sce)$Reverse)] <- ""
    }
    if ("Potential.contaminant" %in% colnames(rowData(sce))) {
        rowData(sce)$Potential.contaminant[is.na(
            rowData(sce)$Potential.contaminant)] <- ""
    }
    if ("Only.identified.by.site" %in% colnames(rowData(sce))) {
        rowData(sce)$Only.identified.by.site[is.na(
            rowData(sce)$Only.identified.by.site)] <- ""
    }

    filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
        dplyr::select(dplyr::any_of(c("Reverse", "Potential.contaminant",
                                      "Only.identified.by.site", "Score",
                                      "Peptides"))) %>%
        dplyr::mutate(across(dplyr::any_of(c("Reverse", "Potential.contaminant",
                                             "Only.identified.by.site")),
                             function(x) as.numeric(x == "+")))
    if ("Score" %in% colnames(filtdf) && !is.null(minScore)) {
        filtdf <- filtdf %>%
            dplyr::mutate(Score = as.numeric((.data$Score < minScore) |
                                                 is.na(.data$Score)))
    } else {
        filtdf$Score <- NULL
    }
    if ("Peptides" %in% colnames(filtdf) && !is.null(minPeptides)) {
        filtdf <- filtdf %>%
            dplyr::mutate(
                Peptides = as.numeric((.data$Peptides < minPeptides) |
                                          is.na(.data$Peptides)))
    } else {
        filtdf$Peptides <- NULL
    }

    keep <- seq_len(nrow(sce))
    if ("Reverse" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Reverse == ""))
    }
    if ("Potential.contaminant" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Potential.contaminant == ""))
    }
    if ("Only.identified.by.site" %in% colnames(rowData(sce))) {
        keep <- intersect(
            keep, which(rowData(sce)$Only.identified.by.site == ""))
    }
    if ("Score" %in% colnames(rowData(sce)) && !is.null(minScore)) {
        keep <- intersect(keep, which(rowData(sce)$Score >= minScore))
    }
    if ("Peptides" %in% colnames(rowData(sce)) && !is.null(minPeptides)) {
        keep <- intersect(keep, which(rowData(sce)$Peptides >= minPeptides))
    }
    exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
    sce <- sce[keep, ]

    if (nrow(filtdf[rowSums(filtdf) == 0, , drop = FALSE]) != nrow(sce)) {
        ## This should not happen
        #nocov start
        stop("Something went wrong in the filtering - filtdf and sce are of ",
             "different sizes")
        #nocov end
    }
    .makeFilterPlot(filtdf = filtdf, plotUpset = plotUpset)

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
    sce
}

#' Filter out features in PD/TMT data
#'
#' If \code{inputLevel} is "Proteins", exclude features with
#' 'Score.Sequest.HT.Sequest.HT' below \code{minScore},
#' 'Number.of.Peptides' below \code{minPeptides}, or identified as
#' 'Contaminant' by \code{ProteomeDiscoverer}.
#' If \code{inputLevel} is "PeptideGroups", exclude features with
#' 'Delta.Score.by.Search.Engine.Sequest.HT' below \code{minDeltaScore},
#' 'Number.of.PSMs' below \code{minPSMs}, or identified as
#' 'Contaminant' by \code{ProteomeDiscoverer}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param inputLevel Either \code{"Proteins"} or \code{"PeptideGroups"},
#'     indicating the type of features in \code{sce}.
#' @param minScore Numeric scalar, the minimum allowed value in the
#'     'Score.Sequest.HT.Sequest.HT' column in order to retain the feature.
#'     Only used if \code{inputLevel} is \code{"Proteins"}.
#' @param minPeptides Numeric scalar, the minimum allowed value in the
#'     'Number.of.Peptides' column in order to retain the feature.
#'     Only used if \code{inputLevel} is \code{"Proteins"}.
#' @param minDeltaScore Numeric scalar, the minimum allowed value in the
#'     'Delta.Score.by.Search.Engine.Sequest.HT' column in order to retain the
#'     feature. Only used if \code{inputLevel} is \code{"PeptideGroups"}.
#' @param minPSMs Numeric scalar, the minimum allowed value in the
#'     'Number.of.PSMs' column in order to retain the feature.
#'     Only used if \code{inputLevel} is \code{"PeptideGroups"}.
#' @param masterProteinsOnly Logical scalar indicating whether only master
#'     proteins (where the \code{Master} column value is
#'     \code{IsMasterProtein}) should be retained.
#' @param modificationsCol Character string pointing to a column containing
#'     modification details. \code{excludeUnmodifiedPeptides} and
#'     \code{keepModifications} will use information from this column. Only
#'     used if \code{inputLevel} is \code{"PeptideGroups"}.
#' @param excludeUnmodifiedPeptides Logical scalar, whether to filter out
#'     peptides without modifications. Only used if \code{inputLevel} is
#'     \code{"PeptideGroups"}.
#' @param keepModifications Character string (or \code{NULL}) indicating
#'     which modifications to retain in the analysis. Can be a regular
#'     expression, which will be matched against the \code{modificationsCol}.
#'     If \code{NULL} (the default), all rows are retained. Only used if
#'     \code{inputLevel} is \code{"PeptideGroups"}.
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @returns A filtered object of the same type as \code{sce}.
#'
#' @examples
#' ## Proteins
#' sce <- importExperiment(
#'     inFile = system.file("extdata", "pdtmt_example",
#'                          "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
#'                          package = "einprot"),
#'     iColPattern = "^Abundance.F.+.Sample.")$sce
#'
#' dim(sce)
#' sce <- filterPDTMT(sce = sce, inputLevel = "Proteins", minScore = 2,
#'                    minPeptides = 2, plotUpset = TRUE)
#' dim(sce)
#'
#' ## PeptideGroups
#' sce <- importExperiment(
#'     inFile = system.file("extdata", "pdtmt_example",
#'                          "Fig2_m23139_RTS_QC_varMods_PeptideGroups.txt",
#'                          package = "einprot"),
#'     iColPattern = "^Abundance.F.+.Sample.")$sce
#'
#' dim(sce)
#' sce <- filterPDTMT(sce = sce, inputLevel = "PeptideGroups",
#'                    minPSMs = 2, plotUpset = TRUE, minDeltaScore = 0.2,
#'                    modificationsCol = "Modifications.in.Master.Proteins",
#'                    excludeUnmodifiedPeptides = TRUE)
#' dim(sce)
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate
#' @importFrom ComplexUpset upset
#'
filterPDTMT <- function(sce, inputLevel, minScore = 0, minPeptides = 0,
                        minDeltaScore = 0, minPSMs = 0,
                        masterProteinsOnly = FALSE,
                        modificationsCol = "Modifications",
                        excludeUnmodifiedPeptides = FALSE,
                        keepModifications = NULL, plotUpset = TRUE,
                        exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = inputLevel, type = "character",
                  validValues = c("Proteins", "PeptideGroups"))
    .assertScalar(x = minScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minDeltaScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPSMs, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = masterProteinsOnly, type = "logical")
    if (inputLevel == "PeptideGroups") {
        .assertScalar(x = modificationsCol, type = "character",
                      validValues = colnames(rowData(sce)), allowNULL = TRUE)
        .assertScalar(x = excludeUnmodifiedPeptides, type = "logical")
        .assertScalar(x = keepModifications, type = "character",
                      allowNULL = TRUE)
    }
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    ## Make sure that the columns used for filtering are character vectors
    if ("Contaminant" %in% colnames(rowData(sce))) {
        rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant)
    }

    if (inputLevel == "Proteins") {
        filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
            dplyr::select(dplyr::any_of(c("Contaminant", "Number.of.Peptides",
                                          "Score.Sequest.HT.Sequest.HT",
                                          "Master"))) %>%
            dplyr::mutate(across(dplyr::any_of(c("Contaminant")),
                                 function(x) as.numeric(x == "True")))
        if ("Number.of.Peptides" %in% colnames(filtdf) &&
            !is.null(minPeptides)) {
            filtdf <- filtdf %>%
                dplyr::mutate(
                    Number.of.Peptides =
                        as.numeric((.data$Number.of.Peptides < minPeptides) |
                                       is.na(.data$Number.of.Peptides)))
        } else {
            filtdf$Number.of.Peptides <- NULL
        }
        if ("Score.Sequest.HT.Sequest.HT" %in% colnames(filtdf) &&
            !is.null(minScore)) {
            filtdf <- filtdf %>%
                dplyr::mutate(
                    Score.Sequest.HT.Sequest.HT =
                        as.numeric((.data$Score.Sequest.HT.Sequest.HT < minScore) |
                                       is.na(.data$Score.Sequest.HT.Sequest.HT)))
        } else {
            filtdf$Score.Sequest.HT.Sequest.HT <- NULL
        }
        if ("Master" %in% colnames(filtdf)) {
            filtdf <- filtdf %>%
                dplyr::mutate(
                    Master = as.numeric((.data$Master != "IsMasterProtein") |
                                            is.na(.data$Master)))
        }

        keep <- seq_len(nrow(sce))
        if ("Contaminant" %in% colnames(rowData(sce))) {
            keep <- intersect(keep, which(rowData(sce)$Contaminant == "False"))
        }
        if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce)) &&
            !is.null(minScore)) {
            keep <- intersect(
                keep,
                which(rowData(sce)$Score.Sequest.HT.Sequest.HT >= minScore)
            )
        }
        if ("Number.of.Peptides" %in% colnames(rowData(sce)) &&
            !is.null(minPeptides)) {
            keep <- intersect(
                keep, which(rowData(sce)$Number.of.Peptides >= minPeptides)
            )
        }
        if (masterProteinsOnly) {
            if ("Master" %in% colnames(rowData(sce))) {
                keep <- intersect(
                    keep, which(rowData(sce)$Master == "IsMasterProtein")
                )
            }
        } else {
            filtdf$Master <- NULL
        }
        exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
        sce <- sce[keep, ]
    } else if (inputLevel == "PeptideGroups") {
        colsToExtract <- c("Contaminant", "Number.of.PSMs",
                           "Delta.Score.by.Search.Engine.Sequest.HT")
        if (!is.null(modificationsCol)) {
            colsToExtract <- c(colsToExtract, modificationsCol)
        }
        filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
            dplyr::select(dplyr::any_of(colsToExtract)) %>%
            dplyr::mutate(across(dplyr::any_of(c("Contaminant")),
                                 function(x) as.numeric(x == "True")))
        if (!is.null(modificationsCol) && excludeUnmodifiedPeptides) {
            filtdf <- filtdf %>%
                dplyr::mutate(Unmodified.peptide =
                                  .data[[modificationsCol]] == "")
        } else {
            filtdf$Unmodified.peptide <- NULL
        }
        if (!is.null(modificationsCol) && !is.null(keepModifications)) {
            filtdf <- filtdf %>%
                dplyr::mutate(No.modification.tokeep =
                                  !grepl(keepModifications,
                                         .data[[modificationsCol]]))
        } else {
            filtdf$No.modification.tokeep <- NULL
        }
        if ("Number.of.PSMs" %in% colnames(filtdf) && !is.null(minPSMs)) {
            filtdf <- filtdf %>%
                dplyr::mutate(Number.of.PSMs =
                                  as.numeric((.data$Number.of.PSMs < minPSMs) |
                                                 is.na(.data$Number.of.PSMs)))
        } else {
            filtdf$Number.of.PSMs <- NULL
        }
        if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(filtdf) &&
            !is.null(minDeltaScore)) {
            filtdf <- filtdf %>%
                dplyr::mutate(
                    Delta.Score.by.Search.Engine.Sequest.HT =
                        as.numeric((.data$Delta.Score.by.Search.Engine.Sequest.HT <
                                        minDeltaScore) |
                                       is.na(.data$Delta.Score.by.Search.Engine.Sequest.HT)))
        } else {
            filtdf$Delta.Score.by.Search.Engine.Sequest.HT <- NULL
        }
        ## Remove the modifications column (not needed anymore)
        if (!is.null(modificationsCol) &&
            modificationsCol %in% colnames(filtdf)) {
            filtdf[[modificationsCol]] <- NULL
        }

        keep <- seq_len(nrow(sce))
        if ("Contaminant" %in% colnames(rowData(sce))) {
            keep <- intersect(keep, which(rowData(sce)$Contaminant == "False"))
        }
        if ("Delta.Score.by.Search.Engine.Sequest.HT" %in%
            colnames(rowData(sce)) &&
            !is.null(minDeltaScore)) {
            keep <- intersect(
                keep,
                which(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= minDeltaScore)
            )
        }
        if ("Number.of.PSMs" %in% colnames(rowData(sce)) && !is.null(minPSMs)) {
            keep <- intersect(
                keep, which(rowData(sce)$Number.of.PSMs >= minPSMs)
            )
        }
        if (!is.null(modificationsCol) && excludeUnmodifiedPeptides) {
            keep <- intersect(
                keep, which(rowData(sce)[[modificationsCol]] != "")
            )
        }
        if (!is.null(modificationsCol) && !is.null(keepModifications)) {
            keep <- intersect(
                keep, which(grepl(keepModifications,
                                  rowData(sce)[[modificationsCol]]))
            )
        }

        exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
        sce <- sce[keep, ]
    }

    if (nrow(filtdf[rowSums(filtdf) == 0, , drop = FALSE]) != nrow(sce)) {
        ## This should not happen - could end up here though if the values in
        ## the Contaminant column are not the expected ones
        stop("Something went wrong in the filtering - filtdf and sce are of ",
             "different sizes")
    }
    .makeFilterPlot(filtdf = filtdf, plotUpset = plotUpset)

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }

    sce
}

#' Filter out features in FragPipe data
#'
#' Exclude features with 'Combined.Total.Peptides' below \code{minPeptides},
#' or identified as either 'Reverse' (Protein name matching
#' \code{revPattern}) or 'Potential.contaminant' (Protein name starting
#' with \code{contam_}) by \code{FragPipe}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param minPeptides Numeric scalar, the minimum allowed value in the
#'     'Combined.Total.Peptides' column in order to retain the feature.
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#' @param revPattern Character scalar providing the pattern (a regular
#'     expression) used to identify decoys (reverse hits). The pattern is
#'     matched against the IDs in the FragPipe \code{Protein} column.
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @returns A filtered object of the same type as \code{sce}.
#'
#' @examples
#' sce <- importExperiment(inFile = system.file("extdata", "fp_example",
#'                                              "combined_protein.tsv",
#'                                              package = "einprot"),
#'                         iColPattern = ".MaxLFQ.Intensity$")$sce
#'
#' dim(sce)
#' sce <- filterFragPipe(sce = sce, minPeptides = 2,
#'                       plotUpset = TRUE,
#'                       revPattern = "^rev_")
#' dim(sce)
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate across
#' @importFrom ComplexUpset upset
#' @importFrom rlang .data
#'
filterFragPipe <- function(sce, minPeptides, plotUpset = TRUE,
                           revPattern = "^rev_", exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = revPattern, type = "character")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    ## Make sure that the columns used for filtering later are character vectors
    rowData(sce)$Potential.contaminant <-
        ifelse(grepl("^contam_", rowData(sce)$Protein), "+", "")
    rowData(sce)$Reverse <- ifelse(grepl(revPattern, rowData(sce)$Protein),
                                   "+", "")

    filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
        dplyr::select(dplyr::any_of(c("Reverse", "Potential.contaminant",
                                      "Combined.Total.Peptides"))) %>%
        dplyr::mutate(across(dplyr::any_of(c("Reverse",
                                             "Potential.contaminant")),
                             function(x) as.numeric(x == "+")))
    if ("Combined.Total.Peptides" %in% colnames(filtdf) &&
        !is.null(minPeptides)) {
        filtdf <- filtdf %>%
            dplyr::mutate(
                Combined.Total.Peptides = as.numeric(
                    (.data$Combined.Total.Peptides < minPeptides) |
                        is.na(.data$Combined.Total.Peptides)))
    } else {
        filtdf$Combined.Total.Peptides <- NULL
    }

    keep <- seq_len(nrow(sce))
    if ("Reverse" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Reverse == ""))
    }
    if ("Potential.contaminant" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Potential.contaminant == ""))
    }
    if ("Combined.Total.Peptides" %in% colnames(rowData(sce)) &&
        !is.null(minPeptides)) {
        keep <- intersect(
            keep, which(rowData(sce)$Combined.Total.Peptides >= minPeptides)
        )
    }
    exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
    sce <- sce[keep, ]

    if (nrow(filtdf[rowSums(filtdf) == 0, , drop = FALSE]) != nrow(sce)) {
        ## This should not happen
        #nocov start
        stop("Something went wrong in the filtering - filtdf and sce are of ",
             "different sizes")
        #nocov end
    }
    .makeFilterPlot(filtdf = filtdf, plotUpset = plotUpset)

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
    sce
}

#' Filter out features in Spectronaut data
#'
#' Exclude features with 'PG.Cscore' below \code{minScore},
#' 'PG.NrOfStrippedSequencesIdentified.Experiment.wide' below
#' \code{minPeptides}, or where the 'PG.ProteinGroups' column contains the
#' specified \code{revPattern}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param minScore Numeric scalar, the minimum allowed value in the 'PG.Cscore'
#'     column in order to retain the feature.
#' @param minPeptides Numeric scalar, the minimum allowed value in the
#'     'PG.NrOfStrippedSequencesIdentified.Experiment.wide' column in order to
#'     retain the feature.
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#' @param revPattern Character scalar providing the pattern (a regular
#'     expression) used to identify decoys (reverse hits). The pattern is
#'     matched against the IDs in the Spectronaut \code{PG.ProteinGroups}
#'     column.
#' @param contamPattern Character scalar providing the pattern (a regular
#'     expression) used to identify contaminants. The pattern is
#'     matched against the IDs in the Spectronaut \code{PG.ProteinGroups}
#'     column.
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @returns A filtered object of the same type as \code{sce}.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate across
#' @importFrom ComplexUpset upset
#' @importFrom rlang .data
#'
filterSpectronaut <- function(sce, minScore, minPeptides, plotUpset = TRUE,
                              revPattern = "_Decoy$",
                              contamPattern = "^contam_", exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = minScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = revPattern, type = "character")
    .assertScalar(x = contamPattern, type = "character")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    ## Make sure that the columns used for filtering later are character vectors
    rowData(sce)$Reverse <- ifelse(grepl(revPattern, rowData(sce)$PG.ProteinGroups),
                                   "+", "")
    rowData(sce)$Contaminant <- ifelse(grepl(contamPattern,
                                             rowData(sce)$PG.ProteinGroups),
                                       "+", "")

    filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
        dplyr::select(dplyr::any_of(c("Reverse", "PG.NrOfStrippedSequencesIdentified.Experiment.wide",
                                      "PG.Cscore", "Contaminant"))) %>%
        dplyr::mutate(across(dplyr::any_of(c("Reverse", "Contaminant")),
                             function(x) as.numeric(x == "+")))
    if ("PG.NrOfStrippedSequencesIdentified.Experiment.wide" %in% colnames(filtdf) &&
        !is.null(minPeptides)) {
        filtdf <- filtdf %>%
            dplyr::mutate(
                PG.NrOfStrippedSequencesIdentified.Experiment.wide = as.numeric(
                    (.data$PG.NrOfStrippedSequencesIdentified.Experiment.wide < minPeptides) |
                        is.na(.data$PG.NrOfStrippedSequencesIdentified.Experiment.wide)))
    } else {
        filtdf$PG.NrOfStrippedSequencesIdentified.Experiment.wide <- NULL
    }
    if ("PG.Cscore" %in% colnames(filtdf) && !is.null(minScore)) {
        filtdf <- filtdf %>%
            dplyr::mutate(PG.Cscore = as.numeric((.data$PG.Cscore < minScore) |
                                                     is.na(.data$PG.Cscore)))
    } else {
        filtdf$PG.Cscore <- NULL
    }

    keep <- seq_len(nrow(sce))
    if ("Reverse" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Reverse == ""))
    }
    if ("Contaminant" %in% colnames(rowData(sce))) {
        keep <- intersect(keep, which(rowData(sce)$Contaminant == ""))
    }
    if ("PG.NrOfStrippedSequencesIdentified.Experiment.wide" %in% colnames(rowData(sce)) &&
        !is.null(minPeptides)) {
        keep <- intersect(
            keep, which(rowData(sce)$PG.NrOfStrippedSequencesIdentified.Experiment.wide >= minPeptides)
        )
    }
    if ("PG.Cscore" %in% colnames(rowData(sce)) && !is.null(minScore)) {
        keep <- intersect(keep, which(rowData(sce)$PG.Cscore >= minScore))
    }
    exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
    sce <- sce[keep, ]

    if (nrow(filtdf[rowSums(filtdf) == 0, , drop = FALSE]) != nrow(sce)) {
        ## This should not happen
        #nocov start
        stop("Something went wrong in the filtering - filtdf and sce are of ",
             "different sizes")
        #nocov end
    }
    .makeFilterPlot(filtdf = filtdf, plotUpset = plotUpset)

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
    sce
}
