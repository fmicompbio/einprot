#' Filter out features in MaxQuant data
#'
#' Exclude features with 'Score' below \code{minScore}, 'Peptides' below
#' \code{minPeptides}, or identified as either 'Reverse',
#' 'Potential.contaminant' or 'Only.identified.by.site' by MaxQuant.
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
#' @return A filtered object of the same type as \code{sce}.
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
        rowData(sce)$Potential.contaminant[is.na(rowData(sce)$Potential.contaminant)] <- ""
    }
    if ("Only.identified.by.site" %in% colnames(rowData(sce))) {
        rowData(sce)$Only.identified.by.site[is.na(rowData(sce)$Only.identified.by.site)] <- ""
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
        keep <- intersect(keep, which(rowData(sce)$Only.identified.by.site == ""))
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
    if (plotUpset && any(rowSums(filtdf) > 0)) {
        print(ComplexUpset::upset(filtdf[rowSums(filtdf) > 0, , drop = FALSE],
                                  intersect = colnames(filtdf)))
    }

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
#' 'Contaminant' by ProteomeDiscoverer.
#' If \code{inputLevel} is "PeptideGroups", exclude features with
#' 'Delta.Score.by.Search.Engine.Sequest.HT' below \code{minDeltaScore},
#' 'Number.of.PSMs' below \code{minPSMs}, or identified as
#' 'Contaminant' by ProteomeDiscoverer.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param inputLevel Either "Proteins" or "PeptideGroups", indicating the type
#'     of features in \code{sce}.
#' @param minScore Numeric scalar, the minimum allowed value in the
#'     'Score.Sequest.HT.Sequest.HT' column in order to retain the feature.
#'     Only used if \code{inputLevel} is "Proteins".
#' @param minPeptides Numeric scalar, the minimum allowed value in the
#'     'Number.of.Peptides' column in order to retain the feature.
#'     Only used if \code{inputLevel} is "Proteins".
#' @param minDeltaScore Numeric scalar, the minimum allowed value in the
#'     'Delta.Score.by.Search.Engine.Sequest.HT' column in order to retain the
#'     feature. Only used if \code{inputLevel} is "PeptideGroups".
#' @param minPSMs Numeric scalar, the minimum allowed value in the
#'     'Number.of.PSMs' column in order to retain the feature.
#'     Only used if \code{inputLevel} is "PeptideGroups".
#' @param masterProteinsOnly Logical scalar indicating whether only master
#'     proteins (where the \code{Master} column value is
#'     \code{IsMasterProtein}) should be retained.
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @return A filtered object of the same type as \code{sce}.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate
#' @importFrom ComplexUpset upset
#'
filterPDTMT <- function(sce, inputLevel, minScore = 0, minPeptides = 0,
                        minDeltaScore = 0, minPSMs = 0,
                        masterProteinsOnly = FALSE, plotUpset = TRUE,
                        exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = inputLevel, type = "character",
                  validValues = c("Proteins", "PeptideGroups"))
    .assertScalar(x = minScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minDeltaScore, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = minPSMs, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = masterProteinsOnly, type = "logical")
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
        if ("Number.of.Peptides" %in% colnames(filtdf) && !is.null(minPeptides)) {
            filtdf <- filtdf %>%
                dplyr::mutate(Number.of.Peptides =
                                  as.numeric((.data$Number.of.Peptides < minPeptides) |
                                                 is.na(.data$Number.of.Peptides)))
        } else {
            filtdf$Number.of.Peptides <- NULL
        }
        if ("Score.Sequest.HT.Sequest.HT" %in% colnames(filtdf) && !is.null(minScore)) {
            filtdf <- filtdf %>%
                dplyr::mutate(Score.Sequest.HT.Sequest.HT =
                                  as.numeric((.data$Score.Sequest.HT.Sequest.HT < minScore) |
                                                 is.na(.data$Score.Sequest.HT.Sequest.HT)))
        } else {
            filtdf$Score.Sequest.HT.Sequest.HT <- NULL
        }
        if ("Master" %in% colnames(filtdf)) {
            filtdf <- filtdf %>%
                dplyr::mutate(Master =
                                  as.numeric((.data$Master != "IsMasterProtein") |
                                                 is.na(.data$Master)))
        }

        keep <- seq_len(nrow(sce))
        if ("Contaminant" %in% colnames(rowData(sce))) {
            keep <- intersect(keep, which(rowData(sce)$Contaminant == "False"))
        }
        if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce)) &&
            !is.null(minScore)) {
            keep <- intersect(keep, which(rowData(sce)$Score.Sequest.HT.Sequest.HT >= minScore))
        }
        if ("Number.of.Peptides" %in% colnames(rowData(sce)) && !is.null(minPeptides)) {
            keep <- intersect(keep, which(rowData(sce)$Number.of.Peptides >= minPeptides))
        }
        if (masterProteinsOnly) {
            if ("Master" %in% colnames(rowData(sce))) {
                keep <- intersect(keep, which(rowData(sce)$Master == "IsMasterProtein"))
            }
        } else {
            filtdf$Master <- NULL
        }
        exclude <- rowData(sce[setdiff(seq_len(nrow(sce)), keep), ])
        sce <- sce[keep, ]
    } else if (inputLevel == "PeptideGroups") {
        filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
            dplyr::select(dplyr::any_of(c("Contaminant", "Number.of.PSMs",
                                          "Delta.Score.by.Search.Engine.Sequest.HT"))) %>%
            dplyr::mutate(across(dplyr::any_of(c("Contaminant")),
                                 function(x) as.numeric(x == "True")))
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
        keep <- seq_len(nrow(sce))
        if ("Contaminant" %in% colnames(rowData(sce))) {
            keep <- intersect(keep, which(rowData(sce)$Contaminant == "False"))
        }
        if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce)) &&
            !is.null(minDeltaScore)) {
            keep <- intersect(keep, which(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= minDeltaScore))
        }
        if ("Number.of.PSMs" %in% colnames(rowData(sce)) && !is.null(minPSMs)) {
            keep <- intersect(keep, which(rowData(sce)$Number.of.PSMs >= minPSMs))
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
    if (plotUpset && any(rowSums(filtdf) > 0)) {
        print(ComplexUpset::upset(filtdf[rowSums(filtdf) > 0, , drop = FALSE],
                                  intersect = colnames(filtdf)))
    }

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }

    sce
}

#' Filter out features in FragPipe data
#'
#' Exclude features with 'Combined.Total.Peptides' below \code{minPeptides},
#' or identified as either 'Reverse' (Protein name starting with rev_) or
#' 'Potential.contaminant' (Protein name starting with contam_) by FragPipe.
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
#' @param exclFile Character scalar, the path to a text file where the
#'     features that are filtered out are written. If \code{NULL} (default),
#'     excluded features are not recorded.
#'
#' @return A filtered object of the same type as \code{sce}.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate across
#' @importFrom ComplexUpset upset
#' @importFrom rlang .data
#'
filterFragPipe <- function(sce, minPeptides, plotUpset = TRUE,
                           exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = minPeptides, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    ## Make sure that the columns used for filtering later are character vectors
    rowData(sce)$Potential.contaminant <- ifelse(grepl("^contam_", rowData(sce)$Protein), "+", "")
    rowData(sce)$Reverse <- ifelse(grepl("^rev_", rowData(sce)$Protein), "+", "")

    filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
        dplyr::select(dplyr::any_of(c("Reverse", "Potential.contaminant",
                                      "Combined.Total.Peptides"))) %>%
        dplyr::mutate(across(dplyr::any_of(c("Reverse", "Potential.contaminant")),
                             function(x) as.numeric(x == "+")))
    if ("Combined.Total.Peptides" %in% colnames(filtdf) && !is.null(minPeptides)) {
        filtdf <- filtdf %>%
            dplyr::mutate(
                Combined.Total.Peptides = as.numeric((.data$Combined.Total.Peptides < minPeptides) |
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
    if ("Combined.Total.Peptides" %in% colnames(rowData(sce)) && !is.null(minPeptides)) {
        keep <- intersect(keep, which(rowData(sce)$Combined.Total.Peptides >= minPeptides))
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
    if (plotUpset && any(rowSums(filtdf) > 0)) {
        print(ComplexUpset::upset(filtdf[rowSums(filtdf) > 0, , drop = FALSE],
                                  intersect = colnames(filtdf)))
    }

    if (!is.null(exclFile)) {
        write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                    row.names = FALSE, col.names = TRUE)
    }
    sce
}
