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
#'
#' @return A filtered object of the same type as \code{sce}.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate across
#' @importFrom ComplexUpset upset
#' @importFrom rlang .data
#'
filterMaxQuant <- function(sce, minScore, minPeptides, plotUpset = TRUE) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = minScore, type = "numeric")
    .assertScalar(x = minPeptides, type = "numeric")
    .assertScalar(x = plotUpset, type = "logical")

    ## Make sure that the columns used for filtering are character vectors
    rowData(sce)$Reverse[is.na(rowData(sce)$Reverse)] <- ""
    rowData(sce)$Potential.contaminant[is.na(rowData(sce)$Potential.contaminant)] <- ""
    rowData(sce)$Only.identified.by.site[is.na(rowData(sce)$Only.identified.by.site)] <- ""

    filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
        dplyr::select("Reverse", "Potential.contaminant",
                      "Only.identified.by.site", "Score", "Peptides") %>%
        dplyr::mutate(across(c("Reverse", "Potential.contaminant",
                               "Only.identified.by.site"),
                             function(x) as.numeric(x == "+"))) %>%
        dplyr::mutate(Score = as.numeric((.data$Score < minScore) |
                                             is.na(.data$Score)),
                      Peptides = as.numeric((.data$Peptides < minPeptides) |
                                                is.na(.data$Peptides)))

    sce <- sce[which(rowData(sce)$Reverse == "" &
                         rowData(sce)$Potential.contaminant == "" &
                         rowData(sce)$Only.identified.by.site == "" &
                         rowData(sce)$Score >= minScore &
                         rowData(sce)$Peptides >= minPeptides), ]

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
#' @param plotUpset Logical scalar, whether to generate an UpSet plot
#'     detailing the reasons for features being filtered out. Only
#'     generated if any feature is in fact filtered out.
#'
#' @return A filtered object of the same type as \code{sce}.
#'
#' @importFrom SummarizedExperiment rowData
#' @importFrom dplyr select mutate
#' @importFrom ComplexUpset upset
#'
filterPDTMT <- function(sce, inputLevel, minScore = 0, minPeptides = 0,
                        minDeltaScore = 0, minPSMs = 0, plotUpset = TRUE) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = inputLevel, type = "character",
                  validValues = c("Proteins", "PeptideGroups"))
    .assertScalar(x = minScore, type = "numeric")
    .assertScalar(x = minPeptides, type = "numeric")
    .assertScalar(x = minDeltaScore, type = "numeric")
    .assertScalar(x = minPSMs, type = "numeric")
    .assertScalar(x = plotUpset, type = "logical")

    ## Make sure that the columns used for filtering are character vectors
    rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant)

    if (inputLevel == "Proteins") {
        filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
            dplyr::select("Contaminant", "Number.of.Peptides",
                          "Score.Sequest.HT.Sequest.HT") %>%
            dplyr::mutate(across(c("Contaminant"),
                                 function(x) as.numeric(x == "True"))) %>%
            dplyr::mutate(Number.of.Peptides =
                              as.numeric((.data$Number.of.Peptides < minPeptides) |
                                             is.na(.data$Number.of.Peptides)),
                          Score.Sequest.HT.Sequest.HT =
                              as.numeric((.data$Score.Sequest.HT.Sequest.HT < minScore) |
                                             is.na(.data$Score.Sequest.HT.Sequest.HT)))

        sce <- sce[which(rowData(sce)$Contaminant == "False" &
                             rowData(sce)$Score.Sequest.HT.Sequest.HT >= minScore &
                             rowData(sce)$Number.of.Peptides >= minPeptides), ]
    } else if (inputLevel == "PeptideGroups") {
        filtdf <- as.data.frame(SummarizedExperiment::rowData(sce)) %>%
            dplyr::select("Contaminant", "Number.of.PSMs",
                          "Delta.Score.by.Search.Engine.Sequest.HT") %>%
            dplyr::mutate(across(c("Contaminant"),
                                 function(x) as.numeric(x == "True"))) %>%
            dplyr::mutate(Number.of.PSMs =
                              as.numeric((.data$Number.of.PSMs < minPSMs) |
                                             is.na(.data$Number.of.PSMs)),
                          Delta.Score.by.Search.Engine.Sequest.HT =
                              as.numeric((.data$Delta.Score.by.Search.Engine.Sequest.HT <
                                              minDeltaScore) |
                                             is.na(.data$Delta.Score.by.Search.Engine.Sequest.HT)))
        sce <- sce[which(rowData(sce)$Contaminant == "False" &
                             rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >=
                             minDeltaScore &
                             rowData(sce)$Number.of.PSMs >= minPSMs), ]
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

    sce
}

