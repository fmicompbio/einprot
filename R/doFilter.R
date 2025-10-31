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

#' Filter features in a SummarizedExperiment object
#'
#' Filter the features (rows) in a SummarizedExperiment object based on a
#' user-defined combination of filters.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @param sce A \code{SummarizedExperiment} object (or a derivative).
#' @param filtersSE A named \code{list}, where each element is a function
#'     that takes a \code{SummarizedExperiment} object as input and returns
#'     a logical vector of the same length as the number of rows in \code{sce},
#'     and where \code{TRUE} implies that the row should be retained. Default
#'     sets of filtering functions for input data from different tools are
#'     provided (see \code{?defaultFiltersSE}).
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
#' # MaxQuant
#' sce <- importExperiment(inFile = system.file("extdata", "mq_example",
#'                                              "1356_proteinGroups.txt",
#'                                              package = "einprot"),
#'                         iColPattern = "^LFQ.intensity.")$sce
#'
#' dim(sce)
#' sce <- filterFeaturesSE(sce = sce, filtersSE = einprotMQfilters,
#'                         plotUpset = TRUE)
#' dim(sce)
#'
#' # ProteomeDiscoverer, Proteins
#' sce <- importExperiment(
#'     inFile = system.file("extdata", "pdtmt_example",
#'                          "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
#'                          package = "einprot"),
#'     iColPattern = "^Abundance.F.+.Sample.")$sce
#'
#' dim(sce)
#' sce <- filterFeaturesSE(sce = sce, filtersSE = einprotPDTMTProteinFilters,
#'                         plotUpset = TRUE)
#' dim(sce)
#'
#' # ProteomeDiscoverer, PeptideGroups
#' sce <- importExperiment(
#'     inFile = system.file("extdata", "pdtmt_example",
#'                          "Fig2_m23139_RTS_QC_varMods_PeptideGroups.txt",
#'                          package = "einprot"),
#'     iColPattern = "^Abundance.F.+.Sample.")$sce
#'
#' dim(sce)
#' sce <- filterFeaturesSE(sce = sce, filtersSE = einprotPDTMTPeptideGroupFilters)
#' dim(sce)
#'
#' # FragPipe
#' sce <- importExperiment(inFile = system.file("extdata", "fp_example",
#'                                              "combined_protein.tsv",
#'                                              package = "einprot"),
#'                         iColPattern = ".MaxLFQ.Intensity$")$sce
#'
#' dim(sce)
#' sce <- filterFeaturesSE(sce = sce, filtersSE = einprotFragPipeFilters,
#'                        plotUpset = TRUE)
#' dim(sce)
#'
filterFeaturesSE <- function(sce, filtersSE = list(), plotUpset = TRUE,
                             exclFile = NULL) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertVector(x = filtersSE, type = "list")
    if (length(filtersSE) > 0) {
        .assertVector(x = names(filtersSE), type = "character")
        for (f in filtersSE) {
            stopifnot(is(f, "function"))
            stopifnot(length(formals(f)) == 1)
        }
    }
    .assertScalar(x = plotUpset, type = "logical")
    .assertScalar(x = exclFile, type = "character", allowNULL = TRUE)

    if (length(filtersSE) > 0) {
        # apply filters and check that all vectors are logical and have
        # no NA values
        tokeep <- lapply(filtersSE, function(f) f(sce))
        for (filts in names(tokeep)) {
            stopifnot(is.logical(tokeep[[filts]]) &&
                          all(!is.na(tokeep[[filts]])))
        }

        filtdf <- do.call(bind_rows, lapply(tokeep, function(x) !x))
        exclude <- rowData(sce)[which(rowSums(filtdf) > 0), , drop = FALSE]
        sce <- sce[which(rowSums(filtdf) == 0), ]

        .makeFilterPlot(filtdf = filtdf, plotUpset = plotUpset)

        if (!is.null(exclFile)) {
            write.table(exclude, file = exclFile, quote = FALSE, sep = "\t",
                        row.names = FALSE, col.names = TRUE)
        }
        sce
    } else {
        sce
    }
}

#' Default filters (SE)
#'
#' einprot provides default filters for data imported from various tools.
#' To see the functions and thresholds used in each of these, print the
#' corresponding object below.
#'
#' @export
#' @name defaultFiltersSE
#' @rdname defaultFiltersSE
einprotMQFilters <- list(
    Score = function(sce) {
        if ("Score" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Score >= 10 & !is.na(rowData(sce)$Score)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Peptides = function(sce) {
        if ("Peptides" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Peptides >= 2 & !is.na(rowData(sce)$Peptides)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Reverse = function(sce) {
        if ("Reverse" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Reverse == "" | is.na(rowData(sce)$Reverse)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Contaminant = function(sce) {
        if ("Potential.contaminant" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Potential.contaminant == "" |
                is.na(rowData(sce)$Potential.contaminant)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    OnlySite = function(sce) {
        if ("Only.identified.by.site" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Only.identified.by.site == "" |
                is.na(rowData(sce)$Only.identified.by.site)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)

#' @export
#' @rdname defaultFiltersSE
einprotPDTMTProteinFilters <- list(
    Score = function(sce) {
        if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Score.Sequest.HT.Sequest.HT >= 2 &
                !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Peptides = function(sce) {
        if ("Number.of.Peptides" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Number.of.Peptides >= 2 & !is.na(rowData(sce)$Number.of.Peptides)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Contaminant = function(sce) {
        if ("Contaminant" %in% colnames(rowData(sce))) {
            rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant)
            kp <- rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)

#' @export
#' @rdname defaultFiltersSE
einprotPDTMTPeptideGroupFilters <- list(
    DeltaScore = function(sce) {
        if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.2 &
                !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    PSMs = function(sce) {
        if ("Number.of.PSMs" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Number.of.PSMs >= 2 & !is.na(rowData(sce)$Number.of.PSMs)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Contaminant = function(sce) {
        if ("Contaminant" %in% colnames(rowData(sce))) {
            rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant)
            kp <- rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)

#' @export
#' @rdname defaultFiltersSE
einprotFragPipeFilters <- list(
    Contaminant = function(sce) {
        if ("Protein" %in% colnames(rowData(sce))) {
            kp <- !grepl("^contam_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Reverse = function(sce) {
        if ("Protein" %in% colnames(rowData(sce))) {
            kp <- !grepl("^rev_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Peptides = function(sce) {
        if ("Combined.Total.Peptides" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$Combined.Total.Peptides >= 2 & !is.na(rowData(sce)$Combined.Total.Peptides)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)

#' @export
#' @rdname defaultFiltersSE
einprotSpectronautFilters <- list(
    Contaminant = function(sce) {
        if ("PG.ProteinGroups" %in% colnames(rowData(sce))) {
            kp <- !grepl("^contam_", rowData(sce)$PG.ProteinGroups) &
                !is.na(rowData(sce)$PG.ProteinGroups)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Reverse = function(sce) {
        if ("PG.ProteinGroups" %in% colnames(rowData(sce))) {
            kp <- !grepl("_Decoy$", rowData(sce)$PG.ProteinGroups) &
                !is.na(rowData(sce)$PG.ProteinGroups)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Peptides = function(sce) {
        if ("PG.NrOfStrippedSequencesIdentified.Experiment.wide" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$PG.NrOfStrippedSequencesIdentified.Experiment.wide >= 2 &
                !is.na(rowData(sce)$PG.NrOfStrippedSequencesIdentified.Experiment.wide)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    },
    Score = function(sce) {
        if ("PG.Cscore" %in% colnames(rowData(sce))) {
            kp <- rowData(sce)$PG.Cscore >= 2 &
                !is.na(rowData(sce)$PG.Cscore)
        } else {
            kp <- rep(TRUE, nrow(sce))
        }
        kp
    }
)
