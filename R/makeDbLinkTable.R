#' @noRd
#' @keywords internal
.makeLinkFromId <- function(id, linktype = "UniProt",
                            removeSuffix = TRUE) {
    .assertScalar(x = id, type = "character")
    .assertScalar(x = linktype, type = "character",
                  validValues = c("UniProt", "AlphaFold", "PomBase",
                                  "WormBase"))

    if (id != "") {
        if (removeSuffix) {
            id <- sub("-[0-9]+$", "", id)
        }
        if (linktype == "UniProt") {
            sprintf('<a href="%s" target="_blank"> %s</a>',
                    paste0("https://www.uniprot.org/uniprot/", id), id)
        } else if (linktype == "AlphaFold") {
            sprintf('<a href="%s" target="_blank"> %s</a>',
                    paste0("https://alphafold.ebi.ac.uk/entry/", id), id)
        } else if (linktype == "PomBase") {
            sprintf('<a href="%s" target="_blank"> %s</a>',
                    paste0("https://www.pombase.org/gene/", id), id)
        } else if (linktype == "WormBase") {
            sprintf('<a href="%s" target="_blank"> %s</a>',
                    paste0("https://wormbase.org/species/c_elegans/gene/", id), id)
        } else {
            ""
        }
    } else {
        ""
    }
}

#' Download and process conversion tables from UniProt IDs to PomBase/WormBase.
#'
#' The PomBase-UniProt ID conversion table is downloaded from
#' https://www.pombase.org/data/names_and_identifiers/PomBase2UniProt.tsv.
#' The WormBase-UniProt ID conversion table is downloaded from
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/CAEEL_6239_idmapping.dat.gz.
#'
#' @param type Character scalar, either "PomBase" or "WormBase"
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A \code{data.frame} with ID conversion information.
#'
#' @examples
#' df <- getConvTable(type = "PomBase")
#' df <- getConvTable(type = "WormBase")
#'
#' @importFrom utils read.delim
#' @importFrom readr read_tsv
#' @importFrom stats setNames
#' @importFrom dplyr filter mutate
#' @importFrom tidyr pivot_wider unnest
#' @importFrom rlang .data
#'
getConvTable <- function(type) {
    .assertScalar(x = type, type = "character",
                  validValues = c("PomBase", "WormBase"))

    if (type == "PomBase") {
        utils::read.delim(
            url("https://www.pombase.org/data/names_and_identifiers/PomBase2UniProt.tsv"),
            header = FALSE
        ) %>%
            stats::setNames(c("PomBaseID", "UniProtID"))
    } else if (type == "WormBase") {
        readr::read_tsv(gzcon(
            url("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/CAEEL_6239_idmapping.dat.gz")),
            col_names = FALSE, col_types = "ccc") %>%
            dplyr::filter(.data$X2 %in% c("UniProtKB-ID", "WormBase")) %>%
            tidyr::pivot_wider(id_cols = .data$X1, names_from = .data$X2,
                               values_from = .data$X3, values_fn = list) %>%
            tidyr::unnest(cols = c("UniProtKB-ID", "WormBase")) %>%
            dplyr::mutate(`UniProtKB-ID` = gsub("_CAEEL$", "",
                                                .data$`UniProtKB-ID`)) %>%
            dplyr::rename(UniProtID = .data$X1, WormBaseID = .data$WormBase,
                          UniProtKB.ID = .data$`UniProtKB-ID`) %>%
            as.data.frame()
    }
}

#' Make a table with database links
#'
#' @param df \code{data.frame} with protein identifiers (UniProt IDs).
#' @param idCol Character scalar giving the column name of the column in
#'     \code{df} containing protein identifiers. Each element of this
#'     column should be a string, which can represent multiple protein IDs
#'     separated by semicolons. See examples for an illustration.
#' @param speciesCommon Character scalar, providing the common species
#'     name (e.g., mouse, roundworm, fission yeast).
#' @param addSpeciesSpecificColumns Logical scalar, indicating whether to
#'     add species-specific columns (whenever applicable).
#' @param convTablePomBase Conversion table between UniProt IDs and
#'     PomBase IDs. Only used if \code{speciesCommon} is
#'     \code{"fission yeast"}. A suitably formatted conversion table can be
#'     downloaded using \code{getConvTable(type = "PomBase")}.
#' @param convTableWormBase Conversion table between UniProt IDs and
#'     WormBase IDs. Only used if \code{speciesCommon} is \code{"roundworm"}.
#'     A suitably formatted conversion table can be
#'     downloaded using \code{getConvTable(type = "WormBase")}.
#' @param removeSuffix Logical scalar indicating whether suffixes of the
#'     form `-[0-9]+` should be removed from the protein ID before
#'     generating the URL. Currently only influencing the AlphaFold URL.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A \code{data.frame} with database links for the proteins.
#'
#' @examples
#' pbconv <- getConvTable(type = "PomBase")
#' makeDbLinkTable(data.frame(id = c("B5BP45", "O13282")), idCol = "id",
#'                            speciesCommon = "fission yeast",
#'                            convTablePomBase = pbconv)
#' wbconv <- getConvTable(type = "WormBase")
#' makeDbLinkTable(data.frame(gid = c("eps-8", "epi-1"),
#'                            pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
#'                 idCol = "pid", speciesCommon = "roundworm",
#'                 convTableWormBase = wbconv)
#'
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
makeDbLinkTable <- function(df, idCol, speciesCommon,
                            addSpeciesSpecificColumns = TRUE,
                            convTablePomBase = NULL,
                            convTableWormBase = NULL,
                            removeSuffix = TRUE) {

    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .assertVector(x = df, type = "data.frame")
    .assertScalar(x = idCol, type = "character", validValues = colnames(df))
    .assertScalar(x = speciesCommon, type = "character",
                  validValues = getSupportedSpecies()$speciesCommon)
    .assertScalar(x = addSpeciesSpecificColumns, type = "logical")
    .assertVector(x = convTablePomBase, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = convTableWormBase, type = "data.frame",
                  allowNULL = TRUE)

    ## --------------------------------------------------------------------- ##
    ## Create UniProt and AlphaFold columns
    ## --------------------------------------------------------------------- ##
    linkTable <- df %>%
        dplyr::mutate(UniProt = vapply(.data[[idCol]], function(mpds) {
            paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                .makeLinkFromId(mpd, linktype = "UniProt",
                                removeSuffix = FALSE)
            }, ""), collapse = ";")
        }, "NA")) %>%
        dplyr::mutate(AlphaFold = vapply(.data[[idCol]], function(mpds) {
            paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                .makeLinkFromId(mpd, linktype = "AlphaFold",
                                removeSuffix = removeSuffix)
            }, ""), collapse = ";")
        }, "NA"))

    if (addSpeciesSpecificColumns &&
        speciesCommon == "fission yeast" && !is.null(convTablePomBase)) {
        ## Add PomBase links
        linkTable <- linkTable %>%
            dplyr::mutate(PomBase = vapply(.data[[idCol]], function(mpds) {
                paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                    pbid <- convTablePomBase$PomBaseID[convTablePomBase$UniProtID == mpd]
                    if (length(pbid) != 0) {
                        paste(vapply(pbid, function(pb) {
                            .makeLinkFromId(pb, linktype = "PomBase",
                                            removeSuffix = FALSE)
                        }, ""), collapse = ";")
                    } else {
                        ""
                    }
                }, ""), collapse = ";")
            }, "NA"))
    }

    if (addSpeciesSpecificColumns &&
        speciesCommon == "roundworm" && !is.null(convTableWormBase)) {
        ## Add WormBase links
        linkTable <- linkTable %>%
            dplyr::mutate(WormBase = vapply(.data[[idCol]], function(mpds) {
                wbids <- unlist(lapply(strsplit(mpds, ";")[[1]], function(mpd) {
                    convTableWormBase$WormBaseID[convTableWormBase$UniProtKB.ID == mpd |
                                                     convTableWormBase$UniProtID == mpd]
                }))
                if (length(wbids) != 0 && length(setdiff(wbids, "")) != 0) {
                    wbids <- setdiff(wbids, "")
                    paste(vapply(wbids, function(wb) {
                        if (!is.na(wb)) {
                            .makeLinkFromId(wb, linktype = "WormBase",
                                            removeSuffix = FALSE)
                        } else {
                            ""
                        }
                    }, ""), collapse = ";")
                } else {
                    ""
                }
            }, "NA"))
    }

    linkTable
}
