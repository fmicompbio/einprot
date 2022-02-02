.makeLinkFromId <- function(id, linktype = "UniProt") {
    .assertScalar(x = id, type = "character")
    .assertScalar(x = linktype, type = "character",
                  validValues = c("UniProt", "AlphaFold", "PomBase",
                                  "WormBase"))

    if (id != "") {
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

#' Download and process conversion tables from UniProt IDs to Pombase/Wormbase
#'
#' @param type Character scalar, either "Pombase" or "Wormbase"
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A \code{data.frame} with conversion information.
#'
#' @examples
#' df <- getConvTable(type = "Pombase")
#' df <- getConvTable(type = "Wormbase")
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
                  validValues = c("Pombase", "Wormbase"))

    if (type == "Pombase") {
        utils::read.delim(
            url("https://www.pombase.org/data/names_and_identifiers/PomBase2UniProt.tsv"),
            header = FALSE
        ) %>%
            stats::setNames(c("PomBaseID", "UniProtID"))
    } else if (type == "Wormbase") {
        readr::read_tsv(gzcon(
            url("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/CAEEL_6239_idmapping.dat.gz")),
            col_names = FALSE, col_types = "ccc") %>%
            dplyr::filter(.data$X2 %in% c("UniProtKB-ID", "WormBase")) %>%
            tidyr::pivot_wider(id_cols = .data$X1, names_from = .data$X2,
                               values_from = .data$X3, values_fn = list) %>%
            tidyr::unnest(cols = c("UniProtKB-ID", "WormBase")) %>%
            dplyr::mutate(`UniProtKB-ID` = gsub("_CAEEL$", "", .data$`UniProtKB-ID`)) %>%
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
#'     separated by semicolon. See examples for an illustration.
#' @param speciesCommon Character scalar, providing the common species
#'     name (e.g., mouse, roundworm, fission yeast).
#' @param addSpeciesSpecificColumns Logical scalar, indicating whether to
#'     add species-specific columns (whenever applicable).
#' @param convTablePombase Conversion table between UniProt IDs and
#'     Pombase IDs. Only used if \code{speciesCommon} is
#'     \code{"fission yeast"}.
#' @param convTableWormbase Conversion table between UniProt IDs and
#'     Wormbase IDs. Only used if \code{speciesCommon} is \code{"roundworm"}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A \code{data.frame} with database links for the proteins.
#'
#' @examples
#' pbconv <- getConvTable(type = "Pombase")
#' makeDbLinkTable(data.frame(id = c("B5BP45", "O13282")), idCol = "id",
#'                            speciesCommon = "fission yeast",
#'                            convTablePombase = pbconv)
#' wbconv <- getConvTable(type = "Wormbase")
#' makeDbLinkTable(data.frame(gid = c("eps-8", "epi-1"),
#'                            pid = c("Q7YTG1;O18250", "C1P641;C1P640")),
#'                 idCol = "pid", speciesCommon = "roundworm",
#'                 convTableWormbase = wbconv)
#'
#' @importFrom dplyr mutate
#' @importFrom rlang .data
#'
makeDbLinkTable <- function(df, idCol, speciesCommon,
                            addSpeciesSpecificColumns = TRUE,
                            convTablePombase = NULL,
                            convTableWormbase = NULL) {

    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .assertVector(x = df, type = "data.frame")
    .assertScalar(x = idCol, type = "character", validValues = colnames(df))
    .assertScalar(x = speciesCommon, type = "character",
                  validValues = getSupportedSpecies()$speciesCommon)
    .assertScalar(x = addSpeciesSpecificColumns, type = "logical")
    .assertVector(x = convTablePombase, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = convTableWormbase, type = "data.frame",
                  allowNULL = TRUE)

    ## --------------------------------------------------------------------- ##
    ## Create UniProt and AlphaFold columns
    ## --------------------------------------------------------------------- ##
    linkTable <- df %>%
        # tibble::rownames_to_column("pid") %>%
        # dplyr::select(.data$pid, .data[[idCol]]) %>%
        dplyr::mutate(UniProt = vapply(.data[[idCol]], function(mpds) {
            paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                .makeLinkFromId(mpd, linktype = "UniProt")
            }, ""), collapse = ";")
        }, "NA")) %>%
        dplyr::mutate(AlphaFold = vapply(.data[[idCol]], function(mpds) {
            paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                .makeLinkFromId(mpd, linktype = "AlphaFold")
            }, ""), collapse = ";")
        }, "NA"))

    if (addSpeciesSpecificColumns &&
        speciesCommon == "fission yeast" && !is.null(convTablePombase)) {
        ## Add PomBase links
        linkTable <- linkTable %>%
            dplyr::mutate(PomBase = vapply(.data[[idCol]], function(mpds) {
                paste(vapply(strsplit(mpds, ";")[[1]], function(mpd) {
                    pbid <- convTablePombase$PomBaseID[convTablePombase$UniProtID == mpd]
                    if (length(pbid) != 0) {
                        paste(vapply(pbid, function(pb) {
                            .makeLinkFromId(pb, linktype = "PomBase")
                        }, ""), collapse = ";")
                    } else {
                        ""
                    }
                }, ""), collapse = ";")
            }, "NA"))
    }

    if (addSpeciesSpecificColumns &&
        speciesCommon == "roundworm" && !is.null(convTableWormbase)) {
        ## Add WormBase links
        linkTable <- linkTable %>%
            dplyr::mutate(WormBase = vapply(.data[[idCol]], function(mpds) {
                wbids <- unlist(lapply(strsplit(mpds, ";")[[1]], function(mpd) {
                    convTableWormbase$WormBaseID[convTableWormbase$UniProtKB.ID == mpd |
                                                     convTableWormbase$UniProtID == mpd]
                }))
                if (length(wbids) != 0 && length(setdiff(wbids, "")) != 0) {
                    wbids <- setdiff(wbids, "")
                    paste(vapply(wbids, function(wb) {
                        if (!is.na(wb)) {
                            .makeLinkFromId(wb, linktype = "WormBase")
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
