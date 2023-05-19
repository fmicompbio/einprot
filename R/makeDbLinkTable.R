#' Format table columns
#'
#' @export
#' @author Charlotte Soneson
#'
#' @param tbl A \code{data.frame} or similar object.
#' @param columns Character vector giving the column names of \code{tbl} that
#'     should be formatted.
#' @param signifDigits Numeric scalar, the number of significant digits to
#'     round numeric columns to. Set to \code{NULL} to skip rounding.
#' @param maxLevels Numeric scalar. If character columns have at most this
#'     number of unique values, they will be encoded as factors.
#'
formatTableColumns <- function(tbl, columns, signifDigits, maxLevels = 10) {
    .assertVector(x = columns, type = "character")
    .assertScalar(x = signifDigits, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = maxLevels, type = "numeric")

    for (nm in intersect(colnames(tbl), columns)) {
        if (is.numeric(tbl[[nm]])) {
            if (all(tbl[[nm]] == round(tbl[[nm]]), na.rm = TRUE)) {
                ## Integers - represent as such
                tbl[[nm]] <- as.integer(tbl[[nm]])
            } else {
                ## Not integers - possibly round
                if (!is.null(signifDigits)) {
                    tbl[[nm]] <- signif(tbl[[nm]], digits = signifDigits)
                }
            }
        } else if (is.logical(tbl[[nm]])) {
            ## Leave logical columns as they are
        } else {
            ## Character - encode as factor if less than or equal to 10 levels
            if (length(unique(tbl[[nm]])) <= maxLevels) {
                tbl[[nm]] <- factor(tbl[[nm]])
            }
        }
    }
    tbl
}

#' @noRd
#' @keywords internal
.makeLinkFromId <- function(id, linktype = "UniProt",
                            removeSuffix = TRUE) {
    .assertScalar(x = id, type = "character")
    .assertScalar(x = linktype, type = "character",
                  validValues = c("UniProt", "AlphaFold", "PomBase",
                                  "WormBase"))

    if (!is.na(id) && id != "") {
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
            tidyr::pivot_wider(id_cols = "X1", names_from = "X2",
                               values_from = "X3", values_fn = list) %>%
            tidyr::unnest(cols = c("UniProtKB-ID", "WormBase")) %>%
            dplyr::mutate(`UniProtKB-ID` = gsub("_CAEEL$", "",
                                                .data$`UniProtKB-ID`)) %>%
            dplyr::rename(UniProtID = "X1", WormBaseID = "WormBase",
                          UniProtKB.ID = "UniProtKB-ID") %>%
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
#'     generated using \code{getConvTable(type = "PomBase")}.
#' @param convTableWormBase Conversion table between UniProt IDs and
#'     WormBase IDs. Only used if \code{speciesCommon} is \code{"roundworm"}.
#'     A suitably formatted conversion table can be
#'     generated using \code{getConvTable(type = "WormBase")}.
#' @param removeSuffix Logical scalar indicating whether suffixes of the
#'     form `-[0-9]+` should be removed from the protein ID before
#'     generating the URL. Currently only influencing the AlphaFold URL.
#' @param signifDigits Numeric scalar giving the number of significant digits
#'     to round numeric columns to. If \code{NULL}, no rounding will be
#'     performed.
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
                            removeSuffix = TRUE, signifDigits = 4) {

    ## -------------------------------------------------------------------------
    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = df, type = "data.frame")
    .assertScalar(x = idCol, type = "character", validValues = colnames(df))
    .assertScalar(x = speciesCommon, type = "character",
                  validValues = getSupportedSpecies()$speciesCommon)
    .assertScalar(x = addSpeciesSpecificColumns, type = "logical")
    .assertVector(x = convTablePomBase, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = convTableWormBase, type = "data.frame",
                  allowNULL = TRUE)
    .assertScalar(x = removeSuffix, type = "logical")
    .assertScalar(x = signifDigits, type = "numeric", allowNULL = TRUE)

    ## -------------------------------------------------------------------------
    ## Create UniProt and AlphaFold columns
    ## -------------------------------------------------------------------------
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

    ## -------------------------------------------------------------------------
    ## Round numeric columns, encode integers as such
    ## Represent non-numeric columns with <= 10 unique values as factors
    ## -------------------------------------------------------------------------
    linkTable <- formatTableColumns(
        tbl = linkTable,
        columns = setdiff(colnames(linkTable), c("UniProt", "AlphaFold",
                                                 "WormBase", "PomBase")),
        signifDigits = signifDigits, maxLevels = 10)

    linkTable
}
