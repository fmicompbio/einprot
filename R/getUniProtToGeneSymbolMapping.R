.getUniProtToIDMappingFile <- function(spi) {
    .assertVector(x = spi, type = "list")
    stopifnot("species" %in% names(spi))

    ## Get file name
    if (spi$species == "Mus musculus") {
        f <- "MOUSE_10090_idmapping.dat.gz"
    } else if (spi$species == "Homo sapiens") {
        f <- "HUMAN_9606_idmapping.dat.gz"
    } else if (spi$species == "Caenorhabditis elegans") {
        f <- "CAEEL_6239_idmapping.dat.gz"
    } else if (spi$species == "Danio rerio") {
        f <- "DANRE_7955_idmapping.dat.gz"
    } else if (spi$species == "Drosophila melanogaster") {
        f <- "DROME_7227_idmapping.dat.gz"
    } else if (spi$species == "Saccharomyces cerevisiae") {
        ## taxId 4932 is not available
        f <- "YEAST_559292_idmapping.dat.gz"
    } else if (spi$species == "Schizosaccharomyces pombe 972h-") {
        f <- "SCHPO_284812_idmapping.dat.gz"
    } else {
        stop("Unsupported species")
    }

    f
}

#' Get mapping from UniProt IDs to another ID type
#'
#' Generate a \code{data.frame} with the mapping between UniProtIDs and
#' another ID type.
#' The mapping is obtained from the UniProt ID mapping files (downloaded from
#' \url{https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/}).
#'
#' @param species Either a taxonomy ID, a species ID or a common species name.
#'     See \code{getSupportedSpecies()} for valid values.
#' @param targetId Character scalar giving the target ID type. Should be one
#'     of the ID types included in the second column in the UniProt ID mapping
#'     files. Examples include "Gene_Name", "Ensembl", "STRING", "RefSeq_NT"
#'     (depending on the species).
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A data.frame with two columns, corresponding to the UniProtID and
#'     the corresponding other ID. Note that both columns can have
#'     duplicated values, if there is not a one-to-one mapping between
#'     UniProtIDs and the other ID type.
#'
#' @examples
#' df <- getUniProtToIDMapping("fission yeast", targetId = "Gene_Name")
#' head(df)
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_wider unnest
#' @importFrom dplyr filter rename all_of
#'
getUniProtToIDMapping <- function(species, targetId = "Gene_Name") {
    .assertScalar(x = targetId, type = "character")

    spi <- getSpeciesInfo(species)

    f <- .getUniProtToIDMappingFile(spi)

    readr::read_tsv(gzcon(
        url(paste0("https://ftp.uniprot.org/pub/databases/uniprot/",
                   "current_release/knowledgebase/idmapping/by_organism/", f))),
        col_names = FALSE, col_types = "ccc") %>%
        dplyr::filter(.data$X2 %in% targetId) %>%
        tidyr::pivot_wider(id_cols = "X1", names_from = "X2",
                           values_from = "X3", values_fn = list) %>%
        tidyr::unnest(cols = dplyr::all_of(targetId)) %>%
        dplyr::rename(UniProtID = "X1") %>%
        as.data.frame()
}
