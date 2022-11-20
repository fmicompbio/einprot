.getUniProtToGeneSymbolMappingFile <- function(spi) {
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

#' Get mapping from UniProt IDs to gene symbol
#'
#' Generate a data.frame with the mapping between UniProtIDs and gene symbols.
#' The mapping is obtained from the UniProt ID mapping files (downloaded from
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/).
#'
#' @param species Either a taxonomy ID, a species ID or a common species name.
#'     See \code{getSupportedSpecies()} for valid values.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A data.frame with two columns, corresponding to the UniProtID and
#'     the corresponding gene symbol. Note that both columns can have
#'     duplicated values, if there is not a one-to-one mapping between
#'     UniProtIDs and gene symbols.
#'
#' @importFrom readr read_tsv
#' @importFrom tidyr pivot_wider unnest
#' @importFrom dplyr filter rename
#'
getUniProtToGeneSymbolMapping <- function(species) {
    spi <- getSpeciesInfo(species)

    f <- .getUniProtToGeneSymbolMappingFile(spi)

    readr::read_tsv(gzcon(
        url(paste0("https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/", f))),
        col_names = FALSE, col_types = "ccc") %>%
        dplyr::filter(.data$X2 %in% c("Gene_Name")) %>%
        tidyr::pivot_wider(id_cols = "X1", names_from = "X2",
                           values_from = "X3", values_fn = list) %>%
        tidyr::unnest(cols = c("Gene_Name")) %>%
        dplyr::rename(UniProtID = "X1", GeneSymbol = "Gene_Name") %>%
        as.data.frame()
}
