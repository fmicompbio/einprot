#' Get a list of species supported by einprot
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A `data.frame` with three columns (`taxId`,
#'     `species` and `speciesCommon`) for each of the species
#'     supported by `einprot`.
#'
#' @examples
#' getSupportedSpecies()
#'
getSupportedSpecies <- function() {
    data.frame(
        taxId = c(10090, 9606, 6239, 7955, 7227, 4932, 284812),
        species = c("Mus musculus", "Homo sapiens", "Caenorhabditis elegans",
                    "Danio rerio", "Drosophila melanogaster",
                    "Saccharomyces cerevisiae",
                    "Schizosaccharomyces pombe 972h-"),
        speciesCommon = c("mouse", "human", "roundworm", "zebrafish",
                          "fruitfly", "baker's yeast", "fission yeast")
    )
}

#' Get species info
#'
#' Get the scientific species name, the common name and the taxonomic ID for
#' any of the species supported by `einprot` (see \code{getSupportedSpecies()}
#' for a list of supported species).
#'
#' @param species Character or numeric scalar, representing either a
#'     scientific species ID, a common species name or a taxonomic ID for
#'     one of the species supported by einprot.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with three elements: the scientific species name, the
#'     common name and the taxonomic ID.
#'
#' @examples
#' getSpeciesInfo("mouse")
#' getSpeciesInfo(6239)
#' getSpeciesInfo("Homo sapiens")
#'
getSpeciesInfo <- function(species) {
    stopifnot(length(species) == 1)
    stopifnot(is.character(species) || is.numeric(species))

    taxTable <- getSupportedSpecies()
    if (tolower(species) %in% tolower(taxTable$speciesCommon)) {
        species_id <- taxTable$species[match(tolower(species),
                                             tolower(taxTable$speciesCommon))]
        species_common <-
            taxTable$speciesCommon[match(tolower(species),
                                         tolower(taxTable$speciesCommon))]
    } else if (tolower(species) %in% tolower(taxTable$species)) {
        species_id <- taxTable$species[match(tolower(species),
                                             tolower(taxTable$species))]
        species_common <-
            taxTable$speciesCommon[match(tolower(species),
                                         tolower(taxTable$species))]
    } else if (species %in% taxTable$taxId) {
        species_id <- taxTable$species[match(species,
                                             taxTable$taxId)]
        species_common <- taxTable$speciesCommon[match(species,
                                                       taxTable$taxId)]
    } else {
        stop("Unknown species ", species)
    }
    tax_id <- taxTable$taxId[match(species_id, taxTable$species)]

    list(species = species_id,
         speciesCommon = species_common,
         taxId = tax_id)
}
