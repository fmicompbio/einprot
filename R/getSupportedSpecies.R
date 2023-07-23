#' List supported species
#'
#' Get a list of species supported by \code{einprot}.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A \code{data.frame} with three columns (\code{taxId},
#'     \code{species} and \code{speciesCommon}) for each of the species
#'     supported by \code{einprot}.
#'
#' @examples
#' getSupportedSpecies()
#'
getSupportedSpecies <- function() {
    data.frame(
        taxId = c(10090, 9606, 6239, 7955, 7227, 4932, 284812, 28377,
                  9913, 9615, 9796, 9685, 9031, 9544, 13616, 9258,
                  9598, 10116, 9823, 8364),
        species = c("Mus musculus", "Homo sapiens", "Caenorhabditis elegans",
                    "Danio rerio", "Drosophila melanogaster",
                    "Saccharomyces cerevisiae",
                    "Schizosaccharomyces pombe 972h-",
                    "Anolis carolinensis", "Bos taurus",
                    "Canis lupus familiaris", "Equus caballus",
                    "Felis catus", "Gallus gallus",
                    "Macaca mulatta", "Monodelphis domestica",
                    "Ornithorhynchus anatinus", "Pan troglodytes",
                    "Rattus norvegicus", "Sus scrofa",
                    "Xenopus tropicalis"),
        speciesCommon = c("mouse", "human", "roundworm", "zebrafish",
                          "fruitfly", "baker's yeast", "fission yeast",
                          "green anole", "bovine", "dog", "horse", "cat",
                          "chicken", "rhesus macaque", "opossum",
                          "platypus", "chimpanzee", "Norway rat",
                          "pig", "tropical clawed frog")
    )
}

#' Get species info
#'
#' Get the scientific species name, the common name and the taxonomic ID for
#' any of the species supported by \code{einprot} (see
#' \code{getSupportedSpecies()} for a list of supported species).
#'
#' @param species Character or numeric scalar, representing either a
#'     scientific species ID, a common species name or a taxonomic ID for
#'     one of the species supported by einprot.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A list with three elements: the scientific species name, the
#'     common name and the taxonomic ID.
#'
#' @examples
#' getSpeciesInfo("mouse")
#' getSpeciesInfo(6239)
#' getSpeciesInfo("Homo sapiens")
#' ## unsupported species
#' getSpeciesInfo("E.coli")
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
        tax_id <- taxTable$taxId[match(species_id, taxTable$species)]
    } else if (tolower(species) %in% tolower(taxTable$species)) {
        species_id <- taxTable$species[match(tolower(species),
                                             tolower(taxTable$species))]
        species_common <-
            taxTable$speciesCommon[match(tolower(species),
                                         tolower(taxTable$species))]
        tax_id <- taxTable$taxId[match(species_id, taxTable$species)]
    } else if (species %in% taxTable$taxId) {
        species_id <- taxTable$species[match(species,
                                             taxTable$taxId)]
        species_common <- taxTable$speciesCommon[match(species,
                                                       taxTable$taxId)]
        tax_id <- taxTable$taxId[match(species_id, taxTable$species)]
    } else {
        warning("Unknown species ", species)
        species_id <- species
        species_common <- ""
        tax_id <- NA_real_
    }

    list(species = species_id,
         speciesCommon = species_common,
         taxId = tax_id)
}
