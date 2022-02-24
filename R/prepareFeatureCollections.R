## Help function to replace every nth ";" by "; "
.createPattern <- function(n) {
    sprintf("(%s[^;]+);", strrep("[^;]+;", n - 1))
}

## Help function to replace IDs in feature collections by other IDs
.replaceIdsInList <- function(chl, dfConv, currentIdCol, newIdCol, pat) {
    ## Replace each complex with the existing row names in qft
    chl <- S4Vectors::endoapply(chl, function(x) {
        dfConv[[newIdCol]][dfConv[[currentIdCol]] %in% x]
    })
    S4Vectors::mcols(chl)$sharedGenes <- vapply(
        chl, function(w) gsub(pat, "\\1; ", paste(w, collapse = ";")), "")
    S4Vectors::mcols(chl)$nSharedGenes <- lengths(chl)
    chl
}

#' Prepare feature collections for testing with camera
#'
#' Prepare feature collections for testing with camera. The function maps
#' the feature IDs in the collections (complexes or GO terms) to the
#' values in the specified \code{idCol} column of \code{rowData(qft[[1]])},
#' and subsequently replaces them with the corresponding row names of the
#' \code{QFeatures} object. Feature sets with too few features (after the
#' matching) are removed.
#'
#' @param qft A \code{QFeatures} object.
#' @param idCol Character scalar, indicating which column in
#'     \code{rowData(qft[[1]])} that contains IDs matching those in the
#'     feature collections (gene symbols).
#' @param includeFeatureCollections Character vector indicating the types
#'     of feature collections to prepare. Should be a subset of
#'     \code{c("complexes", "GO")} or \code{NULL}.
#' @param complexDbPath Character scalar providing the path to the database
#'     of complexes, generated using \code{makeComplexDB()} and serialized
#'     to a .rds file.
#' @param speciesInfo List with at least two entries (\code{species} and
#'     \code{speciesCommon}), providing the species information. Typically
#'     generated using \code{getSpeciesInfo()}.
#' @param complexSpecies Character scalar, either \code{"all"} or
#'     \code{"current"}, indicating whether all complexes should be
#'     tested, or only those defined for the current species.
#' @param customComplexes Named list, for providing any custom complexes
#'     that are not already included in the database provided via
#'     \code{complexDbPath}.
#' @param minSizeToKeep Numeric scalar, indicating the minimum size of a
#'     feature collection to be retained.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A list of \code{CharacterList}s (one for each feature collection).
#'
#' @importFrom S4Vectors mcols DataFrame endoapply
#' @importFrom IRanges CharacterList
#' @importFrom methods as
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr select %>%
#' @importFrom rlang .data
#' @importFrom SummarizedExperiment rowData
#' @importFrom tidyr separate_rows
#'
prepareFeatureCollections <- function(qft, idCol, includeFeatureCollections,
                                      complexDbPath, speciesInfo,
                                      complexSpecies, customComplexes = list(),
                                      minSizeToKeep = 2) {
    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .assertVector(x = qft, type = "QFeatures")
    .assertScalar(x = idCol, type = "character",
                  validValues = colnames(SummarizedExperiment::rowData(qft[[1]])))
    .assertVector(x = includeFeatureCollections, type = "character",
                  validValues = c("complexes", "GO"), allowNULL = TRUE)
    .assertScalar(x = complexDbPath, type = "character", allowNULL = TRUE)
    if (!is.null(complexDbPath)) {
        stopifnot(file.exists(complexDbPath))
    }
    .assertVector(x = speciesInfo, type = "list")
    .assertVector(x = names(speciesInfo), type = "character",
                  validValues = c("species", "speciesCommon", "taxId"))
    .assertScalar(x = complexSpecies, type = "character",
                  validValues = c("current", "all"))
    .assertVector(x = customComplexes, type = "list")
    if (length(customComplexes) > 0) {
        .assertVector(x = names(customComplexes), type = "character")
    }
    .assertScalar(x = minSizeToKeep, type = "numeric", rngIncl = c(0, Inf))

    ## --------------------------------------------------------------------- ##
    ## Initialization
    ## --------------------------------------------------------------------- ##
    pat <- .createPattern(10)
    featureCollections <- list()

    ## --------------------------------------------------------------------- ##
    ## Get matching between rownames and gene names
    ## --------------------------------------------------------------------- ##
    dfGene <- data.frame(
        rowName = rownames(SummarizedExperiment::rowData(qft[[1]])),
        genes = SummarizedExperiment::rowData(qft[[1]])[[idCol]]) %>%
        tidyr::separate_rows(.data$genes, sep = ";")

    ## --------------------------------------------------------------------- ##
    ## Complexes
    ## --------------------------------------------------------------------- ##
    if ("complexes" %in% includeFeatureCollections) {
        complexes <- readRDS(complexDbPath)
        if (speciesInfo$speciesCommon %in% names(complexes)) {
            crl <- complexes[[speciesInfo$speciesCommon]]
        } else if (speciesInfo$species %in% names(complexes)) {
            crl <- complexes[[speciesInfo$species]]
        } else {
            stop("No complex database available for the current species")
        }
        if (complexSpecies == "current") {
            ## Only test complexes defined for the current species
            crl <- crl[S4Vectors::mcols(crl)$Species.common %in%
                           c(speciesInfo$species, speciesInfo$speciesCommon)]
        }
        S4Vectors::mcols(crl)$genes <- vapply(
            crl, function(w) gsub(pat, "\\1; ", paste(w, collapse = ";")), "")
        S4Vectors::mcols(crl)$nGenes <- lengths(crl)
    } else {
        crl <- IRanges::CharacterList()
    }

    ## Custom complexes (to add to the list of complexes above)
    if (length(customComplexes) != 0) {
        tmpcompl <- methods::as(customComplexes, "CharacterList")
        S4Vectors::mcols(tmpcompl) <- S4Vectors::DataFrame(
            Species.common = speciesInfo$speciesCommon,
            Source = "custom",
            PMID = NA_character_,
            All.names = names(tmpcompl),
            genes = sapply(customComplexes,
                           function(w) gsub(pat, "\\1; ",
                                            paste(w, collapse = ";"))),
            nGenes = lengths(tmpcompl)
        )
        if (length(crl) > 0) {
            crl <- c(tmpcompl, crl)
        } else {
            crl <- tmpcompl
        }
    }

    if (length(crl) > 0) {
        crl <- .replaceIdsInList(chl = crl, dfConv = dfGene,
                                 currentIdCol = "genes", newIdCol = "rowName",
                                 pat = pat)
        crl <- crl[lengths(crl) >= minSizeToKeep]
        featureCollections$complexes <- crl
    }

    ## --------------------------------------------------------------------- ##
    ## GO terms
    ## --------------------------------------------------------------------- ##
    if ("GO" %in% includeFeatureCollections) {
        goannots <- msigdbr::msigdbr(species = speciesInfo$species, category = "C5") %>%
            dplyr::select(.data$gs_name, .data$gene_symbol)
        goannots <- methods::as(lapply(split(goannots, f = goannots$gs_name),
                              function(w) unique(w$gene_symbol)),
                       "CharacterList")
        S4Vectors::mcols(goannots)$genes <- sapply(goannots, function(w)
            gsub(pat, "\\1; ", paste(w, collapse = ";"))
        )
        S4Vectors::mcols(goannots)$nGenes <- lengths(goannots)
        goannots <- .replaceIdsInList(chl = goannots, dfConv = dfGene,
                                      currentIdCol = "genes", newIdCol = "rowName",
                                      pat = pat)
        goannots <- goannots[lengths(goannots) >= minSizeToKeep]
        featureCollections$GO <- goannots
    }

    featureCollections
}
