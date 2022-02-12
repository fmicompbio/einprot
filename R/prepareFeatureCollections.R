## Help function to replace every nth ";" by "; "
.createPattern <- function(n) {
    sprintf("(%s[^;]+);", strrep("[^;]+;", n - 1))
}

#' @export
#' @author Charlotte Soneson
#'
#' @importFrom S4Vectors mcols DataFrame
#' @importFrom IRanges CharacterList
#' @importFrom methods as
#' @importFrom msigdbr msigdbr
#' @importFrom dplyr select %>%
#'
prepareFeatureCollections <- function(includeFeatureCollections,
                                      complexDbPath, speciesInfo,
                                      complexSpecies, customComplexes) {
    pat <- .createPattern(10)

    featureCollections <- list()

    ## Complexes
    ## -----------------------------------------------------------------
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
        S4Vectors::mcols(crl)$genes <- sapply(crl,
                                   function(w) gsub(pat, "\\1; ",
                                                    paste(w, collapse = ";")))
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
                                            paste(w, collapse = ";")))
        )
        if (length(crl) > 0) {
            crl <- c(crl, tmpcompl)
        } else {
            crl <- tmpcompl
        }
    }

    if (length(crl) > 0) {
        featureCollections$complexes <- crl
    }

    ## GO terms
    ## -----------------------------------------------------------------
    if ("GO" %in% includeFeatureCollections) {
        goannots <- msigdbr::msigdbr(species = speciesInfo$species, category = "C5") %>%
            dplyr::select(gs_name, gene_symbol)
        goannots <- methods::as(lapply(split(goannots, f = goannots$gs_name),
                              function(w) unique(w$gene_symbol)),
                       "CharacterList")
        S4Vectors::mcols(goannots)$genes <- sapply(goannots, function(w)
            gsub(pat, "\\1; ", paste(w, collapse = ";"))
        )

        featureCollections$GO <- goannots
    }

    featureCollections
}
