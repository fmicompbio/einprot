#' Generate a comprehensive database of complexes
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @param dbDir Path to database directory, where all raw files will be
#'     downloaded and the output will be saved. Will be created if it
#'     doesn't exist.
#' @param customComplexTxt File path to text file with custom complexes
#'     (if any). Should be a tab-delimited text file with four columns:
#'     "Complex.name", "Gene.names", "Species.common", "Source".
#'
#' @return Invisibly, the path to the generated complex database.
#'
#' @importFrom dplyr distinct %>%
#' @importFrom babelgene orthologs
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom IRanges CharacterList
#' @importFrom utils packageVersion read.delim download.file unzip
#' @importFrom methods as
#'
makeComplexDB <- function(dbDir, customComplexTxt = NULL) {
    if (!dir.exists(dbDir)) {
        dir.create(dbDir, recursive = TRUE)
    }

    ## --------------------------------------------------------------------- ##
    ## Download complex DB files
    ## --------------------------------------------------------------------- ##
    ## Yeast (S cerevisiae)
    utils::download.file(
        "http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab",
        destfile = file.path(dbDir, "S_cerevisiae_CYC2008_complex.tab")
    )
    YEAST.in <- read.delim(file.path(dbDir, "S_cerevisiae_CYC2008_complex.tab"))

    ## Mammalian complexes from Corum
    utils::download.file(
        "http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip",
        destfile = file.path(dbDir, "CORUM_allComplexes.txt.zip")
    )
    utils::unzip(file.path(dbDir, "CORUM_allComplexes.txt.zip"), exdir = dbDir)
    stopifnot(file.exists(file.path(dbDir, "allComplexes.txt")))
    file.rename(from = file.path(dbDir, "allComplexes.txt"),
                to = file.path(dbDir, "CORUM_allComplexes.txt"))
    CORUM.in <- read.delim(file.path(dbDir, "CORUM_allComplexes.txt"))

    ## S. pombe complexes
    utils::download.file(
        "https://www.pombase.org/data/annotations/Gene_ontology/GO_complexes/Complex_annotation.tsv",
        destfile = file.path(dbDir, "S_pombe_complex_annotation.tsv")
    )
    SCHPO.in <- read.delim(file.path(dbDir, "S_pombe_complex_annotation.tsv"))

    ### custom complexes source: please provide.
    ## (format: c("Complex.name", "Gene.names", "Species.common", "Source")
    # custom.in <- read.delim(file.path(dbDir, "Complexes_custom.txt"))
    if (!is.null(customComplexTxt)) {
        custom.in <- read.delim(customComplexTxt)
    } else {
        custom.in <- NULL
    }

    ## --------------------------------------------------------------------- ##
    ## Make CharacterLists
    ## --------------------------------------------------------------------- ##
    ## S. cerevisiae
    YEAST.chl <- methods::as(split(YEAST.in$Name,
                                   f = paste("S.cer:", YEAST.in$Complex)),
                             "CharacterList")
    S4Vectors::mcols(YEAST.chl) <- S4Vectors::DataFrame(
        Species.common = "baker's yeast",
        Source = "CYC2008"
    )

    ## CORUM
    CORUM.in <- CORUM.in[, c("Organism", "ComplexName",
                             "subunits.Gene.name.")] %>%
        dplyr::distinct()
    CORUM.list <- split(CORUM.in, CORUM.in$Organism)
    CORUM.list <- lapply(CORUM.list, function(l) {
        l$ComplexName <- make.unique(l$ComplexName, sep = "-variant")
        l
    })
    CORUM.chl <- lapply(CORUM.list, function(l) {
        l0 <- split(l$subunits.Gene.name., f = paste0(tolower(l$Organism),
                                                      ": ", l$ComplexName))
        l0 <- lapply(l0, function(m) strsplit(m, ";")[[1]])
        l0 <- methods::as(l0, "CharacterList")
        S4Vectors::mcols(l0) <- S4Vectors::DataFrame(
            Species.common = tolower(l$Organism[1]),
            Source = "CORUM"
        )
        l0
    })
    CORUM.chl <- c(CORUM.chl$Bovine, CORUM.chl$Dog, CORUM.chl$Human,
                   CORUM.chl$Mouse, CORUM.chl$Pig, CORUM.chl$Rat)

    ## S. pombe
    SCHPO.chl <- methods::as(split(SCHPO.in$symbol,
                                   f = paste("S.pombe:", SCHPO.in$GO_name)),
                             "CharacterList")
    S4Vectors::mcols(SCHPO.chl) <- S4Vectors::DataFrame(
        Species.common = "Schizosaccharomyces pombe 972h-",
        Source = "pombase"
    )

    ## Custom
    if (!is.null(custom.in)) {
        custom.chl <- methods::as(lapply(as.list(custom.in$Gene.names),
                                         function(l) strsplit(l, ";")[[1]]),
                                  "CharacterList")
        names(custom.chl) <- paste0(tolower(custom.in$Organism),
                                    ": ", custom.in$Complex.name)
        S4Vectors::mcols(custom.chl) <- S4Vectors::DataFrame(
            Species.common = tolower(custom.in$Organism),
            Source = custom.in$Source
        )
    }

    ## --------------------------------------------------------------------- ##
    ## Combine
    ## --------------------------------------------------------------------- ##
    all_complexes <- c(YEAST.chl, CORUM.chl, SCHPO.chl, custom.chl)

    saveRDS(all_complexes, file = file.path(
        dbDir, paste0("complexdb_einprot", utils::packageVersion("einprot"), "_",
                       gsub("-", "", Sys.Date()), ".rds")
    ))

    ## --------------------------------------------------------------------- ##
    ## Make databases for individual species (find orthologs)
    ## --------------------------------------------------------------------- ##
    all_orth_complexes <- list()
    for (species_out in c("mouse", "human", "baker's yeast", "roundworm",
                          "Schizosaccharomyces pombe 972h-")) {
        print(species_out)
        orth_complexes <- all_complexes
        for (i in seq_along(all_complexes)) {
            cplx <- names(all_complexes)[i]
            orgfrom <- S4Vectors::mcols(all_complexes)$Species.common[i]
            if (orgfrom == "human") {
                if (species_out == "human") {
                    out <- all_complexes[[cplx]]
                } else {
                    out <- tryCatch({
                        babelgene::orthologs(all_complexes[[cplx]],
                                             species = species_out,
                                             human = TRUE)$symbol},
                        error = function(e) c())
                }
            } else {
                if (species_out == orgfrom) {
                    out <- all_complexes[[cplx]]
                } else {
                    outtmp <- tryCatch({
                        babelgene::orthologs(all_complexes[[cplx]],
                                             species = orgfrom,
                                             human = FALSE)$human_symbol},
                        error = function(e) c())
                    if (species_out == "human") {
                        out <- outtmp
                    } else {
                        out <- tryCatch({
                            babelgene::orthologs(outtmp,
                                                 species = species_out,
                                                 human = TRUE)$symbol},
                            error = function(e) c())
                    }
                }
            }
            orth_complexes[[cplx]] <- out
        }
        all_orth_complexes[[species_out]] <- orth_complexes
    }

    ## --------------------------------------------------------------------- ##
    ## Save
    ## --------------------------------------------------------------------- ##
    orthPath <- file.path(
        dbDir, paste0("complexdb_einprot",
                      utils::packageVersion("einprot"), "_",
                      gsub("-", "", Sys.Date()), "_orthologs.rds")
    )
    saveRDS(all_orth_complexes, file = orthPath)

    invisible(orthPath)
}

