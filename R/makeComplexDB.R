#' Generate a comprehensive database of complexes
#'
#' This function generates a comprehensive cross-species database of
#' protein complexes. It downloads the complex definitions from
#' http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab (S. cerevisiae),
#' http://mips.helmholtz-muenchen.de/corum/download/allComplexes.txt.zip
#' (mammals) and
#' https://www.pombase.org/data/annotations/Gene_ontology/GO_complexes/Complex_annotation.tsv
#' (S. pombe), and next uses the \code{babelgene} package to map the
#' complexes to orthologues in the other species.
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
    .assertScalar(x = dbDir, type = "character")
    .assertScalar(x = customComplexTxt, type = "character", allowNULL = TRUE)
    if (!dir.exists(dbDir)) {
        dir.create(dbDir, recursive = TRUE)
    }
    if (!is.null(customComplexTxt)) {
        stopifnot(file.exists(customComplexTxt))
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
    ## Add PMID
    YEAST.pmid <- YEAST.in[, c("Complex", "PubMed_id")] %>%
        dplyr::filter(!is.na(.data$PubMed_id)) %>%
        dplyr::distinct() %>%
        dplyr::mutate(Complex = paste("S.cer:", .data$Complex))
    stopifnot(!any(duplicated(YEAST.pmid$Complex)))
    S4Vectors::mcols(YEAST.chl)$PMID <-
        YEAST.pmid$PubMed_id[match(names(YEAST.chl),
                                   YEAST.pmid$Complex)]

    ## CORUM
    CORUM.in <- CORUM.in[, c("Organism", "ComplexName",
                             "subunits.Gene.name.", "PubMed.ID")] %>%
        dplyr::group_by(.data$Organism, .data$ComplexName,
                        .data$subunits.Gene.name.) %>%
        dplyr::summarize(PubMed.ID = paste(sort(.data$PubMed.ID),
                                           collapse = ";")) %>%
        dplyr::distinct() %>%
        dplyr::ungroup()
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
        ## Add PMID
        S4Vectors::mcols(l0)$PMID <-
            l$PubMed.ID[match(names(l0), paste0(tolower(l$Organism), ": ",
                                                l$ComplexName))]
        l0
    })

    ## S. pombe
    SCHPO.chl <- methods::as(split(SCHPO.in$symbol,
                                   f = paste("S.pombe:", SCHPO.in$GO_name)),
                             "CharacterList")
    S4Vectors::mcols(SCHPO.chl) <- S4Vectors::DataFrame(
        Species.common = "Schizosaccharomyces pombe 972h-",
        Source = "pombase"
    )
    ## Add PMID
    SCHPO.pmid <- SCHPO.in[, c("GO_name", "source")] %>%
        dplyr::filter(!is.na(source)) %>%
        dplyr::group_by(.data$GO_name) %>%
        dplyr::summarize(source = paste(
            sort(unique(unlist(split(gsub("PMID:", "",
                                          .data$source), ",")[[1]]))),
            collapse = ";")) %>%
        dplyr::distinct() %>%
        dplyr::ungroup() %>%
        dplyr::mutate(GO_name = paste("S.pombe:", .data$GO_name))
    stopifnot(!any(duplicated(SCHPO.pmid$GO_name)))
    S4Vectors::mcols(SCHPO.chl)$PMID <-
        SCHPO.pmid$source[match(names(SCHPO.chl),
                                SCHPO.pmid$GO_name)]

    ## Custom
    if (!is.null(custom.in)) {
        custom.chl <- methods::as(lapply(as.list(custom.in$Gene.names),
                                         function(l) strsplit(l, ";")[[1]]),
                                  "CharacterList")
        names(custom.chl) <- paste0(tolower(custom.in$Organism),
                                    ": ", custom.in$Complex.name)
        S4Vectors::mcols(custom.chl) <- S4Vectors::DataFrame(
            Species.common = tolower(custom.in$Organism),
            Source = custom.in$Source,
            PMID = NA_character_
        )
    } else {
        custom.chl <- NULL
    }

    ## --------------------------------------------------------------------- ##
    ## Combine
    ## --------------------------------------------------------------------- ##
    all_complexes <- c(YEAST.chl, CORUM.chl$Bovine, CORUM.chl$Dog,
                       CORUM.chl$Human, CORUM.chl$Mouse, CORUM.chl$Pig,
                       CORUM.chl$Rat, SCHPO.chl, custom.chl)

    complPath <- file.path(
        dbDir, paste0("complexdb_einprot",
                      utils::packageVersion("einprot"), "_",
                      gsub("-", "", Sys.Date()), ".rds")
    )
    saveRDS(all_complexes, file = complPath)

    ## --------------------------------------------------------------------- ##
    ## Make databases for individual species (find orthologs)
    ## --------------------------------------------------------------------- ##
    all_orth_complexes <- list()
    for (species_out in c("mouse", "human", "baker's yeast", "roundworm",
                          "Schizosaccharomyces pombe 972h-")) {
        print(species_out)

        ## Order complexes depending on the species
        if (species_out == "mouse") {
            all_complexes <- c(CORUM.chl$Mouse, CORUM.chl$Rat,
                               CORUM.chl$Human, CORUM.chl$Bovine, CORUM.chl$Dog,
                               CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "human") {
            all_complexes <- c(CORUM.chl$Human, CORUM.chl$Mouse,
                               CORUM.chl$Rat, CORUM.chl$Bovine, CORUM.chl$Dog,
                               CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "baker's yeast") {
            all_complexes <- c(YEAST.chl, SCHPO.chl,
                               CORUM.chl$Human, CORUM.chl$Mouse,
                               CORUM.chl$Rat, CORUM.chl$Bovine, CORUM.chl$Dog,
                               CORUM.chl$Pig)
        } else if (species_out == "roundworm") {
            all_complexes <- c(CORUM.chl$Human, CORUM.chl$Mouse,
                               CORUM.chl$Rat, CORUM.chl$Bovine, CORUM.chl$Dog,
                               CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "Schizosaccharomyces pombe 972h-") {
            all_complexes <- c(SCHPO.chl, YEAST.chl,
                               CORUM.chl$Human, CORUM.chl$Mouse,
                               CORUM.chl$Rat, CORUM.chl$Bovine, CORUM.chl$Dog,
                               CORUM.chl$Pig)
        } else {
            stop("Unsupported species")
        }

        if (!is.null(custom.chl)) {
            all_complexes <- c(custom.chl, all_complexes)
        }

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
            orth_complexes[[cplx]] <- unique(out)
        }

        ## ----------------------------------------------------------------- ##
        ## Combine duplicates
        ## ----------------------------------------------------------------- ##
        ## Sort all vectors
        orth_compl_sort <- lapply(as.list(orth_complexes),
                                  function(l) sort(l))
        stopifnot(names(orth_complexes) == names(orth_compl_sort))

        ## Uniquify list of complexes
        idx_unique <- which(!duplicated(orth_compl_sort))
        unique_complexes <- orth_compl_sort[idx_unique]

        ## Group unique complexes
        group <- match(orth_compl_sort, unique_complexes)

        ## Modify names of complex list
        S4Vectors::mcols(orth_complexes)$All.names <- names(orth_complexes)
        tmpnames <- names(orth_complexes)
        tmppmid <- S4Vectors::mcols(orth_complexes)$PMID
        for (i in seq_along(orth_complexes)) {
            w <- which(group == group[i])
            if (length(w) > 1) {
                ## Use tmpnames/tmppmid above to avoid combining
                ## already modified names/PMIDs
                S4Vectors::mcols(orth_complexes)$All.names[i] <-
                    paste(tmpnames[w], collapse = ";")
                S4Vectors::mcols(orth_complexes)$PMID[i] <-
                    paste(tmppmid[w], collapse = ";")
                names(orth_complexes)[i] <-
                    paste0(names(orth_complexes)[i],
                           " (+", length(w) - 1,
                           ifelse(length(w) == 2, " alt. ID)", " alt. IDs)"))
            }
        }

        ## Remove duplicated entries
        orth_complexes <- orth_complexes[idx_unique]

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

    invisible(list(complPath = complPath, orthPath = orthPath))
}

