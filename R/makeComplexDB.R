#' @keywords internal
#' @noRd
#' @importFrom utils download.file
.getCyc2008Db <- function(dbDir) {
    #nocov start
    ## Yeast (S cerevisiae)
    tryCatch({
        utils::download.file(
            "http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab",
            destfile = file.path(dbDir, "S_cerevisiae_CYC2008_complex.tab")
        )
        read.delim(file.path(dbDir, "S_cerevisiae_CYC2008_complex.tab"))
    }, error = function(e) {
        warning("CYC2008 data could not be downloaded")
        NULL
    })
    #nocov end
}

#' @keywords internal
#' @noRd
#' @importFrom utils download.file unzip
.getCorumDb <- function(dbDir) {
    #nocov start
    ## Mammalian complexes from Corum
    tryCatch({
        utils::download.file(
            "https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip",
            destfile = file.path(dbDir, "CORUM_allComplexes.txt.zip")
        )
        utils::unzip(file.path(dbDir, "CORUM_allComplexes.txt.zip"),
                     exdir = dbDir)
        stopifnot(file.exists(file.path(dbDir, "allComplexes.txt")))
        file.rename(from = file.path(dbDir, "allComplexes.txt"),
                    to = file.path(dbDir, "CORUM_allComplexes.txt"))
        read.delim(file.path(dbDir, "CORUM_allComplexes.txt"))
    }, error = function(e) {
        warning("CORUM data could not be downloaded")
        NULL
    })
    #nocov end
}

#' @keywords internal
#' @noRd
#' @importFrom utils download.file
.getPombaseDb <- function(dbDir) {
    #nocov start
    ## S. pombe complexes
    tryCatch({
        utils::download.file(
            "https://www.pombase.org/data/annotations/Gene_ontology/GO_complexes/Complex_annotation.tsv",
            destfile = file.path(dbDir, "S_pombe_complex_annotation.tsv")
        )
        read.delim(file.path(dbDir, "S_pombe_complex_annotation.tsv"))
    }, error = function(e) {
        warning("PomBase data could not be downloaded")
        NULL
    })
    #nocov end
}

#' Generate a comprehensive database of complexes
#'
#' This function generates a comprehensive cross-species database of
#' protein complexes. It downloads the complex definitions from
#' http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab (S. cerevisiae),
#' https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip
#' (mammals) and
#' https://www.pombase.org/data/annotations/Gene_ontology/GO_complexes/Complex_annotation.tsv
#' (S. pombe), and next uses the \code{babelgene} package to map the
#' complexes to orthologs in the other species.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @param dbDir Path to database directory, where all raw files will be
#'     downloaded and the output will be saved. Will be created if it
#'     doesn't exist.
#' @param customComplexTxt File path to text file with custom complexes
#'     (if any). Should be a tab-delimited text file with five columns:
#'     "Complex.name", "Gene.names", "Organism", "Source", "PMID".
#' @param Cyc2008Db,CorumDb,PombaseDb data.frames providing annotations from
#'     CYC2008 (S cerevisiae), Corum (mammals) and Pombase (S pombe),
#'     respectively. These arguments are provided mainly to allow testing, and
#'     typically are not specified by the end user, except in cases where the
#'     files have already been downloaded and stored locally. If provided,
#'     it is important that the data frames are obtained by simply reading the
#'     downloaded text files - the function assumes a certain set of columns.
#'     If \code{NULL} (the default), the files will be downloaded from the
#'     paths indicated in the Details.
#'
#' @returns Invisibly, the path to the generated complex database.
#'
#' @importFrom dplyr distinct %>%
#' @importFrom babelgene orthologs
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom IRanges CharacterList
#' @importFrom utils packageVersion read.delim download.file unzip
#' @importFrom methods as
#'
makeComplexDB <- function(dbDir, customComplexTxt = NULL, Cyc2008Db = NULL,
                          CorumDb = NULL, PombaseDb = NULL) {
    .assertScalar(x = dbDir, type = "character")
    .assertScalar(x = customComplexTxt, type = "character", allowNULL = TRUE)
    .assertVector(x = Cyc2008Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CorumDb, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = PombaseDb, type = "data.frame", allowNULL = TRUE)
    if (!dir.exists(dbDir)) {
        dir.create(dbDir, recursive = TRUE)
    }
    if (!is.null(customComplexTxt)) {
        stopifnot(file.exists(customComplexTxt))
    }

    ## -------------------------------------------------------------------------
    ## Download complex DB files
    ## -------------------------------------------------------------------------
    if (!is.null(Cyc2008Db)) {
        YEAST.in <- Cyc2008Db
    } else {
        YEAST.in <- .getCyc2008Db(dbDir = dbDir)
    }

    if (!is.null(CorumDb)) {
        CORUM.in <- CorumDb
    } else {
        CORUM.in <- .getCorumDb(dbDir = dbDir)
    }

    if (!is.null(PombaseDb)) {
        SCHPO.in <- PombaseDb
    } else {
        SCHPO.in <- .getPombaseDb(dbDir = dbDir)
    }

    ### custom complexes source: please provide.
    ## (format: c("Complex.name", "Gene.names", "Species.common", "Source")
    # custom.in <- read.delim(file.path(dbDir, "Complexes_custom.txt"))
    if (!is.null(customComplexTxt)) {
        custom.in <- read.delim(customComplexTxt)
    } else {
        custom.in <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Make CharacterLists
    ## -------------------------------------------------------------------------
    ## S. cerevisiae
    if (!is.null(YEAST.in)) {
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
    } else {
        YEAST.chl <- NULL
    }

    ## CORUM
    if (!is.null(CORUM.in)) {
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
            l0 <- lapply(l0, function(m) strsplit(m, "[ ]*;[ ]*")[[1]])
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
    } else {
        CORUM.chl <- NULL
    }

    ## S. pombe
    if (!is.null(SCHPO.in)) {
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
    } else {
        SCHPO.chl <- NULL
    }

    ## Custom
    if (!is.null(custom.in)) {
        custom.chl <- methods::as(lapply(as.list(custom.in$Gene.names),
                                         function(l) strsplit(l, ";")[[1]]),
                                  "CharacterList")
        names(custom.chl) <- paste0(tolower(custom.in$Organism),
                                    ": ", custom.in$Complex.name)
        if ("PMID" %in% colnames(custom.in)) {
            S4Vectors::mcols(custom.chl) <- S4Vectors::DataFrame(
                Species.common = tolower(custom.in$Organism),
                Source = custom.in$Source,
                PMID = custom.in$PMID
            )
        } else {
            S4Vectors::mcols(custom.chl) <- S4Vectors::DataFrame(
                Species.common = tolower(custom.in$Organism),
                Source = custom.in$Source,
                PMID = NA_character_
            )
        }
    } else {
        custom.chl <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Combine
    ## -------------------------------------------------------------------------
    L <- list(YEAST.chl, CORUM.chl$Bovine, CORUM.chl$Dog,
              CORUM.chl$Human, CORUM.chl$Mouse, CORUM.chl$Pig,
              CORUM.chl$Rat, SCHPO.chl, custom.chl)
    all_complexes <- do.call(c, L[vapply(L, function(x) !is.null(x), FALSE)])

    ## Remove any "" or NA entries
    ## First find all complexes that have any "" or NA entries
    W <- which(any(all_complexes == "" | is.na(all_complexes)))
    for (w in W) {
        tmpset <- setdiff(all_complexes[[w]], "")
        tmpset <- tmpset[!is.na(tmpset)]
        all_complexes[[w]] <- tmpset
    }

    complPath <- file.path(
        dbDir, paste0("complexdb_einprot",
                      utils::packageVersion("einprot"), "_",
                      gsub("-", "", Sys.Date()), ".rds")
    )
    saveRDS(all_complexes, file = complPath)

    ## -------------------------------------------------------------------------
    ## Make databases for individual species (find orthologs)
    ## -------------------------------------------------------------------------
    all_orth_complexes <- list()
    for (species_out in c("mouse", "human", "baker's yeast",
                          "Caenorhabditis elegans",
                          "Schizosaccharomyces pombe 972h-")) {
        message(species_out)

        ## Order complexes depending on the species
        if (species_out == "mouse") {
            all_complexes <- list(CORUM.chl$Mouse, CORUM.chl$Rat,
                                  CORUM.chl$Human, CORUM.chl$Bovine,
                                  CORUM.chl$Dog,
                                  CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "human") {
            all_complexes <- list(CORUM.chl$Human, CORUM.chl$Mouse,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog,
                                  CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "baker's yeast") {
            all_complexes <- list(YEAST.chl, SCHPO.chl,
                                  CORUM.chl$Human, CORUM.chl$Mouse,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog,
                                  CORUM.chl$Pig)
        } else if (species_out == "Caenorhabditis elegans") {
            all_complexes <- list(CORUM.chl$Human, CORUM.chl$Mouse,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog,
                                  CORUM.chl$Pig, YEAST.chl, SCHPO.chl)
        } else if (species_out == "Schizosaccharomyces pombe 972h-") {
            all_complexes <- list(SCHPO.chl, YEAST.chl,
                                  CORUM.chl$Human, CORUM.chl$Mouse,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog,
                                  CORUM.chl$Pig)
        } else {
            #nocov start
            stop("Unsupported species")
            #nocov end
        }

        all_complexes <- do.call(
            c, all_complexes[vapply(all_complexes, function(x) !is.null(x),
                                    FALSE)])
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

        ## ---------------------------------------------------------------------
        ## Combine duplicates
        ## ---------------------------------------------------------------------
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

        ## Remove any "" or NA entries
        W <- which(any(orth_complexes == "" | is.na(orth_complexes)))
        for (w in W) {
            tmpset <- setdiff(orth_complexes[[w]], "")
            tmpset <- tmpset[!is.na(tmpset)]
            orth_complexes[[w]] <- tmpset
        }

        all_orth_complexes[[species_out]] <- orth_complexes
    }

    ## -------------------------------------------------------------------------
    ## Save
    ## -------------------------------------------------------------------------
    orthPath <- file.path(
        dbDir, paste0("complexdb_einprot",
                      utils::packageVersion("einprot"), "_",
                      gsub("-", "", Sys.Date()), "_orthologs.rds")
    )
    saveRDS(all_orth_complexes, file = orthPath)

    invisible(list(complPath = complPath, orthPath = orthPath))
}

