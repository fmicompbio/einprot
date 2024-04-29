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
        warning("CYC2008 data could not be downloaded", immediate. = TRUE)
        NULL
    })
    #nocov end
}

#' @keywords internal
#' @noRd
#' @importFrom utils download.file
.getComplexPortalDb <- function(dbDir, species) {
    #nocov start
    if (!species %in% c("human", "mouse",
                        "Caenorhabditis elegans",
                        "Schizosaccharomyces pombe 972h-")) {
        warning("Unsupported species for ComplexPortal: ", species)
        return(NULL)
    }
    tryCatch({
        taxid <- getSpeciesInfo(species)$taxId
        utils::download.file(
            paste0("https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/",
                   taxid, ".tsv"),
            destfile = file.path(dbDir, paste0("ComplexPortal_", taxid, "_complex.tsv"))
        )
        read.delim(file.path(dbDir, paste0("ComplexPortal_", taxid, "_complex.tsv")))
    }, error = function(e) {
        warning("ComplexPortal ", species, " data could not be downloaded")
        NULL
    })
    #nocov end
}

#' @keywords internal
#' @noRd
#' @importFrom utils download.file
.getHuMAP2Db <- function(dbDir) {
    #nocov start
    tryCatch({
        utils::download.file(
            "http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt",
            destfile = file.path(dbDir, "HuMAP2_complex.txt")
        )
        read.delim(file.path(dbDir, "HuMAP2_complex.txt"), sep = ",")
    }, error = function(e) {
        warning("HuMAP2 data could not be downloaded", immediate. = TRUE)
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
        warning("CORUM data could not be downloaded", immediate. = TRUE)
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
        warning("PomBase data could not be downloaded", immediate. = TRUE)
        NULL
    })
    #nocov end
}

#' Generate a comprehensive database of complexes
#'
#' This function generates a comprehensive cross-species database of
#' protein complexes. It downloads the complex definitions from
#' \url{http://wodaklab.org/cyc2008/resources/CYC2008_complex.tab} (S. cerevisiae),
#' \url{https://mips.helmholtz-muenchen.de/corum/download/releases/current/allComplexes.txt.zip}
#' (mammals),
#' \url{https://www.pombase.org/data/annotations/Gene_ontology/GO_complexes/Complex_annotation.tsv}
#' (S. pombe),
#' \url{http://humap2.proteincomplexes.org/static/downloads/humap2/humap2_complexes_20200809.txt}
#' (human) and
#' \url{https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/}
#' (human, mouse, C. elegans and S. pombe),
#' and next uses the \code{babelgene} package to map the
#' complexes to orthologs in the other species. A pre-generated version of the
#' database is provided with \code{einprot} (see \code{listComplexDBs()}).
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
#'     \code{"Complex.name"}, \code{"Gene.names"} (semi-colon separated),
#'     \code{"Organism"}, \code{"Source"}, \code{"PMID"}.
#' @param Cyc2008Db,CorumDb,PombaseDb,HuMAP2Db,CPortal9606Db,CPortal10090Db,CPortal6239Db,CPortal284812Db data.frames providing
#'     annotations from CYC2008 (S cerevisiae), Corum (mammals),
#'     Pombase (S pombe), HuMAP2 (human) and the Complex Portal (human, mouse,
#'     C. elegans, S. pombe),
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
#' @examples
#' ## Read small subsets of the raw files provided with einprot to make the
#' ## processing faster. Typically, the files would be downloaded as part of
#' ## the process of generating the complex DB.
#' cyc2008db <- read.delim(system.file("extdata", "complexes",
#'                                     "cyc2008_complex_extract.tab",
#'                                     package = "einprot"))
#' corumdb <- read.delim(system.file("extdata", "complexes",
#'                                   "corum_complex_extract.txt",
#'                                   package = "einprot"))
#' pombasedb <- read.delim(system.file("extdata", "complexes",
#'                                     "pombase_complex_extract.tsv",
#'                                     package = "einprot"))
#' humap2db <- read.delim(system.file("extdata", "complexes",
#'                                    "humap2_complex_extract.txt",
#'                                    package = "einprot"), sep = ",")
#' cportal9606db <- read.delim(system.file("extdata", "complexes",
#'                                         "complexportal9606_complex_extract.tsv",
#'                                         package = "einprot"))
#' cportal10090db <- read.delim(system.file("extdata", "complexes",
#'                                          "complexportal10090_complex_extract.tsv",
#'                                          package = "einprot"))
#' cportal6239db <- read.delim(system.file("extdata", "complexes",
#'                                         "complexportal6239_complex_extract.tsv",
#'                                         package = "einprot"))
#' cportal284812db <- read.delim(system.file("extdata", "complexes",
#'                                           "complexportal284812_complex_extract.tsv",
#'                                           package = "einprot"))
#' dbdir <- tempdir()
#' dbs <- makeComplexDB(dbDir = dbdir, Cyc2008Db = cyc2008db,
#'                      CorumDb = corumdb, PombaseDb = pombasedb,
#'                      HuMAP2Db = humap2db, CPortal9606Db = cportal9606db,
#'                      CPortal10090Db = cportal10090db,
#'                      CPortal6239Db = cportal6239db,
#'                      CPortal284812Db = cportal284812db)
#'
#' ## List of complexes
#' compl <- readRDS(dbs$complPath)
#' compl
#'
#' ## Complexes mapped to all species via orthologs
#' orth <- readRDS(dbs$orthPath)
#' orth
#'
#' @importFrom dplyr distinct %>% select
#' @importFrom tidyr unite
#' @importFrom babelgene orthologs
#' @importFrom S4Vectors DataFrame mcols
#' @importFrom IRanges CharacterList
#' @importFrom utils packageVersion read.delim download.file unzip
#' @importFrom methods as
#'
makeComplexDB <- function(dbDir, customComplexTxt = NULL, Cyc2008Db = NULL,
                          CorumDb = NULL, PombaseDb = NULL,
                          HuMAP2Db = NULL, CPortal9606Db = NULL,
                          CPortal10090Db = NULL, CPortal6239Db = NULL,
                          CPortal284812Db = NULL) {
    .assertScalar(x = dbDir, type = "character")
    .assertScalar(x = customComplexTxt, type = "character", allowNULL = TRUE)
    .assertVector(x = Cyc2008Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CorumDb, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = PombaseDb, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = HuMAP2Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CPortal9606Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CPortal10090Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CPortal6239Db, type = "data.frame", allowNULL = TRUE)
    .assertVector(x = CPortal284812Db, type = "data.frame", allowNULL = TRUE)
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

    if (!is.null(HuMAP2Db)) {
        HUMAP2.in <- HuMAP2Db
    } else {
        HUMAP2.in <- .getHuMAP2Db(dbDir = dbDir)
    }

    if (!is.null(CPortal9606Db)) {
        CPortal9606.in <- CPortal9606Db
    } else {
        CPortal9606.in <- .getComplexPortalDb(dbDir = dbDir,
                                              species = getSpeciesInfo(9606)$speciesCommon)
    }

    if (!is.null(CPortal10090Db)) {
        CPortal10090.in <- CPortal10090Db
    } else {
        CPortal10090.in <- .getComplexPortalDb(dbDir = dbDir,
                                               species = getSpeciesInfo(10090)$speciesCommon)
    }

    if (!is.null(CPortal6239Db)) {
        CPortal6239.in <- CPortal6239Db
    } else {
        CPortal6239.in <- .getComplexPortalDb(dbDir = dbDir,
                                              species = getSpeciesInfo(6239)$speciesCommon)
    }

    if (!is.null(CPortal284812Db)) {
        CPortal284812.in <- CPortal284812Db
    } else {
        CPortal284812.in <- .getComplexPortalDb(dbDir = dbDir,
                                                species = getSpeciesInfo(284812)$speciesCommon)
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

    ## HuMAP2
    if (!is.null(HUMAP2.in)) {
        HUMAP2.in$HuMAP2_ID <- paste0(HUMAP2.in$HuMAP2_ID, "_conf",
                                      HUMAP2.in$Confidence)
        HUMAP2.chl <- split(HUMAP2.in$genenames,
                            f = paste("human:", HUMAP2.in$HuMAP2_ID))
        HUMAP2.chl <- lapply(HUMAP2.chl, function(m) strsplit(m, " ")[[1]])
        HUMAP2.chl <- methods::as(HUMAP2.chl, "CharacterList")
        S4Vectors::mcols(HUMAP2.chl) <- S4Vectors::DataFrame(
            Species.common = tolower("human"),
            Source = "HuMAP2",
            PMID = ""
        )
    } else {
        HUMAP2.chl <- NULL
    }

    ## ComplexPortal
    if (!is.null(CPortal9606.in)) {
        mapdf <- getUniProtToIDMapping(9606, targetId = "Gene_Name")
        tmp <- CPortal9606.in %>%
            dplyr::select(X.Complex.ac, Recommended.name, Expanded.participant.list) %>%
            tidyr::unite(X.Complex.ac, Recommended.name, col = "Name", sep = "_")
        spc <- getSpeciesInfo(9606)$speciesCommon
        CPortal9606.chl <- split(tmp$Expanded.participant.list,
                                 f = paste0(spc, ": ", tmp$Name))
        CPortal9606.chl <- lapply(CPortal9606.chl, function(m) {
            mapdf$Gene_Name[match(sub("\\([0-9]+\\)", "", strsplit(m, "\\|")[[1]]),
                                  mapdf$UniProtID)]
        })
        CPortal9606.chl <- methods::as(CPortal9606.chl, "CharacterList")
        S4Vectors::mcols(CPortal9606.chl) <- S4Vectors::DataFrame(
            Species.common = spc,
            Source = "ComplexPortal",
            PMID = ""
        )
    } else {
        CPortal9606.chl <- NULL
    }

    if (!is.null(CPortal10090.in)) {
        mapdf <- getUniProtToIDMapping(10090, targetId = "Gene_Name")
        tmp <- CPortal10090.in %>%
            dplyr::select(X.Complex.ac, Recommended.name, Expanded.participant.list) %>%
            tidyr::unite(X.Complex.ac, Recommended.name, col = "Name", sep = "_")
        spc <- getSpeciesInfo(10090)$speciesCommon
        CPortal10090.chl <- split(tmp$Expanded.participant.list,
                                  f = paste0(spc, ": ", tmp$Name))
        CPortal10090.chl <- lapply(CPortal10090.chl, function(m) {
            mapdf$Gene_Name[match(sub("\\([0-9]+\\)", "", strsplit(m, "\\|")[[1]]),
                                  mapdf$UniProtID)]
        })
        CPortal10090.chl <- methods::as(CPortal10090.chl, "CharacterList")
        S4Vectors::mcols(CPortal10090.chl) <- S4Vectors::DataFrame(
            Species.common = spc,
            Source = "ComplexPortal",
            PMID = ""
        )
    } else {
        CPortal10090.chl <- NULL
    }

    if (!is.null(CPortal284812.in)) {
        mapdf <- getUniProtToIDMapping(284812, targetId = "Gene_Name")
        tmp <- CPortal284812.in %>%
            dplyr::select(X.Complex.ac, Recommended.name, Expanded.participant.list) %>%
            tidyr::unite(X.Complex.ac, Recommended.name, col = "Name", sep = "_")
        spc <- getSpeciesInfo(284812)$species
        CPortal284812.chl <- split(tmp$Expanded.participant.list,
                                   f = paste0("S.pombe: ", tmp$Name))
        CPortal284812.chl <- lapply(CPortal284812.chl, function(m) {
            mapdf$Gene_Name[match(sub("\\([0-9]+\\)", "", strsplit(m, "\\|")[[1]]),
                                  mapdf$UniProtID)]
        })
        CPortal284812.chl <- methods::as(CPortal284812.chl, "CharacterList")
        S4Vectors::mcols(CPortal284812.chl) <- S4Vectors::DataFrame(
            Species.common = spc,
            Source = "ComplexPortal",
            PMID = ""
        )
    } else {
        CPortal284812.chl <- NULL
    }

    if (!is.null(CPortal6239.in)) {
        mapdf <- getUniProtToIDMapping(6239, targetId = "Gene_Name")
        tmp <- CPortal6239.in %>%
            dplyr::select(X.Complex.ac, Recommended.name, Expanded.participant.list) %>%
            tidyr::unite(X.Complex.ac, Recommended.name, col = "Name", sep = "_")
        spc <- getSpeciesInfo(6239)$species
        CPortal6239.chl <- split(tmp$Expanded.participant.list,
                                 f = paste0("C.elegans: ", tmp$Name))
        CPortal6239.chl <- lapply(CPortal6239.chl, function(m) {
            mapdf$Gene_Name[match(sub("\\([0-9]+\\)", "", strsplit(m, "\\|")[[1]]),
                                  mapdf$UniProtID)]
        })
        CPortal6239.chl <- methods::as(CPortal6239.chl, "CharacterList")
        S4Vectors::mcols(CPortal6239.chl) <- S4Vectors::DataFrame(
            Species.common = spc,
            Source = "ComplexPortal",
            PMID = ""
        )
    } else {
        CPortal6239.chl <- NULL
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
              CORUM.chl$Rat, SCHPO.chl, HUMAP2.chl, CPortal6239.chl,
              CPortal284812.chl, CPortal10090.chl, CPortal9606.chl, custom.chl)
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
            all_complexes <- list(CORUM.chl$Mouse, CPortal10090.chl,
                                  CORUM.chl$Rat, CORUM.chl$Human,
                                  CPortal9606.chl, HUMAP2.chl, CORUM.chl$Bovine,
                                  CORUM.chl$Dog, CORUM.chl$Pig, CPortal6239.chl,
                                  YEAST.chl, SCHPO.chl, CPortal284812.chl)
        } else if (species_out == "human") {
            all_complexes <- list(CORUM.chl$Human, CPortal9606.chl, HUMAP2.chl,
                                  CORUM.chl$Mouse, CPortal10090.chl,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog, CORUM.chl$Pig, CPortal6239.chl,
                                  YEAST.chl, SCHPO.chl, CPortal284812.chl)
        } else if (species_out == "baker's yeast") {
            all_complexes <- list(YEAST.chl, SCHPO.chl, CPortal284812.chl,
                                  CORUM.chl$Human, CPortal9606.chl, HUMAP2.chl,
                                  CORUM.chl$Mouse, CPortal10090.chl,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog, CORUM.chl$Pig, CPortal6239.chl)
        } else if (species_out == "Caenorhabditis elegans") {
            all_complexes <- list(CPortal6239.chl, CORUM.chl$Human,
                                  CPortal9606.chl, HUMAP2.chl, CORUM.chl$Mouse,
                                  CPortal10090.chl,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog, CORUM.chl$Pig,
                                  YEAST.chl, SCHPO.chl, CPortal284812.chl)
        } else if (species_out == "Schizosaccharomyces pombe 972h-") {
            all_complexes <- list(SCHPO.chl, CPortal284812.chl, YEAST.chl,
                                  CORUM.chl$Human, CPortal9606.chl, HUMAP2.chl,
                                  CORUM.chl$Mouse, CPortal10090.chl,
                                  CORUM.chl$Rat, CORUM.chl$Bovine,
                                  CORUM.chl$Dog, CORUM.chl$Pig, CPortal6239.chl)
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

