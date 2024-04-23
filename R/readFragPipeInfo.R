#' Extract information from FragPipe logs
#'
#' Read FragPipe config/log files and extract information about the run.
#'
#' @param fragpipeDir Character scalar, the path to a FragPipe output
#'     directory. Should contain files <fragpipeDir>/fragpipe_*.config (or
#'     <fragpipeDir>/fragpipe.workflow) and
#'     <fragpipeDir>/log_*.txt (if not, the corresponding fields will
#'     be empty in the output).
#'
#' @author Charlotte Soneson, Jan Seebacher
#' @export
#'
#' @examples
#' readFragPipeInfo(system.file("extdata", "fp_example",
#'                              package = "einprot"))
#'
#' @returns A list with extracted information about the \code{FragPipe} run.
#'
#' @importFrom tidyr separate_wider_delim separate
#' @importFrom dplyr mutate filter bind_rows
#'
readFragPipeInfo <- function(fragpipeDir) {
    .assertScalar(x = fragpipeDir, type = "character")

    ## -------------------------------------------------------------------------
    ## Get file names from FragPipe directory
    ## -------------------------------------------------------------------------
    fpConfigFile <- list.files(fragpipeDir, pattern = "^fragpipe.+.config$",
                               full.names = TRUE)
    ## No config file - try workflow file
    if (length(fpConfigFile) == 0) {
        fpConfigFile <- list.files(fragpipeDir, pattern = "^fragpipe.workflow$",
                                   full.names = TRUE)
    }
    stopifnot(length(fpConfigFile) <= 1)
    fpLogFile <- list.files(fragpipeDir, pattern = "^log_.+.txt$",
                            full.names = TRUE)
    stopifnot(length(fpLogFile) <= 1)

    ## -------------------------------------------------------------------------
    ## Read config file and log file
    ## -------------------------------------------------------------------------
    if (length(fpConfigFile) == 1) {
        configDf <- read.delim(fpConfigFile, header = FALSE, sep = "=") %>%
            setNames(c("parameter", "value"))
    } else {
        configDf <- NULL
    }
    if (length(fpLogFile) == 1) {
        logDf <- read.delim(fpLogFile, blank.lines.skip = TRUE,
                            header = FALSE) %>%
            setNames("info")
    } else {
        logDf <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Version info
    ## -------------------------------------------------------------------------
    if (!is.null(configDf)) {
        fpVersion <- gsub(".+\\((.+)\\).+", "\\1", configDf$parameter[1])
    } else if (!is.null(logDf)) {
        fpVersion <- paste(gsub("# (.+)ui state cache", "\\1",
                                logDf$info[grep("# FragPipe v", logDf$info)]),
                           collapse = ", ")
    } else {
        fpVersion <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Get info about raw files/experiments
    ## -------------------------------------------------------------------------
    if (!is.null(logDf)) {
        fpRawTxt <- gsub("(.+)\\:.+", "\\1",
                         gsub("\\\\", "/",
                              unlist(logDf[grep("\\.raw: Scans = ",
                                                logDf$info), ])))
        fpRawDirs <- paste(unique(dirname(fpRawTxt)), collapse = ", ")
        fpRawFiles <- paste(basename(fpRawTxt), collapse = ", ")
    } else {
        fpRawTxt <- fpRawDirs <- fpRawFiles <- NULL
    }

    ## Experiments
    if (!is.null(logDf)) {
        fpExperiments <- paste(
            gsub("  Experiment/Group: ", "",
                 grep("Experiment/Group: ", logDf$info, value = TRUE)),
            collapse = ", ")
    } else {
        fpExperiments <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Simplify logDf
    ## -------------------------------------------------------------------------
    if (!is.null(logDf)) {
        logDfConf <- data.frame(
            info = logDf$info[(grep("fragpipe.config", logDf$info)[1] + 4):
                                  grep("workflow.threads=", logDf$info)[1]]) %>%
            tidyr::separate_wider_delim("info", names = c("parameter", "value"),
                                        delim = "=", too_many = "merge")
    } else {
        logDfConf <- NULL
    }

    ## Use configDf if it exists, otherwise logDfConf
    if (!is.null(configDf)) {
        dataDf <- configDf
    } else if (!is.null(logDfConf)) {
        dataDf <- logDfConf
    } else {
        dataDf <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Fixed modifications
    ## -------------------------------------------------------------------------
    if (!is.null(dataDf)) {
        fpFixedMods <- dataDf$value[which(dataDf$parameter ==
                                              "msfragger.table.fix-mods")]
        fpFixedModsDf <- as.data.frame(unlist(strsplit(fpFixedMods, ";"))) %>%
            setNames("orig") %>%
            tidyr::separate(.data$orig,
                            into = c("modMass", "modSymbol", "logi", "num"),
                            sep = ",", remove = FALSE) %>%
            dplyr::mutate(modSymbol = gsub(" $", "", .data$modSymbol),
                          modMass = round(as.numeric(.data$modMass), 4),
                          modAA = vapply(strsplit(.data$modSymbol, " "),
                                         .subset, 1,
                                         FUN.VALUE = ""),
                          logi = as.logical(.data$logi),
                          modType = "fixed")
    } else {
        fpFixedModsDf <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Variable modifications
    ## -------------------------------------------------------------------------
    if (!is.null(dataDf)) {
        fpVarMods <- dataDf$value[which(dataDf$parameter ==
                                            "msfragger.table.var-mods")]
        fpVarModsDf <- as.data.frame(unlist(strsplit(fpVarMods, ";"))) %>%
            setNames("orig") %>%
            tidyr::separate(.data$orig,
                            into = c("modMass", "modSymbol", "logi", "num"),
                            sep = ",", remove = FALSE) %>%
            dplyr::mutate(modSymbol = gsub(" $", "", .data$modSymbol),
                          modMass = round(as.numeric(.data$modMass), 4),
                          modAA = .data$modSymbol,
                          logi = as.logical(.data$logi),
                          modType = "variable")
    } else {
        fpVarModsDf <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Put modifications together
    ## -------------------------------------------------------------------------
    if (!is.null(fpFixedModsDf) && !is.null(fpVarModsDf)) {
        fpAllModsDf <- dplyr::bind_rows(
            fpFixedModsDf %>% dplyr::filter(.data$modMass > 0),
            fpVarModsDf %>% dplyr::filter(.data$logi)
        ) %>%
            dplyr::mutate(modName = paste0(.data$modAA, "(",
                                           .data$modMass, ")")) %>%
            dplyr::mutate(modName = gsub("\\[\\^", "N-term", .data$modName)) %>%
            dplyr::mutate(modNameRegex = gsub("\\(", "\\\\(",
                                              gsub("\\)", "\\\\)",
                                                   .data$modName))) %>%
            dplyr::mutate(modAA = gsub("\\[\\^", "n", .data$modAA))
    } else {
        fpAllModsDf <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Search parameters
    ## -------------------------------------------------------------------------
    if (!is.null(configDf)) {
        fpSearchEngine <- paste(gsub(
            ".+(MSFragger-.+).jar", "\\1",
            configDf$value[configDf$parameter ==
                               "fragpipe-config.bin-msfragger"]),
            collapse = ", ")
        fpFastaFiles <- paste(gsub("\\\\", "/",
                                   configDf$value[configDf$parameter ==
                                                      "database.db-path"]),
                              collapse = ", ")
        # fpContaminants <- "cRAP" ## cRAP in FASTA format can be obtained from
        # the GPM FTP site, using the URL ftp://ftp.thegpm.org/fasta/cRAP.
    } else if (!is.null(logDfConf)) {
        fpSearchEngine <- paste(gsub(
            ".+(MSFragger-.+).jar", "\\1",
            logDfConf$value[logDfConf$parameter ==
                                "fragpipe-config.bin-msfragger"]),
            collapse = ", ")
        fpFastaFiles <- paste(gsub(
            "/:/", ":/",
            gsub("//", "/",
                 gsub("\\\\", "/",
                      logDfConf$value[logDfConf$parameter ==
                                          "database.db-path"]))),
            collapse = ", ")
        # fpContaminants <- "cRAP"
    } else {
        fpSearchEngine <- fpFastaFiles <- NULL
        # fpContaminants <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Masses
    ## -------------------------------------------------------------------------
    if (!is.null(logDfConf)) {
        fpPepLengthL <-
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.digest_min_length"]
        fpPepLengthU <-
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.digest_max_length"]
        fpPepMassLo <-
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.misc.fragger.digest-mass-lo"]
        fpPepMassHi <-
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.misc.fragger.digest-mass-hi"]
        fpPepSel <- paste(paste0("length: ", fpPepLengthL, "-", fpPepLengthU,
                                 " AA; mass: ", fpPepMassLo, "-", fpPepMassHi, " Da"),
                          collapse = ", ")

        fpPmassTolUnits <- ifelse(
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.precursor_mass_units"] == 0,
            "Da", "ppm")
        fpPmassTolL <- logDfConf$value[logDfConf$parameter ==
                                           "msfragger.precursor_mass_lower"]
        fpPmassTolU <- logDfConf$value[logDfConf$parameter ==
                                           "msfragger.precursor_mass_upper"]
        fpPmassTol <- paste0(paste(fpPmassTolL, fpPmassTolU, sep = "-"),
                             " [", fpPmassTolUnits, "]")

        fpFmassTolUnits <- ifelse(
            logDfConf$value[logDfConf$parameter ==
                                "msfragger.fragment_mass_units"] == 0,
            "Da", "ppm")
        fpFmassTolU <- logDfConf$value[logDfConf$parameter ==
                                           "msfragger.fragment_mass_tolerance"]
        fpFmassTol <- paste0(fpFmassTolU, " [", fpFmassTolUnits, "]")
        fpFmassTolO <- gsub("New fragment_mass_tolerance = ", "",
                            logDf$info[grep("New fragment_mass_tolerance = ",
                                            logDf$info)])
        fpMassTol <- paste(paste0("precursor:", fpPmassTol, "; fragment:",
                                  fpFmassTol, " (after optimization:",
                                  fpFmassTolO, ")"), collapse = ", ")
    } else {
        fpPepSel <- fpPmassTol <- fpMassTol <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Quantification parameters
    ## -------------------------------------------------------------------------
    if (!is.null(dataDf)) {
        ### for IonQuant
        fpQuantMethods <- paste(
            paste("IonQuant:",
                  dataDf$value[dataDf$parameter == "ionquant.run-ionquant"] ==
                      "true"),
            paste("Calculate MaxLFQ intensity:",
                  dataDf$value[dataDf$parameter == "ionquant.maxlfq"] == 1),
            paste("Normalization:",
                  dataDf$value[dataDf$parameter == "ionquant.normalization"] ==
                      1),
            paste("match-between runs (MBR):",
                  dataDf$value[dataDf$parameter == "ionquant.mbr"] == 1),
            paste("min. ions:",
                  dataDf$value[dataDf$parameter == "ionquant.minions"]),
            paste("Top N ions:",
                  dataDf$value[dataDf$parameter == "ionquant.tp"]),
            sep = ", ")
        fpEnzymes <- paste(paste0(
            dataDf$value[dataDf$parameter == "msfragger.search_enzyme_name_1"],
            paste0("[", paste(
                dataDf$value[dataDf$parameter ==
                                 "msfragger.search_enzyme_cut_1"],
                paste0(dataDf$value[dataDf$parameter ==
                                        "msfragger.search_enzyme_sense_1"],
                       "-terminal"),
                paste0(dataDf$value[dataDf$parameter ==
                                        "msfragger.allowed_missed_cleavage_2"],
                       " missed cleavages", "]"), sep = ", "))),
            collapse = ", ")

        if ("msfragger.search_enzyme_name_2" %in% dataDf$parameter &&
            dataDf$value[dataDf$parameter ==
                         "msfragger.search_enzyme_name_2"] != "null") {
            fpEnzymes <- paste(
                fpEnzymes,
                dataDf$value[dataDf$parameter ==
                                 "msfragger.search_enzyme_name_2"],
                paste0("[", paste(
                    dataDf$value[dataDf$parameter ==
                                     "msfragger.search_enzyme_cut_2"],
                    paste0(dataDf$value[dataDf$parameter ==
                                            "msfragger.search_enzyme_sense_2"],
                           "-terminal"),
                    paste0(
                        dataDf$value[dataDf$parameter ==
                                         "msfragger.allowed_missed_cleavage_2"],
                        " missed cleavages", "]"), sep = ", ")), sep = "; ")
        }

        fpVariableModifications <- paste(
            fpAllModsDf$modName[fpAllModsDf$modType == "variable"],
            collapse = ", ")
        fpFixedModifications <- paste(
            fpAllModsDf$modName[fpAllModsDf$modType == "fixed"],
            collapse = ", ")
    } else {
        fpQuantMethods <- fpEnzymes <- fpVariableModifications <-
            fpFixedModifications <- NULL
    }

    ## -------------------------------------------------------------------------
    ## Database decoy tag
    ## -------------------------------------------------------------------------
    if (!is.null(dataDf)) {
        fpDecoyTag <- paste(dataDf$value[dataDf$parameter ==
                                             "database.decoy-tag"],
                            collapse = ";")
    } else {
        fpDecoyTag <- "rev_"
    }

    l <- list("FragPipe version" = fpVersion,
              "FragPipe parameter file" = fpConfigFile,
              "FragPipe log file" = fpLogFile,
              "Search engine" = fpSearchEngine,
              "Raw file location" = fpRawDirs,
              "Raw files" = fpRawFiles,
              "Sample names" = fpExperiments,
              "Databases" = fpFastaFiles,
              # "Contaminants" = fpContaminants,
              "Peptides (ranges)" = fpPepSel,
              "Mass error tolerances" = fpMassTol,
              "Quantification settings (LFQ)" = fpQuantMethods,
              "Enzymes" = fpEnzymes,
              "Variable modifications" = fpVariableModifications,
              "Fixed modifications" = fpFixedModifications,
              "Database decoy tag" = fpDecoyTag)
    l[vapply(l, length, 0) > 0]
}
