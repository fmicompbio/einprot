#' Read FragPipe config/log files and extract information
#'
#' @param fragpipeDir Character scalar, the path to a FragPipe output
#'     directory. Should contain files <fragpipeDir>/fragpipe_*.config and
#'     <fragpipeDir>/log_*.txt (if not, the corresponding fields will
#'     be empty in the output).
#'
#' @author Charlotte Soneson, Jan Seebacher
#' @export
#'
#' @return A list with extracted information about the FragPipe run.
#'
readFragPipeInfo <- function(fragpipeDir) {
    .assertScalar(x = fragpipeDir, type = "character")

    ## ---------------------------------------------------------------------- ##
    ## Get file names from FragPipe directory
    ## ---------------------------------------------------------------------- ##
    fpConfigFile <- list.files(fragpipeDir, pattern = "^fragpipe.+.config$",
                               full.names = TRUE)
    stopifnot(length(fpConfigFile) <= 1)
    fpLogFile <- list.files(fragpipeDir, pattern = "^log_.+.txt$",
                            full.names = TRUE)
    stopifnot(length(fpLogFile) <= 1)

    ## ---------------------------------------------------------------------- ##
    ## Read config file and log file
    ## ---------------------------------------------------------------------- ##
    if (length(fpConfigFile) == 1) {
        configDf <- read.delim(fpConfigFile, header = FALSE, sep = "=") %>%
            setNames(c("parameter", "value"))
    } else {
        configDf <- NULL
    }
    if (length(fpLogFile) == 1) {
        logDf <- read.delim(fpLogFile, blank.lines.skip = TRUE, header = FALSE) %>%
            setNames("info")
    } else {
        logDf <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Get info about raw files
    ## ---------------------------------------------------------------------- ##
    if (!is.null(logDf)) {
        fpRawTxt <- gsub("(.+)\\:.+", "\\1",
                         gsub("\\\\", "/",
                              unlist(logDf[grep("\\.raw: Scans = ", logDf$info), ])))
        fpRawDirs <- paste(unique(dirname(fpRawTxt)), collapse = ", ")
        fpRawFiles <- paste(basename(fpRawTxt), collapse = ", ")
    } else {
        fpRawTxt <- fpRawDirs <- fpRawFiles <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Fixed modifications
    ## ---------------------------------------------------------------------- ##
    if (!is.null(configDf)) {
        fpFixedMods <- configDf$value[which(configDf$parameter ==
                                                "msfragger.table.fix-mods")]
        fpFixedModsDf <- as.data.frame(unlist(strsplit(fpFixedMods, ";"))) %>%
            setNames("orig") %>%
            tidyr::separate(.data$orig,
                            into = c("modMass", "modSymbol", "logi", "num"),
                            sep = ",", remove = FALSE) %>%
            dplyr::mutate(modSymbol = gsub(" $", "", .data$modSymbol),
                          modMass = round(as.numeric(.data$modMass), 4),
                          modAA = vapply(strsplit(.data$modSymbol, " "), .subset, 1,
                                         FUN.VALUE = ""),
                          logi = as.logical(.data$logi),
                          modType = "fixed")
    } else {
        fpFixedModsDf <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Variable modifications
    ## ---------------------------------------------------------------------- ##
    if (!is.null(configDf)) {
        fpVarMods <- configDf$value[which(configDf$parameter ==
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

    ## ---------------------------------------------------------------------- ##
    ## Put modifications together
    ## ---------------------------------------------------------------------- ##
    if (!is.null(fpFixedModsDf) && !is.null(fpVarModsDf)) {
        fpAllModsDf <- dplyr::bind_rows(
            fpFixedModsDf %>% dplyr::filter(.data$modMass > 0),
            fpVarModsDf %>% dplyr::filter(.data$logi)
        ) %>%
            dplyr::mutate(modName = paste0(.data$modAA, "(", .data$modMass, ")")) %>%
            dplyr::mutate(modName = gsub("\\[\\^", "N-term", .data$modName)) %>%
            dplyr::mutate(modNameRegex = gsub("\\(", "\\\\(",  gsub("\\)", "\\\\)",
                                                                    .data$modName))) %>%
            dplyr::mutate(modAA = gsub("\\[\\^", "n", .data$modAA))
    } else {
        fpAllModsDf <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Version info
    ## ---------------------------------------------------------------------- ##
    if (!is.null(configDf)) {
        fpVersion <- gsub(".+\\((.+)\\).+", "\\1", configDf$parameter[1])
    } else {
        fpVersion <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Search parameters
    ## ---------------------------------------------------------------------- ##
    if (!is.null(configDf)) {
        fpSearchEngine <- gsub(
            ".+(MSFragger-.+).jar", "\\1",
            configDf$value[configDf$parameter == "fragpipe-config.bin-msfragger"])
        fpFastaFiles <- gsub("\\\\", "/",
                             configDf$value[configDf$parameter == "database.db-path"])
        fpContaminants <- "cRAP" ## cRAP in FASTA format can be obtained from the
        ## GPM FTP site, using the URL ftp://ftp.thegpm.org/fasta/cRAP.
    } else {
        fpSearchEngine <- fpFastaFiles <- fpContaminants <- NULL
    }

    if (!is.null(logDf)) {
        fpExperiments <- paste(
            gsub("  Experiment/Group: ", "",
                 grep("Experiment/Group: ", logDf$orig, value = TRUE)),
            collapse = ", ")
    } else {
        fpExperiments <- NULL
    }

    ## ---------------------------------------------------------------------- ##
    ## Quantification parameters
    ## ---------------------------------------------------------------------- ##
    if (!is.null(configDf)) {
        ### for IonQuant
        fpQuantMethods <- paste(
            paste("IonQuant:",
                  configDf$value[configDf$parameter == "ionquant.run-ionquant"] == "true"),
            paste("Calculate MaxLFQ intensity:",
                  configDf$value[configDf$parameter == "ionquant.maxlfq"] == 1),
            paste("Normalization:",
                  configDf$value[configDf$parameter == "ionquant.normalization"] == 1),
            paste("match-between runs (MBR):",
                  configDf$value[configDf$parameter == "ionquant.mbr"] == 1),
            paste("min. ions:",
                  configDf$value[configDf$parameter == "ionquant.minions"]), sep = ", ")
        fpEnzymes <- paste0(
            configDf$value[configDf$parameter == "msfragger.search_enzyme_name_1"],
            paste0("[", paste(
                configDf$value[configDf$parameter == "msfragger.search_enzyme_cut_1"],
                paste0(configDf$value[configDf$parameter ==
                                          "msfragger.search_enzyme_sense_1"], "-terminal"),
                paste0(configDf$value[configDf$parameter ==
                                          "msfragger.allowed_missed_cleavage_2"],
                       " missed cleavages", "]"), sep = ", ")))

        if (configDf$value[configDf$parameter == "msfragger.search_enzyme_name_2"] != "null") {
            fpEnzymes <- paste(
                fpEnzymes,
                configDf$value[configDf$parameter == "msfragger.search_enzyme_name_2"],
                paste0("[", paste(
                    configDf$value[configDf$parameter == "msfragger.search_enzyme_cut_2"],
                    paste0(configDf$value[configDf$parameter ==
                                              "msfragger.search_enzyme_sense_2"], "-terminal"),
                    paste0(configDf$value[configDf$parameter == "msfragger.allowed_missed_cleavage_2"],
                           " missed cleavages", "]"), sep = ", ")), sep = "; ")
        }

        fpVariableModifications <- paste(fpAllModsDf$modName[fpAllModsDf$modType == "variable"],
                                         collapse = ", ")
        fpFixedModifications <- paste(fpAllModsDf$modName[fpAllModsDf$modType == "fixed"],
                                      collapse = ", ")
    } else {
        fpQuantMethods <- fpEnzymes <- fpVariableModifications <-
            fpFixedModifications <- NULL
    }

    l <- list("FragPipe version" = fpVersion,
              "FragPipe parameter file" = fpConfigFile,
              "FragPipe log file" = fpLogFile,
              "Search engine" = fpSearchEngine,
              "Raw file location" = fpRawDirs,
              "Raw files" = fpRawFiles,
              "Sample names" = fpExperiments,
              "Databases" = fpFastaFiles,
              "Contaminants" = fpContaminants,
              "Quantification settings (LFQ)" = fpQuantMethods,
              "Enzymes" = fpEnzymes,
              "Variable modifications" = fpVariableModifications,
              "Fixed modifications" = fpFixedModifications)
    l[vapply(l, length, 0) > 0]
}
