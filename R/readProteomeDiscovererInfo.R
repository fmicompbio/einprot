#' Read Proteome Discoverer metadata
#'
#' Read metadata from Proteome Discoverer files. If some of the files are
#' missing, return the information from the ones that are available.
#'
#' @param pdOutputFolder Character string pointing to the PD/TMT output folder.
#'     Should contain the files \code{pdResultName_InputFiles.txt} and
#'     \code{pdResultName_StudyInformation.txt}.
#' @param pdResultName Character string providing the base name for the
#'     files in the \code{pdOutputFolder}.
#' @param pdAnalysisFile Path to a .pdAnalysis file from Proteome Discoverer.
#'     Can be \code{NULL}.
#'
#' @seealso querypdAnalysis
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A \code{list} with settings read from the Proteome Discoverer files.
#'
#' @examples
#' pdi <- readProteomeDiscovererInfo(
#'     pdOutputFolder = system.file("extdata", "pdtmt_example",
#'                                  package = "einprot"),
#'     pdResultName = "Fig2_m23139_RTS_QC_varMods",
#'     pdAnalysisFile = system.file("extdata", "pdtmt_example",
#'                                  "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
#'                                  package = "einprot"))
#'
#' @importFrom utils read.delim
readProteomeDiscovererInfo <- function(pdOutputFolder, pdResultName,
                                       pdAnalysisFile) {
    .assertScalar(x = pdOutputFolder, type = "character")
    .assertScalar(x = pdResultName, type = "character")
    .assertScalar(x = pdAnalysisFile, type = "character", allowNULL = TRUE)

    reqFiles <- c(
        structure(file.path(pdOutputFolder, paste0(
            pdResultName, c("_InputFiles.txt", "_StudyInformation.txt"))),
            names = c("input", "studyinf")),
        pdAnalysisFile)
    msg <- !file.exists(reqFiles)
    if (any(msg)) {
        message("Missing files: ", paste(reqFiles[msg], collapse = ", "))
    }
    if (all(msg)) {
        ## All files missing - return empty list
        return(list())
    }

    ## Data from InputFiles file
    if (file.exists(reqFiles["input"])) {
        pd_InputFiles <- utils::read.delim(reqFiles["input"], sep = "\t")

        pd_version <- sub("Created with Discoverer version: ", "",
                          pd_InputFiles$Software.Revision[1])
        pd_instruments <- paste(setdiff(unique(pd_InputFiles$Instrument.Name[-1]), ""),
                                collapse = ", ")
        pd_raw_dirs <- ifelse(
            length(unique(gsub(
                "(.+)\\\\.*.raw", "\\1",
                pd_InputFiles$File.Name[grep(".raw", pd_InputFiles$File.Name)],
                perl = TRUE))) == 1,
            gsub("\\\\", "/", unique(gsub(
                "(.+\\\\).*.raw", "\\1",
                pd_InputFiles$File.Name[grep(".raw", pd_InputFiles$File.Name)],
                perl = TRUE))),
            "multiple")
        pd_raw_files <- paste(gsub(
            ".+\\\\(.+)", "\\1",
            pd_InputFiles$File.Name[grep(".raw", pd_InputFiles$File.Name)]),
            collapse = ", ")
        pd_search_results <- paste(gsub(
            ".+\\\\(.+)", "\\1",
            pd_InputFiles$File.Name[grep(".msf", pd_InputFiles$File.Name)]),
            collapse = ", ")
    } else {
        pd_version <- pd_instruments <- pd_raw_dirs <- pd_raw_files <-
            pd_search_results <- NULL
    }

    ## Data from StudyInformation file
    if (file.exists(reqFiles["studyinf"])) {
        pd_StudyInformation <-
            utils::read.delim(reqFiles["studyinf"], sep = "\t")

        pd_experiments <- paste(unique(pd_StudyInformation$Sample.Group),
                                collapse = ", ")
    } else {
        pd_experiments <- NULL
    }

    ## Data from pdAnalysis file
    if (!is.null(pdAnalysisFile) && file.exists(pdAnalysisFile)) {
        pd_db_info <- getContaminantsDatabaseFrompdAnalysis(pdAnalysisFile)
        pd_contaminants <- paste(pd_db_info$contaminantDb, collapse = ", ")
        pd_ProtMarker <- paste(pd_db_info$protMarkers, collapse = ", ")

        pd_search_parameters <- getSearchParametersFrompdAnalysis(pdAnalysisFile)
        pd_fasta_files <- paste(pd_search_parameters$fasta_db, collapse = ", ")
        pd_search_engine <- paste(pd_search_parameters$search_engine, collapse = ", ")
        pd_enzymes <- paste(pd_search_parameters$enzymes, collapse = ", ")
        pd_fixed_modifications <-
            paste(pd_search_parameters$staticModifications, collapse = ", ")
        pd_variable_modifications <-
            paste(pd_search_parameters$dynamicModifications, collapse = ", ")

        pd_quant_info <- getQuantInfoFrompdAnalysis(pdAnalysisFile)
        pd_quant_mode <- paste(pd_quant_info$quant_mode, collapse = ", ")
        pd_abundance <- paste(pd_quant_info$abundance_based_on, collapse = ", ")
        pd_quanvaluecorrection <- paste(pd_quant_info$quan_value_corr, collapse = ", ")
        pd_peptides_for_quantification <- paste(pd_quant_info$peptides_to_use, collapse = ", ")
        pd_CoIsolationThr <- paste(pd_quant_info$co_isolation_thr, collapse = ", ")
        pd_avReporterSNThr <- paste(pd_quant_info$ave_reporter_sn_thr, collapse = ", ")
        pd_SPSMMpct <- paste(pd_quant_info$sps_mm_pct_thr, collapse = ", ")
        pd_normMode <- paste(pd_quant_info$norm_mode, collapse = ", ")
        pd_ImputationMode <- paste(pd_quant_info$imputation_mode, collapse = ", ")

        pd_max_missed_cleavages <- paste(getMaxMissedCleavagesFrompdAnalysis(pdAnalysisFile),
                                         collapse = ", ")

        pd_validation <- getValidationInfoFrompdAnalysis(pdAnalysisFile)
        pd_confidence_threshold <-
            paste0("strict: ", paste(pd_validation$targetFDRstrictPSM, collapse = ", "),
                   ", relaxed: ", paste(pd_validation$targetFDRrelaxedPSM, collapse = ", "))
        pd_validation_based_on <- pd_validation$validationBasedOn
        pd_validation_method <- paste(pd_validation$validationMethod, collapse = ", ")

        pd_templates <- getTemplateNamesFrompdAnalysis(pdAnalysisFile)
        pd_PSM_validation <- paste(getPSMValidationInfoFrompdAnalysis(pdAnalysisFile),
                                   collapse = ", ")
        pd_calibration <- getCalibrationFrompdAnalysis(pdAnalysisFile)
        pd_quant_order <- paste(getQuantOrderFrompdAnalysis(pdAnalysisFile), collapse = ", ")

        pd_quant_methods <- paste0(
            "Peptides used:", pd_peptides_for_quantification,
            ", quan. method: ", pd_quant_mode,
            ", quan. MS order: ", pd_quant_order,
            ", abundance type: ", pd_abundance,
            ", quan. correction: ", pd_quanvaluecorrection,
            ", MS1 co-isolation threshold: ", pd_CoIsolationThr,
            ", av. reporter SN threshold: ", pd_avReporterSNThr,
            ", PSP mass matches [%]: ", pd_SPSMMpct,
            ", norm. method: ", pd_normMode,
            ", PD imputation: ", pd_ImputationMode
        )
    } else {
        pd_contaminants <- pd_ProtMarker <- pd_fasta_files <-
            pd_search_engine <- pd_enzymes <- pd_fixed_modifications <-
            pd_variable_modifications <- pd_quant_mode <- pd_abundance <-
            pd_quanvaluecorrection <- pd_peptides_for_quantification <-
            pd_CoIsolationThr <- pd_avReporterSNThr <- pd_SPSMMpct <-
            pd_normMode <- pd_ImputationMode <- pd_confidence_threshold <-
            pd_validation_method <- pd_validation_based_on <- pd_templates <-
            pd_PSM_validation <- pd_calibration <- pd_quant_order <-
            pd_quant_methods <- pd_max_missed_cleavages <- NULL
    }

    ## Return value
    L <- list("PD version" = pd_version,
              "PD output folder" = pdOutputFolder,
              "PD result name" = pdResultName,
              "PD analysis file" = pdAnalysisFile,
              "PD Processing WF" = pd_templates$Processing,
              "PD Consensus WF" = pd_templates$Consensus,
              "Search engine" = pd_search_engine,
              "Instruments" = pd_instruments,
              "Raw file location" = pd_raw_dirs,
              "Raw files" = pd_raw_files,
              "Sample names" = pd_experiments,
              "Databases" = pd_fasta_files,
              "Contaminants" = pd_contaminants,
              "Quantification settings (LFQ)" = pd_quant_methods,
              "Enzymes" = pd_enzymes,
              "Variable modifications" = pd_variable_modifications,
              "Fixed modifications" = pd_fixed_modifications,
              "Validation method" = pd_validation_method,
              "Validation based on" = pd_validation_based_on,
              "Confidence thresholds" = pd_confidence_threshold,
              "Max missed cleavages" = pd_max_missed_cleavages
    )
    L[!vapply(L, is.null, TRUE)]
}
