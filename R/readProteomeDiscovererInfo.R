#' Read Proteome Discoverer metadata
#'
#' @param pdOutputFolder Character string pointing to the PD/TMT output folder.
#'     Should contain the files \code{pdResultName_InputFiles.txt},
#'     \code{pdResultName_StudyInformation.txt} and
#'     \code{pdResultName_Proteins.txt}.
#' @param pdResultName Character string providing the base name for the
#'     files in the \code{pdOutputFolder}.
#' @param pdAnalysisFile Path to a .pdAnalysis file from Proteome Discoverer
#'
#' @seealso querypdAnalysis
#'
#' @export
#' @author Charlotte Soneson
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
    .assertScalar(x = pdAnalysisFile, type = "character")
    reqFiles <- c(file.path(pdOutputFolder, paste0(
        pdResultName, c("_InputFiles.txt", "_StudyInformation.txt",
                        "_Proteins.txt"))), pdAnalysisFile)
    msg <- !file.exists(reqFiles)
    if (any(msg)) {
        stop("Missing files: ", paste(reqFiles[msg], collapse = ", "))
    }

    pd_InputFiles <- utils::read.delim(
        file.path(pdOutputFolder,
                  paste0(pdResultName,  "_InputFiles.txt")), sep = "\t")
    pd_StudyInformation <- utils::read.delim(
        file.path(pdOutputFolder,
                  paste0(pdResultName, "_StudyInformation.txt")), sep = "\t")

    pd_version <- sub("Created with Discoverer version: ", "",
                      pd_InputFiles$Software.Revision[1])
    pd_instruments <- unique(pd_InputFiles$Instrument.Name[-1])
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
    pd_experiments <- paste(unique(pd_StudyInformation$Sample.Group),
                            collapse = ", ")

    pd_db_info <- getContaminantsDatabaseFrompdAnalysis(pdAnalysisFile)
    pd_contaminants <- paste(pd_db_info$contaminantDb, collapse = ", ")
    pd_ProtMarker <- paste(pd_db_info$protMarkers, collapse = ", ")

    pd_search_parameters <- getSearchParametersFrompdAnalysis(pdAnalysisFile)
    pd_fasta_files <- pd_search_parameters$fasta_db
    pd_search_engine <- pd_search_parameters$search_engine
    pd_enzymes <- pd_search_parameters$enzymes
    pd_fixed_modifications <-
        paste(pd_search_parameters$staticModifications, collapse = ", ")
    pd_variable_modifications <-
        paste(pd_search_parameters$dynamicModifications, collapse = ", ")

    pd_quant_info <- getQuantInfoFrompdAnalysis(pdAnalysisFile)
    pd_quant_mode <- pd_quant_info$quant_mode
    pd_abundance <- pd_quant_info$abundance_based_on
    pd_quanvaluecorrection <- pd_quant_info$quan_value_corr
    pd_peptides_for_quantification <- pd_quant_info$peptides_to_use
    pd_CoIsolationThr <- pd_quant_info$co_isolation_thr
    pd_avReporterSNThr <- pd_quant_info$ave_reporter_sn_thr
    pd_SPSMMpct <- pd_quant_info$sps_mm_pct_thr
    pd_normMode <- pd_quant_info$norm_mode
    pd_ImputationMode <- pd_quant_info$imputation_mode

    pd_validation <- getValidationInfoFrompdAnalysis(pdAnalysisFile)
    pd_confidence_threshold <-
        paste0("strict: ", pd_validation$targetFDRstrictPSM,
               ", relaxed: ", pd_validation$targetFDRrelaxedPSM)
    pd_validation_based_on <- pd_validation$validationBasedOn
    pd_validation_method <- pd_validation$validationMethod

    pd_templates <- getTemplateNamesFrompdAnalysis(pdAnalysisFile)
    pd_PSM_validation <- getPSMValidationInfoFrompdAnalysis(pdAnalysisFile)
    pd_calibration <- getCalibrationFrompdAnalysis(pdAnalysisFile)
    pd_quant_order <- getQuantOrderFrompdAnalysis(pdAnalysisFile)

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

    list("PD version" = pd_version,
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
         "Confidence thresholds" = pd_confidence_threshold
    )
}
