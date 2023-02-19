#' Help functions for extracting information from a pdAnalysis file
#'
#' @param pdAnalysisFile Path to a .pdAnalysis file from Proteome Discoverer
#'
#' @author Charlotte Soneson
#'
#' @name querypdAnalysis
#'
#' @return The extracted information from the .pdAnalysis file
#'
#' @importFrom xml2 read_xml xml_find_all xml_attr xml_text
NULL

#' @rdname querypdAnalysis
#' @export
getContaminantsDatabaseFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") ==
                       "ContaminantsDetectorNode"]
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    ## Contaminant database
    nodescont <- nodes[xml2::xml_attr(nodes, "Name") == "ContaminantsDatabase"]
    contdb <- setdiff(xml2::xml_attr(nodescont, "DisplayValue"), "")

    ## Additional databases
    nodes <- nodes[grep("AdditionalDatabase",
                        xml2::xml_attr(nodes, "Name"))]
    nodes <- nodes[xml2::xml_attr(nodes, "IsValueSet") == "True"]
    adbs <- xml2::xml_attr(nodes, "DisplayValue")
    protMarkers <- c(if (length(contdb) > 0) "Contaminants", adbs)

    list(contaminantDb = contdb,
         protMarkers = protMarkers)
}

#' @rdname querypdAnalysis
#' @export
getSearchParametersFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") ==
                       "IseNode"]
    search_engine <- xml2::xml_attr(nodes, "FriendlyName")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    fastadb <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "FastaDatabase"],
        "DisplayValue")
    enzymes <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "Enzyme"],
        "DisplayValue")

    ## Dynamic modifications
    nodesdyn <- nodes[grep("Dynamic Modifications", xml2::xml_attr(nodes, "Category"))]
    maxEqualDynamicModificationsPerPeptide <-
        xml2::xml_attr(nodesdyn[xml2::xml_attr(nodesdyn, "Name") == "MaxEqualModificationsPerPeptide"],
                       "DisplayValue")
    nodesdyn <- nodesdyn[grep("DynamicModification", xml2::xml_attr(nodesdyn, "Name"))]
    nodesdyn <- nodesdyn[xml2::xml_attr(nodesdyn, "IsValueSet") == "True"]
    dynamicModifications <- xml2::xml_attr(nodesdyn, "DisplayValue")

    ## Static modifications
    nodesstat <- nodes[grep("Static Modifications", xml2::xml_attr(nodes, "Category"))]
    nodesstat <- nodesstat[grep("StaticModification", xml2::xml_attr(nodesstat, "Name"))]
    nodesstat <- nodesstat[xml2::xml_attr(nodesstat, "IsValueSet") == "True"]
    staticModifications <- xml2::xml_attr(nodesstat, "DisplayValue")

    list(search_engine = unique(search_engine),
         fasta_db = unique(setdiff(fastadb, "")),
         enzymes = unique(enzymes),
         dynamicModifications = unique(dynamicModifications),
         staticModifications = unique(staticModifications))
}

#' @rdname querypdAnalysis
#' @export
getMaxMissedCleavagesFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") ==
                       "IseNode"]
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    nodes <- nodes[xml2::xml_attr(nodes, "Name") == "MaxMissedCleavages"]
    max_missed_cleavages <- xml2::xml_attr(nodes, "DisplayValue")
    unique(max_missed_cleavages)
}

#' @rdname querypdAnalysis
#' @export
getQuantOrderFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") ==
                       "ReporterIonQuantifierNode"]
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    nodes <- nodes[xml2::xml_attr(nodes, "Name") == "MSOrderFilter"]
    ms_order <- xml2::xml_attr(nodes, "DisplayValue")
    unique(ms_order)
}

#' @rdname querypdAnalysis
#' @export
getTemplateNamesFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowInfo")
    templates <- list()
    for (n in nodes) {
        w <- xml_find_all(n, ".//WorkflowType")
        w <- xml2::xml_text(w)
        d <- xml_find_all(n, ".//WorkflowDefinition")
        d <- xml_find_all(d, ".//Workflow")
        d <- xml_attr(d, "TemplateName")
        templates[[w]] <- d
    }
    templates
}

#' @rdname querypdAnalysis
#' @export
getValidationInfoFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") ==
                       "ReportPeptideValidatorNode"]
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    targetfdrstrictpsm <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "TargetFPRHigh"],
        "DisplayValue")
    targetfdrrelaxedpsm <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "TargetFPRMiddle"],
        "DisplayValue")
    targetfdrstrictpept <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "TargetPeptideFPRHigh"],
        "DisplayValue")
    targetfdrrelaxedpept <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "TargetPeptideFPRMiddle"],
        "DisplayValue")

    validation_based_on <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "ValidationBasedOn"],
        "DisplayValue")
    validation_purpose_details <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "ValidationBasedOn"],
        "PurposeDetails")

    target_decoy_selection <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "TargetDecoySelection"],
        "IsValueSet")

    list(targetFDRstrictPSM = targetfdrstrictpsm,
         targetFDRrelaxedPSM = targetfdrrelaxedpsm,
         targetFDRstrictPeptide = targetfdrstrictpept,
         targetFDRrelaxedPeptide = targetfdrrelaxedpept,
         validationMethod = strsplit(validation_purpose_details, "\\.")[[1]][2],
         validationBasedOn = paste(c(if (target_decoy_selection == "True") "Target/Decoy",
                                     validation_based_on), collapse = ", "))
}

#' @rdname querypdAnalysis
#' @export
getPSMValidationInfoFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "Category") == "PSM Validation"]
    xml2::xml_attr(nodes, "FriendlyName")
}

#' @rdname querypdAnalysis
#' @export
getCalibrationFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "ProcessingNodeName") == "SpectrumRecalibrationNode"]
    length(nodes) > 0  ## if the node is there, return TRUE
}

#' @rdname querypdAnalysis
#' @export
getQuantInfoFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "Category") == "Quantification"]

    quant_mode <- xml2::xml_attr(nodes, "FriendlyName")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameters")
    nodes <- xml2::xml_find_all(nodes, ".//ProcessingNodeParameter")
    peptides_to_use <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "PeptidesToUse"],
        "DisplayValue"
    )
    abundance_based_on <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "ReporterAbundanceBasedOn"],
        "DisplayValue"
    )
    quan_value_corr <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "ApplyQuanValueCorrections"],
        "DisplayValue"
    )
    co_isolation_thr <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "IsolationInterferenceThreshold"],
        "DisplayValue"
    )
    ave_reporter_sn_thr <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "AverageReporterSignalToNoiseThreshold"],
        "DisplayValue"
    )
    sps_mm_pct_thr <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "SpsMassMatchesPercentageThreshold"],
        "DisplayValue"
    )
    norm_mode <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "NormalizationMode"],
        "DisplayValue"
    )
    imputation_mode <- xml2::xml_attr(
        nodes[xml2::xml_attr(nodes, "Name") == "MissingValueImputationMode"],
        "DisplayValue"
    )

    list(quant_mode = unique(quant_mode),
         peptides_to_use = unique(peptides_to_use),
         abundance_based_on = unique(abundance_based_on),
         quan_value_corr = unique(quan_value_corr),
         co_isolation_thr = unique(co_isolation_thr),
         ave_reporter_sn_thr = unique(ave_reporter_sn_thr),
         sps_mm_pct_thr = unique(sps_mm_pct_thr),
         norm_mode = unique(norm_mode),
         imputation_mode = unique(imputation_mode))
}

