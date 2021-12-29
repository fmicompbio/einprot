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
    nodes <- nodes[xml2::xml_attr(nodes, "Name") == "ContaminantsDatabase"]
    # nodes <- nodes[xml2::xml_attr(nodes, "Name") %in%
    #                    c("ContaminantsDatabase",
    #                      paste0("AdditionalDatabase", seq_len(maxAdditional)))]
    setdiff(xml_attr(nodes, "DisplayValue"), "")
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
    list(search_engine = search_engine,
         fasta_db = setdiff(fastadb, ""),
         enzymes = enzymes)
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
getPSMValidationInfoFrompdAnalysis <- function(pdAnalysisFile) {
    pda <- xml2::read_xml(pdAnalysisFile)
    nodes <- xml2::xml_find_all(pda, ".//WorkflowNode")
    nodes <- nodes[xml2::xml_attr(nodes, "Category") == "PSM Validation"]
    xml2::xml_attr(nodes, "FriendlyName")
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

    list(quant_mode = quant_mode,
         peptides_to_use = peptides_to_use,
         abundance_based_on = abundance_based_on,
         quan_value_corr = quan_value_corr,
         co_isolation_thr = co_isolation_thr,
         ave_reporter_sn_thr = ave_reporter_sn_thr,
         sps_mm_pct_thr = sps_mm_pct_thr,
         norm_mode = norm_mode,
         imputation_mode = imputation_mode)
}

