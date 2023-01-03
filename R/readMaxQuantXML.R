#' Read MaxQuant (mqpar.xml) file and extract information
#'
#' @param mqParameterFile Character scalar, the path to a MaxQuant parameter
#'     file. Can be \code{NULL} (in this case, an empty list is returned).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A list with extracted information about the MaxQuant run.
#'
#' @examples
#' mq <- readMaxQuantXML(system.file("extdata", "mq_example",
#'                                   "1356_mqpar.xml",
#'                                   package = "einprot"))
#'
#' @importFrom XML xmlToList xmlParse
#'
readMaxQuantXML <- function(mqParameterFile) {
    .assertScalar(x = mqParameterFile, type = "character", allowNULL = TRUE)

    if (is.null(mqParameterFile)) {
        return(list())
    } else {
        if (!file.exists(mqParameterFile)) {
            stop(mqParameterFile, " doesn't exist.")
        }

        mq_pars <- XML::xmlToList(XML::xmlParse(mqParameterFile))

        mq_version <- mq_pars$maxQuantVersion
        mq_contaminants <- paste0("*/MaxQuant_", mq_version,
                                  "/MaxQuant/bin/conf/contaminants.fasta")
        mq_search_engine <- "Andromeda"
        mq_raw_dirs <- ifelse(
            length(unique(gsub("(.+)\\\\.*.raw", "\\1", mq_pars$filePaths,
                               perl = TRUE))) == 1,
            gsub("\\\\", "/", unique(gsub("(.+\\\\).*.raw", "\\1",
                                          mq_pars$filePaths, perl = TRUE))),
            "multiple")
        mq_raw_files <- ifelse(
            mq_raw_dirs == "multiple", paste(
                gsub("\\\\", "/", unlist(mq_pars$filePaths)), collapse = ", "),
            gsub(mq_raw_dirs, "", paste(
                gsub("\\\\", "/", unlist(mq_pars$filePaths)), collapse = ", ")))
        mq_experiments <- paste(unlist(mq_pars$experiments), collapse = ", ")
        mq_fasta_files <- paste(gsub(
            "\\\\", "/", unlist(mq_pars$fastaFiles)), collapse = ", ")
        mq_quant_mode <- mq_pars$quantMode
        mq_LFQ_min_ratio_counts <- mq_pars$minRatioCount
        mq_fast_LFQ <- mq_pars$parameterGroups$parameterGroup$fastLfq
        mq_MBR <- mq_pars$matchBetweenRuns
        mq_peptides_for_quantification <- mq_pars$minRazorPeptides
        mq_iBAQ <- mq_pars$ibaq
        mq_requantify <- mq_pars$parameterGroups$parameterGroup$reQuantify
        mq_quant_methods <- paste0(
            "LFQ min. ratio count: ", mq_LFQ_min_ratio_counts,
            ", fastLFQ: ", mq_fast_LFQ,
            ", match-between runs (MBR): ", mq_MBR,
            ", Intensity based absolute quantification (iBAQ):", mq_iBAQ)
        mq_enzymes <- paste(
            unlist(mq_pars$parameterGroups$parameterGroup$enzymes), collapse = ", "
        )
        mq_fixed_modifications <- paste(unlist(mq_pars$fixedModifications),
                                        collapse = ", ")
        mq_variable_modifications <- paste(
            unlist(mq_pars$parameterGroups$parameterGroup$variableModifications),
            collapse = ", "
        )

        return(list(
            "MaxQuant version" = mq_version,
            "Parameter file" = mqParameterFile,
            "Search engine" = mq_search_engine,
            "Raw file location" = mq_raw_dirs,
            "Raw files" = mq_raw_files,
            "Sample names" = mq_experiments,
            "Databases" = mq_fasta_files,
            "Contaminants" = mq_contaminants,
            "Quantification mode" = mq_quant_mode,
            "Quantification settings (LFQ)" = mq_quant_methods,
            "Min. razor peptides" = mq_peptides_for_quantification,
            "Requantify" = mq_requantify,
            "Enzymes" = mq_enzymes,
            "Variable modifications" = mq_variable_modifications,
            "Fixed modifications" = mq_fixed_modifications
        ))
    }
}
