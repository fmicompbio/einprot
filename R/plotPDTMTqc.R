# Help function to rename modifications
#' @keywords internal
#' @noRd
.renameModifications <- function(mod, isTMT = TRUE, labelTMT = "") {
    ## Rename modification  names, to make them shorter for barplotting
    ## rename mods names that contain numbers
    mod <- gsub("DSSO:H2O", "DSSO-OH", mod)
    ## rename mods names that contain numbers
    mod <- gsub("DSSO:Tris", "DSSO-Tris", mod)
    mod <- gsub("Delta\\:H\\(2\\)C\\(3\\)O\\(1\\)", "DSSO-C=C", mod)
    mod <- gsub("DSSO sulfenic acid", "DSSO-SOH", mod)
    mod <- gsub("fragment ", "", mod)
    mod <- gsub("Carbamidomethyl", "CAM", mod)
    mod <- gsub("Oxidation", "Ox", mod)
    mod <- gsub("N-Term\\(Prot\\)", "NT", mod) #protein N-term
    mod <- gsub("N-Term", "n-", mod) # peptide Nterm
    mod <- gsub("Met\\-loss", "\\-M", mod) #protein Nterm Met removal
    mod <- gsub("Acetyl", "ac", mod)
    mod <- gsub("Phospho", "phos", mod)
    mod <- gsub("Cys\\-alk\\ NHS\\-iST", "alk\\ NHS\\-iST", mod)

    ## Remove mod positions from mods, N-Term will be unaffected,
    ## but TMT6plex will become TMTplex
    mod <- gsub("[0-9]", "", mod) ## remove numbers from AA pos. "M6(Oxidation)

    ## fix TMT label name
    if (isTMT) {
        mod <- gsub("TMTplex", labelTMT, mod)
    }
    mod
}

#' Generate summary QC plot for TMT data quantified with Proteome Discoverer
#'
#' Based mainly on the information from the PSMs.txt file (but also using
#' Proteins.txt and PeptideGroups.txt), this function creates a summary
#' QC plot for a proteomics experiment quantified with Proteome Discoverer.
#'
#' @param pdOutputFolder Character scalar pointing to the output folder from
#'     Proteome Discoverer. The folder must contain the files
#'     \code{pdResultName_Proteins.txt}, \code{pdResultName_PeptideGroups.txt},
#'     \code{pdResultName_PSMs.txt}, \code{pdResultName_MSMSSpectrumInfo.txt}
#'     and \code{pdResultName_QuanSpectra.txt}.
#' @param pdResultName Character scalar giving the shared prefix for the
#'     Proteome Discoverer output files (see \code{pdOutputFolder}).
#' @param masterOnly Logical scalar, whether to only consider Master proteins.
#' @param poiText Character scalar, pattern of interest to search for in the
#'     protein description.
#' @param doPlot Logical scalar, whether to render the plot or not.
#' @param textSize Numeric scalar, the font size of the text shown in the
#'     plots.
#'
#' @returns Invisibly, a list of ggplot objects.
#'
#' @author Charlotte Soneson, Jan Seebacher
#'
#' @export
#'
#' @importFrom utils read.delim
#' @importFrom dplyr %>% filter group_by tally mutate bind_rows summarize
#'     rename select arrange slice_head desc
#' @importFrom stats quantile median
#' @importFrom stringr str_split
#' @importFrom rlang .data
#' @importFrom scales percent
#' @importFrom ggplot2 ggplot geom_histogram geom_vline aes labs annotate
#'     theme_minimal geom_bar geom_text stat_density2d scale_fill_continuous
#'     theme coord_cartesian geom_hline element_text element_blank layer_data
#'     geom_rect after_stat
#' @importFrom cowplot plot_grid
#'
#' @examples
#' plots <- plotPDTMTqc(
#'     pdOutputFolder = system.file("extdata", "pdtmt_example",
#'                                  package = "einprot"),
#'     pdResultName = "Fig2_m23139_RTS_QC_varMods",
#'     doPlot = TRUE)
#'
plotPDTMTqc <- function(pdOutputFolder, pdResultName, masterOnly = FALSE,
                        poiText = "", doPlot = TRUE, textSize = 4) {
    ## -------------------------------------------------------------------------
    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertScalar(x = pdOutputFolder, type = "character")
    .assertScalar(x = pdResultName, type = "character")
    .assertScalar(x = masterOnly, type = "logical")
    .assertScalar(x = poiText, type = "character")
    .assertScalar(x = doPlot, type = "logical")
    .assertScalar(x = textSize, type = "numeric")
    reqFiles <- file.path(pdOutputFolder, paste0(
        pdResultName, c("_Proteins.txt", "_PSMs.txt",
                        "_PeptideGroups.txt", "_MSMSSpectrumInfo.txt",
                        "_QuanSpectra.txt")))
    msg <- !file.exists(reqFiles)
    if (any(msg)) {
        stop("Missing files: ", paste(reqFiles[msg], collapse = ", "))
    }

    ## -------------------------------------------------------------------------
    ## Read data
    ## -------------------------------------------------------------------------
    ## Read protein, peptide and peptide-to-spectrum matches files
    proteins <- utils::read.delim(
        file = file.path(pdOutputFolder,
                         paste0(pdResultName, "_Proteins.txt")))
    psms <- utils::read.delim(
        file = file.path(pdOutputFolder,
                         paste0(pdResultName, "_PSMs.txt")))
    pepgroups <- utils::read.delim(
        file = file.path(pdOutputFolder,
                         paste0(pdResultName, "_PeptideGroups.txt")))

    ## Read MSMS data, to count MS2 spectra, to compute ID rate later
    nmsms <- nrow(utils::read.delim(
        file = file.path(pdOutputFolder, paste0(pdResultName,
                                                "_MSMSSpectrumInfo.txt"))))

    ## Load Quan spectra
    quanspec <- utils::read.delim(
        file = file.path(pdOutputFolder, paste0(pdResultName,
                                                "_QuanSpectra.txt")))

    ## -------------------------------------------------------------------------
    ## Filter protein data
    ## -------------------------------------------------------------------------
    ## Remove proteins with a single peptide
    proteinssub <- proteins %>%
        dplyr::filter(.data$Number.of.Peptides > 1)

    if (masterOnly) {
        proteinssub <- proteinssub %>%
            dplyr::filter(!(.data$Master %in%
                                c("IsMasterProteinCandidateRejected",
                                  "IsMasterProteinRejected", "None")))
    }

    ## -------------------------------------------------------------------------
    ## Complete PSM file if necessary
    ## -------------------------------------------------------------------------
    ## Check if PSM file has Sequence column, otherwise generate it
    if ("Sequence" %in% colnames(psms)) {
        psms$Sequence <- as.character(psms$Sequence)
    } else if (all(grepl("\\[", psms$Annotated.Sequence))) {
        psms$Sequence <- toupper(gsub("\\[.+\\]\\.(.+)\\.\\[.+\\]", "\\1",
                                      psms$Annotated.Sequence))
    } else {
        psms$Sequence <- toupper(gsub("(.+)\\..+", "",
                                      as.character(psms$Annotated.Sequence)))
    }

    ## -------------------------------------------------------------------------
    ## Get values to use for later
    ## -------------------------------------------------------------------------
    ## Get minimum search score to indicate in plot later
    scoreNames <- c("XCorr", "Ions.Score")  ## possible values of score column
    scoreName <- scoreNames[min(which(scoreNames %in% colnames(psms)))]
    minScore <- stats::quantile(psms[[scoreName]], probs = 0.25)

    ## Check modifications
    mods <- unlist(stringr::str_split(psms$Modifications, pattern = "; "))
    mods <- mods[nchar(mods) > 0]
    modsUnique <- unique(gsub(".+\\((.+)\\)", "\\1", mods))
    if (any(grepl("TMT", modsUnique))) {
        ## TMT data
        mods <- .renameModifications(mods, isTMT = TRUE,
                                     labelTMT = grep("TMT", modsUnique,
                                                     value = TRUE))
    } else {
        mods <- .renameModifications(mods, isTMT = FALSE, labelTMT = "")
    }

    ## -------------------------------------------------------------------------
    ## Generate plots
    ## -------------------------------------------------------------------------
    plots <- list()

    ## Peptide lengths
    if ("Sequence" %in% colnames(psms)) {
        df <- data.frame(pepLength = nchar(psms$Sequence))
        mpl <- round(mean(df$pepLength), 1)
        plots[[1]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$pepLength)) +
            ggplot2::geom_histogram(bins = 30, fill = "lightgrey",
                                    color = "grey20") +
            ggplot2::geom_vline(xintercept = mpl,
                                linetype = 2, color = "red") +
            ggplot2::labs(x = "Peptide length (amino acids)", y = "Frequency") +
            ggplot2::annotate("text", mpl, Inf, hjust = -0.1, vjust = 1,
                              label = paste0("average:\n", mpl, "AAs"),
                              color = "red", size = textSize) +
            ggplot2::theme_minimal()
    } else {
        plots[[1]] <- NULL
    }

    ## Missed cleavages
    if ("Number.of.Missed.Cleavages" %in% colnames(psms)) {
        plots[[2]] <- ggplot2::ggplot(
            psms, ggplot2::aes(x = .data$Number.of.Missed.Cleavages,
                               label = scales::percent(
                                   prop.table(ggplot2::after_stat(.data$count)),
                                   accuracy = 0.1))) +
            ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
            ggplot2::geom_text(
                stat = "count", size = textSize, color = "red",
                angle = 90, y = max(table(psms$Number.of.Missed.Cleavages))/2) +
            ggplot2::labs(x = "PSMs - # of missed cleavages", y = "") +
            ggplot2::theme_minimal()
    } else {
        plots[[2]] <- NULL
    }

    ## Retention time distribution
    if (all(c("RT.in.min", "Intensity") %in% colnames(psms))) {
        plots[[3]] <- ggplot2::ggplot(
            psms, ggplot2::aes(x = .data$RT.in.min,
                               y = log10(.data$Intensity))) +
            ggplot2::stat_density2d(
                ggplot2::aes(fill = ggplot2::after_stat(density)^0.25),
                geom = "tile", contour = FALSE, n = 200) +
            ggplot2::scale_fill_continuous(low = "white", high = "darkblue") +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = "Retention time (min)",
                          y = "log10(MS1 Intensity)")
    } else {
        plots[[3]] <- NULL
    }

    ## Charge z distribution
    if ("Charge" %in% colnames(psms)) {
        plots[[4]] <- ggplot2::ggplot(
            psms, ggplot2::aes(x = .data$Charge,
                               label = scales::percent(
                                   prop.table(ggplot2::after_stat(.data$count)),
                                   accuracy = 0.1))) +
            ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
            ggplot2::geom_text(stat = "count", size = textSize, color = "red",
                               angle = 90, y = max(table(psms$Charge))/2) +
            ggplot2::labs(x = "PSMs - charge (z)", y = "") +
            ggplot2::theme_minimal()
    } else {
        plots[[4]] <- NULL
    }

    ## ppm vs score
    ppmName <- "Delta.M.in.ppm"
    if (all(c(scoreName, ppmName) %in% colnames(psms))) {
        massTol <- 20
        maxScore <- max(psms[[scoreName]])
        ppm <- stats::median(psms[[ppmName]])
        plots[[5]] <- ggplot2::ggplot(
            psms, ggplot2::aes(x = .data[[ppmName]], y = .data[[scoreName]])) +
            ggplot2::stat_density2d(
                ggplot2::aes(fill = ggplot2::after_stat(density)^0.25),
                geom = "tile",
                contour = FALSE, n = 200) +
            ggplot2::scale_fill_continuous(low = "white", high = "darkblue") +
            ggplot2::theme_minimal() +
            ggplot2::theme(legend.position = "none") +
            ggplot2::labs(x = "Mass deviation [ppm]",
                          y = paste0("Score (", scoreName, ")")) +
            ggplot2::coord_cartesian(xlim = c(-massTol, massTol)) +
            ggplot2::geom_vline(xintercept = ppm, linetype = 2, color = "red") +
            ggplot2::geom_hline(yintercept = minScore, linetype = 2,
                                color = "red") +
            ggplot2::annotate("text", x = ppm + 5, y = maxScore,
                              label = paste(ppm, "ppm"),
                              size = textSize, color = "red") +
            ggplot2::annotate("text", x = 0.9 * massTol,
                              y = minScore + (maxScore - minScore) / 6,
                              label = "Q1(score)", size = textSize,
                              color = "red", angle = 90)
    } else {
        plots[[5]] <- NULL
    }

    ## Ion injection times
    if ("Ion.Inject.Time.in.ms" %in% colnames(psms)) {
        gg <- ggplot2::ggplot(psms,
                              ggplot2::aes(x = .data$Ion.Inject.Time.in.ms)) +
            ggplot2::geom_histogram(bins = 30, fill = "lightgrey",
                                    color = "grey20") +
            ggplot2::labs(x = "Ion Inject Time [ms]", y = "Frequency") +
            ggplot2::theme_minimal()
        plotdf <- ggplot2::layer_data(gg, i = 1L)
        nr <- nrow(plotdf)
        gg <- gg +
            ggplot2::geom_rect(xmin = plotdf$xmin[nr], xmax = plotdf$xmax[nr],
                               ymin = plotdf$ymin[nr], ymax = plotdf$ymax[nr],
                               fill = "red", alpha = 0.1) +
            ggplot2::annotate("text", x = plotdf$x[nr], y = 0.8 * max(plotdf$y),
                              label = scales::percent(
                                  plotdf$y[nr]/sum(plotdf$y),
                                  accuracy = 0.01),
                              angle = 90, color = "red", size = textSize)
        plots[[6]] <- gg
    } else {
        plots[[6]] <- NULL
    }

    ## Summary mods and yields
    if (all(c("Annotated.Sequence", "Sequence") %in% colnames(psms))) {
        getTotal <- function(mod) {
            if (mod == "n-") nrow(psms)
            else if (mod == "NT") sum(grepl("\\[-\\]\\.",
                                            psms$Annotated.Sequence))
            else sum(stringr::str_count(toupper(as.character(psms$Sequence)),
                                        mod))
        }
        df <- data.frame(mod = mods) %>%
            dplyr::group_by(.data$mod) %>% dplyr::tally() %>%
            dplyr::mutate(aaCount = vapply(gsub("(.+)\\(.+", "\\1", .data$mod),
                                           getTotal, FUN.VALUE = 1))
        plots[[7]] <- ggplot2::ggplot(
            df, ggplot2::aes(x = .data$mod, y = .data$n,
                             label = scales::percent(.data$n/.data$aaCount,
                                                     accuracy = 0.01))) +
            ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                              color = "grey20") +
            ggplot2::geom_text(y = max(df$n)/2, size = textSize,
                               angle = 90, color = "red") +
            ggplot2::labs(x = "", y = "PSMs") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                               hjust = 1,
                                                               vjust = 0.5))
    } else {
        plots[[7]] <- NULL
    }

    ## Number of proteins, peptide groups, PSMs
    nQuan <- sum(quanspec$Quan.Info == "")
    df <- data.frame(
        cls = c("Master proteins", "Proteins", "Peptide groups",
                "PSMs", "MSMS quantified"),
        n = c(sum(proteins$Master == "IsMasterProtein"), nrow(proteins),
              nrow(pepgroups), nrow(psms), nQuan),
        lab1 = c(sum(proteins$Master == "IsMasterProtein"), nrow(proteins),
                 nrow(pepgroups), nrow(psms), nQuan),
        lab2 = c("", "", "",
                 paste0(scales::percent(nrow(psms)/nmsms, accuracy = 0.01),
                        " ident."),
                 paste0(scales::percent(nQuan/nmsms, accuracy = 0.01),
                        " quant.")))
    df$cls <- factor(df$cls, levels = df$cls)
    plots[[8]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$cls,
                                                   y = .data$n)) +
        ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                          color = "grey20") +
        ggplot2::geom_text(aes(label = .data$lab1), y = 0.2 * max(df$n),
                           size = textSize, angle = 90, color = "red") +
        ggplot2::geom_text(aes(label = .data$lab2), y = 0.7 * max(df$n),
                           size = textSize, angle = 90, color = "red") +
        ggplot2::labs(x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust = 1,
                                                           vjust = 0.5))

    ## Proteins per DB
    if (all(c("Marked.as", "Contaminant") %in% colnames(psms))) {
        df <- psms %>%
            dplyr::group_by(.data$Marked.as, .data$Contaminant) %>%
            dplyr::tally()
        dfplot <- dplyr::bind_rows(
            df %>% dplyr::group_by(.data$Marked.as) %>%
                dplyr::summarize(n = sum(.data$n)) %>%
                dplyr::rename(Label = "Marked.as") %>%
                dplyr::select("Label", "n"),
            df %>% dplyr::filter(.data$Contaminant == "True") %>%
                dplyr::group_by(.data$Contaminant) %>%
                dplyr::summarize(n = sum(.data$n)) %>%
                dplyr::mutate(Label = "CON") %>%
                dplyr::select("Label", "n")
        ) %>% dplyr::mutate(frac = .data$n/sum(df$n))
        plots[[9]] <- ggplot2::ggplot(
            dfplot, ggplot2::aes(x = .data$Label, y = .data$n,
                                 label = scales::percent(.data$frac,
                                                         accuracy = 0.1))) +
            ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                              color = "grey20") +
            ggplot2::geom_text(size = textSize, angle = 90,
                               y = max(dfplot$n)/2, color = "red") +
            ggplot2::labs(x = "PSMs - # of proteins by DB", y = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                               hjust = 1,
                                                               vjust = 0.5))
    } else {
        plots[[9]] <- NULL
    }

    ## Master proteins
    if ("Master" %in% colnames(proteins)) {
        df <- proteins %>%
            dplyr::group_by(.data$Master) %>%
            dplyr::tally()
        plots[[10]] <- ggplot2::ggplot(df, ggplot2::aes(x = .data$Master,
                                                        y = .data$n)) +
            ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                              color = "grey20") +
            ggplot2::geom_text(ggplot2::aes(label = .data$n), size = textSize,
                               angle = 90, y = 0.1 * max(df$n),
                               color = "red") +
            ggplot2::geom_text(ggplot2::aes(label = .data$Master),
                               size = textSize,
                               angle = 90, y = 0.6 * max(df$n), color = "red") +
            ggplot2::labs(x = "Master proteins", y = "") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank())
    } else {
        plots[[10]] <- NULL
    }

    ## Sequence coverage for top proteins
    if ("Coverage.in.Percent" %in% colnames(proteinssub)) {
        df <- proteinssub %>%
            dplyr::arrange(dplyr::desc(.data$Coverage.in.Percent)) %>%
            dplyr::slice_head(n = 50) %>%
            dplyr::mutate(idx = seq_along(.data$Coverage.in.Percent))
        mcov <- round(mean(proteinssub$Coverage.in.Percent), 1)
        plots[[11]] <- ggplot2::ggplot(
            df, ggplot2::aes(x = .data$idx, y = .data$Coverage.in.Percent)) +
            ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                              color = "grey20") +
            ggplot2::geom_hline(yintercept = mcov,
                                color = "red", linetype = 2) +
            ggplot2::annotate("text", x = 25, y = mcov * 1.2,
                              label = paste0("mean: ", mcov),
                              color = "red", size = textSize) +
            ggplot2::labs(x = "Top 50 proteins",
                          y = "% Protein sequence coverage") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank())
    } else {
        plots[[11]] <- NULL
    }

    ## POI - Proteins of interest
    if ("Protein.Descriptions" %in% colnames(psms)) {
        df <- data.frame(poiFound = grepl(poiText, psms$Protein.Descriptions))
        plots[[12]] <- ggplot2::ggplot(
            df, ggplot2::aes(x = .data$poiFound,
                             label = scales::percent(
                                 prop.table(ggplot2::after_stat(.data$count)),
                                 accuracy = 0.1))) +
            ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
            ggplot2::geom_text(stat = "count", size = textSize, color = "red",
                               angle = 90, y = nrow(df)/2) +
            ggplot2::labs(x = paste0("PSMs - ", poiText), y = "") +
            ggplot2::theme_minimal()
    } else {
        plots[[12]] <- NULL
    }

    if (doPlot) {
        print(cowplot::plot_grid(plotlist = plots, nrow = 3,
                                 align = "vh", axis = "tl"))
    }
    return(invisible(plots))
}
