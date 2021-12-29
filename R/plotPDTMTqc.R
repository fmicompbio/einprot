# Help function to rename modifications
.renameModifications <- function(mod, isTMT = TRUE, labelTMT = "") {
    ## Rename modification  names, to make them shorter for barplotting
    mod <- gsub("DSSO:H2O", "DSSO-OH", mod) ## rename mods names that contain numbers
    mod <- gsub("DSSO:Tris", "DSSO-Tris", mod) ## rename mods names that contain numbers
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
#' @return Invisibly, a list of ggplot objects.
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
#'     theme coord_cartesian geom_hline element_text element_blank
#' @importFrom cowplot plot_grid
#'
plotPDTMTqc <- function(pdOutputFolder, pdResultName, masterOnly = TRUE,
                        poiText = "", doPlot = TRUE, textSize = 4) {
    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .assertScalar(pdOutputFolder, type = "character")
    .assertScalar(pdResultName, type = "character")
    .assertScalar(masterOnly, type = "logical")
    .assertScalar(poiText, type = "character")
    .assertScalar(doPlot, type = "logical")
    .assertScalar(textSize, type = "numeric")

    ## --------------------------------------------------------------------- ##
    ## Read data
    ## --------------------------------------------------------------------- ##
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
        file = file.path(pdOutputFolder, paste0(pdResultName, "_QuanSpectra.txt")))

    ## --------------------------------------------------------------------- ##
    ## Filter protein data
    ## --------------------------------------------------------------------- ##
    ## Remove proteins with a single peptide
    proteinssub <- proteins %>%
        dplyr::filter(Number.of.Peptides > 1)

    if (masterOnly) {
        proteinssub <- proteinssub %>%
            dplyr::filter(!(Master %in% c("IsMasterProteinCandidateRejected",
                                          "IsMasterProteinRejected", "None")))
    }

    ## --------------------------------------------------------------------- ##
    ## Complete PSM file if necessary
    ## --------------------------------------------------------------------- ##
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

    ## --------------------------------------------------------------------- ##
    ## Get values to use for later
    ## --------------------------------------------------------------------- ##
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

    ## --------------------------------------------------------------------- ##
    ## Generate plots
    ## --------------------------------------------------------------------- ##
    plots <- list()

    ## Peptide lengths
    df <- data.frame(pepLength = nchar(psms$Sequence))
    mpl <- round(mean(df$pepLength), 1)
    plots[[1]] <- ggplot2::ggplot(df, ggplot2::aes(x = pepLength)) +
        ggplot2::geom_histogram(bins = 30, fill = "lightgrey",
                                color = "grey20") +
        ggplot2::geom_vline(xintercept = mpl,
                            linetype = 2, color = "red") +
        ggplot2::labs(x = "Peptide length (amino acids)", y = "Frequency") +
        ggplot2::annotate("text", mpl, Inf, hjust = -0.1, vjust = 1,
                          label = paste0("average:\n", mpl, "AAs"),
                          color = "red", size = textSize) +
        ggplot2::theme_minimal()

    ## Missed cleavages
    plots[[2]] <- ggplot2::ggplot(
        psms, ggplot2::aes(x = Number.of.Missed.Cleavages,
                           label = scales::percent(prop.table(stat(count)),
                                                   accuracy = 0.1))) +
        ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
        ggplot2::geom_text(
            stat = "count", size = textSize, color = "red",
            angle = 90, y = max(table(psms$Number.of.Missed.Cleavages))/2) +
        ggplot2::labs(x = "PSMs - # of missed cleavages", y = "") +
        ggplot2::theme_minimal()

    ## Retention time distribution
    plots[[3]] <- ggplot2::ggplot(
        psms, ggplot2::aes(x = RT.in.min, y = log10(Intensity))) +
        ggplot2::stat_density2d(ggplot2::aes(fill = ..density..^0.25),
                                geom = "tile", contour = FALSE, n = 200) +
        ggplot2::scale_fill_continuous(low = "white", high = "darkblue") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Retention time (min)", y = "log10(MS1 Intensity)")

    ## Charge z distribution
    plots[[4]] <- ggplot2::ggplot(
        psms, ggplot2::aes(x = Charge,
                           label = scales::percent(prop.table(stat(count)),
                                                   accuracy = 0.1))) +
        ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
        ggplot2::geom_text(stat = "count", size = textSize, color = "red",
                           angle = 90, y = max(table(psms$Charge))/2) +
        ggplot2::labs(x = "PSMs - charge (z)", y = "") +
        ggplot2::theme_minimal()

    ## ppm vs score
    massTol <- 20
    maxScore <- max(psms[[scoreName]])
    ppmName <- "Delta.M.in.ppm"
    ppm <- stats::median(psms[[ppmName]])
    plots[[5]] <- ggplot2::ggplot(
        psms, ggplot2::aes(x = .data[[ppmName]], y = .data[[scoreName]])) +
        ggplot2::stat_density2d(
            ggplot2::aes(fill = ..density..^0.25), geom = "tile",
            contour = FALSE, n = 200) +
        ggplot2::scale_fill_continuous(low = "white", high = "darkblue") +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none") +
        ggplot2::labs(x = "Mass deviation [ppm]",
                      y = paste0("Score (", scoreName, ")")) +
        ggplot2::coord_cartesian(xlim = c(-massTol, massTol)) +
        ggplot2::geom_vline(xintercept = ppm, linetype = 2, color = "red") +
        ggplot2::geom_hline(yintercept = minScore, linetype = 2, color = "red") +
        ggplot2::annotate("text", x = ppm + 5, y = maxScore,
                          label = paste(ppm, "ppm"),
                          size = textSize, color = "red") +
        ggplot2::annotate("text", x = 0.9 * massTol,
                          y = minScore + (maxScore - minScore) / 6,
                          label = "Q1(score)", size = textSize,
                          color = "red", angle = 90)

    ## Ion injection times
    df <- data.frame(injectTime = psms$Ion.Inject.Time.in.ms)
    plots[[6]] <- ggplot2::ggplot(df, ggplot2::aes(x = injectTime)) +
        ggplot2::geom_histogram(bins = 30, fill = "lightgrey",
                                color = "grey20") +
        ggplot2::labs(x = "Ion Inject Time [ms]", y = "Frequency") +
        ggplot2::theme_minimal()

    # ### 4. Ion injection times
    # ## Ion Inject Time [ms]
    # inject.t.col <- grep("INJECT", toupper(colnames(PSMs.r)))
    # ## histogram of inj times
    # h <- hist(PSMs.r[,inject.t.col], xlab = "Ion Inject Time [ms]", main = "")
    # rect(h$mids[length(h$mids)]-h$breaks[2]/2, 0, h$mids[length(h$mids)]+h$breaks[2]/2, h$counts[length(h$mids)], col = "red")
    # text(x= h$mids[length(h$mids)], y = 0.9*max(h$counts), paste(round(100*h$counts[length(h$mids)]/nrow(PSMs.r), 1), "%"), srt = 90, cex = 0.8)

    ## Summary mods and yields
    getTotal <- function(mod) {
        if (mod == "n-") nrow(psms)
        else if (mod == "NT") sum(grepl("\\[-\\]\\.", psms$Annotated.Sequence))
        else sum(stringr::str_count(toupper(as.character(psms$Sequence)), mod))
    }
    df <- data.frame(mod = mods) %>%
        dplyr::group_by(mod) %>% dplyr::tally() %>%
        dplyr::mutate(aaCount = vapply(gsub("(.+)\\(.+", "\\1", mod),
                                       getTotal, FUN.VALUE = 1))
    plots[[7]] <- ggplot2::ggplot(
        df, ggplot2::aes(x = mod, y = n,
                         label = scales::percent(n/aaCount,
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

    ## Number of proteins, peptide groups, PSMs
    nQuan <- sum(quanspec$Quan.Info == "")
    df <- data.frame(cls = c("Master proteins", "Proteins", "Peptide groups", "PSMs",
                             "MSMS quantified"),
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
    plots[[8]] <- ggplot2::ggplot(df, ggplot2::aes(x = cls, y = n)) +
        ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                          color = "grey20") +
        ggplot2::geom_text(aes(label = lab1), y = max(df$n)/4,
                           size = textSize, angle = 90, color = "red") +
        ggplot2::geom_text(aes(label = lab2), y = 0.6 * max(df$n),
                           size = textSize, angle = 90, color = "red") +
        ggplot2::labs(x = "", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust = 1,
                                                           vjust = 0.5))

    ## Proteins per DB
    df <- psms %>%
        dplyr::group_by(Marked.as, Contaminant) %>%
        dplyr::tally()
    dfplot <- dplyr::bind_rows(
        df %>% dplyr::group_by(Marked.as) %>%
            dplyr::summarize(n = sum(n)) %>%
            dplyr::rename(Label = Marked.as) %>%
            dplyr::select(Label, n),
        df %>% dplyr::filter(Contaminant == "True") %>%
            dplyr::group_by(Contaminant) %>%
            dplyr::summarize(n = sum(n)) %>%
            dplyr::mutate(Label = "CON") %>%
            dplyr::select(Label, n)
    ) %>% dplyr::mutate(frac = n/sum(df$n))
    plots[[9]] <- ggplot2::ggplot(
        dfplot, ggplot2::aes(x = Label, y = n,
                             label = scales::percent(frac,
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

    ## Master proteins
    df <- proteins %>%
        dplyr::group_by(Master) %>%
        dplyr::tally()
    plots[[10]] <- ggplot2::ggplot(df, ggplot2::aes(x = Master, y = n)) +
        ggplot2::geom_bar(stat = "identity", fill = "lightgrey",
                          color = "grey20") +
        ggplot2::geom_text(ggplot2::aes(label = n), size = textSize,
                           angle = 90, y = 0.1 * max(df$n),
                           color = "red") +
        ggplot2::geom_text(ggplot2::aes(label = Master), size = textSize,
                           angle = 90, y = 0.6 * max(df$n), color = "red") +
        ggplot2::labs(x = "Master proteins", y = "") +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_blank())

    ## Sequence coverage for top proteins
    df <- proteinssub %>%
        dplyr::arrange(dplyr::desc(Coverage.in.Percent)) %>%
        dplyr::slice_head(n = 50) %>%
        dplyr::mutate(idx = seq_along(Coverage.in.Percent))
    mcov <- round(mean(proteinssub$Coverage.in.Percent), 1)
    plots[[11]] <- ggplot2::ggplot(
        df, ggplot2::aes(x = idx, y = Coverage.in.Percent)) +
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

    ## POI - Proteins of interest
    df <- data.frame(poiFound = grepl(poiText, psms$Protein.Descriptions))
    plots[[12]] <- ggplot2::ggplot(
        df, ggplot2::aes(x = poiFound,
                         label = scales::percent(prop.table(stat(count)),
                                                 accuracy = 0.1))) +
        ggplot2::geom_bar(fill = "lightgrey", color = "grey20") +
        ggplot2::geom_text(stat = "count", size = textSize, color = "red",
                           angle = 90, y = nrow(df)/2) +
        ggplot2::labs(x = paste0("PSMs - ", poiText), y = "") +
        ggplot2::theme_minimal()

    if (doPlot) {
        print(cowplot::plot_grid(plotlist = plots, nrow = 3,
                                 align = "vh", axis = "tl"))
    }
    return(invisible(plots))
}
