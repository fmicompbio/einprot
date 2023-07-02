#' Run PCA and generate plots
#'
#' Apply PCA to the assay defined by \code{assayName} and extract the top
#' \code{ncomponent} components. For each pair defined in \code{plotpairs},
#' generate a collection of plots.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assayName Character scalar defining the assay in \code{sce} to
#'     use for the PCA.
#' @param ncomponents Numeric scalar, the (maximal) number of components to
#'     extract. The actual number can be lower if the number of samples is
#'     too small.
#' @param ntop Number of features (with highest variance) to use to
#'     generate the PCA. Will be passed on to \code{scater::runPCA}.
#' @param plotpairs A list of numeric vectors of length 2, indicating which
#'     pairs of PCs to generate plots for.
#' @param maxNGroups Numeric scalar, the maximum number of groups to display
#'     in the legend in the scatter plot in the combined plot. If there are
#'     more than \code{maxNGroups} groups, the legend is suppressed.
#' @param maxTextWidthBarplot Numeric scalar giving the maximum allowed width
#'     for text labels in the bar plot of log-fold changes. If not \code{NULL},
#'     the size of the labels will be scaled down in an attempt to keep the
#'     labels inside the canvas. Typically set to half the width of the
#'     plot device (in inches).
#' @param colourBy Character scalar indicating the name of the column of
#'     \code{colData(sce)} that will be used to colour the points.
#' @param subset_row Vector specifying the subset of features to use for
#'     dimensionality reduction. Can be a character vector of row names, an
#'     integer vector of row indices or a logical vector. Will be passed to
#'     \code{scater::runPCA}.
#' @param scale Logical scalar indicating whether the values should be scaled
#'     before the PCA is applied. Will be passed to \code{scater::runPCA}.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A list with the following components:
#' \itemize{
#'  \item{sce}{the input sce, expanded with the calculated PCs, in addition the
#'  feature coefficients will be added to the `rowData`}
#'  \item{plotcoord}{a list of `ggplot` objects containing coordinate plots for
#'  the desired pairs of components}
#'  \item{plotcombined}{a list of `ggplot` objects containing combined
#'  coordinate, scree and coefficient plots for the desired pairs of components}
#'  \item{plotpairs}{a `ggpairs` plot with all extracted components}
#' }
#'
#' @importFrom scater runPCA
#' @importFrom scuttle makePerCellDF
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr bind_rows filter arrange mutate
#' @importFrom scales percent_format
#' @importFrom ggalt geom_encircle
#' @importFrom ggplot2 ggplot aes geom_point ggtitle coord_fixed labs geom_col
#'     scale_fill_manual scale_y_continuous scale_x_continuous theme coord_flip
#'     geom_text element_blank element_line scale_fill_gradient2 theme_bw
#' @importFrom cowplot theme_cowplot plot_grid
#' @importFrom BiocSingular ExactParam
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom GGally ggpairs
#' @importFrom utils head tail
#'
doPCA <- function(sce, assayName, ncomponents = 10, ntop = Inf,
                  plotpairs = list(c(1, 2)), maxNGroups = 10,
                  maxTextWidthBarplot = NULL, colourBy = "group",
                  subset_row = NULL, scale = FALSE) {

    ## -------------------------------------------------------------------------
    ## Check arguments
    ## -------------------------------------------------------------------------
    .assertVector(x = sce, type = "SingleCellExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = ncomponents, type = "numeric")
    .assertScalar(x = ntop, type = "numeric", rngIncl = c(1, Inf))
    .assertVector(x = plotpairs, type = "list")
    .assertScalar(x = maxNGroups, type = "numeric")
    .assertScalar(x = maxTextWidthBarplot, type = "numeric", allowNULL = TRUE)
    .assertScalar(x = colourBy, type = "character", allowNULL = TRUE,
                  validValues = colnames(SummarizedExperiment::colData(sce)))
    .assertScalar(x = scale, type = "logical")
    for (elm in plotpairs) {
        .assertVector(x = elm, type = "numeric", len = 2)
    }
    if (ncomponents > (ncol(sce) - 1)) {
        message("Not enough samples - only ", (ncol(sce) - 1),
                " components will be extracted")
        ncomponents <- min(ncomponents, ncol(sce) - 1)
    }
    if (any(unlist(plotpairs) > ncomponents)) {
        stop("'plotpairs' requests components that will not be extracted")
    }

    ## -------------------------------------------------------------------------
    ## Run PCA and add to sce
    ## -------------------------------------------------------------------------
    sce <- scater::runPCA(sce, exprs_values = assayName,
                          ncomponents = ncomponents, ntop = ntop,
                          BSPARAM = BiocSingular::ExactParam(),
                          name = paste0("PCA_", assayName),
                          subset_row = subset_row, scale = scale)

    ## -------------------------------------------------------------------------
    ## Add coefficients to rowData (replace any existing columns with the
    ## same names)
    ## -------------------------------------------------------------------------
    cf <- attr(SingleCellExperiment::reducedDim(sce, paste0("PCA_", assayName)),
               "rotation")
    colnames(cf) <- paste0("PCA_", assayName, "_", colnames(cf))
    stopifnot(rownames(cf) %in% rownames(sce))
    cf <- cf[match(rownames(sce), rownames(cf)), , drop = FALSE]
    rownames(cf) <- rownames(sce)
    for (nm in colnames(cf)) {
        if (nm %in% colnames(SummarizedExperiment::rowData(sce))) {
            SummarizedExperiment::rowData(sce)[, nm] <- cf[, nm]
        } else {
            SummarizedExperiment::rowData(sce) <- cbind(
                SummarizedExperiment::rowData(sce),
                cf[, nm, drop = FALSE]
            )
        }
    }

    ## -------------------------------------------------------------------------
    ## Sample plot
    ## -------------------------------------------------------------------------
    percvar <- attr(SingleCellExperiment::reducedDim(sce,
                                                     paste0("PCA_", assayName)),
                    "percentVar")
    plist <- list()
    for (pr in plotpairs) {
        pcadf <- data.frame(
            sampleLabel = colnames(sce),
            SingleCellExperiment::reducedDim(sce, paste0("PCA_",
                                                         assayName))[, pr])
        if (!is.null(colourBy)) {
            pcadf[[colourBy]] <- SummarizedExperiment::colData(sce)[[colourBy]]
        }
        p <- ggplot2::ggplot(pcadf,
                             ggplot2::aes(x = .data[[paste0("PC", pr[1])]],
                                          y = .data[[paste0("PC", pr[2])]],
                                          label = .data$sampleLabel))
        if (!is.null(colourBy)) {
            p <- p + ggplot2::aes(color = .data[[colourBy]])
        }
        p <- p + ggplot2::geom_point(alpha = 0.5, size = 4) +
            ggplot2::ggtitle(paste0("Assay: ", assayName)) +
            ggplot2::coord_fixed() +
            ggplot2::labs(
                x = paste0("PC", pr[1], " (", round(percvar[pr[1]], 2), "%)"),
                y = paste0("PC", pr[2], " (", round(percvar[pr[2]], 2), "%)")) +
            cowplot::theme_cowplot()
        plist[[paste0("PC", pr[1], "_", pr[2])]] <- p
    }

    ## -------------------------------------------------------------------------
    ## Combined plot
    ## -------------------------------------------------------------------------
    plistcomb <- list()
    for (pr in plotpairs) {
        ## Scree plot
        pscree <- ggplot2::ggplot(data.frame(comp = seq_along(percvar),
                                             percvar = percvar),
                                  ggplot2::aes(x = .data$comp, y = percvar,
                                               fill = .data$comp %in% pr)) +
            ggplot2::geom_col() +
            ggplot2::scale_fill_manual(values = c(`TRUE` = "red",
                                                  `FALSE` = "grey")) +
            ggplot2::labs(x = "Component", y = "Percent\nexplained\nvariance") +
            ggplot2::scale_y_continuous(
                labels = scales::percent_format(scale = 1)
            ) +
            ggplot2::scale_x_continuous(breaks = seq_along(percvar)) +
            cowplot::theme_cowplot() +
            ggplot2::theme(legend.position = "none")

        ## Coefficient plot
        pcoef <- cowplot::plot_grid(
            plotlist = lapply(pr, function(i) {
                pdt <- dplyr::bind_rows(
                    data.frame(coef = cf[, paste0("PCA_", assayName,
                                                  "_PC", i)]) %>%
                        dplyr::filter(!is.na(.data$coef) & .data$coef > 0) %>%
                        dplyr::arrange(dplyr::desc(.data$coef)) %>%
                        utils::head(n = 10) %>%
                        tibble::rownames_to_column("pid"),
                    data.frame(coef = cf[, paste0("PCA_", assayName,
                                                  "_PC", i)]) %>%
                        dplyr::filter(!is.na(.data$coef) & .data$coef < 0) %>%
                        dplyr::arrange(dplyr::desc(.data$coef)) %>%
                        utils::tail(n = 10) %>%
                        tibble::rownames_to_column("pid")
                )
                if (!is.null(maxTextWidthBarplot)) {
                    ## graphics::strwidth or systemfonts::string_width
                    ## can be used to get the width of a string,
                    ## but create a new plot. We make a crude approximation by
                    ## the number of characters in the word, and use a
                    ## manually estimated 'average' conversion factor to
                    ## get the width in inches
                    maxLabelLength <- max(nchar(pdt$pid) / 10, na.rm = TRUE)
                    text_size <- min(4 * maxTextWidthBarplot / maxLabelLength,
                                     4)
                } else {
                    text_size <- 4
                }
                rng <- c(-max(abs(pdt$coef)), max(abs(pdt$coef)))
                pdt <- pdt %>%
                    dplyr::mutate(label_y = ifelse(.data$coef < 0, rng[2]/20,
                                                   rng[1]/20),
                                  label_hjust = ifelse(.data$coef < 0, 0, 1))
                ggplot2::ggplot(pdt, ggplot2::aes(
                    x = forcats::fct_reorder(.data$pid, .data$coef),
                    y = .data$coef, fill = sign(.data$coef))) +
                    ggplot2::geom_col() +
                    ggplot2::coord_flip() +
                    ggplot2::geom_text(ggplot2::aes(label = .data$pid,
                                                    y = .data$label_y,
                                                    hjust = .data$label_hjust),
                                       size = text_size) +
                    cowplot::theme_cowplot() +
                    ggplot2::theme(
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        legend.position = "none",
                        panel.grid.major.x = ggplot2::element_line(
                            colour = "grey80", linetype = "dashed"),
                        axis.line = element_blank()) +
                    ggplot2::scale_y_continuous(limits = rng) +
                    ggplot2::scale_fill_gradient2(low = scales::muted("red"),
                                                  high = scales::muted("blue"),
                                                  limits = c(-1, 1)) +
                    ggplot2::labs(title = paste0("PC", i), y = "coefficient")
            }),
            nrow = 1
        )

        ## Combine
        pscat <- plist[[paste0("PC", pr[1], "_", pr[2])]]
        if (!is.null(colourBy)) {
            pscat <- pscat +
                ggalt::geom_encircle(
                    ggplot2::aes(fill = .data[[colourBy]],
                                 group = .data[[colourBy]]),
                    alpha = 0.5, show.legend = FALSE, na.rm = TRUE,
                    s_shape = 0.5, expand = 0.05, spread = 0.1)
            if (length(unique(pcadf[[colourBy]])) > maxNGroups) {
                pscat <- pscat +
                    theme(legend.position = "none")
            }
        }
        plistcomb[[paste0("PC", pr[1], "_", pr[2])]] <-
            cowplot::plot_grid(
                pscat,
                pscree,
                pcoef,
                ncol = 1,
                rel_heights = c(1.25, 0.5, 1))
    }

    ## -------------------------------------------------------------------------
    ## Pairs plot
    ## -------------------------------------------------------------------------
    if (!is.null(colourBy)) {
        pairsdf <- scuttle::makePerCellDF(sce, use.dimred = paste0("PCA_",
                                                                   assayName),
                                          use.coldata = colourBy)
        pcidx <- which(colnames(pairsdf) != colourBy)
        colnames(pairsdf) <- sub(paste0("PCA_", assayName, "."), "PC",
                                 colnames(pairsdf))
        ppairs <- GGally::ggpairs(pairsdf, columns = pcidx,
                                  ggplot2::aes(colour = .data[[colourBy]]),
                                  upper = "blank")
    } else {
        pairsdf <- scuttle::makePerCellDF(sce, use.dimred = paste0("PCA_",
                                                                   assayName),
                                          use.coldata = FALSE)
        pcidx <- seq_len(ncol(pairsdf))
        colnames(pairsdf) <- sub(paste0("PCA_", assayName, "."), "PC",
                                 colnames(pairsdf))
        ppairs <- GGally::ggpairs(pairsdf, columns = pcidx,
                                  upper = "blank")
    }
    ppairs <- ppairs +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0("Assay: ", assayName))

    ## -------------------------------------------------------------------------
    ## Return values
    ## -------------------------------------------------------------------------
    return(list(sce = sce, plotcoord = plist, plotcombined = plistcomb,
                plotpairs = ppairs))
}
