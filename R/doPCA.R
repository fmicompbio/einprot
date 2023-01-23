#' Run PCA and generate plots
#'
#' Apply PCA to the assay defined by \code{assayName} and extract the top
#' \code{ncomponent} components. For each pair defined in \code{plotpairs},
#' generate a collection of plots.
#'
#' @param sce A \code{SingleCellExperiment} object.
#' @param assayName A character scalar defining the assay in \code{sce} to
#'     use for the PCA.
#' @param ncomponents Numeric scalar, the (maximal) number of components to
#'     extract. The actual number can be lower if the number of samples is
#'     too small.
#' @param ntop Number of features (with highest variance) to use to
#'     generate the PCA. Will be passed on to \code{scater::runPCA}.
#' @param plotpairs A list of numeric vectors of length 2, indicating which
#'     pairs of PCs to generate plots for.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A list with the following components: sce (the input sce,
#' expanded with the calculated PCAs, in addition the feature coefficients
#' will be added to the rowData), plotcoord (a list of ggplot objects
#' containing coordinate plots for the desired pairs of components),
#' plotcombined (a list of ggplot objects containing combined coordinate,
#' scree and coefficient plots for the desired pairs of components),
#' plotpairs (a ggpairs plot with all extracted components).
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
#' @importFrom SummarizedExperiment rowData
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom GGally ggpairs
#' @importFrom utils head tail
#'
doPCA <- function(sce, assayName, ncomponents = 10, ntop = Inf,
                  plotpairs = list(c(1, 2))) {

    ## ---------------------------------------------------------------------- ##
    ## Check arguments
    ## ---------------------------------------------------------------------- ##
    .assertVector(x = sce, type = "SingleCellExperiment")
    .assertScalar(x = assayName, type = "character",
                  validValues = SummarizedExperiment::assayNames(sce))
    .assertScalar(x = ncomponents, type = "numeric")
    .assertScalar(x = ntop, type = "numeric", rngIncl = c(1, Inf))
    .assertVector(x = plotpairs, type = "list")
    for (elm in plotpairs) {
        .assertVector(x = elm, type = "numeric", len = 2)
    }
    ncomponents <- min(ncomponents, ncol(sce) - 1)
    if (any(unlist(plotpairs) > ncomponents)) {
        stop("'plotpairs' requests components that will not be extracted")
    }

    ## ---------------------------------------------------------------------- ##
    ## Run PCA and add to sce
    ## ---------------------------------------------------------------------- ##
    sce <- scater::runPCA(sce, exprs_values = assayName,
                          ncomponents = ncomponents, ntop = ntop,
                          BSPARAM = BiocSingular::ExactParam(),
                          name = paste0("PCA_", assayName))

    ## ---------------------------------------------------------------------- ##
    ## Add coefficients to rowData (replace any existing columns with the
    ## same names)
    ## ---------------------------------------------------------------------- ##
    cf <- attr(SingleCellExperiment::reducedDim(sce, paste0("PCA_", assayName)),
               "rotation")
    colnames(cf) <- paste0("PCA_", assayName, "_", colnames(cf))
    stopifnot(rownames(sce) %in% rownames(cf))
    cf <- cf[match(rownames(sce), rownames(cf)), , drop = FALSE]
    stopifnot(rownames(sce) == rownames(cf))
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

    ## ---------------------------------------------------------------------- ##
    ## Sample plot
    ## ---------------------------------------------------------------------- ##
    percvar <- attr(SingleCellExperiment::reducedDim(sce, paste0("PCA_", assayName)),
                    "percentVar")
    plist <- list()
    for (pr in plotpairs) {
        pcadf <- data.frame(
            sampleLabel = colnames(sce),
            SingleCellExperiment::reducedDim(sce, paste0("PCA_", assayName))[, pr],
            group = sce$group)
        p <- ggplot2::ggplot(pcadf,
                             ggplot2::aes(x = .data[[paste0("PC", pr[1])]],
                                          y = .data[[paste0("PC", pr[2])]],
                                          color = .data$group,
                                          label = .data$sampleLabel)) +
            ggplot2::geom_point(alpha = 0.5, size = 4) +
            ggplot2::ggtitle(paste0("Assay: ", assayName)) +
            ggplot2::coord_fixed() +
            ggplot2::labs(
                x = paste0("PC", pr[1], " (", round(percvar[pr[1]], 2), "%)"),
                y = paste0("PC", pr[2], " (", round(percvar[pr[2]], 2), "%)")) +
            cowplot::theme_cowplot()
        plist[[paste0("PC", pr[1], "_", pr[2])]] <- p
    }

    ## ---------------------------------------------------------------------- ##
    ## Combined plot
    ## ---------------------------------------------------------------------- ##
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
            ggplot2::scale_y_continuous(labels = scales::percent_format(scale = 1)) +
            ggplot2::scale_x_continuous(breaks = seq_along(percvar)) +
            cowplot::theme_cowplot() +
            ggplot2::theme(legend.position = "none")

        ## Coefficient plot
        pcoef <- cowplot::plot_grid(
            plotlist = lapply(pr, function(i) {
                pdt <- dplyr::bind_rows(
                    data.frame(coef = cf[, paste0("PCA_", assayName, "_PC", i)]) %>%
                        dplyr::filter(.data$coef > 0) %>%
                        dplyr::arrange(dplyr::desc(.data$coef)) %>%
                        utils::head(n = 10) %>%
                        tibble::rownames_to_column("pid"),
                    data.frame(coef = cf[, paste0("PCA_", assayName, "_PC", i)]) %>%
                        dplyr::filter(.data$coef < 0) %>%
                        dplyr::arrange(dplyr::desc(.data$coef)) %>%
                        utils::tail(n = 10) %>%
                        tibble::rownames_to_column("pid")
                )
                rng <- c(-max(abs(pdt$coef)), max(abs(pdt$coef)))
                pdt <- pdt %>%
                    dplyr::mutate(label_y = ifelse(.data$coef < 0, rng[2]/20, rng[1]/20),
                                  label_hjust = ifelse(.data$coef < 0, 0, 1))
                ggplot2::ggplot(pdt, ggplot2::aes(
                    x = forcats::fct_reorder(.data$pid, .data$coef),
                    y = .data$coef, fill = sign(.data$coef))) +
                    ggplot2::geom_col() +
                    ggplot2::coord_flip() +
                    ggplot2::geom_text(ggplot2::aes(label = .data$pid, y = .data$label_y,
                                                    hjust = .data$label_hjust)) +
                    cowplot::theme_cowplot() +
                    ggplot2::theme(
                        axis.text.y = ggplot2::element_blank(),
                        axis.ticks.y = ggplot2::element_blank(),
                        axis.title.y = ggplot2::element_blank(),
                        legend.position = "none",
                        panel.grid.major.x = ggplot2::element_line(colour = "grey80",
                                                                   linetype = "dashed"),
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
        plistcomb[[paste0("PC", pr[1], "_", pr[2])]] <-
            cowplot::plot_grid(
                plist[[paste0("PC", pr[1], "_", pr[2])]] +
                    ggalt::geom_encircle(
                        ggplot2::aes(fill = .data$group, group = .data$group),
                        alpha = 0.5, show.legend = FALSE, na.rm = TRUE,
                        s_shape = 0.5, expand = 0.05, spread = 0.1),
                pscree,
                pcoef,
                ncol = 1,
                rel_heights = c(1.25, 0.5, 1))
    }

    ## ---------------------------------------------------------------------- ##
    ## Pairs plot
    ## ---------------------------------------------------------------------- ##
    pairsdf <- scuttle::makePerCellDF(sce, use.dimred = paste0("PCA_", assayName),
                                      use.coldata = "group")
    pcidx <- which(colnames(pairsdf) != "group")
    colnames(pairsdf) <- sub(paste0("PCA_", assayName, "."), "PC",
                             colnames(pairsdf))
    ppairs <- GGally::ggpairs(pairsdf, columns = pcidx,
                              ggplot2::aes(colour = .data$group),
                              upper = "blank") +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(paste0("Assay: ", assayName))

    ## ---------------------------------------------------------------------- ##
    ## Return values
    ## ---------------------------------------------------------------------- ##
    return(list(sce = sce, plotcoord = plist, plotcombined = plistcomb,
                plotpairs = ppairs))
}
