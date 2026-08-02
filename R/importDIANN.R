#' Default filters (DF)
#'
#' einprot provides default filters for data imported from various tools.
#' To see the functions and thresholds used in each of these, print the
#' corresponding object below.
#'
#' @export
#' @name defaultFiltersDF
#' @rdname defaultFiltersDF
einprotDIANNFiltersDF <- list(Filter = function(df) {
    if ("Proteotypic" %in% colnames(df)) {
        df <- df |> dplyr::filter(Proteotypic == 1)
    }
    if ("PG.MaxLFQ.Quality" %in% colnames(df)) {
        df <- df |> dplyr::filter(PG.MaxLFQ.Quality >= 0.7)
    }
    if ("Global.PG.Q.Value" %in% colnames(df)) {
        df <- df |> dplyr::filter(Global.PG.Q.Value <= 0.01)
    }
    if ("Protein.Q.Value" %in% colnames(df)) {
        df <- df |> dplyr::filter(Protein.Q.Value <= 0.005)
    }
    # df <- dplyr::collect(dplyr::compute(df))
    if (all(c("Protein.Ids", "Precursor.Id") %in% colnames(df))) {
        df <- df |> dplyr::group_by(Protein.Ids) |>
            dplyr::mutate(N.Proteotypic.Sequences = dplyr::n_distinct(Precursor.Id)) |>
            dplyr::ungroup() |>
            dplyr::filter(N.Proteotypic.Sequences >= 2)
    }
    df
})

#' Import data from DIA-NN into a SingleCellExperiment object
#'
#' Import data from a DIA-NN quantification file into a
#' \code{SingleCellExperiment} object. Note that DIA-NN support in einprot is
#' currently experimental - please be aware that the interface may change, and
#' interpret results with caution.
#'
#' @param inFile Path to a tab-delimited input text file from DIA-NN; either
#'     \code{*.parquet}, \code{pg_matrix.tsv}, \code{pr_matrix.tsv} or
#'     \code{report.tsv}.
#' @param fileType Character scalar indicating the type of input file; either
#'     \code{"parquet"}, \code{"pg_matrix"}, \code{"pr_matrix"} or
#'     \code{"main_report"}.
#' @param outLevel Character scalar indicating the desired output level;
#'     either \code{"pg"} (protein group) or \code{"pr"} (precursor).
#' @param includeOnlySamples,excludeSamples Character vectors defining
#'     regular expressions to match against the extracted columns to
#'     retain or exclude samples.
#' @param stopIfEmpty Logical scalar, whether to raise an error if no
#'     columns matching the patterns are found.
#' @param aName Character scalar giving the name of the main assay (if
#'     \code{fileType} is \code{"pg_matrix"} or \code{"pr_matrix"}), or the
#'     column from which to get the values for the main assay (if
#'     \code{fileType} is \code{"main_report"} or \code{"parquet"}).
#' @param filtersDF Named list where each element is a filtering function
#'     that takes a \code{data.frame} (or \code{tibble}) as
#'     input and returns a filtered object of the same class. No other
#'     arguments are allowed. The filtering functions will be applied in the
#'     order they are provided in the list, after reading the long-format text
#'     file or parquet file, and before creating the wide-format assay. Only
#'     applies if \code{fileType} is \code{"main_report"} or \code{"parquet"}.
#' @param ... Additional arguments that will be passed on to
#'     \code{QFeatures::readSummarizedExperiment} (e.g., the number of rows
#'     to import).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A list with two elements: a \code{SingleCellExperiment} object and
#' a character scalar with the main assay name.
#'
#' @examples
#' sceL <- importDIANN(system.file("extdata", "diann_example",
#'                                 "PXD028735.pg_matrix.tsv",
#'                                 package = "einprot"),
#'                     fileType = "pg_matrix", outLevel = "pg",
#'                     aName = "MaxLFQ")
#' sceL
#'
#' @importFrom tibble tibble
#' @importFrom dplyr select filter inner_join collect compute all_of group_by
#' @importFrom dplyr summarize count distinct
#' @importFrom tidyr pivot_wider
#' @importFrom QFeatures readSummarizedExperiment
#' @importFrom SummarizedExperiment assayNames colnames
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom S4Vectors metadata DataFrame
#' @importFrom tools file_path_sans_ext
#' @importFrom methods as
#' @importFrom utils read.delim
#' @importFrom rlang .data
importDIANN <- function(inFile, fileType = "pg_matrix", outLevel = "pg",
                        includeOnlySamples = "",
                        excludeSamples = "", stopIfEmpty = FALSE,
                        aName = "PG.MaxLFQ", filtersDF = list(), ...) {
    ## Check input arguments
    .assertScalar(x = inFile, type = "character")
    stopifnot(file.exists(inFile))
    .assertScalar(x = fileType, type = "character",
                  validValues = c("pg_matrix", "pr_matrix", "main_report",
                                  "parquet"))
    .assertScalar(x = outLevel, type = "character",
                  validValues = c("pg", "pr"))
    .assertVector(x = includeOnlySamples, type = "character")
    .assertVector(x = excludeSamples, type = "character")
    .assertScalar(x = stopIfEmpty, type = "logical")
    .assertScalar(x = aName, type = "character")
    .assertVector(x = filtersDF, type = "list")
    if (length(filtersDF) > 0) {
        .assertVector(x = names(filtersDF), type = "character")
        for (f in filtersDF) {
            stopifnot(is(f, "function"))
            stopifnot(length(formals(f)) == 1)
        }
    }

    if (fileType == "pg_matrix" && outLevel == "pr") {
        stop("To obtain precursor output, please provide either the ",
             "pr_matrix or the main report.")
    }
    if (fileType == "pr_matrix" && outLevel == "pg") {
        stop("To obtain protein group output, please provide either the ",
             "pg_matrix or the main report.")
    }

    if (fileType %in% c("pg_matrix", "pr_matrix")) {
        cn <- getColumnNames(inFile, check.names = FALSE)
        iColsAll <- grep("\\.mzML$|\\.wiff$|\\.dia$|\\.d$|\\.raw$|\\|/", cn, value = TRUE)
        iCols <- .getiCols(iColsAll = iColsAll,
                           includeOnlySamples = includeOnlySamples,
                           excludeSamples = excludeSamples,
                           stopIfEmpty = stopIfEmpty)

        sce <- QFeatures::readSummarizedExperiment(
            inFile, quantCols = iCols, sep = "\t", check.names = FALSE, ...
        )

        SummarizedExperiment::assayNames(sce) <- aName

        ## Add list of columns to metadata
        S4Vectors::metadata(sce)$colList <- list()
        S4Vectors::metadata(sce)$colList[[aName]] <- iCols

        ## Remove directory name and extension from colnames
        colnames(sce) <- tools::file_path_sans_ext(basename(
            gsub("\\\\", "/", colnames(sce))))

        sce <- methods::as(sce, "SingleCellExperiment")
        return(list(sce = sce, aName = aName))

    } else if (fileType %in% c("main_report", "parquet")) {
        if (fileType == "main_report") {
            tmp <- read.delim(inFile, header = TRUE, sep = "\t")
        } else {
            # read data as tibble
            .assertPackagesAvailable("arrow")
            tmp <- arrow::read_parquet(inFile, as_data_frame = TRUE)
        }
        iColsAll <- tmp |> dplyr::select("Run") |> distinct() |> pull()
        iCols <- .getiCols(iColsAll = iColsAll,
                           includeOnlySamples = includeOnlySamples,
                           excludeSamples = excludeSamples,
                           stopIfEmpty = stopIfEmpty)
        tmp <- tmp |> dplyr::filter(.data$Run %in% iCols)
        stopifnot(aName %in% colnames(tmp))

        # Filter
        for (f in filtersDF) {
            tmp <- f(tmp)
        }

        # tmp <- dplyr::collect(dplyr::compute(tmp))
        if (outLevel == "pg") {
            # if Genes and/or Protein.Names columns are not present, add NA
            # columns
            if (!"Genes" %in% colnames(tmp)) {
                tmp$Genes <- NA_character_
            }
            if (!"Protein.Names" %in% colnames(tmp)) {
                tmp$Protein.Names <- NA_character_
            }
            tmp <- tmp |> dplyr::select(all_of(c("Run", "Protein.Group",
                                                 "Genes", "Protein.Names",
                                                 "Protein.Ids", aName,
                                                 grep("^PG\\.", colnames(tmp), value = TRUE))))
            tmp <- tmp |>
                dplyr::group_by(dplyr::across(-c("Protein.Ids", "Genes", "Protein.Names"))) %>%
                dplyr::summarize(
                    Protein.Ids =
                        paste(unique(unlist(strsplit(.data$Protein.Ids, ";"))), collapse = ";"),
                    Protein.Names =
                        paste(unique(unlist(strsplit(.data$Protein.Names, ";"))), collapse = ";"),
                    Genes =
                        paste(unique(unlist(strsplit(.data$Genes, ";"))), collapse = ";"),
                    .groups = "drop") |>
                dplyr::distinct()
            rd <- DataFrame(Protein.Group = unique(tmp$Protein.Group))
            aL <- list()
            for (nm in setdiff(colnames(tmp), c("Run", "Protein.Group"))) {
                ## Check if it's a sample-specific value or not
                tmpsub <- tmp %>%
                    dplyr::select(dplyr::all_of(c("Protein.Group", nm)))
                tmpcount <- tmpsub %>%
                    dplyr::count(.data$Protein.Group)
                if (all(tmpcount$n == 1) && length(iCols) != 1) {
                    ## One value per protein group -> annotation
                    tmpsub <- tmpsub %>%
                        dplyr::distinct()
                    rd[[nm]] <- tmpsub[[nm]][match(rd$Protein.Group, tmpsub$Protein.Group)]
                } else {
                    ## One value per protein group/sample -> assay
                    fv <- ifelse(is.numeric(tmp[[nm]]), 0, "0")
                    tmpsub <- tmp %>%
                        dplyr::select(dplyr::all_of(c("Run", "Protein.Group", nm))) %>%
                        tidyr::pivot_wider(names_from = "Run",
                                           values_from = dplyr::all_of(nm),
                                           values_fill = fv) %>%
                        as.data.frame()
                    rownames(tmpsub) <- tmpsub$Protein.Group
                    tmpsub$Protein.Group <- NULL
                    tmpsub <- as.matrix(tmpsub)
                    aL[[nm]] <- tmpsub[match(rd$Protein.Group, rownames(tmpsub)),
                                       match(iCols, colnames(tmpsub)), drop = FALSE]
                }
            }
            ## Add annotation columns corresponding to specific assays
            for (aNm in c("Genes", "Protein.Names", "Protein.Ids")) {
                if (aNm %in% names(aL)) {
                    ## Aggregate across columns for each row
                    rd[[aNm]] <- apply(aL[[aNm]], 1, function(x) {
                        paste(setdiff(unique(unlist(strsplit(x, ";"))), "0"), collapse = ";")
                    })
                }
            }
            sce <- SingleCellExperiment::SingleCellExperiment(
                assays = aL,
                rowData = rd,
                colData = DataFrame(sampleName = iCols)
            )

            ## Add list of columns to metadata
            S4Vectors::metadata(sce)$colList <- list()
            return(list(sce = sce, aName = aName))
        } else {
            stop("Not yet implemented")
        }
    }
}
