% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotMissingValues.R
\name{plotMissingValuesHeatmap}
\alias{plotMissingValuesHeatmap}
\title{Plot heatmap of missing values}
\usage{
plotMissingValuesHeatmap(sce, assayMissing)
}
\arguments{
\item{sce}{A \code{SummarizedExperiment} object.}

\item{assayMissing}{Character scalar indicating the name of a
logical assay of \code{sce} representing the missingness pattern.
\code{"FALSE"} entries should represent observed values, while
\code{"TRUE"} entries represent missing values.}
}
\value{
A \code{ComplexHeatmap} object.
}
\description{
Create a heatmap of the data matrix, indicating which values are missing
and observed, respectively.
}
\examples{
sce <- importExperiment(system.file("extdata", "mq_example",
                                    "1356_proteinGroups.txt",
                                    package = "einprot"),
                        iColPattern = "^iBAQ\\\\.")$sce
SummarizedExperiment::assay(sce, "iBAQ")[
    SummarizedExperiment::assay(sce, "iBAQ") == 0] <- NA
SummarizedExperiment::assay(sce, "missing") <-
    is.na(SummarizedExperiment::assay(sce, "iBAQ"))
plotMissingValuesHeatmap(sce, "missing")

}
\author{
Charlotte Soneson
}
