% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryPlots.R
\name{makeIntensityBoxplots}
\alias{makeIntensityBoxplots}
\title{Make intensity boxplot}
\usage{
makeIntensityBoxplots(sce, assayName, doLog, ylab, maxNGroups = 25)
}
\arguments{
\item{sce}{A \code{SummarizedExperiment} object (or a derivative). The
\code{colData} must have columns named \code{"sample"} and
\code{"group"}, for grouping and coloring the values, respectively.}

\item{assayName}{Character scalar, the name of the assay of \code{sce}
to use for the plots.}

\item{doLog}{Logical scalar, whether to log-transform the y-axis.}

\item{ylab}{Character scalar, the label to use for the y-axis.}

\item{maxNGroups}{Numeric scalar, the maximum number of groups to display
in the legend. If there are more than \code{maxNGroups} groups, the
legend is suppressed.}
}
\value{
A \code{ggplot} object.
}
\description{
Generate boxplots of the intensity distributions for the samples in a
\code{SummarizedExperiment} object.
}
\examples{
mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                      package = "einprot")
samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
             "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
             "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\\\.",
                        includeOnlySamples = samples)
sampleAnnot <- data.frame(sample = samples,
                          group = gsub("_IP.*", "", samples))
sce <- addSampleAnnots(out$sce, sampleAnnot = sampleAnnot)
makeIntensityBoxplots(sce, assayName = "iBAQ", doLog = TRUE,
                      ylab = "log intensity")

}
\author{
Charlotte Soneson
}
