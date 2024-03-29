% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addSampleAnnots.R
\name{addSampleAnnots}
\alias{addSampleAnnots}
\title{Add sample annotations to SummarizedExperiment object}
\usage{
addSampleAnnots(sce, sampleAnnot)
}
\arguments{
\item{sce}{A \code{SummarizedExperiment} object (or a derivative).}

\item{sampleAnnot}{A \code{data.frame} with sample annotations. Must
have at least columns named \code{sample} (which must contain
all the column names of \code{sce}) and \code{group} (which contains
the group assignment for each sample).}
}
\value{
An object of the same type as \code{sce} with additional sample
    annotations.
}
\description{
Add sample annotations from an external annotation table to an existing
\code{SummarizedExperiment} object.
}
\examples{
## Import example data
mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                      package = "einprot")
samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
             "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
             "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\\\.",
                        includeOnlySamples = samples)

## Define sample annotations
sampleAnnot <- data.frame(sample = samples,
                          group = gsub("_IP.*", "", samples),
                          random = sample(seq_len(9), 9))

## Add sample annotations to SCE
sce <- addSampleAnnots(out$sce, sampleAnnot = sampleAnnot)
SummarizedExperiment::colData(sce)  ## group information added to sce

}
\author{
Charlotte Soneson
}
