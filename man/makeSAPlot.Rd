% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summaryPlots.R
\name{makeSAPlot}
\alias{makeSAPlot}
\title{Construct SA plot from limma results}
\usage{
makeSAPlot(testList)
}
\arguments{
\item{testList}{List of test results, typically generated using
\code{runTest()}.}
}
\value{
A \code{cowplot} object
}
\description{
Given a list of \code{data.frames} with \code{limma} test results, create
an SA plot for each contrast.
}
\examples{
sce <- readRDS(system.file("extdata", "mq_example", "1356_sce.rds",
                           package = "einprot"))
tres <- runTest(sce, comparisons = list(c("RBC_ctrl", "Adnp"),
                                        c("RBC_ctrl", "Chd4BF")),
                testType = "limma", assayForTests = "log2_LFQ.intensity",
                assayImputation = "imputed_LFQ.intensity")
makeSAPlot(tres$tests)

}
\author{
Charlotte Soneson
}
