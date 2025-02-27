% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/textSnippets.R
\name{textSnippets}
\alias{textSnippets}
\alias{emptySampleText}
\alias{testText}
\alias{normText}
\alias{saText}
\alias{expDesignText}
\alias{introText}
\alias{filterByModText}
\alias{inputText}
\alias{featureCollectionText}
\title{Text snippets for use in analysis reports}
\usage{
emptySampleText(sce, assayName)

testText(testType, minlFC = 0, samSignificance = TRUE)

normText(normMethod)

saText(testType)

expDesignText(testType)

introText(expType)

filterByModText(excludeUnmodifiedPeptides, keepModifications)

inputText(expTypeLevel)

featureCollectionText(featureCollections)
}
\arguments{
\item{sce}{A SummarizedExperiment object.}

\item{assayName}{Character scalar representing the assay that will be used
to check whether there are samples with missing (NA) values for all
features.}

\item{testType}{Character scalar giving the statistical test, either
\code{"limma"}, \code{"ttest"}, or \code{"proDA"}.}

\item{minlFC}{Numeric scalar giving the minimum logFC threshold.}

\item{samSignificance}{Logical scalar indicating whether the SAM statistic
should be used to determine significance for the t-test.}

\item{normMethod}{Character scalar giving the normalization method.}

\item{expType}{The quantification/identification tool used to generate the
data, one of "MaxQuant", "ProteomeDiscoverer", "FragPipe", "DIANN"
or "Spectronaut".}

\item{expTypeLevel}{The quantification/identification tool used to generate
the data (including the data level for ProteomeDiscoverer), one of
"MaxQuant", "ProteomeDiscovererProteins",
"ProteomeDiscovererPeptideGroups", "FragPipe", "DIANN"
or "Spectronaut".}
}
\value{
A character string.
}
\description{
These text snippets are pasted in the analysis reports and describe the
analysis that was performed. The functions are not intended to be called
directly by the user.
}
\examples{
testText(testType = "limma")
testText(testType = "ttest")
normText(normMethod = "none")
expDesignText(testType = "limma")
}
\author{
Charlotte Soneson
}
