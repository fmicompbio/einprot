% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makeDbLinkTable.R
\name{formatTableColumns}
\alias{formatTableColumns}
\title{Format table columns}
\usage{
formatTableColumns(tbl, columns, signifDigits, maxLevels = 10)
}
\arguments{
\item{tbl}{A \code{data.frame} or similar object.}

\item{columns}{Character vector giving the column names of \code{tbl} that
should be formatted.}

\item{signifDigits}{Numeric scalar, the number of significant digits to
round numeric columns to. Set to \code{NULL} to skip rounding.}

\item{maxLevels}{Numeric scalar. If character columns have at most this
number of unique values, they will be encoded as factors.}
}
\value{
The input \code{data.frame} with formatted columns.
}
\description{
Format columns in a \code{data.frame}.
}
\examples{
df <- formatTableColumns(data.frame(x = rnorm(10),
                                    y = sample(LETTERS[seq_len(3)], 10,
                                               replace = TRUE),
                                    z = rnorm(10)),
                         columns = c("x", "y"),
                         signifDigits = 2, maxLevels = 10)
df
summary(df)

df <- formatTableColumns(data.frame(x = rnorm(10),
                                    y = sample(LETTERS[seq_len(3)], 10,
                                               replace = TRUE),
                                    z = rnorm(10)),
                         columns = c("x", "y", "z"),
                         signifDigits = 2, maxLevels = 2)
df
summary(df)

}
\author{
Charlotte Soneson
}
