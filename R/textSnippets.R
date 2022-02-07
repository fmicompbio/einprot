#' Text snippets for use in analysis reports
#'
#' These text snippets are pasted in the analysis reports and describes the
#' analysis that was performed. They are not intended to be called
#' directly by the user.
#'
#' @author Charlotte Soneson
#' @name textSnippets
NULL

#' @rdname textSnippets
#' @export
limmaText <- function() {
    cat(paste0("For this, we use the treat function from the ",
               "[limma](https://bioconductor.org/packages/limma/) R/Bioconductor ",
               "package [@McCarthy2009treat; @Ritchie2015limma; @Phipson2016robust]. ",
               "For more information about the df.prior, representing the amount of ",
               "extra information that is borrowed from the full set of features ",
               "in order to improve the inference for each feature, see section 13.2 in the ",
               "[limma user guide](https://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf). ",
               "In addition to the feature-wise tests, we apply the ",
               "camera method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistic returned from limma."))
}

#' @rdname textSnippets
#' @export
ttestText <- function() {
    cat(paste0("For this, we use a Student's t-test. To determine which features ",
               "show significant changes, we calculate the SAM statistic [@Tusher2001sam], ",
               "and estimate the false discovery rate at different thresholds ",
               "using permutations, mimicking the approach used by Perseus [@Tyanova2016perseus].",
               "In addition to the feature-wise tests, we apply the ",
               "camera method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "SAM statistic calculated from the t-statistic and the specified S0."))
}
