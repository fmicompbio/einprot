#' Text snippets for use in analysis reports
#'
#' These text snippets are pasted in the analysis reports and describe the
#' analysis that was performed. They are not intended to be called
#' directly by the user.
#'
#' @param testType Character scalar giving the statistical test, either
#'     "limma" or "ttest".
#' @param minlFC Numeric scalar giving the minimum logFC threshold.
#' @param normMethod Character scalar giving the normalization method.
#'
#' @author Charlotte Soneson
#' @name textSnippets
#'
#' @return A character string.
#' @examples
#' testText(testType = "limma")
#' testText(testType = "ttest")
NULL

#' @rdname textSnippets
#' @export
testText <- function(testType, minlFC = 0, samSignificance = TRUE) {
    .assertScalar(x = testType, type = "character",
                  validValues = c("ttest", "limma", "proDA"))
    .assertScalar(x = minlFC, type = "numeric")
    .assertScalar(x = samSignificance, type = "logical")
    if (testType == "limma" && minlFC != 0) {
        paste0("For this, we use the treat function from the ",
               "[limma](https://bioconductor.org/packages/limma/) ",
               "R/Bioconductor package [@McCarthy2009treat; ",
               "@Ritchie2015limma; @Phipson2016robust]. For more ",
               "information about the df.prior, representing the amount of ",
               "extra information that is borrowed from the full set of ",
               "features in order to improve the inference for each ",
               "feature, see section 13.2 in the [limma user guide]",
               "(https://www.bioconductor.org/packages/devel/bioc/vignettes",
               "/limma/inst/doc/usersguide.pdf). ",
               "In addition to the feature-wise tests, we apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistic returned from limma.")
    } else if (testType == "limma" && minlFC == 0) {
        paste0("For this, we use the ",
               "[limma](https://bioconductor.org/packages/limma/) ",
               "R/Bioconductor package [",
               "@Ritchie2015limma; @Phipson2016robust]. For more ",
               "information about the df.prior, representing the amount of ",
               "extra information that is borrowed from the full set of ",
               "features in order to improve the inference for each ",
               "feature, see section 13.2 in the [limma user guide]",
               "(https://www.bioconductor.org/packages/devel/bioc/vignettes",
               "/limma/inst/doc/usersguide.pdf). ",
               "In addition to the feature-wise tests, we apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistic returned from limma.")
    } else if (testType == "ttest" && samSignificance) {
        paste0("For this, we use a Student's t-test. To determine which ",
               "features show significant changes, we calculate the SAM ",
               "statistic [@Tusher2001sam], and estimate the false ",
               "discovery rate at different thresholds using permutations, ",
               "mimicking the approach used by Perseus [@Tyanova2016perseus].",
               "In addition to the feature-wise tests, we apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "SAM statistic calculated from the t-statistic and the ",
               "specified S0.")
    } else if (testType == "ttest" && !samSignificance) {
        paste0("For this, we use a Student's t-test. To determine which ",
               "features show significant changes, we calculate ",
               "adjusted p-values using the Benjamini-Hochberg method ",
               "[@BenjaminiHochberg1995fdr].",
               "In addition to the feature-wise tests, we apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistic.")
    } else if (testType == "proDA") {
        paste0("For this, we use the ",
               "[proDA](https://bioconductor.org/packages/proDA/) ",
               "R/Bioconductor package [@AhlmannEltze2020proda]. ",
               "In addition to the feature-wise tests, we apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistic returned from proDA.")
    }
}

#' @rdname textSnippets
#' @export
normText <- function(normMethod) {
    .assertScalar(x = normMethod, type = "character")

    if (normMethod == "none") {
        paste0("The log2 intensities are not normalized further across ",
               "samples since 'normMethod' is set to 'none'.")
    } else {
        paste0("The log2 intensities are next normalized across samples ",
               "using the ", normMethod, " method.")
    }
}

#' @rdname textSnippets
#' @export
saText <- function(testType) {
    .assertScalar(x = testType, type = "character",
                  validValues = c("limma", "ttest", "proDA"))

    if (testType == "limma") {
        paste0("We first show a diagnostic plot for each comparison. These ",
               "plots display the square root of the residual standard deviation ",
               "(y-axis) versus the mean abundance (across all the groups used to ",
               "perform the model fit, x-axis). The curve indicated in the plots show ",
               "the mean-variance trend inferred by `limma`.")
    } else {
        ""
    }
}

#' @rdname textSnippets
#' @export
expDesignText <- function(testType) {
    .assertScalar(x = testType, type =  "character",
                  validValues = c("limma", "ttest", "proDA"))

    if (testType == "limma") {
        paste0("The plots below illustrate the experimental design used ",
               "for the linear model(s) and contrasts by `limma`. The plot ",
               "to the right shows the number of samples for each combination ",
               "of factor levels across the predictors, and is useful for ",
               "detecting imbalances between group sizes for different ",
               "conditions. The plot to the left summarizes the expected ",
               "response value for each combination of predictor levels, ",
               "expressed in terms of the linear model coefficients. ",
               "For more details on how to interpret the plots, we refer to ",
               "@Soneson2020emm or @Law2020design. Clicking on the arrow ",
               "below the plots will reveal the design matrix used by ",
               "limma, as well as the contrasts that were fit for each ",
               "comparison.")
    } else {
        ""
    }
}
