#' Text snippets for use in analysis reports
#'
#' These text snippets are pasted in the analysis reports and describe the
#' analysis that was performed. The functions are not intended to be called
#' directly by the user.
#'
#' @param testType Character scalar giving the statistical test, either
#'     \code{"limma"}, \code{"ttest"}, or \code{"proDA"}.
#' @param minlFC Numeric scalar giving the minimum logFC threshold.
#' @param samSignificance Logical scalar indicating whether the SAM statistic
#'     should be used to determine significance for the t-test.
#' @param normMethod Character scalar giving the normalization method.
#' @param sce A SummarizedExperiment object.
#' @param assayName Character scalar representing the assay that will be used
#'     to check whether there are samples with missing (NA) values for all
#'     features.
#' @param expType The quantification/identification tool used to generate the
#'     data, one of "MaxQuant", "ProteomeDiscoverer", "FragPipe", "DIANN"
#'     or "Spectronaut".
#' @param expTypeLevel The quantification/identification tool used to generate
#'     the data (including the data level for ProteomeDiscoverer), one of
#'     "MaxQuant", "ProteomeDiscovererProteins",
#'     "ProteomeDiscovererPeptideGroups", "FragPipe", "DIANN"
#'     or "Spectronaut".
#'
#' @author Charlotte Soneson
#' @name textSnippets
#'
#' @returns A character string.
#' @examples
#' testText(testType = "limma")
#' testText(testType = "ttest")
#' normText(normMethod = "none")
#' expDesignText(testType = "limma")
NULL

#' @rdname textSnippets
#' @export
#' @importFrom SummarizedExperiment assayNames assay
emptySampleText <- function(sce, assayName) {
    .assertVector(x = sce, type = "SummarizedExperiment")
    .assertScalar(x = assayName, type = "character")
    stopifnot(assayName %in% SummarizedExperiment::assayNames(sce))

    no_feats <- which(colSums(!is.na(
        SummarizedExperiment::assay(sce, assayName))) == 0)
    if (length(no_feats) > 0) {
        paste0("The following sample(s) do not have any detected features ",
               "and will be removed from further analysis: ",
               paste(colnames(sce)[no_feats], collapse = ", "), ".")
    } else {
        ""
    }
}

#' @rdname textSnippets
#' @export
testText <- function(testType, minlFC = 0, samSignificance = TRUE) {
    .assertScalar(x = testType, type = "character",
                  validValues = c("ttest", "limma", "proDA", "none"))
    .assertScalar(x = minlFC, type = "numeric")
    .assertScalar(x = samSignificance, type = "logical")

    if (testType == "limma" && minlFC != 0) {
        paste0("For each feature, we then compare the (possibly imputed) ",
               "log2 intensities between groups. ",
               "For this, we use the treat function from the ",
               "[limma](https://bioconductor.org/packages/limma/) ",
               "R/Bioconductor package [@McCarthy2009treat; ",
               "@Ritchie2015limma; @Phipson2016robust]. For more ",
               "information about the df.prior, representing the amount of ",
               "extra information that is borrowed from the full set of ",
               "features in order to improve the inference for each ",
               "feature, see section 13.2 in the [limma user guide]",
               "(https://www.bioconductor.org/packages/devel/bioc/vignettes",
               "/limma/inst/doc/usersguide.pdf). ",
               "If requested, in addition to the feature-wise tests, we ",
               "apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistics returned from limma.")
    } else if (testType == "limma" && minlFC == 0) {
        paste0("For each feature, we then compare the (possibly imputed) ",
               "log2 intensities between groups. ",
               "For this, we use the ",
               "[limma](https://bioconductor.org/packages/limma/) ",
               "R/Bioconductor package [",
               "@Ritchie2015limma; @Phipson2016robust]. For more ",
               "information about the df.prior, representing the amount of ",
               "extra information that is borrowed from the full set of ",
               "features in order to improve the inference for each ",
               "feature, see section 13.2 in the [limma user guide]",
               "(https://www.bioconductor.org/packages/devel/bioc/vignettes",
               "/limma/inst/doc/usersguide.pdf). ",
               "If requested, in addition to the feature-wise tests, we ",
               "apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistics returned from limma.")
    } else if (testType == "ttest" && samSignificance) {
        paste0("For each feature, we then compare the (possibly imputed) ",
               "log2 intensities between groups. ",
               "For this, we use a Student's t-test. To determine which ",
               "features show significant changes, we calculate the SAM ",
               "statistic [@Tusher2001sam], and estimate the false ",
               "discovery rate at different thresholds using permutations, ",
               "mimicking the approach used by Perseus [@Tyanova2016perseus]. ",
               "If requested, in addition to the feature-wise tests, we ",
               "apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "SAM statistics calculated from the t-statistics and the ",
               "specified S0.")
    } else if (testType == "ttest" && !samSignificance) {
        paste0("For each feature, we then compare the (possibly imputed) ",
               "log2 intensities between groups. ",
               "For this, we use a Student's t-test. To determine which ",
               "features show significant changes, we calculate ",
               "adjusted p-values using the Benjamini-Hochberg method ",
               "[@BenjaminiHochberg1995fdr]. ",
               "If requested, in addition to the feature-wise tests, we ",
               "apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistics.")
    } else if (testType == "proDA") {
        paste0("For each feature, we then compare the (possibly imputed) ",
               "log2 intensities between groups. ",
               "For this, we use the ",
               "[proDA](https://bioconductor.org/packages/proDA/) ",
               "R/Bioconductor package [@AhlmannEltze2020proda]. ",
               "If requested, in addition to the feature-wise tests, we ",
               "apply the camera ",
               "method [@Wu2012camera] to test for significance of each ",
               "included feature collection. These tests are based on the ",
               "t-statistics returned from proDA.")
    } else if (testType == "none") {
        paste0("Since testType = 'none', no statistical testing will be",
               "performed.")
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
                  validValues = c("limma", "ttest", "proDA", "none"))

    if (testType == "limma") {
        paste0("We first show a diagnostic plot for each comparison. These ",
               "plots display the square root of the residual standard ",
               "deviation (y-axis) versus the mean abundance (across all the ",
               "groups used to perform the model fit, x-axis). The curve ",
               "indicated in the plots show the mean-variance trend ",
               "inferred by `limma`.")
    } else {
        ""
    }
}

#' @rdname textSnippets
#' @export
expDesignText <- function(testType) {
    .assertScalar(x = testType, type =  "character",
                  validValues = c("limma", "ttest", "proDA", "none"))

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

#' @rdname textSnippets
#' @export
introText <- function(expType) {
    .assertScalar(x = expType, type = "character",
                  validValues = c("MaxQuant", "ProteomeDiscoverer",
                                  "FragPipe", "DIANN", "Spectronaut"))

    if (expType == "MaxQuant") {
        paste0("This report describes a reproducible end-to-end analysis of ",
               "a proteomics dataset quantified with ",
               "[MaxQuant](https://www.maxquant.org/) [@Cox2008maxquant]. ")
    } else if (expType == "FragPipe") {
        paste0("This report describes a reproducible end-to-end analysis of ",
               "a proteomics dataset quantified with ",
               "[FragPipe](https://fragpipe.nesvilab.org/). ")
    } else if (expType == "ProteomeDiscoverer") {
        paste0("This report describes a reproducible end-to-end analysis of ",
               "a proteomics dataset quantified with ",
               "[Proteome Discoverer](https://www.thermofisher.com/ch/en/home/",
               "industrial/mass-spectrometry/liquid-chromatography-mass-",
               "spectrometry-lc-ms/lc-ms-software/multi-omics-data-analysis/",
               "proteome-discoverer-software.html) [@Orsburn2021pd]. ")
    } else if (expType == "DIANN") {
        paste0("This report describes a reproducible end-to-end analysis of ",
               "a proteomics dataset quantified with ",
               "[DIA-NN](https://github.com/vdemichev/DiaNN) [@Demichev2020diann]. ")
    } else if (expType == "Spectronaut") {
        paste0("This report describes a reproducible end-to-end analysis of ",
               "a proteomics dataset quantified with ",
               "[Spectronaut](https://biognosys.com/resources/spectronaut-",
               "a-groundbreaking-increase-in-identifications/).")
    }
}

#' @rdname textSnippets
#' @export
filterByModText <- function(excludeUnmodifiedPeptides, keepModifications) {
    if (excludeUnmodifiedPeptides && !is.null(keepModifications)) {
        paste0("Next, we filter out unmodified peptides and peptides ",
               "without any of the requested modifications ",
               "(", paste(keepModifications, collapse = ", "), ").")
    } else if (excludeUnmodifiedPeptides) {
        paste0("Next, we filter out unmodified peptides.")
    } else if (!is.null(keepModifications)) {
        paste0("Next, we filter out peptides ",
               "without any of the requested modifications ",
               "(", paste(keepModifications, collapse = ", "), ").")
    } else {
        ""
    }
}

#' @rdname textSnippets
#' @export
inputText <- function(expTypeLevel) {
    .assertScalar(x = expTypeLevel, type = "character",
                  validValues = c("MaxQuant", "ProteomeDiscovererProteins",
                                  "ProteomeDiscovererPeptideGroups",
                                  "FragPipe", "DIANN", "Spectronaut"))

    if (expTypeLevel == "MaxQuant") {
        paste0("The input to this workflow is a `proteinGroups.txt` file ",
               "from MaxQuant (see path in the table above). We read the MaxQuant ",
               "intensities into `R` and store them in a ",
               "[SingleCellExperiment]",
               "(https://bioconductor.org/packages/SingleCellExperiment/) object. ")
    } else if (expTypeLevel == "ProteomeDiscovererProteins") {
        paste0("The input to this workflow is a `Proteins.txt` file from ",
               "Proteome Discoverer (see path in the table above). We read the PD ",
               "intensities into `R` and store them in a ",
               "[SingleCellExperiment]" ,
               "(https://bioconductor.org/packages/SingleCellExperiment/) object. ")
    } else if (expTypeLevel == "ProteomeDiscovererPeptideGroups") {
        paste0("The input to this workflow is a `PeptideGroups.txt` file from ",
               "Proteome Discoverer (see path in the table above). We read the PD ",
               "intensities into `R` and store them in a ",
               "[SingleCellExperiment]" ,
               "(https://bioconductor.org/packages/SingleCellExperiment/) object. ")
    } else if (expTypeLevel == "FragPipe") {
        paste0("The input to this workflow is a `combined_protein` file ",
               "from FragPipe (see path in the table above). We read the FragPipe " ,
               "intensities into `R` and store them in a ",
               "[SingleCellExperiment]",
               "(https://bioconductor.org/packages/SingleCellExperiment/) object. ")
    } else if (expTypeLevel == "DIANN") {
        paste0("The input to this workflow is a `pg_matrix.tsv`, ",
               "`pr_matrix.tsv` or `report.tsv` file from DIA-NN (see path in the ",
               "table above). We read the DIA-NN intensities into `R` and store ",
               "them in a [SingleCellExperiment]",
               "(https://bioconductor.org/packages/SingleCellExperiment/) object. ")
    } else if (expTypeLevel == "Spectronaut") {
        paste0("The input to this workflow is a `Report.tsv` file from ",
               "Spectronaut (see path in the table above). We read the ",
               "Spectronaut intensities into `R` and store them in a ",
               "[SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/) ",
               "object. ")
    }
}

#' @rdname textSnippets
#' @export
featureCollectionText <- function(featureCollections) {
    if (length(featureCollections) == 0) {
        out <- paste0("No feature collections were tested.")
    } else {
        out <- ""
        if ("complexes" %in% featureCollections) {
            out <- paste0(out, "Complexes are obtained from Corum, CYC2008, ",
                          "PomBase, HuMAP2 and the Complex Portal (see ",
                          "the `makeComplexDB` function for more details). ")
        }
        if ("GO" %in% featureCollections) {
            out <- paste0(out, "GO terms are obtained from MSigDB. ")
        }
        if ("pathways" %in% featureCollections) {
            out <- paste0(out, "Pathway information is obtained from ",
                          "BIOCARTA, KEGG, PID, REACTOME and WIKIPATHWAYS, ",
                          "via MSigDB (see the `prepareFeatureCollections` ",
                          "function for more details). ")
        }
    }
    out
}
