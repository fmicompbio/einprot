#' Compile a list of comparisons to make
#'
#' @param allGroups Character vector containing all group labels in the
#'     dataset.
#' @param comparisons List of character vectors defining comparisons to
#'     perform. The first element of each vector represents the
#'     denominator of the comparison. If not empty, \code{ctrlGroup} and
#'     \code{allPairwiseComparisons} are ignored.
#' @param ctrlGroup Character scalar defining the sample group to use as
#'     control group in comparisons.
#' @param allPairwiseComparisons Logical scalar, should all pairwise
#'     comparisons be performed?
#'
#' @export
#' @author Charlotte Soneson
#'
#' @return A list of comparisons to perform.
#'
#' @examples
#' ## Perform all pairwise comparisons
#' makeListOfComparisons(allGroups = c("g1", "g2", "g3"),
#'                       comparisons = list(),
#'                       allPairwiseComparisons = TRUE,
#'                       ctrlGroup = "g1")
#'
#' ## Pre-specify the comparisons
#' makeListOfComparisons(allGroups = c("g1", "g2", "g3"),
#'                       comparisons = list(c("g1", "g3")),
#'                       allPairwiseComparisons = TRUE,
#'                       ctrlGroup = "g1")
#'
#' @importFrom utils combn
#'
makeListOfComparisons <- function(allGroups, comparisons,
                                  allPairwiseComparisons, ctrlGroup) {
    .assertVector(x = allGroups, type = "character")
    .assertVector(x = comparisons, type = "list")
    .assertScalar(x = allPairwiseComparisons, type = "logical",
                  allowNULL = TRUE)
    .assertScalar(x = ctrlGroup, type = "character", allowNULL = TRUE)

    if (length(comparisons) > 0) {
        ## Pre-specified comparisons - just check that they are valid
        for (cmp in comparisons) {
            if (!all(cmp %in% allGroups)) {
                stop("Misspecified 'comparisons', the following group(s) ",
                     "is/are not present: ", paste(setdiff(cmp, allGroups),
                                                   collapse = ", "))
            }
        }
    } else {
        ## Need to generate the comparisons list
        if ((allPairwiseComparisons &&
             !all(setdiff(ctrlGroup, "") %in% allGroups)) ||
            (!allPairwiseComparisons && !all(ctrlGroup %in% allGroups))) {
            stop("Misspecified 'ctrlGroup'")
        }

        if (allPairwiseComparisons) {
            comparisons <- utils::combn(allGroups, 2, simplify = FALSE)
            comparisons <- lapply(comparisons, function(cmp) {
                if (ctrlGroup %in% cmp && ctrlGroup != cmp[1]) {
                    rev(cmp)
                } else {
                    cmp
                }
            })
        } else {
            ## Get groups to compare to the control group
            other_groups <- setdiff(allGroups, ctrlGroup)
            comparisons <- lapply(as.list(other_groups),
                                  function(x) c(ctrlGroup, x))
        }
    }

    comparisons
}
