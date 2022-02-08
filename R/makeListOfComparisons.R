#' @export
makeListOfComparisons <- function(allGroups, comparisons,
                                  allPairwiseComparisons, ctrlGroup) {
    ## Get list of comparisons to make
    if (length(comparisons) == 0) {
        ## No manually defined comparisons
        if (allPairwiseComparisons) {
            comparisons <- combn(allGroups, 2, simplify = FALSE)
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
            comparisons <- lapply(as.list(other_groups), function(x) c(ctrlGroup, x))
        }
    }
    comparisons
}
