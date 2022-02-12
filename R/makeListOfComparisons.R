#' @export
makeListOfComparisons <- function(allGroups, comparisons,
                                  allPairwiseComparisons, ctrlGroup) {
    if (length(comparisons) == 0) {
        if ((allPairwiseComparisons && !all(setdiff(ctrlGroup, "") %in% allGroups)) ||
            (!allPairwiseComparisons && !all(ctrlGroup %in% allGroups))) {
            stop("Misspecified 'ctrlGroup'")
        }
    } else {
        ## Check that all comparisons specify valid group names
        for (cmp in comparisons) {
            if (!all(cmp %in% allGroups)) {
                stop("Misspecified 'comparisons', the following group(s) is/are ",
                     "not present: ", paste(setdiff(cmp, allGroups), collapse = ", "))
            }
        }
    }

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
