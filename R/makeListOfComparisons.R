#' @author Charlotte Soneson
#' @noRd
#' @keywords internal
#'
.assignNamesToComparisons <- function(comparisons) {
    for (i in seq_along(comparisons)) {
        if (is.null(names(comparisons)) ||
            (!is.null(names(comparisons)) &&
             (is.na(names(comparisons)[i]) || names(comparisons)[i] == ""))) {
            names(comparisons)[i] <- paste0(comparisons[[i]][2], "_vs_",
                                            comparisons[[i]][1])
        }
    }
    if (any(duplicated(names(comparisons)))) {
        stop("Duplicated comparison names not allowed: ",
             paste(names(comparisons)[duplicated(names(comparisons))],
                   collapse = ", "))
    }
    comparisons
}

#' Compile a list of comparisons
#'
#' Construct a list defining the pairwise group comparisons that will be
#' performed.
#'
#' @param allGroups Character vector containing all group labels in the
#'     dataset.
#' @param comparisons List of length-2 character vectors defining comparisons
#'     to perform. The first element of each vector represents the
#'     denominator of the comparison, the second the numerator.
#'     If \code{comparisons} is non-empty, \code{ctrlGroup} and
#'     \code{allPairwiseComparisons} are ignored. If one of the two elements
#'     in the vector is "complement", it will be created as the complement
#'     of the other group.
#' @param mergeGroups Named list defining groups to merge. Each entry of the
#'     list corresponds to a new group, and consists of a vector with the
#'     original group names (a subset of \code{allGroups}). Only comparisons
#'     contrasting non-overlapping sets of groups will be performed.
#' @param ctrlGroup Character scalar defining the sample group to use as
#'     control group in comparisons.
#' @param allPairwiseComparisons Logical scalar, whether all pairwise
#'     comparisons shall be performed.
#' @param discardGroup Character vector. Any comparison including any of
#'     these groups will be discarded.
#'
#' @export
#' @author Charlotte Soneson
#'
#' @returns A list with two elements:
#' \itemize{
#'  \item{\code{comparisons}}{ - a list of vectors of length 2, indicating
#'  the pairwise group comparisons.}
#'  \item{\code{groupComposition}}{ - a list specifying the composition of
#'  each group label used in \code{comparisons}, in terms of the original
#'  group names. This is used to keep track of the composition of merged
#'  groups.}
#' }
#'
#' @examples
#' ## Perform all pairwise comparisons
#' makeListOfComparisons(allGroups = c("g1", "g2", "g3"),
#'                       comparisons = list(), mergeGroups = list(),
#'                       allPairwiseComparisons = TRUE,
#'                       ctrlGroup = "g1")
#'
#' ## Pre-specify the comparisons
#' makeListOfComparisons(allGroups = c("g1", "g2", "g3"),
#'                       comparisons = list(c("g1", "g3")),
#'                       mergeGroups = list(),
#'                       allPairwiseComparisons = TRUE,
#'                       ctrlGroup = "g1")
#'
#' ## Compare each group to its complement
#' makeListOfComparisons(allGroups = c("g1", "g2", "g3"),
#'                       comparisons = lapply(c("g1", "g2", "g3"),
#'                                            function(g) c("complement", g)),
#'                       mergeGroups = list(),
#'                       allPairwiseComparisons = TRUE,
#'                       ctrlGroup = "")
#'
#'
#' @importFrom utils combn
#'
makeListOfComparisons <- function(allGroups, comparisons, mergeGroups = list(),
                                  allPairwiseComparisons, ctrlGroup,
                                  discardGroup = NULL) {
    .assertVector(x = allGroups, type = "character")
    .assertVector(x = comparisons, type = "list")
    if (length(comparisons) > 0) {
        .assertVector(x = vapply(comparisons, length, 0), type = "numeric",
                      validValues = 2)
    }
    .assertVector(x = mergeGroups, type = "list")
    if (length(mergeGroups) > 0) {
        .assertVector(x = names(mergeGroups), type = "character")
        if (is.null(names(mergeGroups)) || any(names(mergeGroups) == "") ||
            any(duplicated(names(mergeGroups)))) {
            stop("'mergeGroups' must be a named list, without duplicated names")
        }
    }
    .assertScalar(x = allPairwiseComparisons, type = "logical",
                  allowNULL = TRUE)
    .assertScalar(x = ctrlGroup, type = "character", allowNULL = FALSE)
    .assertVector(x = discardGroup, type = "character", allowNULL = TRUE)

    ## Check that mergeGroups only contains values from allGroups (and
    ## optionally 'complement')
    stdf <- setdiff(unlist(mergeGroups), c(allGroups, "complement"))
    if (length(stdf) > 0) {
        stop("Unknown group included in mergeGroups: ",
             paste(stdf, collapse = ", "))
    }

    ## Expand mergeGroups to have also single-group entries
    for (g in allGroups) {
        if (!(g %in% names(mergeGroups))) {
            mergeGroups[g] <- g
        }
    }

    if (length(comparisons) > 0) {
        ## Pre-specified comparisons - just check that they are valid
        for (j in seq_along(comparisons)) {
            cmp <- comparisons[[j]]
            ## First expand any "complement"
            if ("complement" %in% cmp) {
                i <- which(cmp == "complement")
                specGrName <- setdiff(cmp, "complement")
                newComplName <- paste0(specGrName, "_complement")
                cmp[i] <- newComplName
                mergeGroups[[newComplName]] <- setdiff(
                    allGroups,  mergeGroups[[specGrName]])
            }
            ## Check that the comparisons only use valid (merged) groups
            if (!all(cmp %in% names(mergeGroups))) {
                stop("Misspecified 'comparisons', the following group(s) ",
                     "is/are not present: ",
                     paste(setdiff(cmp, names(mergeGroups)),
                           collapse = ", "))
            }
            ## Check that the compared groups don't overlap
            if (length(intersect(mergeGroups[[cmp[1]]],
                                 mergeGroups[[cmp[2]]])) != 0) {
                stop("Invalid comparison, groups overlap: ",
                     cmp[[1]], " vs ", cmp[[2]])
            }
            comparisons[[j]] <- cmp
        }
    } else {
        ## Need to generate the comparisons list
        if ((allPairwiseComparisons &&
             !all(setdiff(ctrlGroup, "") %in% names(mergeGroups))) ||
            (!allPairwiseComparisons && !all(ctrlGroup %in%
                                             names(mergeGroups)))) {
            stop("Misspecified 'ctrlGroup'")
        }

        if (allPairwiseComparisons) {
            ## All pairwise comparisons, except overlapping groups
            comparisons <- utils::combn(names(mergeGroups), 2, simplify = FALSE)
            comparisons <- lapply(comparisons, function(cmp) {
                if (length(intersect(mergeGroups[[cmp[1]]],
                                     mergeGroups[[cmp[2]]])) != 0) {
                    NULL
                } else {
                    if (ctrlGroup %in% cmp && ctrlGroup != cmp[1]) {
                        rev(cmp)
                    } else {
                        cmp
                    }
                }
            })
        } else {
            ## Each group vs the control group
            ## Get groups to compare to the control group
            other_groups <- setdiff(names(mergeGroups), ctrlGroup)
            comparisons <- lapply(as.list(other_groups), function(x) {
                if (length(intersect(mergeGroups[[ctrlGroup]],
                                     mergeGroups[[x]])) != 0) {
                    NULL
                } else {
                    c(ctrlGroup, x)
                }
            })
        }
    }

    ## Remove NULL entries
    comparisons <- comparisons[vapply(comparisons, function(cps) {
        !is.null(cps)
    }, FALSE)]

    ## Remove comparisons containing the groups to discard
    comparisons <- comparisons[vapply(comparisons, function(cps) {
        !any(discardGroup %in% unlist(mergeGroups[cps]))
    }, FALSE)]

    comparisons <- .assignNamesToComparisons(comparisons)

    list(comparisons = comparisons, groupComposition = mergeGroups)
}
