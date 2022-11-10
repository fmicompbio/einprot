test_that("makeListOfComparisons processes input argument correctly", {
    expect_error(makeListOfComparisons(
        allGroups = 1, comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "'allGroups' must be of class 'character'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = c("g1", "g2"),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "'comparisons' must be of class 'list'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2", "g1")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "All values in 'vapplycomparisonslength0' must be one of")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = NULL,
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "'comparisons' must not be NULL")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = 1,
        allPairwiseComparisons = FALSE, ctrlGroup = "g2", discardGroup = NULL),
        "'mergeGroups' must be of class 'list'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = NULL,
        allPairwiseComparisons = FALSE, ctrlGroup = "g2", discardGroup = NULL),
        "'mergeGroups' must not be NULL")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(1),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "'namesmergeGroups' must not be NULL")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(g12 = "g1", g12 = "g2"),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "'mergeGroups' must be a named list")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(g12 = c("g1", "missing")),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = NULL),
        "Unknown group included in mergeGroups")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = 1, ctrlGroup = "g2", discardGroup = NULL),
        "'allPairwiseComparisons' must be of class 'logical'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = c(TRUE, FALSE), ctrlGroup = "g2",
        discardGroup = NULL),
        "'allPairwiseComparisons' must have length 1")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = TRUE, discardGroup = NULL),
        "'ctrlGroup' must be of class 'character'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = NULL, discardGroup = NULL),
        "'ctrlGroup' must not be NULL")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2", discardGroup = 1),
        "'discardGroup' must be of class 'character'")

    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g4", discardGroup = NULL),
        "Misspecified 'ctrlGroup'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g3", "g2")),
        mergeGroups = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g4", discardGroup = NULL),
        "Misspecified 'comparisons'")
})

test_that("makeListOfComparisons works as expected", {
    allGroups <- c("g1", "g2", "g3", "g4")
    c1out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "g1", discardGroup = NULL)
    expect_type(c1out, "list")
    expect_equal(length(c1out), 2)
    c1 <- c1out$comparisons
    c1gc <- c1out$groupComposition
    expect_type(c1, "list")
    expect_type(c1gc, "list")
    expect_equal(length(c1), 6)
    expect_equal(length(c1gc), 4)
    expect_named(c1, c("g2_vs_g1", "g3_vs_g1", "g4_vs_g1",
                       "g3_vs_g2", "g4_vs_g2", "g4_vs_g3"))
    expect_named(c1gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c1[[1]], c("g1", "g2"))
    expect_equal(c1[[2]], c("g1", "g3"))
    expect_equal(c1[[3]], c("g1", "g4"))
    expect_equal(c1[[4]], c("g2", "g3"))
    expect_equal(c1[[5]], c("g2", "g4"))
    expect_equal(c1[[6]], c("g3", "g4"))
    expect_equal(c1gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Change control group
    c2out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "g3", discardGroup = NULL)
    expect_type(c2out, "list")
    expect_equal(length(c2out), 2)
    c2 <- c2out$comparisons
    c2gc <- c2out$groupComposition
    expect_type(c2, "list")
    expect_type(c2gc, "list")
    expect_equal(length(c2), 6)
    expect_equal(length(c2gc), 4)
    expect_named(c2, c("g2_vs_g1", "g1_vs_g3", "g4_vs_g1",
                       "g2_vs_g3", "g4_vs_g2", "g4_vs_g3"))
    expect_named(c2gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c2[[1]], c("g1", "g2"))
    expect_equal(c2[[2]], c("g3", "g1"))
    expect_equal(c2[[3]], c("g1", "g4"))
    expect_equal(c2[[4]], c("g3", "g2"))
    expect_equal(c2[[5]], c("g2", "g4"))
    expect_equal(c2[[6]], c("g3", "g4"))
    expect_equal(c2gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Only comparisons to control group
    c3out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = FALSE,
                                   ctrlGroup = "g3", discardGroup = NULL)
    expect_type(c3out, "list")
    expect_equal(length(c3out), 2)
    c3 <- c3out$comparisons
    c3gc <- c3out$groupComposition
    expect_type(c3, "list")
    expect_type(c3gc, "list")
    expect_equal(length(c3), 3)
    expect_equal(length(c3gc), 4)
    expect_named(c3, c("g1_vs_g3", "g2_vs_g3", "g4_vs_g3"))
    expect_named(c3gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c3[[1]], c("g3", "g1"))
    expect_equal(c3[[2]], c("g3", "g2"))
    expect_equal(c3[[3]], c("g3", "g4"))
    expect_equal(c3gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Specify comparisons
    c4out <- makeListOfComparisons(allGroups = allGroups,
                                   comparisons = list(c("g1", "g2"),
                                                      c("g4", "g3"),
                                                      c("g3", "g1")),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = FALSE,
                                   ctrlGroup = "g3", discardGroup = NULL)
    expect_type(c4out, "list")
    expect_equal(length(c4out), 2)
    c4 <- c4out$comparisons
    c4gc <- c4out$groupComposition
    expect_type(c4, "list")
    expect_type(c4gc, "list")
    expect_equal(length(c4), 3)
    expect_equal(length(c4gc), 4)
    expect_named(c4, c("g2_vs_g1", "g3_vs_g4", "g1_vs_g3"))
    expect_named(c4gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c4[[1]], c("g1", "g2"))
    expect_equal(c4[[2]], c("g4", "g3"))
    expect_equal(c4[[3]], c("g3", "g1"))
    expect_equal(c4gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Specify comparisons, and name some of them
    c4bout <- makeListOfComparisons(allGroups = allGroups,
                                    comparisons = list(firstcomp = c("g1", "g2"),
                                                       c("g4", "g3"),
                                                       thirdcomp = c("g3", "g1")),
                                    mergeGroups = list(),
                                    allPairwiseComparisons = FALSE,
                                    ctrlGroup = "g3", discardGroup = NULL)
    expect_type(c4bout, "list")
    expect_equal(length(c4bout), 2)
    c4b <- c4bout$comparisons
    c4bgc <- c4bout$groupComposition
    expect_type(c4b, "list")
    expect_type(c4bgc, "list")
    expect_equal(length(c4b), 3)
    expect_equal(length(c4bgc), 4)
    expect_named(c4b, c("firstcomp", "g3_vs_g4", "thirdcomp"))
    expect_named(c4bgc, c("g1", "g2", "g3", "g4"))
    expect_equal(c4b[[1]], c("g1", "g2"))
    expect_equal(c4b[[2]], c("g4", "g3"))
    expect_equal(c4b[[3]], c("g3", "g1"))
    expect_equal(c4bgc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Don't specify control group - fine since allPairwiseComparisons = TRUE
    c5out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "", discardGroup = NULL)
    expect_type(c5out, "list")
    expect_equal(length(c5out), 2)
    c5 <- c5out$comparisons
    c5gc <- c5out$groupComposition
    expect_type(c5, "list")
    expect_type(c5gc, "list")
    expect_equal(length(c5), 6)
    expect_equal(length(c5gc), 4)
    expect_named(c5, c("g2_vs_g1", "g3_vs_g1", "g4_vs_g1",
                       "g3_vs_g2", "g4_vs_g2", "g4_vs_g3"))
    expect_named(c5gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c5[[1]], c("g1", "g2"))
    expect_equal(c5[[2]], c("g1", "g3"))
    expect_equal(c5[[3]], c("g1", "g4"))
    expect_equal(c5[[4]], c("g2", "g3"))
    expect_equal(c5[[5]], c("g2", "g4"))
    expect_equal(c5[[6]], c("g3", "g4"))
    expect_equal(c5gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Discard comparisons with group g2
    c6out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "", discardGroup = "g2")
    expect_type(c6out, "list")
    expect_equal(length(c6out), 2)
    c6 <- c6out$comparisons
    c6gc <- c6out$groupComposition
    expect_type(c6, "list")
    expect_type(c6gc, "list")
    expect_equal(length(c6), 3)
    expect_equal(length(c6gc), 4)
    expect_named(c6, c("g3_vs_g1", "g4_vs_g1", "g4_vs_g3"))
    expect_named(c6gc, c("g1", "g2", "g3", "g4"))
    expect_equal(c6[[1]], c("g1", "g3"))
    expect_equal(c6[[2]], c("g1", "g4"))
    expect_equal(c6[[3]], c("g3", "g4"))
    expect_equal(c6gc, list(g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Don't specify control group - error since allPairwiseComparisons = FALSE
    expect_error(makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                       mergeGroups = list(),
                                       allPairwiseComparisons = FALSE,
                                       ctrlGroup = ""),
                 "Misspecified 'ctrlGroup'")

    ## Merge groups
    ## Don't specify control group - fine since allPairwiseComparisons = TRUE
    c7out <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                   mergeGroups = list(gr12 = c("g1", "g2")),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "", discardGroup = NULL)
    expect_type(c7out, "list")
    expect_equal(length(c7out), 2)
    c7 <- c7out$comparisons
    c7gc <- c7out$groupComposition
    expect_type(c7, "list")
    expect_type(c7gc, "list")
    expect_equal(length(c7), 8)
    expect_equal(length(c7gc), 5)
    expect_named(c7, c("g3_vs_gr12", "g4_vs_gr12", "g2_vs_g1", "g3_vs_g1",
                       "g4_vs_g1", "g3_vs_g2", "g4_vs_g2", "g4_vs_g3"))
    expect_named(c7gc, c("gr12", "g1", "g2", "g3", "g4"))
    expect_equal(c7[[1]], c("gr12", "g3"))
    expect_equal(c7[[2]], c("gr12", "g4"))
    expect_equal(c7[[3]], c("g1", "g2"))
    expect_equal(c7[[4]], c("g1", "g3"))
    expect_equal(c7[[5]], c("g1", "g4"))
    expect_equal(c7[[6]], c("g2", "g3"))
    expect_equal(c7[[7]], c("g2", "g4"))
    expect_equal(c7[[8]], c("g3", "g4"))
    expect_equal(c7gc, list(gr12 = c("g1", "g2"),
                            g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4"))

    ## Merge groups, with complement
    ## Don't specify control group - fine since allPairwiseComparisons = TRUE
    c8out <- makeListOfComparisons(allGroups = allGroups,
                                   comparisons = list(c("gr12", "complement")),
                                   mergeGroups = list(gr12 = c("g1", "g2")),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "", discardGroup = NULL)
    expect_type(c8out, "list")
    expect_equal(length(c8out), 2)
    c8 <- c8out$comparisons
    c8gc <- c8out$groupComposition
    expect_type(c8, "list")
    expect_type(c8gc, "list")
    expect_equal(length(c8), 1)
    expect_equal(length(c8gc), 6)
    expect_named(c8, c("gr12_complement_vs_gr12"))
    expect_named(c8gc, c("gr12", "g1", "g2", "g3", "g4",
                         "gr12_complement"))
    expect_equal(c8[[1]], c("gr12", "gr12_complement"))
    expect_equal(c8gc, list(gr12 = c("g1", "g2"),
                            g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4",
                            gr12_complement = c("g3", "g4")))

    ## Merge groups, with complement, discard group
    ## Don't specify control group - fine since allPairwiseComparisons = TRUE
    c9out <- makeListOfComparisons(allGroups = allGroups,
                                   comparisons = list(c("gr12", "complement"),
                                                      c("gr13", "complement"),
                                                      c("gr13", "g4")),
                                   mergeGroups = list(gr12 = c("g1", "g2"),
                                                      gr13 = c("g3", "g1")),
                                   allPairwiseComparisons = TRUE,
                                   ctrlGroup = "", discardGroup = "g2")
    expect_type(c9out, "list")
    expect_equal(length(c9out), 2)
    c9 <- c9out$comparisons
    c9gc <- c9out$groupComposition
    expect_type(c9, "list")
    expect_type(c9gc, "list")
    expect_equal(length(c9), 1)
    expect_equal(length(c9gc), 8)
    expect_named(c9, c("g4_vs_gr13"))
    expect_named(c9gc, c("gr12", "gr13", "g1", "g2", "g3", "g4",
                         "gr12_complement", "gr13_complement"))
    expect_equal(c9[[1]], c("gr13", "g4"))
    expect_equal(c9gc, list(gr12 = c("g1", "g2"), gr13 = c("g3", "g1"),
                            g1 = "g1", g2 = "g2", g3 = "g3", g4 = "g4",
                            gr12_complement = c("g3", "g4"),
                            gr13_complement = c("g2", "g4")))

})
