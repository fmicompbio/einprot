test_that("makeListOfComparisons processes input argument correctly", {
    expect_error(makeListOfComparisons(
        allGroups = 1, comparisons = list(c("g1", "g2")),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2"),
        "'allGroups' must be of class 'character'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = c("g1", "g2"),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2"),
        "'comparisons' must be of class 'list'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2", "g1")),
        allPairwiseComparisons = TRUE, ctrlGroup = "g2"),
        "All values in 'vapplycomparisonslength0' must be one of")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = NULL,
        allPairwiseComparisons = TRUE, ctrlGroup = "g2"),
        "'comparisons' must not be NULL")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        allPairwiseComparisons = 1, ctrlGroup = "g2"),
        "'allPairwiseComparisons' must be of class 'logical'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        allPairwiseComparisons = c(TRUE, FALSE), ctrlGroup = "g2"),
        "'allPairwiseComparisons' must have length 1")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g1", "g2")),
        allPairwiseComparisons = TRUE, ctrlGroup = TRUE),
        "'ctrlGroup' must be of class 'character'")

    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(),
        allPairwiseComparisons = TRUE, ctrlGroup = "g4"),
        "Misspecified 'ctrlGroup'")
    expect_error(makeListOfComparisons(
        allGroups = c("g1", "g2"), comparisons = list(c("g3", "g2")),
        allPairwiseComparisons = TRUE, ctrlGroup = "g4"),
        "Misspecified 'comparisons'")
})

test_that("makeListOfComparisons works as expected", {
    allGroups <- c("g1", "g2", "g3", "g4")
    c1 <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                allPairwiseComparisons = TRUE,
                                ctrlGroup = "g1")
    expect_type(c1, "list")
    expect_equal(length(c1), 6)
    expect_equal(c1[[1]], c("g1", "g2"))
    expect_equal(c1[[2]], c("g1", "g3"))
    expect_equal(c1[[3]], c("g1", "g4"))
    expect_equal(c1[[4]], c("g2", "g3"))
    expect_equal(c1[[5]], c("g2", "g4"))
    expect_equal(c1[[6]], c("g3", "g4"))

    c2 <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                allPairwiseComparisons = TRUE,
                                ctrlGroup = "g3")
    expect_type(c2, "list")
    expect_equal(length(c2), 6)
    expect_equal(c2[[1]], c("g1", "g2"))
    expect_equal(c2[[2]], c("g3", "g1"))
    expect_equal(c2[[3]], c("g1", "g4"))
    expect_equal(c2[[4]], c("g3", "g2"))
    expect_equal(c2[[5]], c("g2", "g4"))
    expect_equal(c2[[6]], c("g3", "g4"))

    c3 <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                allPairwiseComparisons = FALSE,
                                ctrlGroup = "g3")
    expect_type(c3, "list")
    expect_equal(length(c3), 3)
    expect_equal(c3[[1]], c("g3", "g1"))
    expect_equal(c3[[2]], c("g3", "g2"))
    expect_equal(c3[[3]], c("g3", "g4"))

    c4 <- makeListOfComparisons(allGroups = allGroups,
                                comparisons = list(c("g1", "g2"),
                                                   c("g4", "g3"),
                                                   c("g3", "g1")),
                                allPairwiseComparisons = FALSE,
                                ctrlGroup = "g3")
    expect_type(c4, "list")
    expect_equal(length(c4), 3)
    expect_equal(c4[[1]], c("g1", "g2"))
    expect_equal(c4[[2]], c("g4", "g3"))
    expect_equal(c4[[3]], c("g3", "g1"))

    c5 <- makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                allPairwiseComparisons = TRUE,
                                ctrlGroup = "")
    expect_type(c5, "list")
    expect_equal(length(c5), 6)
    expect_equal(c5[[1]], c("g1", "g2"))
    expect_equal(c5[[2]], c("g1", "g3"))
    expect_equal(c5[[3]], c("g1", "g4"))
    expect_equal(c5[[4]], c("g2", "g3"))
    expect_equal(c5[[5]], c("g2", "g4"))
    expect_equal(c5[[6]], c("g3", "g4"))

    expect_error(makeListOfComparisons(allGroups = allGroups, comparisons = list(),
                                       allPairwiseComparisons = FALSE,
                                       ctrlGroup = ""),
                 "Misspecified 'ctrlGroup'")
})
