test_that("text snippet generation works", {
    ## testText
    expect_error(testText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(testText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(testText(testType = "missing"),
                 "All values in 'testType' must be one of")
    expect_error(testText(testType = "limma", minlFC = "1"),
                 "'minlFC' must be of class 'numeric'")
    expect_error(testText(testType = "limma", minlFC = c(1, 2)),
                 "'minlFC' must have length 1")
    expect_error(testText(testType = "limma", samSignificance = "1"),
                 "'samSignificance' must be of class 'logical'")
    expect_error(testText(testType = "limma", samSignificance = c(TRUE, FALSE)),
                 "'samSignificance' must have length 1")

    expect_type(testText(testType = "limma", minlFC = 1), "character")
    expect_equal(length(testText(testType = "limma")), 1)
    expect_true(grepl("the treat function", testText(testType = "limma", minlFC = 1)))

    expect_type(testText(testType = "limma", minlFC = 0), "character")
    expect_equal(length(testText(testType = "limma")), 1)
    expect_false(grepl("the treat function", testText(testType = "limma", minlFC = 0)))
    expect_true(grepl("limma", testText(testType = "limma", minlFC = 0)))

    expect_type(testText(testType = "ttest", samSignificance = TRUE), "character")
    expect_equal(length(testText(testType = "ttest", samSignificance = TRUE)), 1)
    expect_true(grepl("a Student's t-test", testText(testType = "ttest",
                                                     samSignificance = TRUE)))
    expect_true(grepl("Tusher", testText(testType = "ttest",
                                         samSignificance = TRUE)))
    expect_type(testText(testType = "ttest", samSignificance = FALSE), "character")
    expect_equal(length(testText(testType = "ttest", samSignificance = FALSE)), 1)
    expect_true(grepl("a Student's t-test", testText(testType = "ttest",
                                                     samSignificance = FALSE)))
    expect_false(grepl("Tusher", testText(testType = "ttest",
                                          samSignificance = FALSE)))

    expect_type(testText(testType = "proDA"), "character")
    expect_equal(length(testText(testType = "proDA")), 1)
    expect_true(grepl("proDA", testText(testType = "proDA")))

    expect_type(testText(testType = "none"), "character")
    expect_equal(length(testText(testType = "none")), 1)
    expect_true(grepl("no statistical testing", testText(testType = "none")))

    ## normText
    expect_error(normText(normMethod = 1),
                 "'normMethod' must be of class 'character'")
    expect_error(normText(normMethod = c("limma", "ttest")),
                 "'normMethod' must have length 1")

    expect_type(normText(normMethod = "none"), "character")
    expect_equal(length(normText(normMethod = "none")), 1)
    expect_true(grepl("are not normalized", normText(normMethod = "none")))
    expect_type(normText(normMethod = "center.median"), "character")
    expect_equal(length(normText(normMethod = "center.median")), 1)
    expect_true(grepl("using the center.median method",
                      normText(normMethod = "center.median")))

    ## saText
    expect_error(saText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(saText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(saText(testType = "missing"),
                 "All values in 'testType' must be one of")

    expect_type(saText(testType = "limma"), "character")
    expect_equal(length(saText(testType = "limma")), 1)
    expect_true(grepl("square root of the residual standard",
                      saText(testType = "limma")))

    expect_type(saText(testType = "ttest"), "character")
    expect_equal(length(saText(testType = "ttest")), 1)
    expect_equal(saText(testType = "ttest"), "")

    expect_type(saText(testType = "proDA"), "character")
    expect_equal(length(saText(testType = "proDA")), 1)
    expect_equal(saText(testType = "proDA"), "")

    expect_type(saText(testType = "none"), "character")
    expect_equal(length(saText(testType = "none")), 1)
    expect_equal(saText(testType = "none"), "")

    ## expDesignText
    expect_error(expDesignText(testType = 1),
                 "'testType' must be of class 'character'")
    expect_error(expDesignText(testType = c("limma", "ttest")),
                 "'testType' must have length 1")
    expect_error(expDesignText(testType = "missing"),
                 "All values in 'testType' must be one of")

    expect_type(expDesignText(testType = "limma"), "character")
    expect_equal(length(expDesignText(testType = "limma")), 1)
    expect_true(grepl("experimental design", expDesignText(testType = "limma")))

    expect_type(expDesignText(testType = "ttest"), "character")
    expect_equal(length(expDesignText(testType = "ttest")), 1)
    expect_equal(expDesignText(testType = "ttest"), "")

    expect_type(expDesignText(testType = "proDA"), "character")
    expect_equal(length(expDesignText(testType = "proDA")), 1)
    expect_equal(expDesignText(testType = "proDA"), "")

    expect_type(expDesignText(testType = "none"), "character")
    expect_equal(length(expDesignText(testType = "none")), 1)
    expect_equal(expDesignText(testType = "none"), "")
})
