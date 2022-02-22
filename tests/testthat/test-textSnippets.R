test_that("text snippet generation works", {
    expect_type(limmaText(), "character")
    expect_equal(length(limmaText()), 1)
    expect_type(ttestText(), "character")
    expect_equal(length(ttestText()), 1)
})
