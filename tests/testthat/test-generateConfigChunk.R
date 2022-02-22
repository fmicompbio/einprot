test_that("config chunk generation works", {
    expect_error(.generateConfigChunk(1), "must be of class 'list'")
    expect_error(.generateConfigChunk(list(1, 2)), "must not be NULL")

    chk <- .generateConfigChunk(list(aa = 1, bb = 1:2, cc = list(3, 4)))
    expect_type(chk, "character")
    expect_equal(length(chk), 1)
    expect_equal(chk, "```{r config, eval = TRUE}\n## The following variables were specified as input arguments when calling the rendering function.\n## They will be used in the workflow below.\n\naa <- 1\nbb <- 1:2\ncc <- list(3, 4)\n```")
})
