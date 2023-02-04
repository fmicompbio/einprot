test_that("config chunk generation works", {
    expect_error(.generateConfigChunk(1),
                 "'configlist' must be of class 'list'")
    expect_error(.generateConfigChunk(list(1, 2)),
                 "'namesconfiglist' must not be NULL")
    expect_error(.generateConfigChunk(list(a = 1), fcn = 1),
                 "'fcn' must be of class 'character'")
    expect_error(.generateConfigChunk(list(a = 1), c("dump", "deparse")),
                 "'fcn' must have length 1")
    expect_error(.generateConfigChunk(list(a = 1), fcn = "missing"),
                 "All values in 'fcn' must be one of")

    ## deparse
    chk <- .generateConfigChunk(list(aa = 1, bb = 1:2, cc = list(3, 4),
                                     dd = data.frame(a = 1:3, k = 4:6)),
                                fcn = "deparse")
    expect_type(chk, "character")
    expect_equal(length(chk), 1)
    expect_equal(chk, "```{r config, eval = TRUE}\n## The following variables were specified as input arguments when calling the rendering function.\n## They will be used in the workflow below.\n\naa <- 1\nbb <- 1:2\ncc <- list(3, 4)\ndd <- structure(list(a = 1:3, k = 4:6), class = \"data.frame\", row.names = c(NA, -3L))\n```")

    ## dump
    chk <- .generateConfigChunk(list(aa = 1, bb = 1:2, cc = list(3, 4),
                                     dd = data.frame(a = 1:3, k = 4:6)),
                                fcn = "dump")
    expect_type(chk, "character")
    expect_equal(length(chk), 1)
    expect_equal(chk, "```{r config, eval = TRUE}\n## The following variables were specified as input arguments when calling the rendering function.\n## They will be used in the workflow below.\n\n\naa <-\n1\nbb <-\n1:2\ncc <-\nlist(3, 4)\ndd <-\nstructure(list(a = 1:3, k = 4:6), class = \"data.frame\", row.names = c(NA, \n-3L))\n```")

})
