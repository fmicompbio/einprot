test_that("makeTableFromList works", {
    expect_error(makeTableFromList(l = 1),
                 "'l' must be of class 'list'")
    expect_error(makeTableFromList(l = "list"),
                 "'l' must be of class 'list'")
    expect_error(makeTableFromList(l = list(1, 2)),
                 "'namesl' must not be NULL")
    expect_error(makeTableFromList(l = list(x = c(1, 2), y = 3)),
                 "All values in 'vapplyllength0' must be one of")

    expect_s3_class(makeTableFromList(list(first = "one",
                                           second = "two", third = 3)),
                    "kableExtra")
    expect_s3_class(makeTableFromList(list()),
                    "kableExtra")

    a <- makeTableFromList(list(first = "abc$; def$; efg"))
    expect_s3_class(a, "kableExtra")
})
