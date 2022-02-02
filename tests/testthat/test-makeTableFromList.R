test_that("makeTableFromList works", {
    expect_error(makeTableFromList(l = 1))
    expect_error(makeTableFromList(l = "list"))
    expect_error(makeTableFromList(l = list(1, 2)))
    expect_error(makeTableFromList(l = list(x = c(1, 2), y = 3)))

    expect_s3_class(makeTableFromList(list(first = "one",
                                           second = "two", third = 3)),
                    "kableExtra")
})
