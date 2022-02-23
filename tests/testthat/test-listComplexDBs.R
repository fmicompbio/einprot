test_that("listing complex dbs works", {
    expect_error(listComplexDBs(dbDir = 1),
                 "'dbDir' must be of class 'character'")
    expect_error(listComplexDBs(dbDir = c("path1", "path2")),
                 "'dbDir' must have length 1")

    ldb <- listComplexDBs()
    expect_s3_class(ldb, "data.frame")
    expect_named(ldb, c("complexDbPath", "genDate"))
    expect_true(all(substr(basename(ldb$complexDbPath), start = 1,
                           stop = 10) == "complexdb_"))
})
