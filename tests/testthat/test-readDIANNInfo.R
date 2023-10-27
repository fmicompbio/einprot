test_that("reading DIA-NN log file works", {
    ## Mis-specified arguments
    expect_error(readDIANNInfo(1),
                 "'diannLog' must be of class 'character'")
    expect_error(readDIANNInfo("missing"), "doesn't exist")

    ## Read proper file and check that output is a
    ## correctly named list of length 17
    dn <- readDIANNInfo(
        system.file("extdata", "diann_example", "diann-output.log.txt",
                    package = "einprot"))
    expect_type(dn, "list")
    expect_equal(length(dn), 8)
    expect_named(dn, c("DIA-NN version", "DIA-NN log file", "MS2 Mass accuracy (ppm)",
                       "MS1 Mass accuracy (ppm)", "Input files", "Spectral library",
                       "Sequence databases", "DIA-NN command"))

    ## All entries should be scalar values
    expect_true(all(vapply(dn, length, 0) == 1))

    ## Check individual values
    expect_equal(dn$`DIA-NN version`, "DIA-NN 1.8.2 beta 8 ")
    expect_equal(dn$`DIA-NN log file`,
                 system.file("extdata", "diann_example", "diann-output.log.txt",
                             package = "einprot"))
    expect_equal(dn$`MS2 Mass accuracy (ppm)`, "15 ")
    expect_equal(dn$`MS1 Mass accuracy (ppm)`, "15 ")
    expect_equal(dn$`Spectral library`,
                 "/scratch/cpanse/PXD028735/fasta/uniprotkb_proteome_UP000000625_UP000002311_UP000005640_iRTkit_2023_07_04.2and3plus.predicted.speclib")
    expect_equal(dn$`Sequence databases`,
                 "uniprotkb_proteome_UP000000625_2023_07_04.fasta; uniprotkb_proteome_UP000005640_2023_07_04.fasta; uniprotkb_proteome_UP000002311_2023_07_04.fasta; iRTkit.fasta")

    ## Missing file
    dn <- readDIANNInfo(diannLog = NULL)
    expect_type(dn, "list")
    expect_equal(length(dn), 0)
})
