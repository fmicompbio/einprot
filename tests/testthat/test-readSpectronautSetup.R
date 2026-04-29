test_that("reading Spectronaut setup file works", {
    ## Mis-specified arguments
    expect_error(readSpectronautSetup(1),
                 "'spectronautSetupFile' must be of class 'character'")
    expect_error(readSpectronautSetup("missing"), "doesn't exist")

    ## Read proper file and check that output is a
    ## correctly named list of the right length
    setupfile <- system.file("extdata", "spectronaut_example",
                             "PG.QValRUN_005_PGPEP_075_Minlog2PrecInt_0_reduced_pivot_filtered.setup.txt",
                             package = "einprot")
    sn <- readSpectronautSetup(setupfile)
    expect_type(sn, "list")
    expect_equal(length(sn), 4)
    expect_named(sn, c("Spectronaut version", "Setup file", "Raw files",
                       "Databases"))

    ## All entries should be scalar values
    expect_true(all(vapply(sn, length, 0) == 1))

    ## Check individual values
    expect_equal(sn$`Spectronaut version`, "Spectronaut 20.4.260109.92449")
    expect_equal(sn$`Setup file`, setupfile)
    expect_equal(sn$`Raw files`,
                 "LFQ_24min_500ng_2Th_3ms_A_01.raw, LFQ_24min_500ng_2Th_3ms_A_02.raw, LFQ_24min_500ng_2Th_3ms_A_03.raw, LFQ_24min_500ng_2Th_3ms_B_01.raw, LFQ_24min_500ng_2Th_3ms_B_02.raw, LFQ_24min_500ng_2Th_3ms_B_03.raw")
    expect_equal(sn$`Databases`,
                 "ECOLI__250312.fasta, HUMAN__230321.fasta, YEAST__230321.fasta")

    ## Missing file
    sn <- readSpectronautSetup(spectronautSetupFile = NULL)
    expect_type(sn, "list")
    expect_equal(length(sn), 0)
})
