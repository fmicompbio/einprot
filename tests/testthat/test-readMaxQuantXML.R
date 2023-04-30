test_that("reading MQ XML file works", {
    ## Mis-specified arguments
    expect_error(readMaxQuantXML(1),
                 "'mqParameterFile' must be of class 'character'")
    expect_error(readMaxQuantXML("missing"), "doesn't exist")
    expect_error(readMaxQuantXML(system.file("extdata", "mq_example",
                                             "1356_proteinGroups.txt",
                                             package = "einprot")),
                 "'<' not found")

    ## Read proper file and check that output is a
    ## correctly named list of length 17
    mq <- readMaxQuantXML(
        system.file("extdata", "mq_example", "1356_mqpar.xml",
                    package = "einprot"))
    expect_type(mq, "list")
    expect_equal(length(mq), 17)
    expect_named(mq, c("MaxQuant version", "Parameter file", "Search engine",
                       "Raw file location", "Raw files", "Sample names",
                       "Databases", "Contaminants", "Quantification mode",
                       "Quantification settings (LFQ)", "Min. razor peptides",
                       "Requantify", "Enzymes", "Variable modifications",
                       "Fixed modifications", "Max peptide mass",
                       "Min peptide length"))

    ## All entries should be scalar values
    expect_true(all(vapply(mq, length, 0) == 1))

    ## Check individual values
    expect_equal(mq$`MaxQuant version`, "1.5.3.8")
    expect_equal(mq$`Sample names`, "Adnp_IP04, Adnp_IP05, Adnp_IP06, Chd4BF_IP07, Chd4BF_IP08, Chd4BF_IP09, RBC_ctrl_IP01, RBC_ctrl_IP02, RBC_ctrl_IP03")
    expect_equal(mq$`Quantification mode`, "1")
    expect_equal(mq$`Quantification settings (LFQ)`,
                 "LFQ min. ratio count: 1, fastLFQ: false, match-between runs (MBR): true, Intensity based absolute quantification (iBAQ): true")
    expect_equal(mq$`Min. razor peptides`, "1")
    expect_equal(mq$Requantify, "false")
    expect_equal(mq$Enzymes, "Trypsin/P")
    expect_equal(mq$`Variable modifications`, "Oxidation (M), Acetyl (Protein N-term)")
    expect_equal(mq$`Fixed modifications`, "")
    expect_equal(mq$`Max peptide mass`, "8000")
    expect_equal(mq$`Min peptide length`, "7")

    ## Missing file
    mq <- readMaxQuantXML(mqParameterFile = NULL)
    expect_type(mq, "list")
    expect_equal(length(mq), 0)
})
