test_that("readFragPipeInfo works", {
    ## Mis-specified arguments
    expect_error(readFragPipeInfo(1),
                 "'fragpipeDir' must be of class 'character'")

    ## Missing file -> only decoy tag returned
    expect_equal(readFragPipeInfo("missing"),
                 list(`Database decoy tag` = "rev_"))

    ## Read proper file and check that output is a
    ## correctly named list of length 15
    fp <- readFragPipeInfo(
        system.file("extdata", "fp_example",
                    package = "einprot"))
    expect_type(fp, "list")
    expect_equal(length(fp), 16)
    expect_named(fp, c("FragPipe version", "FragPipe parameter file",
                       "FragPipe log file", "Search engine",
                       "Raw file location", "Raw files", "Sample names",
                       "Databases", "Contaminants", "Peptides (ranges)",
                       "Mass error tolerances", "Quantification settings (LFQ)",
                       "Enzymes", "Variable modifications",
                       "Fixed modifications", "Database decoy tag"))

    ## All entries should be scalar values
    expect_true(all(vapply(fp, length, 0) == 1))

    ## Check individual values
    expect_equal(fp$`FragPipe version`, "19.1")
    expect_equal(basename(fp$`FragPipe parameter file`), "fragpipe.workflow")
    expect_equal(basename(fp$`FragPipe log file`), "log_2023-04-12_20-12-46.txt")
    expect_equal(fp$`Search engine`, "MSFragger-3.7")
    expect_equal(fp$`Raw file location`, "D:/Data/FUSION")
    expect_equal(fp$`Raw files`, "F_160817_AdnpFB_IP06.raw, F_160817_AdnpFB_IP05.raw, F_160817_RBC_ctrl_IP02.raw, F_160817_AdnpFB_IP04.raw, F_160817_RBC_ctrl_IP01.raw, F_160817_Chd4BF_IP09.raw, F_160817_Chd4BF_IP07.raw, F_160817_Chd4BF_IP08.raw, F_160817_RBC_ctrl_IP03.raw")
    expect_equal(fp$`Sample names`, "Adnp_IP04, Adnp_IP05, Adnp_IP06, Chd4BF_IP07, Chd4BF_IP08, Chd4BF_IP09, RBC_ctrl_IP01, RBC_ctrl_IP02, RBC_ctrl_IP03")
    expect_equal(fp$Databases, "D/://Data//FASTA//2023-04-12-decoys-contam_MOUSE__190410.fasta.fas")
    expect_equal(fp$Contaminants, "cRAP")
    expect_equal(fp$`Peptides (ranges)`, "length: 7-50 AA; mass: 500-5000 Da")
    expect_equal(fp$`Mass error tolerances`, "precursor:-20-20 [ppm]; fragment:0.7 [Da] (after optimization:200 PPM)")
    expect_equal(fp$`Quantification settings (LFQ)`, "IonQuant: TRUE, Calculate MaxLFQ intensity: TRUE, Normalization: TRUE, match-between runs (MBR): FALSE, min. ions: 2")
    expect_equal(fp$Enzymes, "stricttrypsin[KR, C-terminal, 2 missed cleavages]")
    expect_equal(fp$`Variable modifications`, "M(15.9949), N-term(42.0106)")
    expect_equal(fp$`Fixed modifications`, "C(57.0215)")
    expect_equal(fp$`Database decoy tag`, "rev_")

    ## Workflow file missing -> get values from log file
    ## -------------------------------------------------------------------------
    dir.create(file.path(tempdir(), "tempfp"))
    file.copy(from = system.file("extdata", "fp_example", "log_2023-04-12_20-12-46.txt",
                                 package = "einprot"),
              to = file.path(tempdir(), "tempfp"))
    fp <- readFragPipeInfo(file.path(tempdir(), "tempfp"))
    expect_type(fp, "list")
    expect_equal(length(fp), 15)
    expect_named(fp, c("FragPipe version",
                       "FragPipe log file", "Search engine",
                       "Raw file location", "Raw files", "Sample names",
                       "Databases", "Contaminants", "Peptides (ranges)",
                       "Mass error tolerances", "Quantification settings (LFQ)",
                       "Enzymes", "Variable modifications",
                       "Fixed modifications", "Database decoy tag"))

    ## All entries should be scalar values
    expect_true(all(vapply(fp, length, 0) == 1))

    ## Check individual values
    expect_equal(fp$`FragPipe version`, "FragPipe v19.1")
    expect_equal(basename(fp$`FragPipe log file`), "log_2023-04-12_20-12-46.txt")
    expect_equal(fp$`Search engine`, "MSFragger-3.7")
    expect_equal(fp$`Raw file location`, "D:/Data/FUSION")
    expect_equal(fp$`Raw files`, "F_160817_AdnpFB_IP06.raw, F_160817_AdnpFB_IP05.raw, F_160817_RBC_ctrl_IP02.raw, F_160817_AdnpFB_IP04.raw, F_160817_RBC_ctrl_IP01.raw, F_160817_Chd4BF_IP09.raw, F_160817_Chd4BF_IP07.raw, F_160817_Chd4BF_IP08.raw, F_160817_RBC_ctrl_IP03.raw")
    expect_equal(fp$`Sample names`, "Adnp_IP04, Adnp_IP05, Adnp_IP06, Chd4BF_IP07, Chd4BF_IP08, Chd4BF_IP09, RBC_ctrl_IP01, RBC_ctrl_IP02, RBC_ctrl_IP03")
    expect_equal(fp$Databases, "D:/Data//FASTA/2023-04-12-decoys-contam_MOUSE__190410.fasta.fas")
    expect_equal(fp$Contaminants, "cRAP")
    expect_equal(fp$`Peptides (ranges)`, "length: 7-50 AA; mass: 500-5000 Da")
    expect_equal(fp$`Mass error tolerances`, "precursor:-20-20 [ppm]; fragment:0.7 [Da] (after optimization:200 PPM)")
    expect_equal(fp$`Quantification settings (LFQ)`, "IonQuant: TRUE, Calculate MaxLFQ intensity: TRUE, Normalization: TRUE, match-between runs (MBR): FALSE, min. ions: 2")
    expect_equal(fp$Enzymes, "stricttrypsin[KR, C-terminal, 2 missed cleavages]")
    expect_equal(fp$`Variable modifications`, "M(15.9949), N-term(42.0106)")
    expect_equal(fp$`Fixed modifications`, "C(57.0215)")
    expect_equal(fp$`Database decoy tag`, "rev_")

    ## -------------------------------------------------------------------------
    ## Create file where msfragger.search_enzyme_name_2 is not null
    ## Copy also workflow file to folder above
    file.copy(from = system.file("extdata", "fp_example", "fragpipe.workflow",
                                 package = "einprot"),
              to = file.path(tempdir(), "tempfp"))
    tmp <- readLines(file.path(tempdir(), "tempfp", "fragpipe.workflow"))
    i <- grep("msfragger.search_enzyme_name_2", tmp)
    expect_equal(i, 142L)
    tmp[i] <- "msfragger.search_enzyme_name_2=trypsin"
    writeLines(tmp, file.path(tempdir(), "tempfp", "fragpipe.workflow"))
    fp <- readFragPipeInfo(file.path(tempdir(), "tempfp"))
    expect_type(fp, "list")
    expect_equal(length(fp), 16)
    expect_named(fp, c("FragPipe version", "FragPipe parameter file",
                       "FragPipe log file", "Search engine",
                       "Raw file location", "Raw files", "Sample names",
                       "Databases", "Contaminants", "Peptides (ranges)",
                       "Mass error tolerances", "Quantification settings (LFQ)",
                       "Enzymes", "Variable modifications",
                       "Fixed modifications", "Database decoy tag"))

    ## All entries should be scalar values
    expect_true(all(vapply(fp, length, 0) == 1))

    ## Check individual values
    expect_equal(fp$`FragPipe version`, "19.1")
    expect_equal(basename(fp$`FragPipe parameter file`), "fragpipe.workflow")
    expect_equal(basename(fp$`FragPipe log file`), "log_2023-04-12_20-12-46.txt")
    expect_equal(fp$`Search engine`, "MSFragger-3.7")
    expect_equal(fp$`Raw file location`, "D:/Data/FUSION")
    expect_equal(fp$`Raw files`, "F_160817_AdnpFB_IP06.raw, F_160817_AdnpFB_IP05.raw, F_160817_RBC_ctrl_IP02.raw, F_160817_AdnpFB_IP04.raw, F_160817_RBC_ctrl_IP01.raw, F_160817_Chd4BF_IP09.raw, F_160817_Chd4BF_IP07.raw, F_160817_Chd4BF_IP08.raw, F_160817_RBC_ctrl_IP03.raw")
    expect_equal(fp$`Sample names`, "Adnp_IP04, Adnp_IP05, Adnp_IP06, Chd4BF_IP07, Chd4BF_IP08, Chd4BF_IP09, RBC_ctrl_IP01, RBC_ctrl_IP02, RBC_ctrl_IP03")
    expect_equal(fp$Databases, "D/://Data//FASTA//2023-04-12-decoys-contam_MOUSE__190410.fasta.fas")
    expect_equal(fp$Contaminants, "cRAP")
    expect_equal(fp$`Peptides (ranges)`, "length: 7-50 AA; mass: 500-5000 Da")
    expect_equal(fp$`Mass error tolerances`, "precursor:-20-20 [ppm]; fragment:0.7 [Da] (after optimization:200 PPM)")
    expect_equal(fp$`Quantification settings (LFQ)`, "IonQuant: TRUE, Calculate MaxLFQ intensity: TRUE, Normalization: TRUE, match-between runs (MBR): FALSE, min. ions: 2")
    expect_equal(fp$Enzymes, "stricttrypsin[KR, C-terminal, 2 missed cleavages]; trypsin; [, C-terminal, 2 missed cleavages]")
    expect_equal(fp$`Variable modifications`, "M(15.9949), N-term(42.0106)")
    expect_equal(fp$`Fixed modifications`, "C(57.0215)")
    expect_equal(fp$`Database decoy tag`, "rev_")

})
