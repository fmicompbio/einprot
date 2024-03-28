test_that("makeComplexDB works", {
    expect_error(makeComplexDB(dbDir = 1, customComplexTxt = NULL),
                 "'dbDir' must be of class 'character'")
    expect_error(makeComplexDB(dbDir = c("dir1", "dir2"),
                               customComplexTxt = NULL),
                 "'dbDir' must have length 1")
    expect_error(makeComplexDB(dbDir = "subdir", customComplexTxt = 1),
                 "'customComplexTxt' must be of class 'character'")
    expect_error(makeComplexDB(dbDir = "subdir",
                               customComplexTxt = c("txt1", "txt2")),
                 "'customComplexTxt' must have length 1")
    expect_error(makeComplexDB(dbDir = "subdir", customComplexTxt = "missing"),
                 "file.exists(customComplexTxt) is not TRUE", fixed = TRUE)

    ## Check that species exist in babelgene
    ## As of May 2022, neither C elegans nor S pombe has a 'common name'
    sps <- babelgene::species()
    expect_true(any(grepl("mouse", sps$common_name)))
    expect_true(any(grepl("baker's yeast", sps$common_name)))
    expect_true(any(grepl("Caenorhabditis elegans", sps$scientific_name)))
    expect_true(any(grepl("Schizosaccharomyces pombe 972h-", sps$scientific_name)))

    ## -------------------------------------------------------------------------
    ## Create complexDb from small example files
    ## -------------------------------------------------------------------------
    cyc2008db <- read.delim(system.file("extdata", "complexes",
                                        "cyc2008_complex_extract.tab",
                                        package = "einprot"))
    corumdb <- read.delim(system.file("extdata", "complexes",
                                      "corum_complex_extract.txt",
                                      package = "einprot"))
    pombasedb <- read.delim(system.file("extdata", "complexes",
                                        "pombase_complex_extract.tsv",
                                        package = "einprot"))
    humap2db <- read.delim(system.file("extdata", "complexes",
                                       "humap2_complex_extract.txt",
                                       package = "einprot"), sep = ",")
    cportal9606db <- read.delim(system.file("extdata", "complexes",
                                            "complexportal9606_complex_extract.tsv",
                                            package = "einprot"))
    cportal10090db <- read.delim(system.file("extdata", "complexes",
                                             "complexportal10090_complex_extract.tsv",
                                             package = "einprot"))
    cportal6239db <- read.delim(system.file("extdata", "complexes",
                                            "complexportal6239_complex_extract.tsv",
                                            package = "einprot"))
    cportal284812db <- read.delim(system.file("extdata", "complexes",
                                              "complexportal284812_complex_extract.tsv",
                                              package = "einprot"))

    dbdir <- tempdir()
    custompath <- file.path(dbdir, "custom_complex.txt")
    custom <- data.frame(Complex.name = "myComplex",
                         Gene.names = "HDAC4;HDAC5",
                         Organism = "Human",
                         Source = "Custom",
                         PMID = "123456")
    write.table(custom, file = custompath,
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    dbs <- makeComplexDB(dbDir = dbdir, customComplexTxt = custompath,
                         Cyc2008Db = cyc2008db, CorumDb = corumdb,
                         PombaseDb = pombasedb, HuMAP2Db = humap2db,
                         CPortal9606Db = cportal9606db,
                         CPortal10090Db = cportal10090db,
                         CPortal6239Db = cportal6239db,
                         CPortal284812Db = cportal284812db)
    compl <- readRDS(dbs$complPath)
    orth <- readRDS(dbs$orthPath)

    ## compl - no ortholog mapping
    ## -------------------------------------------------------------------------
    expect_s4_class(compl, "CharacterList")
    expect_equal(length(compl), 26L)
    expect_named(compl, c("human: BCL6-HDAC4 complex", "human: BCL6-HDAC5 complex",
                          "mouse: BLOC-2 (biogenesis of lysosome-related organelles complex 2)",
                          "rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                          "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                          "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                          "S.pombe: nucleotide-excision repair factor 1 complex",
                          "S.pombe: nucleotide-excision repair factor 2 complex",
                          "S.pombe: nucleotide-excision repair factor 3 complex",
                          "human: myComplex",
                          "human: HuMAP2_00000_conf3", "human: HuMAP2_00001_conf4",
                          "human: HuMAP2_00002_conf5",
                          "C.elegans: CPX-365_GARP tethering complex",
                          "C.elegans: CPX-368_Polycomb Repressive Complex 2",
                          "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                          "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                          "S.pombe: CPX-546_DNA replication factor C complex",
                          "S.pombe: CPX-547_PCNA homotrimer",
                          "mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex",
                          "mouse: CPX-6_bZIP transcription factor complex, Atf4-Creb1",
                          "mouse: CPX-7_bZIP transcription factor complex, Atf1-Atf4",
                          "human: CPX-1_SMAD2-SMAD3-SMAD4 complex",
                          "human: CPX-8_bZIP transcription factor complex, ATF4-CREB1",
                          "human: CPX-9_bZIP transcription factor complex, ATF1-ATF4"),
                 ignore.order = TRUE)
    expect_equal(compl$`human: BCL6-HDAC4 complex`, c("BCL6", "HDAC4"))
    expect_equal(mcols(compl)["human: BCL6-HDAC4 complex", "Species.common"],
                 "human")
    expect_equal(mcols(compl)["human: BCL6-HDAC4 complex", "Source"], "CORUM")
    expect_equal(mcols(compl)["human: BCL6-HDAC4 complex", "PMID"], "11929873")

    expect_equal(compl$`human: HuMAP2_00000_conf3`, c("FYCO1", "TRUB2"))
    expect_equal(mcols(compl)["human: HuMAP2_00000_conf3", "Species.common"],
                 "human")
    expect_equal(mcols(compl)["human: HuMAP2_00000_conf3", "Source"], "HuMAP2")
    expect_equal(mcols(compl)["human: HuMAP2_00000_conf3", "PMID"], "")

    expect_equal(compl$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`,
                 c("Dnm1l", "Clta", "Bcl2l1", "Mff"))
    expect_equal(mcols(compl)["rat: Bcl2l1-Dnm1l-Mff-Clta complex", "Species.common"],
                 "rat")
    expect_equal(mcols(compl)["rat: Bcl2l1-Dnm1l-Mff-Clta complex", "Source"], "CORUM")
    expect_equal(mcols(compl)["rat: Bcl2l1-Dnm1l-Mff-Clta complex", "PMID"], "23792689")

    expect_equal(compl$`human: myComplex`, c("HDAC4", "HDAC5"))

    expect_equal(compl$`S.cer: histone deacetylase complex`, c("HDA1", "HDA2", "HDA3"))
    expect_equal(mcols(compl)["S.cer: histone deacetylase complex", "Species.common"],
                 "baker's yeast")
    expect_equal(mcols(compl)["S.cer: histone deacetylase complex", "Source"], "CYC2008")
    expect_equal(mcols(compl)["S.cer: histone deacetylase complex", "PMID"], "11287668")

    expect_equal(compl$`S.pombe: nucleotide-excision repair factor 3 complex`,
                 c("tfb1", "ptr8", "rad15", "tfb2", "pmh1"))
    expect_equal(mcols(compl)["S.pombe: nucleotide-excision repair factor 3 complex",
                              "Species.common"], "Schizosaccharomyces pombe 972h-")
    expect_equal(mcols(compl)["S.pombe: nucleotide-excision repair factor 3 complex",
                              "Source"], "pombase")
    expect_equal(mcols(compl)["S.pombe: nucleotide-excision repair factor 3 complex",
                              "PMID"], "1534406;GO_REF:0000024")

    expect_equal(compl$`C.elegans: CPX-365_GARP tethering complex`,
                 c("vps-52", "vps-51", "vps-53", "vps-54"))
    expect_equal(mcols(compl)["C.elegans: CPX-365_GARP tethering complex",
                              "Species.common"], "Caenorhabditis elegans")
    expect_equal(mcols(compl)["C.elegans: CPX-365_GARP tethering complex",
                              "Source"], "ComplexPortal")
    expect_equal(mcols(compl)["C.elegans: CPX-365_GARP tethering complex",
                              "PMID"], "")

    expect_equal(compl$`human: CPX-8_bZIP transcription factor complex, ATF4-CREB1`,
                 c("CREB1", "ATF4"))
    expect_equal(mcols(compl)["human: CPX-8_bZIP transcription factor complex, ATF4-CREB1",
                              "Species.common"], "human")
    expect_equal(mcols(compl)["human: CPX-8_bZIP transcription factor complex, ATF4-CREB1",
                              "Source"], "ComplexPortal")
    expect_equal(mcols(compl)["human: CPX-8_bZIP transcription factor complex, ATF4-CREB1",
                              "PMID"], "")

    expect_equal(compl$`mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex`,
                 c("Smad4", "Smad2", "Smad3"))
    expect_equal(mcols(compl)["mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex",
                              "Species.common"], "mouse")
    expect_equal(mcols(compl)["mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex",
                              "Source"], "ComplexPortal")
    expect_equal(mcols(compl)["mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex",
                              "PMID"], "")

    expect_equal(compl$`S.pombe: CPX-546_DNA replication factor C complex`,
                 c("rfc3", "rfc1", "rfc4", "rfc5", "rfc2"))
    expect_equal(mcols(compl)["S.pombe: CPX-546_DNA replication factor C complex",
                              "Species.common"], "Schizosaccharomyces pombe 972h-")
    expect_equal(mcols(compl)["S.pombe: CPX-546_DNA replication factor C complex",
                              "Source"], "ComplexPortal")
    expect_equal(mcols(compl)["S.pombe: CPX-546_DNA replication factor C complex",
                              "PMID"], "")

    ## orth - ortholog mapping to all species
    ## -------------------------------------------------------------------------
    expect_type(orth, "list")
    expect_named(orth, c("mouse", "human", "baker's yeast",
                         "Caenorhabditis elegans",
                         "Schizosaccharomyces pombe 972h-"))

    expect_equal(length(orth$mouse), 23L)
    expect_equal(length(orth$human), 23L)
    expect_equal(length(orth$`baker's yeast`), 13L)
    expect_equal(length(orth$`Caenorhabditis elegans`), 18L)
    expect_equal(length(orth$`Schizosaccharomyces pombe 972h-`), 14L)

    expect_named(orth$mouse,
                 c("human: BCL6-HDAC4 complex", "human: BCL6-HDAC5 complex",
                   "mouse: BLOC-2 (biogenesis of lysosome-related organelles complex 2)",
                   "rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                   "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                   "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                   "S.pombe: nucleotide-excision repair factor 1 complex",
                   "S.pombe: nucleotide-excision repair factor 2 complex",
                   "S.pombe: nucleotide-excision repair factor 3 complex",
                   "human: myComplex",
                   "human: HuMAP2_00000_conf3", "human: HuMAP2_00001_conf4",
                   "human: HuMAP2_00002_conf5",
                   "C.elegans: CPX-365_GARP tethering complex",
                   "C.elegans: CPX-368_Polycomb Repressive Complex 2",
                   "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                   "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                   "S.pombe: CPX-546_DNA replication factor C complex",
                   "S.pombe: CPX-547_PCNA homotrimer",
                   "mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)",
                   "mouse: CPX-6_bZIP transcription factor complex, Atf4-Creb1 (+1 alt. ID)",
                   "mouse: CPX-7_bZIP transcription factor complex, Atf1-Atf4 (+1 alt. ID)"),
                 ignore.order = TRUE)
    expect_named(orth$human,
                 c("human: BCL6-HDAC4 complex", "human: BCL6-HDAC5 complex",
                   "mouse: BLOC-2 (biogenesis of lysosome-related organelles complex 2)",
                   "rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                   "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                   "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                   "S.pombe: nucleotide-excision repair factor 1 complex",
                   "S.pombe: nucleotide-excision repair factor 2 complex",
                   "S.pombe: nucleotide-excision repair factor 3 complex",
                   "human: myComplex",
                   "human: HuMAP2_00000_conf3", "human: HuMAP2_00001_conf4",
                   "human: HuMAP2_00002_conf5",
                   "C.elegans: CPX-365_GARP tethering complex",
                   "C.elegans: CPX-368_Polycomb Repressive Complex 2",
                   "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                   "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                   "S.pombe: CPX-546_DNA replication factor C complex",
                   "S.pombe: CPX-547_PCNA homotrimer",
                   "human: CPX-1_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)",
                   "human: CPX-8_bZIP transcription factor complex, ATF4-CREB1 (+1 alt. ID)",
                   "human: CPX-9_bZIP transcription factor complex, ATF1-ATF4 (+1 alt. ID)"),
                 ignore.order = TRUE)
    expect_named(orth$`baker's yeast`,
                 c("rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                   "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                   "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                   "S.pombe: nucleotide-excision repair factor 1 complex",
                   "S.pombe: nucleotide-excision repair factor 2 complex",
                   "S.pombe: nucleotide-excision repair factor 3 complex",
                   "human: HuMAP2_00001_conf4",
                   "C.elegans: CPX-365_GARP tethering complex",
                   "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                   "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                   "S.pombe: CPX-546_DNA replication factor C complex",
                   "S.pombe: CPX-547_PCNA homotrimer"),
                 ignore.order = TRUE)
    expect_named(orth$`Caenorhabditis elegans`,
                 c("rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                   "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                   "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                   "S.pombe: nucleotide-excision repair factor 1 complex",
                   "S.pombe: nucleotide-excision repair factor 2 complex",
                   "S.pombe: nucleotide-excision repair factor 3 complex",
                   "human: HuMAP2_00000_conf3", "human: HuMAP2_00001_conf4",
                   "human: HuMAP2_00002_conf5",
                   "C.elegans: CPX-365_GARP tethering complex",
                   "C.elegans: CPX-368_Polycomb Repressive Complex 2",
                   "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                   "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                   "S.pombe: CPX-546_DNA replication factor C complex",
                   "S.pombe: CPX-547_PCNA homotrimer",
                   "human: CPX-1_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)",
                   "human: CPX-8_bZIP transcription factor complex, ATF4-CREB1 (+3 alt. IDs)"),
                 ignore.order = TRUE)
    expect_named(orth$`Schizosaccharomyces pombe 972h-`,
                 c("rat: Bcl2l1-Dnm1l-Mff-Clta complex",
                   "S.cer: TRAPP complex", "S.cer: histone deacetylase complex",
                   "S.cer: Ada2p/Gcn5p/Ada3 transcription activator complex",
                   "S.pombe: nucleotide-excision repair factor 1 complex",
                   "S.pombe: nucleotide-excision repair factor 2 complex",
                   "S.pombe: nucleotide-excision repair factor 3 complex",
                   "human: HuMAP2_00001_conf4", "human: HuMAP2_00002_conf5",
                   "C.elegans: CPX-365_GARP tethering complex",
                   "C.elegans: CPX-372_Zyg-9/Tac-1 complex",
                   "S.pombe: CPX-540_Mitochondrial inner membrane pre-sequence translocase complex",
                   "S.pombe: CPX-546_DNA replication factor C complex",
                   "S.pombe: CPX-547_PCNA homotrimer"),
                 ignore.order = TRUE)

    expect_equal(orth$human$`human: BCL6-HDAC4 complex`, c("BCL6", "HDAC4"))
    expect_equal(orth$mouse$`human: BCL6-HDAC4 complex`, c("Bcl6", "Hdac4"))

    expect_equal(orth$human$`human: BCL6-HDAC5 complex`, c("BCL6", "HDAC5"))
    expect_equal(orth$mouse$`human: BCL6-HDAC5 complex`, c("Bcl6", "Hdac5"))

    expect_equal(orth$human$`mouse: BLOC-2 (biogenesis of lysosome-related organelles complex 2)`, c("HPS3", "HPS5", "HPS6"))
    expect_equal(orth$mouse$`mouse: BLOC-2 (biogenesis of lysosome-related organelles complex 2)`, c("Hps5", "Hps6", "Hps3"))

    expect_equal(orth$human$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`, c("BCL2L1", "CLTA", "DNM1L", "MFF"))
    expect_equal(orth$mouse$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`, c("Bcl2l1", "Clta", "Dnm1l", "Mff"))
    expect_equal(orth$`baker's yeast`$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`, c("DNM1"))
    expect_equal(orth$`Caenorhabditis elegans`$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`, c("clic-1", "drp-1", "mff-2"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`rat: Bcl2l1-Dnm1l-Mff-Clta complex`, c("clc1", "dnm1"))

    expect_equal(orth$human$`S.cer: TRAPP complex`, c("TRAPPC1", "TRAPPC2", "TRAPPC3", "TRAPPC4", "TRAPPC5", "TRAPPC6B"))
    expect_equal(orth$mouse$`S.cer: TRAPP complex`, c("Trappc1", "Trappc2", "Trappc3", "Trappc4", "Trappc5", "Trappc6b"))
    expect_equal(orth$`baker's yeast`$`S.cer: TRAPP complex`, c("BET3", "BET5", "GSG1", "KRE11", "TRS120", "TRS130", "TRS20", "TRS23", "TRS31", "TRS33"))
    expect_equal(orth$`Caenorhabditis elegans`$`S.cer: TRAPP complex`, c("trpp-1", "sedl-1", "trpp-3", "trpp-4", "trpp-5", "trpp-6"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`S.cer: TRAPP complex`, c("bet5", "trs20", "bet3", "trs23", "trs31", "trs33"))

    expect_equal(orth$human$`S.pombe: nucleotide-excision repair factor 3 complex`, c("ERCC2", "ERCC3", "GTF2H1", "GTF2H4", "MNAT1"))
    expect_equal(orth$mouse$`S.pombe: nucleotide-excision repair factor 3 complex`, c("Ercc2", "Ercc3", "Gtf2h1", "Gtf2h4", "Mnat1"))
    expect_equal(orth$`baker's yeast`$`S.pombe: nucleotide-excision repair factor 3 complex`, c("RAD3", "SSL2", "TFB1", "TFB2", "TFB3"))
    expect_equal(orth$`Caenorhabditis elegans`$`S.pombe: nucleotide-excision repair factor 3 complex`, c("xpd-1", "xpb-1", "gtf-2H1", "gtf-2H4", "mnat-1"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`S.pombe: nucleotide-excision repair factor 3 complex`, c("tfb1", "ptr8", "rad15", "tfb2", "pmh1"))

    expect_equal(orth$human$`human: HuMAP2_00001_conf4`, c("ANXA6", "SMAD1", "TBCB", "EHD1", "PAFAH1B2", "PAFAH1B3"))
    expect_equal(orth$mouse$`human: HuMAP2_00001_conf4`, c("Anxa6", "Ehd1", "Pafah1b2", "Pafah1b3", "Smad1", "Tbcb"))
    expect_equal(orth$`baker's yeast`$`human: HuMAP2_00001_conf4`, c("ALF1"))
    expect_equal(orth$`Caenorhabditis elegans`$`human: HuMAP2_00001_conf4`, c("rme-1", "sma-2", "tbcb-1"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`human: HuMAP2_00001_conf4`, c("alp11"))

    expect_equal(orth$human$`human: CPX-1_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)`, c("SMAD3", "SMAD4", "SMAD2"))
    expect_equal(orth$mouse$`mouse: CPX-10_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)`, c("Smad4", "Smad2", "Smad3"))
    expect_equal(orth$`Caenorhabditis elegans`$`human: CPX-1_SMAD2-SMAD3-SMAD4 complex (+1 alt. ID)`, c("sma-2", "sma-4"))

    expect_equal(orth$human$`C.elegans: CPX-365_GARP tethering complex`, c("VPS51", "VPS52", "VPS53", "VPS54"))
    expect_equal(orth$mouse$`C.elegans: CPX-365_GARP tethering complex`, c("Vps51", "Vps52", "Vps53", "Vps54"))
    expect_equal(orth$`Caenorhabditis elegans`$`C.elegans: CPX-365_GARP tethering complex`, c("vps-52", "vps-51", "vps-53", "vps-54"))
    expect_equal(orth$`baker's yeast`$`C.elegans: CPX-365_GARP tethering complex`, c("VPS52", "VPS53", "VPS54"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`C.elegans: CPX-365_GARP tethering complex`, c("vps51", "vps52", "vps53", "vps54"))

    expect_equal(orth$human$`S.pombe: CPX-547_PCNA homotrimer`, c("PCNA"))
    expect_equal(orth$mouse$`S.pombe: CPX-547_PCNA homotrimer`, c("Pcna"))
    expect_equal(orth$`Caenorhabditis elegans`$`S.pombe: CPX-547_PCNA homotrimer`, c("pcn-1"))
    expect_equal(orth$`baker's yeast`$`S.pombe: CPX-547_PCNA homotrimer`, c("POL30"))
    expect_equal(orth$`Schizosaccharomyces pombe 972h-`$`S.pombe: CPX-547_PCNA homotrimer`, c("pcn1"))

    ## Without PMIDs for custom complex
    ## Custom complex identical to an annotated one
    ## -------------------------------------------------------------------------
    custompath2 <- file.path(dbdir, "custom_complex_2.txt")
    custom <- data.frame(Complex.name = "myComplex",
                         Gene.names = "BCL6;HDAC5",
                         Organism = "Human",
                         Source = "Custom")
    write.table(custom, file = custompath2,
                col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
    dbs <- makeComplexDB(dbDir = dbdir, customComplexTxt = custompath2,
                         Cyc2008Db = cyc2008db, CorumDb = corumdb,
                         PombaseDb = pombasedb, HuMAP2Db = humap2db,
                         CPortal9606Db = cportal9606db,
                         CPortal10090Db = cportal10090db,
                         CPortal6239Db = cportal6239db,
                         CPortal284812Db = cportal284812db)
    compl <- readRDS(dbs$complPath)
    orth <- readRDS(dbs$orthPath)
    expect_s4_class(compl, "CharacterList")
    expect_equal(length(compl), 26L)
    expect_equal(compl$`human: myComplex`, c("BCL6", "HDAC5"))
    expect_type(orth, "list")
    expect_named(orth, c("mouse", "human", "baker's yeast",
                         "Caenorhabditis elegans",
                         "Schizosaccharomyces pombe 972h-"))
    expect_equal(length(orth$mouse), 22L)
    expect_equal(length(orth$human), 22L)
    expect_equal(length(orth$`baker's yeast`), 13L)
    expect_equal(length(orth$`Caenorhabditis elegans`), 18L)
    expect_equal(length(orth$`Schizosaccharomyces pombe 972h-`), 14L)
    expect_equal(mcols(orth$mouse)["human: myComplex (+1 alt. ID)",
                                   "All.names"],
                 "human: myComplex;human: BCL6-HDAC5 complex")
    expect_equal(mcols(orth$mouse)["human: myComplex (+1 alt. ID)",
                                   "PMID"],
                 "NA;11929873")

    ## Without custom complex
    ## -------------------------------------------------------------------------
    dbs <- makeComplexDB(dbDir = dbdir, customComplexTxt = NULL,
                         Cyc2008Db = cyc2008db, CorumDb = corumdb,
                         PombaseDb = pombasedb, HuMAP2Db = humap2db,
                         CPortal9606Db = cportal9606db,
                         CPortal10090Db = cportal10090db,
                         CPortal6239Db = cportal6239db,
                         CPortal284812Db = cportal284812db)
    compl <- readRDS(dbs$complPath)
    orth <- readRDS(dbs$orthPath)
    expect_s4_class(compl, "CharacterList")
    expect_equal(length(compl), 25L)
    expect_type(orth, "list")
    expect_named(orth, c("mouse", "human", "baker's yeast",
                         "Caenorhabditis elegans",
                         "Schizosaccharomyces pombe 972h-"))
    expect_equal(length(orth$mouse), 22L)
    expect_equal(length(orth$human), 22L)
    expect_equal(length(orth$`baker's yeast`), 13L)
    expect_equal(length(orth$`Caenorhabditis elegans`), 18L)
    expect_equal(length(orth$`Schizosaccharomyces pombe 972h-`), 14L)
})
