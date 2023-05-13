test_that("fixing feature IDs works", {
    mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                          package = "einprot")
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce

    ## -------------------------------------------------------------------------
    ## getFirstId
    ## -------------------------------------------------------------------------
    rd <- as.data.frame(SummarizedExperiment::rowData(sce))
    expect_error(getFirstId(df = rd, colName = 1, separator = ";"),
                 "'colName' must be of class 'character'")
    expect_error(getFirstId(df = rd, colName = c("Protein.IDs", "Protein.IDs"),
                            separator = ";"),
                 "'colName' must have length 1")
    expect_error(getFirstId(df = rd, colName = "Missing", separator = ";"),
                 "All values in 'colName' must be one of")
    expect_equal(getFirstId(df = rd, colName = "Protein.IDs", separator = ";")[1],
                 "A0A023T672")
    expect_equal(getFirstId(df = rd, colName = "Protein.IDs", separator = ",")[1],
                 "A0A023T672;Q9CWZ3-2;Q9CWZ3")
    expect_equal(getFirstId(df = rd, colName = "Oxidation.M.site.IDs",
                            separator = ";")[1], NA_character_)  ## empty value

    ## -------------------------------------------------------------------------
    ## getNthId
    ## -------------------------------------------------------------------------
    rd <- as.data.frame(SummarizedExperiment::rowData(sce))
    expect_error(getNthId(df = rd, colName = 1, N = 1, separator = ";"),
                 "'colName' must be of class 'character'")
    expect_error(getNthId(df = rd, colName = c("Protein.IDs", "Protein.IDs"),
                          N = 1, separator = ";"),
                 "'colName' must have length 1")
    expect_error(getNthId(df = rd, colName = "Missing", N = 1, separator = ";"),
                 "All values in 'colName' must be one of")
    expect_error(getNthId(df = rd, colName = "Protein.IDs", N = "1",
                          separator = ";"),
                 "'N' must be of class 'numeric'")
    expect_error(getNthId(df = rd, colName = "Protein.IDs", N = c(1, 2),
                          separator = ";"),
                 "'N' must have length 1")
    expect_error(getNthId(df = rd, colName = "Protein.IDs", N = 0,
                          separator = ";"),
                 "'N' must be within (0,Inf)", fixed = TRUE)
    expect_equal(getNthId(df = rd, colName = "Protein.IDs", N = 1,
                          separator = ";")[1],
                 "A0A023T672")
    expect_equal(getNthId(df = rd, colName = "Protein.IDs", N = 2,
                          separator = ";")[1],
                 "Q9CWZ3-2")
    expect_equal(getNthId(df = rd, colName = "Protein.IDs", N = 3,
                          separator = ";")[1],
                 "Q9CWZ3")
    expect_equal(getNthId(df = rd, colName = "Protein.IDs", N = 1,
                          separator = ",")[1],
                 "A0A023T672;Q9CWZ3-2;Q9CWZ3")
    expect_equal(getNthId(df = rd, colName = "Oxidation.M.site.IDs",
                          N = 1, separator = ";")[1], NA_character_)  ## empty value
    expect_equal(getNthId(df = rd, colName = "Oxidation.M.site.IDs",
                          N = 2, separator = ";")[1], NA_character_)  ## empty value

    ## -------------------------------------------------------------------------
    ## combineIds
    ## -------------------------------------------------------------------------
    rd <- as.data.frame(SummarizedExperiment::rowData(sce))
    expect_error(combineIds(df = rd, combineCols = 1, combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 "'combineCols' must be of class 'character'")
    expect_error(combineIds(df = rd, combineCols = "Missing",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 "All values in 'combineCols' must be one of")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = 1,
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 "'combineWhen' must be of class 'character'")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "Nonsense",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 "All values in 'combineWhen' must be one of")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = 1, joinSeparator = ".",
                            makeUnique = TRUE),
                 "'splitSeparator' must be of class 'character'")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = c(";", "."), joinSeparator = ".",
                            makeUnique = TRUE),
                 "length(splitSeparator) == length(combineCols) is not TRUE",
                 fixed = TRUE)
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = 1,
                            makeUnique = TRUE),
                 "'joinSeparator' must be of class 'character'")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = c(".", ";"),
                            makeUnique = TRUE),
                 "'joinSeparator' must have length 1")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = 1),
                 "'makeUnique' must be of class 'logical'")
    expect_error(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = c(TRUE, FALSE)),
                 "'makeUnique' must have length 1")

    ## Check that combinations work
    expect_equal(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE)[1],
                 "A0A023T672")
    expect_equal(combineIds(df = rd, combineCols = "Protein.IDs",
                            combineWhen = "nonunique",
                            splitSeparator = ",", joinSeparator = ".",
                            makeUnique = TRUE)[1],
                 "A0A023T672;Q9CWZ3-2;Q9CWZ3")
    expect_equal(combineIds(df = rd, combineCols = c("Protein.IDs", "Majority.protein.IDs"),
                            combineWhen = "always",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE)[1],
                 "A0A023T672.A0A023T672")
    expect_equal(combineIds(df = rd, combineCols = c("Protein.IDs", "Majority.protein.IDs"),
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE)[1],
                 "A0A023T672")

    ## Check combinations with some missing values
    sce1 <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                             nrows = 45)$sce[36:45, ]
    rd1 <- as.data.frame(SummarizedExperiment::rowData(sce1))
    expect_equal(
        combineIds(df = rd1, combineCols = c("Gene.names", "Majority.protein.IDs"),
                   combineWhen = "nonunique", splitSeparator = ";",
                   joinSeparator = ".", makeUnique = FALSE),
        combineIds(df = rd1, combineCols = c("Gene.names", "Majority.protein.IDs"),
                   combineWhen = "missing", splitSeparator = ";",
                   joinSeparator = ".", makeUnique = FALSE)
    )

    ## Check combinations with small synthetic data set
    rd2 <- data.frame(col1 = c("A;P", "B;L", "A", "C", "A;M"),
                      col2 = c("a,q", "b", "c,o", "d", "a,r"),
                      col3 = c("n1", "n2", "n3", "n4", "n5"),
                      col4 = c("", "A", NA, "B", "A"))
    expect_equal(combineIds(df = rd2, combineCols = c("col1", "col2"),
                            combineWhen = "always",
                            splitSeparator = c(";", ","),
                            joinSeparator = ".", makeUnique = FALSE),
                 c("A.a", "B.b", "A.c", "C.d", "A.a"))
    expect_equal(combineIds(df = rd2, combineCols = c("col1", "col2"),
                            combineWhen = "always",
                            splitSeparator = c(";", ","),
                            joinSeparator = ".", makeUnique = TRUE),
                 c("A.a", "B.b", "A.c", "C.d", "A.a.1"))
    expect_equal(combineIds(df = rd2, combineCols = c("col1", "col2"),
                            combineWhen = "nonunique",
                            splitSeparator = c(";", ","),
                            joinSeparator = ".", makeUnique = TRUE),
                 c("A.a", "B", "A.c", "C", "A.a.1"))
    expect_equal(combineIds(df = rd2, combineCols = c("col2", "col1"),
                            combineWhen = "nonunique",
                            splitSeparator = c(",", ";"),
                            joinSeparator = ".", makeUnique = TRUE),
                 c("a.A", "b", "c", "d", "a.A.1"))
    expect_equal(combineIds(df = rd2, combineCols = c("col1", "col2"),
                            combineWhen = "nonunique",
                            splitSeparator = c(",", ";"),
                            joinSeparator = ".", makeUnique = TRUE),
                 c("A;P", "B;L", "A", "C", "A;M"))
    expect_equal(combineIds(df = rd2, combineCols = c("col1", "col2", "col3"),
                            combineWhen = "nonunique",
                            splitSeparator = c(";", ",", ";"),
                            joinSeparator = ".", makeUnique = TRUE),
                 c("A.a.n1", "B", "A.c", "C", "A.a.n5"))
    expect_equal(combineIds(df = rd2, combineCols = c("col4", "col1"),
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = FALSE),
                 c("A", "A.B", "A", "B", "A.A"))
    expect_equal(combineIds(df = rd2, combineCols = c("col4", "col1"),
                            combineWhen = "nonunique",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 c("A", "A.B", "A.1", "B", "A.A"))
    expect_equal(combineIds(df = rd2, combineCols = c("col4", "col1"),
                            combineWhen = "missing",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = FALSE),
                 c("A", "A", "A", "B", "A"))
    expect_equal(combineIds(df = rd2, combineCols = c("col4", "col1"),
                            combineWhen = "missing",
                            splitSeparator = ";", joinSeparator = ".",
                            makeUnique = TRUE),
                 c("A", "A.1", "A.2", "B", "A.3"))

    ## -------------------------------------------------------------------------
    ## fixFeatureIds
    ## -------------------------------------------------------------------------
    ## Fail with wrong arguments
    expect_error(fixFeatureIds(sce = 1,
                               colDefs = list(einprotId = "Gene.names",
                                              einprotLabel = "Gene.names",
                                              einprotGene = "Gene.names",
                                              einprotProtein = "Majority.protein.IDs",
                                              IDsForSTRING = "Gene.names")),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(fixFeatureIds(sce = sce,
                               colDefs = list(einprotId = 1,
                                              einprotLabel = "Gene.names",
                                              einprotGene = "Gene.names",
                                              einprotProtein = "Majority.protein.IDs",
                                              IDsForSTRING = "Gene.names")),
                 "must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce,
                               colDefs = list(einprotId = "Gene.names",
                                              einprotLabel = 1,
                                              einprotGene = "Gene.names",
                                              einprotProtein = "Majority.protein.IDs",
                                              IDsForSTRING = "Gene.names")),
                 "must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce,
                               colDefs = list(einprotId = "Gene.names",
                                              einprotLabel = "Gene.names",
                                              einprotGene = 1,
                                              einprotProtein = "Majority.protein.IDs",
                                              IDsForSTRING = "Gene.names")),
                 "must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce,
                               colDefs = list(einprotId = "Gene.names",
                                              einprotLabel = "Gene.names",
                                              einprotGene = "Gene.names",
                                              einprotProtein = 1,
                                              IDsForSTRING = "Gene.names")),
                 "must be of class 'character'")
    expect_error(fixFeatureIds(sce = sce,
                               colDefs = list(einprotId = "Gene.names",
                                              einprotLabel = "Gene.names",
                                              einprotGene = "Gene.names",
                                              einprotProtein = "Majority.protein.IDs",
                                              IDsForSTRING = 1)),
                 "must be of class 'character'")

    ## Test that it does the right thing
    ## -------------------------------------------------------------------------
    ## All gene and protein names are unique
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 25)$sce
    gns0 <- SummarizedExperiment::rowData(sce)$Gene.names
    pns0 <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    gns <- vapply(strsplit(gns0, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns0, ";"), .subset, 1, FUN.VALUE = "")

    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                       einprotLabel = "Majority.protein.IDs",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                             combineWhen = "missing", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns)

    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = c("Gene.names", "Majority.protein.IDs"),
                       einprotLabel = "Majority.protein.IDs",
                       einprotGene = "Gene.names",
                       einprotProtein = "Majority.protein.IDs",
                       IDsForSTRING = "Gene.names"))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, paste(gns0, pns0, sep = "."))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns0)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns0)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, gns0)

    ## Missing gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce[36:45, ]
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs

    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    midx <- which(gns == "")  ## indices for missing gene IDs
    eidx <- which(gns != "")  ## indices for present gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[eidx], gns[eidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[eidx], gns[eidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## Duplicated gene names
    sce <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                            nrows = 45)$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.names
    pns <- SummarizedExperiment::rowData(sce)$Majority.protein.IDs
    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs")),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    midx <- which(gns == "")  ## indices for missing gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    didx <- setdiff(which(gns %in% gns[duplicated(gns)]), midx)  ## indices for duplicated gene IDs
    uidx <- setdiff(seq_along(gns), c(midx, didx))  ## indices for unique gene IDs
    expect_length(uidx, 38)
    expect_length(midx, 4)
    expect_length(didx, 3)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[didx], paste0(gns[didx], ".", pns[didx]))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names")),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names"),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## Change the separator
    scesep <- sce
    rowData(scesep)$Majority.protein.IDs <- gsub(";", ",",
                                                rowData(scesep)$Majority.protein.IDs)
    rowData(scesep)$Gene.names <- gsub(";", ",",
                                       rowData(scesep)$Gene.names)
    sce1 <- fixFeatureIds(
        sce = scesep,
        colDefs = list(einprotId = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names"),
                                                           splitSeparator = ","),
                       einprotLabel = "Gene.names",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.names",
                                                             separator = ","),
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs",
                                                                separator = ","),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.names", "Majority.protein.IDs"),
                                                              combineWhen = "missing",
                                                              splitSeparator = ",", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)

    ## Don't extract gene and STRING IDs
    sce1 <- fixFeatureIds(
        sce = scesep,
        colDefs = list(einprotId = function(df) combineIds(df, c("Majority.protein.IDs", "Gene.names"),
                                                           splitSeparator = ","),
                       einprotLabel = "Gene.names",
                       einprotGene = NULL,
                       einprotProtein = function(df) getFirstId(df, colName = "Majority.protein.IDs",
                                                                separator = ","),
                       IDsForSTRING = NULL))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, rep(NA_character_, nrow(sce1)))
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING, rep(NA_character_, nrow(sce1)))

    ## --------------------------------------------------------------------- ##
    ## PD data
    sce <- importExperiment(
        inFile = system.file("extdata", "pdtmt_example",
                             "Fig2_m23139_RTS_QC_varMods_Proteins.txt",
                             package = "einprot"),
        iColPattern = "^Abundance\\.F.+\\.Sample\\.")$sce
    gns <- SummarizedExperiment::rowData(sce)$Gene.Symbol
    pns <- SummarizedExperiment::rowData(sce)$Accession
    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Gene.Symbol", "Accession")),
                       einprotLabel = "Gene.Symbol",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.Symbol"),
                       einprotProtein = function(df) getFirstId(df, colName = "Accession"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.Symbol", "Accession"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    midx <- which(gns == "")  ## indices for missing gene IDs
    gns <- vapply(strsplit(gns, ";"), .subset, 1, FUN.VALUE = "")
    pns <- vapply(strsplit(pns, ";"), .subset, 1, FUN.VALUE = "")
    didx <- setdiff(which(gns %in% gns[duplicated(gns)]), midx)  ## indices for duplicated gene IDs
    uidx <- setdiff(seq_along(gns), c(midx, didx))  ## indices for unique gene IDs
    expect_length(uidx, 1620)
    expect_length(midx, 75)
    expect_length(didx, 28)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId[didx], paste0(gns[didx], ".", pns[didx]))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[midx], pns[midx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[uidx], gns[uidx])
    expect_equal(SummarizedExperiment::rowData(sce1)$IDsForSTRING[didx], gns[didx])

    ## Protein names are still there and unique
    sce1 <- fixFeatureIds(
        sce = sce,
        colDefs = list(einprotId = function(df) combineIds(df, c("Accession", "Gene.Symbol")),
                       einprotLabel = "Gene.Symbol",
                       einprotGene = function(df) getFirstId(df, colName = "Gene.Symbol"),
                       einprotProtein = function(df) getFirstId(df, colName = "Accession"),
                       IDsForSTRING = function(df) combineIds(df, c("Gene.Symbol", "Accession"),
                                                              combineWhen = "missing", makeUnique = FALSE)))
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotId, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotProtein, pns)
    expect_equal(SummarizedExperiment::rowData(sce1)$einprotGene, gns)
})
