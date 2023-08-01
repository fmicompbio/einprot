test_that("filtering works (MaxQuant)", {
    ## Fails with wrong argument specification
    expect_error(filterMaxQuant(sce = 1, minScore = 10, minPeptides = 2,
                                plotUpset = TRUE, exclFile = NULL),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = "10",
                                minPeptides = 2, plotUpset = TRUE,
                                exclFile = NULL),
                 "'minScore' must be of class 'numeric'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = c(1, 2),
                                minPeptides = 2, plotUpset = TRUE,
                                exclFile = NULL),
                 "'minScore' must have length 1")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = "2", plotUpset = TRUE,
                                exclFile = NULL),
                 "'minPeptides' must be of class 'numeric'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = c(1, 2), plotUpset = TRUE,
                                exclFile = NULL),
                 "'minPeptides' must have length 1")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = 1,
                                exclFile = NULL),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = c(TRUE, FALSE),
                                exclFile = NULL),
                 "'plotUpset' must have length 1")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = TRUE,
                                exclFile = 1),
                 "'exclFile' must be of class 'character'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = TRUE,
                                exclFile = c(tempfile(), tempfile())),
                 "'exclFile' must have length 1")

    ## Works with correct argument specification
    tfl <- tempfile(fileext = ".txt")
    out <- filterMaxQuant(sce_mq_final, minScore = 10, minPeptides = 3,
                          plotUpset = FALSE, exclFile = tfl)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Score >= 10 &
            rowData(sce_mq_final)$Peptides >= 3 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 53L)  ## same test as above, just with precomputed answer
    expect_true(file.exists(tfl))
    tmpin <- read.delim(tfl)
    expect_equal(nrow(tmpin), 97L)

    out <- filterMaxQuant(sce_mq_final, minScore = 2, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Score >= 2 &
            rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 97L)

    ## Don't filter on score
    out <- filterMaxQuant(sce_mq_final, minScore = NULL, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 97L)

    ## Don't filter on minPeptides
    out <- filterMaxQuant(sce_mq_final, minScore = 5, minPeptides = NULL,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Score >= 5 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 83L)

    ## Missing columns - Score
    tmp <- sce_mq_final
    rowData(tmp)$Score <- NULL
    out <- filterMaxQuant(tmp, minScore = 7, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 97L)

    ## Missing columns - Score and Only.identified.by.site
    tmp <- sce_mq_final
    rowData(tmp)$Score <- NULL
    rowData(tmp)$Only.identified.by.site <- NULL
    out <- filterMaxQuant(tmp, minScore = 7, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant))
    )))
    expect_equal(nrow(out), 112L)

    ## Only one column present
    tmp <- sce_mq_final
    rowData(tmp)$Score <- NULL
    rowData(tmp)$Only.identified.by.site <- NULL
    rowData(tmp)$Peptides <- NULL
    rowData(tmp)$Reverse <- NULL
    out <- filterMaxQuant(tmp, minScore = 7, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant))
    )))
    expect_equal(nrow(out), 112L)

    ## Missing columns - Potential.contaminant
    tmp <- sce_mq_final
    rowData(tmp)$Potential.contaminant <- NULL
    out <- filterMaxQuant(tmp, minScore = 2, minPeptides = 1,
                          plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Score >= 2 &
            rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Only.identified.by.site == "" |
                 is.na(rowData(sce_mq_final)$Only.identified.by.site))
    )))
    expect_equal(nrow(out), 130L)
})

test_that("filtering works (PD/TMT - proteins)", {
    ## Fails with wrong argument specification
    expect_error(filterPDTMT(sce = 1, inputLevel = "Proteins", minScore = 10,
                             minPeptides = 2, minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = 1,
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'inputLevel' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = c("Proteins", "PeptideGroups"),
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'inputLevel' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Nonsense",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "All values in 'inputLevel' must be one of")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = "10", minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minScore' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = c(1, 2), minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minScore' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = "2",
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minPeptides' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = c(1, 2),
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minPeptides' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = "TRUE", plotUpset = TRUE,
                             exclFile = NULL),
                 "'masterProteinsOnly' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = c(TRUE, FALSE), plotUpset = TRUE,
                             exclFile = NULL),
                 "'masterProteinsOnly' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = 1),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = c(TRUE, FALSE),
                             exclFile = NULL),
                 "'plotUpset' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = 1),
                 "'exclFile' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE, plotUpset = TRUE,
                             exclFile = c(tempfile(), tempfile())),
                 "'exclFile' must have length 1")

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 5, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 23L)  ## same test as above, just with precomputed answer

    ## Don't filter on minPeptides
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = 10,
                       minPeptides = NULL, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Don't filter on minScore
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = NULL,
                       minPeptides = 5, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 23L)  ## same test as above, just with precomputed answer

    ## Don't filter on minScore, filter on master proteins
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = NULL,
                       minPeptides = 5, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = TRUE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False" &
            rowData(sce_pd_final)$Master == "IsMasterProtein"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer

    ## Only master proteins
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = TRUE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Number.of.Peptides >= 3 &
            rowData(sce_pd_final)$Contaminant == "False" &
            rowData(sce_pd_final)$Master == "IsMasterProtein"
    )))
    expect_equal(nrow(out), 27L)  ## same test as above, just with precomputed answer

    ## Missing columns - Score
    tfl <- tempfile(fileext = ".txt")
    tmp <- sce_pd_final
    rowData(tmp)$Score.Sequest.HT.Sequest.HT <- NULL
    out <- filterPDTMT(tmp, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = FALSE,
                       exclFile = tfl)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Number.of.Peptides >= 3 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 35L)  ## same test as above, just with precomputed answer
    expect_true(file.exists(tfl))
    tmpin <- read.delim(tfl)
    expect_equal(nrow(tmpin), 35L)

    ## Missing columns - Number.of.Peptides
    tmp <- sce_pd_final
    rowData(tmp)$Number.of.Peptides <- NULL
    out <- filterPDTMT(tmp, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Only one column present
    tmp <- sce_pd_final
    rowData(tmp)$Number.of.Peptides <- NULL
    rowData(tmp)$Contaminant <- NULL
    out <- filterPDTMT(tmp, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Only one column present, but no features excluded (should not plot)
    tmp <- sce_pd_final
    rowData(tmp)$Number.of.Peptides <- NULL
    rowData(tmp)$Score.Sequest.HT.Sequest.HT <- NULL
    out <- filterPDTMT(tmp, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = FALSE, plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), nrow(tmp))
    expect_equal(nrow(out), 70L)  ## same test as above, just with precomputed answer

    ## Missing columns - Master
    tmp <- sce_pd_final
    rowData(tmp)$Master <- NULL
    out <- filterPDTMT(tmp, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       masterProteinsOnly = TRUE, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Number.of.Peptides >= 3 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer
})

test_that("filtering works (PD/TMT - peptidegroups)", {
    ## Fails with wrong argument specification
    expect_error(filterPDTMT(sce = 1, inputLevel = "PeptideGroups", minScore = 10,
                             minPeptides = 2, minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = 1,
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'inputLevel' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = c("Proteins", "PeptideGroups"),
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'inputLevel' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "Nonsense",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "All values in 'inputLevel' must be one of")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = "10", minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minDeltaScore' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = c(1, 2), minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minDeltaScore' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = "2",
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minPSMs' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = c(1, 2),
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'minPSMs' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = 1,
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'modificationsCol' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = c("Modifications", "Modifications"),
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'modificationsCol' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Nonsense",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "All values in 'modificationsCol' must be one of")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = 1,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'excludeUnmodifiedPeptides' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = c(TRUE, FALSE),
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "'excludeUnmodifiedPeptides' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = 1, plotUpset = TRUE,
                             exclFile = NULL),
                 "'keepModifications' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = 1,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = c("mod1", "mod2"), plotUpset = TRUE,
                             exclFile = NULL),
                 "'keepModifications' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = 1,
                             exclFile = NULL),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = c(TRUE, FALSE),
                             exclFile = NULL),
                 "'plotUpset' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = 1),
                 "'exclFile' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2,
                             masterProteinsOnly = FALSE,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = c(tempfile(), tempfile())),
                 "'exclFile' must have length 1")

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = 2,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 4L)  ## same test as above, just with precomputed answer

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.5, minPSMs = 1,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.5 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer

    ## Works with correct argument specification - with modificationsCol = NULL
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.5, minPSMs = 1,
                       masterProteinsOnly = FALSE, modificationsCol = NULL,
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.5 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer

    ## Don't filter on Number.of.PSMs
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = NULL,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 62L)  ## same test as above, just with precomputed answer

    ## Don't filter on delta score
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = NULL, minPSMs = 2,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 4L)  ## same test as above, just with precomputed answer

    ## Missing column - Number.of.PSMs
    tfl <- tempfile(fileext = ".txt")
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Number.of.PSMs <- NULL
    out <- filterPDTMT(tmp, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = 2,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = FALSE,
                       exclFile = tfl)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 62L)  ## same test as above, just with precomputed answer
    expect_true(file.exists(tfl))
    tmpin <- read.delim(tfl)
    expect_equal(nrow(tmpin), 8L)

    ## Missing column - Contaminant
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Contaminant <- NULL
    out <- filterPDTMT(tmp, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = 2,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = FALSE,
                       keepModifications = NULL, plotUpset = FALSE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2
    )))
    expect_equal(nrow(out), 5L)  ## same test as above, just with precomputed answer

    ## Exclude unmodified peptides
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Modifications[c(1, 5, 10, 20)] <- ""
    out <- filterPDTMT(tmp, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = NULL,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = TRUE,
                       keepModifications = NULL, plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(tmp)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(tmp)$Contaminant == "False" &
            rowData(tmp)$Modifications != ""
    )))
    expect_equal(nrow(out), 59L)  ## same test as above, just with precomputed answer

    ## Exclude unmodified peptides, keep only Carbamidomethyl modifications
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Modifications[c(1, 5, 10, 20)] <- ""
    out <- filterPDTMT(tmp, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = NULL,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = TRUE,
                       keepModifications = "Carbamidomethyl", plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(tmp)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(tmp)$Contaminant == "False" &
            rowData(tmp)$Modifications != "" &
            grepl("Carbamidomethyl", rowData(tmp)$Modifications)
    )))
    expect_equal(nrow(out), 9L)  ## same test as above, just with precomputed answer

    ## Exclude unmodified peptides, keep only Carbamidomethyl or 1xTMTpro [K12] modifications
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Modifications[c(1, 5, 10, 20)] <- ""
    out <- filterPDTMT(tmp, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = NULL,
                       masterProteinsOnly = FALSE, modificationsCol = "Modifications",
                       excludeUnmodifiedPeptides = TRUE,
                       keepModifications = "Carbamidomethyl|1xTMTpro \\[K12\\]", plotUpset = TRUE,
                       exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(tmp)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(tmp)$Contaminant == "False" &
            rowData(tmp)$Modifications != "" &
            (grepl("Carbamidomethyl", rowData(tmp)$Modifications) |
                 grepl("1xTMTpro \\[K12\\]", rowData(tmp)$Modifications))
    )))
    expect_equal(nrow(out), 13L)  ## same test as above, just with precomputed answer

    ## Fails if the Contaminant column has the wrong type of values (not True/False)
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Contaminant[rowData(tmp)$Contaminant == "True"] <- "Trrue"
    expect_error(filterPDTMT(tmp, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 0, minDeltaScore = 0.5, minPSMs = 1,
                             modificationsCol = "Modifications",
                             excludeUnmodifiedPeptides = FALSE,
                             keepModifications = NULL, plotUpset = TRUE,
                             exclFile = NULL),
                 "Something went wrong in the filtering")
})

test_that("filtering works (FragPipe)", {
    ## Fails with wrong argument specification
    expect_error(filterFragPipe(sce = 1, minPeptides = 2,
                                plotUpset = TRUE,
                                revPattern = "^rev_", exclFile = NULL),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = "2", plotUpset = TRUE,
                                revPattern = "^rev_", exclFile = NULL),
                 "'minPeptides' must be of class 'numeric'")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = c(1, 2), plotUpset = TRUE,
                                revPattern = "^rev_", exclFile = NULL),
                 "'minPeptides' must have length 1")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = 1,
                                revPattern = "^rev_", exclFile = NULL),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = c(TRUE, FALSE),
                                revPattern = "^rev_", exclFile = NULL),
                 "'plotUpset' must have length 1")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = TRUE,
                                revPattern = 1, exclFile = NULL),
                 "'revPattern' must be of class 'character'")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = TRUE,
                                revPattern = c("p1", "p2"), exclFile = NULL),
                 "'revPattern' must have length 1")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = TRUE,
                                revPattern = "^rev_", exclFile = 1),
                 "'exclFile' must be of class 'character'")
    expect_error(filterFragPipe(sce = sce_fp_final,
                                minPeptides = 2, plotUpset = TRUE,
                                revPattern = "^rev_",
                                exclFile = c(tempfile(), tempfile())),
                 "'exclFile' must have length 1")

    ## Works with correct argument specification
    tfl <- tempfile(fileext = ".txt")
    out <- filterFragPipe(sce_fp_final, minPeptides = 3,
                          plotUpset = FALSE, revPattern = "^rev_",
                          exclFile = tfl)
    expect_equal(nrow(out), length(which(
        rowData(sce_fp_final)$Combined.Total.Peptides >= 3 &
            !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 87L)  ## same test as above, just with precomputed answer
    expect_true(file.exists(tfl))
    tmpin <- read.delim(tfl)
    expect_equal(nrow(tmpin), 150L - 87L)

    out <- filterFragPipe(sce_fp_final, minPeptides = 1,
                          plotUpset = TRUE, revPattern = "^rev_",
                          exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_fp_final)$Combined.Total.Peptides >= 1 &
            !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)

    ## Don't filter on minPeptides
    out <- filterFragPipe(sce_fp_final, minPeptides = NULL,
                          plotUpset = TRUE, revPattern = "^rev_",
                          exclFile = NULL)
    expect_equal(nrow(out), length(which(
        !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)

    ## Missing columns - Combined.Total.Peptides
    tmp <- sce_fp_final
    rowData(tmp)$Combined.Total.Peptides <- NULL
    out <- filterFragPipe(tmp, minPeptides = 3,
                          plotUpset = TRUE, revPattern = "^rev_",
                          exclFile = NULL)
    expect_equal(nrow(out), length(which(
        !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)
})
