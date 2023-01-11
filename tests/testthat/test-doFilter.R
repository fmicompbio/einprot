test_that("filtering works (MaxQuant)", {
    ## Fails with wrong argument specification
    expect_error(filterMaxQuant(sce = 1, minScore = 10, minPeptides = 2,
                                plotUpset = TRUE),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = "10",
                                minPeptides = 2, plotUpset = TRUE),
                 "'minScore' must be of class 'numeric'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = c(1, 2),
                                minPeptides = 2, plotUpset = TRUE),
                 "'minScore' must have length 1")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = "2", plotUpset = TRUE),
                 "'minPeptides' must be of class 'numeric'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = c(1, 2), plotUpset = TRUE),
                 "'minPeptides' must have length 1")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = 1),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterMaxQuant(sce = sce_mq_final, minScore = 10,
                                minPeptides = 2, plotUpset = c(TRUE, FALSE)),
                 "'plotUpset' must have length 1")

    ## Works with correct argument specification
    out <- filterMaxQuant(sce_mq_final, minScore = 10, minPeptides = 3,
                          plotUpset = FALSE)
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

    out <- filterMaxQuant(sce_mq_final, minScore = 2, minPeptides = 1,
                          plotUpset = TRUE)
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
})

test_that("filtering works (PD/TMT - proteins)", {
    ## Fails with wrong argument specification
    expect_error(filterPDTMT(sce = 1, inputLevel = "Proteins", minScore = 10,
                             minPeptides = 2, minDeltaScore = 0, minPSMs = 1,
                             plotUpset = TRUE),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = 1,
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'inputLevel' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = c("Proteins", "PeptideGroups"),
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'inputLevel' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Nonsense",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "All values in 'inputLevel' must be one of")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = "10", minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'minScore' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = c(1, 2), minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'minScore' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = "2",
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'minPeptides' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = c(1, 2),
                             minDeltaScore = 0, minPSMs = 1, plotUpset = TRUE),
                 "'minPeptides' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = 1),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_final, inputLevel = "Proteins",
                             minScore = 10, minPeptides = 2,
                             minDeltaScore = 0, minPSMs = 1, plotUpset = c(TRUE, FALSE)),
                 "'plotUpset' must have length 1")

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_final, inputLevel = "Proteins", minScore = 10,
                       minPeptides = 3, minDeltaScore = 0, minPSMs = 1,
                       plotUpset = FALSE)
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
                             plotUpset = TRUE),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = 1,
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2, plotUpset = TRUE),
                 "'inputLevel' must be of class 'character'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = c("Proteins", "PeptideGroups"),
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2, plotUpset = TRUE),
                 "'inputLevel' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "Nonsense",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2, plotUpset = TRUE),
                 "All values in 'inputLevel' must be one of")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = "10", minPSMs = 2, plotUpset = TRUE),
                 "'minDeltaScore' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = c(1, 2), minPSMs = 1, plotUpset = TRUE),
                 "'minDeltaScore' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = "2", plotUpset = TRUE),
                 "'minPSMs' must be of class 'numeric'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 0, minPSMs = c(1, 2), plotUpset = TRUE),
                 "'minPSMs' must have length 1")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2, plotUpset = 1),
                 "'plotUpset' must be of class 'logical'")
    expect_error(filterPDTMT(sce = sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                             minScore = 0, minPeptides = 1,
                             minDeltaScore = 10, minPSMs = 2, plotUpset = c(TRUE, FALSE)),
                 "'plotUpset' must have length 1")

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.1, minPSMs = 2,
                       plotUpset = FALSE)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 4L)  ## same test as above, just with precomputed answer

    ## Works with correct argument specification
    out <- filterPDTMT(sce_pd_peptide_initial, inputLevel = "PeptideGroups",
                       minScore = 0, minPeptides = 0, minDeltaScore = 0.5, minPSMs = 1,
                       plotUpset = TRUE)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.5 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer
})

