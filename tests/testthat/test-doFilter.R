test_that("filtering works (MaxQuant)", {
    ## Fails with wrong argument specification
    expect_error(filterFeaturesSE(
        sce = 1,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = TRUE, exclFile = NULL),
        "'sce' must be of class 'SummarizedExperiment'")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = 1,
        plotUpset = TRUE, exclFile = NULL),
        "'filtersSE' must be of class 'list'")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(function(se) rowData(se)$Score >= 10,
                         function(se) rowData(se)$Peptides >= 2),
        plotUpset = TRUE, exclFile = NULL),
        "'namesfiltersSE' must not be NULL")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(Score = function(se, minScore) rowData(se)$Score >= minScore,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = TRUE, exclFile = NULL),
        "is not TRUE")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = 1, exclFile = NULL),
        "'plotUpset' must be of class 'logical'")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = c(TRUE, FALSE), exclFile = NULL),
        "'plotUpset' must have length 1")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = TRUE, exclFile = 1),
        "'exclFile' must be of class 'character'")
    expect_error(filterFeaturesSE(
        sce = sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 2),
        plotUpset = TRUE, exclFile = c(tempfile(), tempfile())),
        "'exclFile' must have length 1")

    ## Works with correct argument specification
    tfl <- tempfile(fileext = ".txt")
    out <- filterFeaturesSE(
        sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 10,
                         Peptides = function(se) rowData(se)$Peptides >= 3,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)),
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

    out <- filterFeaturesSE(
        sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 2,
                         Peptides = function(se) rowData(se)$Peptides >= 1,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)),
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

    ## A single filter, with plot
    out <- filterFeaturesSE(
        sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 2),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Score >= 2
    )))
    expect_equal(nrow(out), 130L)

    ## Don't filter on score
    out <- filterFeaturesSE(
        sce_mq_final,
        filtersSE = list(Peptides = function(se) rowData(se)$Peptides >= 1,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)),
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
    out <- filterFeaturesSE(
        sce_mq_final,
        filtersSE = list(Score = function(se) rowData(se)$Score >= 5,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)),
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(Score = function(se) if ("Score" %in% colnames(rowData(se))) rowData(se)$Score >= 2 else rep(TRUE, nrow(se)),
                         Peptides = function(se) rowData(se)$Peptides >= 1,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)),
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(Score = function(se) if ("Score" %in% colnames(rowData(se))) rowData(se)$Score >= 7 else rep(TRUE, nrow(se)),
                         Peptides = function(se) rowData(se)$Peptides >= 1,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant),
                         OnlySite = function(se) if ("Only.identified.by.site" %in% colnames(rowData(se))) (as.character(rowData(se)$Only.identified.by.site) == "" |
                             is.na(rowData(se)$Only.identified.by.site)) else rep(TRUE, nrow(se))),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_mq_final)$Peptides >= 1 &
            (rowData(sce_mq_final)$Reverse == "" |
                 is.na(rowData(sce_mq_final)$Reverse)) &
            (rowData(sce_mq_final)$Potential.contaminant == "" |
                 is.na(rowData(sce_mq_final)$Potential.contaminant))
    )))
    expect_equal(nrow(out), 112L)

    ## Missing columns - Potential.contaminant
    tmp <- sce_mq_final
    rowData(tmp)$Potential.contaminant <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(Score = function(se) if ("Score" %in% colnames(rowData(se))) rowData(se)$Score >= 2 else rep(TRUE, nrow(se)),
                         Peptides = function(se) rowData(se)$Peptides >= 1,
                         Reverse = function(se) as.character(rowData(se)$Reverse) == "" |
                             is.na(rowData(se)$Reverse),
                         Contaminant = function(se) if ("Potential.contaminant" %in% colnames(rowData(se))) (as.character(rowData(se)$Potential.contaminant) == "" |
                             is.na(rowData(se)$Potential.contaminant)) else rep(TRUE, nrow(se)),
                         OnlySite = function(se) if ("Only.identified.by.site" %in% colnames(rowData(se))) (as.character(rowData(se)$Only.identified.by.site) == "" |
                                                                                                                is.na(rowData(se)$Only.identified.by.site)) else rep(TRUE, nrow(se))),
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
    out <- filterFeaturesSE(
        sce_pd_final,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 5 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 23L)  ## same test as above, just with precomputed answer

    ## Don't filter on minPeptides
    out <- filterFeaturesSE(
        sce_pd_final,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Don't filter on minScore
    out <- filterFeaturesSE(
        sce_pd_final,
        filtersSE = list(
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 5 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 23L)  ## same test as above, just with precomputed answer

    ## Don't filter on minScore, filter on master proteins
    out <- filterFeaturesSE(
        sce_pd_final,
        filtersSE = list(
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 5 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            Master = function(sce)  if ("Master" %in% colnames(rowData(sce))) rowData(sce)$Master == "IsMasterProtein" & !is.na(rowData(sce)$Master) else rep(TRUE, nrow(sce))
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Number.of.Peptides >= 5 &
            rowData(sce_pd_final)$Contaminant == "False" &
            rowData(sce_pd_final)$Master == "IsMasterProtein"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer

    ## Only master proteins
    out <- filterFeaturesSE(
        sce_pd_final,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 3 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            Master = function(sce)  if ("Master" %in% colnames(rowData(sce))) rowData(sce)$Master == "IsMasterProtein" & !is.na(rowData(sce)$Master) else rep(TRUE, nrow(sce))
        ),
        plotUpset = FALSE, exclFile = NULL)
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 3 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = tfl)
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 3 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Only one column present
    tmp <- sce_pd_final
    rowData(tmp)$Number.of.Peptides <- NULL
    rowData(tmp)$Contaminant <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 5 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer

    ## Only one column present, but no features excluded (should not plot)
    tmp <- sce_pd_final
    rowData(tmp)$Number.of.Peptides <- NULL
    rowData(tmp)$Score.Sequest.HT.Sequest.HT <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 5 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), nrow(tmp))
    expect_equal(nrow(out), 70L)  ## same test as above, just with precomputed answer

    ## Missing columns - Master
    tmp <- sce_pd_final
    rowData(tmp)$Master <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Score = function(sce) if ("Score.Sequest.HT.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Score.Sequest.HT.Sequest.HT >= 10 & !is.na(rowData(sce)$Score.Sequest.HT.Sequest.HT) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Number.of.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Number.of.Peptides >= 3 & !is.na(rowData(sce)$Number.of.Peptides) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            Master = function(sce)  if ("Master" %in% colnames(rowData(sce))) rowData(sce)$Master == "IsMasterProtein" & !is.na(rowData(sce)$Master) else rep(TRUE, nrow(sce))
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_final)$Score.Sequest.HT.Sequest.HT >= 10 &
            rowData(sce_pd_final)$Number.of.Peptides >= 3 &
            rowData(sce_pd_final)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 30L)  ## same test as above, just with precomputed answer
})

test_that("filtering works (PD/TMT - peptidegroups)", {
    out <- filterFeaturesSE(
        sce_pd_peptide_initial,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            PSMs = function(sce) if ("Number.of.PSMs" %in% colnames(rowData(sce))) rowData(sce)$Number.of.PSMs >= 2 & !is.na(rowData(sce)$Number.of.PSMs) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 4L)  ## same test as above, just with precomputed answer

    ## Works with correct argument specification
    out <- filterFeaturesSE(
        sce_pd_peptide_initial,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.5 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            PSMs = function(sce) if ("Number.of.PSMs" %in% colnames(rowData(sce))) rowData(sce)$Number.of.PSMs >= 1 & !is.na(rowData(sce)$Number.of.PSMs) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.5 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 20L)  ## same test as above, just with precomputed answer

    ## Don't filter on Number.of.PSMs
    out <- filterFeaturesSE(
        sce_pd_peptide_initial,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 62L)  ## same test as above, just with precomputed answer

    ## Don't filter on delta score
    out <- filterFeaturesSE(
        sce_pd_peptide_initial,
        filtersSE = list(
            PSMs = function(sce) if ("Number.of.PSMs" %in% colnames(rowData(sce))) rowData(sce)$Number.of.PSMs >= 2 & !is.na(rowData(sce)$Number.of.PSMs) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2 &
            rowData(sce_pd_peptide_initial)$Contaminant == "False"
    )))
    expect_equal(nrow(out), 4L)  ## same test as above, just with precomputed answer

    ## Missing column - Number.of.PSMs
    tfl <- tempfile(fileext = ".txt")
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Number.of.PSMs <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            PSMs = function(sce) if ("Number.of.PSMs" %in% colnames(rowData(sce))) rowData(sce)$Number.of.PSMs >= 2 & !is.na(rowData(sce)$Number.of.PSMs) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = tfl)
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            PSMs = function(sce) if ("Number.of.PSMs" %in% colnames(rowData(sce))) rowData(sce)$Number.of.PSMs >= 2 & !is.na(rowData(sce)$Number.of.PSMs) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))}
        ),
        plotUpset = FALSE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_pd_peptide_initial)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(sce_pd_peptide_initial)$Number.of.PSMs >= 2
    )))
    expect_equal(nrow(out), 5L)  ## same test as above, just with precomputed answer

    ## Exclude unmodified peptides
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Modifications[c(1, 5, 10, 20)] <- ""
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            ExclUnmodified = function(sce) if ("Modifications" %in% colnames(rowData(sce))) rowData(sce)$Modifications != "" & !is.na(rowData(sce)$Modifications) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(tmp)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(tmp)$Contaminant == "False" &
            rowData(tmp)$Modifications != ""
    )))
    expect_equal(nrow(out), 59L)  ## same test as above, just with precomputed answer

    ## Exclude unmodified peptides, keep only Carbamidomethyl modifications
    tmp <- sce_pd_peptide_initial
    rowData(tmp)$Modifications[c(1, 5, 10, 20)] <- ""
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            ExclUnmodified = function(sce) if ("Modifications" %in% colnames(rowData(sce))) rowData(sce)$Modifications != "" & !is.na(rowData(sce)$Modifications) else rep(TRUE, nrow(sce)),
            KeepMods = function(sce) if ("Modifications" %in% colnames(rowData(sce))) grepl("Carbamidomethyl", rowData(sce)$Modifications) & !is.na(rowData(sce)$Modifications) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
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
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            DeltaScore = function(sce) if ("Delta.Score.by.Search.Engine.Sequest.HT" %in% colnames(rowData(sce))) rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 & !is.na(rowData(sce)$Delta.Score.by.Search.Engine.Sequest.HT) else rep(TRUE, nrow(sce)),
            Contaminant = function(sce) { if ("Contaminant" %in% colnames(rowData(sce))) { rowData(sce)$Contaminant <- as.character(rowData(sce)$Contaminant);rowData(sce)$Contaminant == "False" & !is.na(rowData(sce)$Contaminant) } else rep(TRUE, nrow(sce))},
            ExclUnmodified = function(sce) if ("Modifications" %in% colnames(rowData(sce))) rowData(sce)$Modifications != "" & !is.na(rowData(sce)$Modifications) else rep(TRUE, nrow(sce)),
            KeepMods = function(sce) if ("Modifications" %in% colnames(rowData(sce))) grepl("Carbamidomethyl|1xTMTpro \\[K12\\]", rowData(sce)$Modifications) & !is.na(rowData(sce)$Modifications) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(tmp)$Delta.Score.by.Search.Engine.Sequest.HT >= 0.1 &
            rowData(tmp)$Contaminant == "False" &
            rowData(tmp)$Modifications != "" &
            (grepl("Carbamidomethyl", rowData(tmp)$Modifications) |
                 grepl("1xTMTpro \\[K12\\]", rowData(tmp)$Modifications))
    )))
    expect_equal(nrow(out), 13L)  ## same test as above, just with precomputed answer
})

test_that("filtering works (FragPipe)", {
    ## Works with correct argument specification
    tfl <- tempfile(fileext = ".txt")
    out <- filterFeaturesSE(
        sce_fp_final,
        filtersSE = list(
            Contaminant = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^contam_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Reverse = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^rev_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Combined.Total.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Combined.Total.Peptides >= 3 & !is.na(rowData(sce)$Combined.Total.Peptides) else rep(TRUE, nrow(sce))
        ),
        plotUpset = FALSE, exclFile = tfl)
    expect_equal(nrow(out), length(which(
        rowData(sce_fp_final)$Combined.Total.Peptides >= 3 &
            !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 87L)  ## same test as above, just with precomputed answer
    expect_true(file.exists(tfl))
    tmpin <- read.delim(tfl)
    expect_equal(nrow(tmpin), 150L - 87L)

    out <- filterFeaturesSE(
        sce_fp_final,
        filtersSE = list(
            Contaminant = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^contam_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Reverse = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^rev_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Combined.Total.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Combined.Total.Peptides >= 1 & !is.na(rowData(sce)$Combined.Total.Peptides) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        rowData(sce_fp_final)$Combined.Total.Peptides >= 1 &
            !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)

    ## Don't filter on minPeptides
    out <- filterFeaturesSE(
        sce_fp_final,
        filtersSE = list(
            Contaminant = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^contam_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Reverse = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^rev_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)

    ## Missing columns - Combined.Total.Peptides
    tmp <- sce_fp_final
    rowData(tmp)$Combined.Total.Peptides <- NULL
    out <- filterFeaturesSE(
        tmp,
        filtersSE = list(
            Contaminant = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^contam_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Reverse = function(sce) if ("Protein" %in% colnames(rowData(sce))) !grepl("^rev_", rowData(sce)$Protein) & !is.na(rowData(sce)$Protein) else rep(TRUE, nrow(sce)),
            Peptides = function(sce) if ("Combined.Total.Peptides" %in% colnames(rowData(sce))) rowData(sce)$Combined.Total.Peptides >= 3 & !is.na(rowData(sce)$Combined.Total.Peptides) else rep(TRUE, nrow(sce))
        ),
        plotUpset = TRUE, exclFile = NULL)
    expect_equal(nrow(out), length(which(
        !grepl("^rev_", rownames(sce_fp_final)) &
            !grepl("^contam_", rownames(sce_fp_final))
    )))
    expect_equal(nrow(out), 113L)

    ## Empty filter list -> no filtering
    out <- filterFeaturesSE(sce_fp_final, filtersSE = list(), plotUpset = TRUE,
                            exclFile = NULL)
    expect_identical(sce_fp_final, out)
})
