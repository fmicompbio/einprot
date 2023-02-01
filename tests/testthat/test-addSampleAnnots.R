test_that("adding sample annotations works", {

    ## --------------------------------------------------------------------- ##
    ## Fail with wrong arguments
    expect_error(addSampleAnnots(1,
                                 sampleAnnot = mqSampleAnnot),
                 "'sce' must be of class 'SummarizedExperiment'")
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = 1),
                 "'sampleAnnot' must be of class 'data.frame'")
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = as.matrix(mqSampleAnnot)),
                 "'sampleAnnot' must be of class 'data.frame'")
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = S4Vectors::DataFrame(mqSampleAnnot)),
                 "'sampleAnnot' must be of class 'data.frame'")
    sampleAnnot1 <- mqSampleAnnot
    colnames(sampleAnnot1) <- c("sample", "group1")
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = sampleAnnot1),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    sampleAnnot1 <- mqSampleAnnot
    colnames(sampleAnnot1) <- c("sample1", "group")
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = sampleAnnot1),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = rbind(mqSampleAnnot, mqSampleAnnot)),
                 "all(!duplicated(sampleAnnot$sample)) is not TRUE", fixed = TRUE)

    ## --------------------------------------------------------------------- ##
    ## Add sample annotations
    sampleAnnot <- data.frame(sample = mqSamples,
                              group = gsub("_IP.*", "", mqSamples))
    sce1 <- addSampleAnnots(sce_mq_initial,
                            sampleAnnot = sampleAnnot)
    cdt1 <- SummarizedExperiment::colData(sce1)
    expect_equal(colnames(sce1), colnames(sce_mq_initial))
    expect_s4_class(sce1, "SummarizedExperiment")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group"))
    expect_equal(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$group, sampleAnnot$group)

    ## Try to add again - sample and group columns already exist
    expect_error(addSampleAnnots(sce1, sampleAnnot = sampleAnnot),
                 "'sce' already have column(s) named", fixed = TRUE)

    ## Provide sample annotations in different order
    sampleAnnot <- data.frame(sample = mqSamples,
                              group = gsub("_IP.*", "", mqSamples))
    set.seed(123)
    sampleAnnot <- sampleAnnot[sample(seq_len(nrow(sampleAnnot)),
                                      nrow(sampleAnnot)), ]
    sce1 <- addSampleAnnots(sce_mq_initial,
                            sampleAnnot = sampleAnnot)
    cdt1 <- SummarizedExperiment::colData(sce1)
    expect_equal(colnames(sce1), colnames(sce_mq_initial))
    expect_s4_class(sce1, "SummarizedExperiment")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group"))
    ordr <- match(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$sample, sampleAnnot$sample[ordr])
    expect_equal(cdt1$group, sampleAnnot$group[ordr])

    ## Extra samples in annotation
    sampleAnnot <- data.frame(sample = c(mqSamples, paste0(mqSamples, "_rep2")),
                              group = rep(gsub("_IP.*", "", mqSamples), 2))
    sce1 <- addSampleAnnots(sce_mq_initial,
                            sampleAnnot = sampleAnnot)
    cdt1 <- SummarizedExperiment::colData(sce1)
    expect_equal(colnames(sce1), colnames(sce_mq_initial))
    expect_s4_class(sce1, "SummarizedExperiment")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group"))
    expect_equal(cdt1$sample, mqSamples)
    expect_equal(cdt1$group, sampleAnnot$group[seq_len(nrow(cdt1))])

    ## Additional columns in sampleAnnot
    set.seed(123)
    sampleAnnot <- data.frame(sample = mqSamples,
                              group = gsub("_IP.*", "", mqSamples),
                              batch = rep(c("b1", "b2"), c(4, 5)),
                              age = 100 * runif(length(mqSamples)))
    sce1 <- addSampleAnnots(sce_mq_initial,
                            sampleAnnot = sampleAnnot)
    cdt1 <- SummarizedExperiment::colData(sce1)
    expect_equal(colnames(sce1), colnames(sce_mq_initial))
    expect_s4_class(sce1, "SummarizedExperiment")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group", "batch", "age"))
    expect_equal(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$group, sampleAnnot$group)
    expect_equal(cdt1$age, sampleAnnot$age)
    expect_equal(cdt1$batch, sampleAnnot$batch)

    ## Try to add existing column
    SummarizedExperiment::colData(sce1) <-
        SummarizedExperiment::colData(sce1)[, c("age"), drop = FALSE]
    expect_error(addSampleAnnots(sce1, sampleAnnot = sampleAnnot),
                 "Column already exists in SummarizedExperiment: age",
                 fixed = TRUE)

    ## Missing samples in annotation - should fail
    sampleAnnot <- data.frame(sample = mqSamples[1:7],
                              group = gsub("_IP.*", "", mqSamples[1:7]))
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = sampleAnnot),
                 "Some samples are missing")

    ## Include iColPattern in sample names - should fail
    sampleAnnot <- data.frame(sample = paste0("iBAQ.", mqSamples),
                              group = gsub("_IP.*", "", mqSamples))
    expect_error(addSampleAnnots(sce_mq_initial,
                                 sampleAnnot = sampleAnnot),
                 "Some samples are missing")
})
