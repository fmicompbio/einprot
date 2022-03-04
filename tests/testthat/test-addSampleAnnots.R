test_that("adding sample annotations works", {

    ## --------------------------------------------------------------------- ##
    ## Fail with wrong arguments
    expect_error(addSampleAnnots(1, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list()),
                 "'qft' must be of class 'QFeatures'")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = 1,
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list()),
                 "'iColPattern' must be of class 'character'")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = c("a", "b"),
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list()),
                 "'iColPattern' must have length 1")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = 1,
                                 mergeGroups = list()),
                 "'sampleAnnot' must be of class 'data.frame'")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = as.matrix(sampleAnnot),
                                 mergeGroups = list()),
                 "'sampleAnnot' must be of class 'data.frame'")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = S4Vectors::DataFrame(sampleAnnot),
                                 mergeGroups = list()),
                 "'sampleAnnot' must be of class 'data.frame'")
    sampleAnnot1 <- sampleAnnot
    colnames(sampleAnnot1) <- c("sample", "group1")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot1,
                                 mergeGroups = list()),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    sampleAnnot1 <- sampleAnnot
    colnames(sampleAnnot1) <- c("sample1", "group")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot1,
                                 mergeGroups = list()),
                 "colnames(sampleAnnot)) is not TRUE", fixed = TRUE)
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = rbind(sampleAnnot, sampleAnnot),
                                 mergeGroups = list()),
                 "all(!duplicated(sampleAnnot$sample)) is not TRUE", fixed = TRUE)
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = c("Adnp", "RBC_ctrl")),
                 "'mergeGroups' must be of class 'list'")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list(c("Adnp", "RBC_ctrl"))),
                 "'namesmergeGroups' must not be NULL")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list(g1 = "Adnp",
                                                    g1 = "RBC_ctrl")),
                 "'mergeGroups' must be a named list, without duplicated")
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot,
                                 mergeGroups = list(g1 = "Adnp",
                                                    g12 = c("Adnp", "RBC_ctrl"))),
                 "A given name can just be part of one merged group")

    ## --------------------------------------------------------------------- ##
    ## Add sample annotations
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    qft1 <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                            sampleAnnot = sampleAnnot, mergeGroups = list())
    cdt1 <- SummarizedExperiment::colData(qft1)
    expect_s4_class(qft1, "QFeatures")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group_orig", "group"))
    expect_equal(cdt1$group_orig, cdt1$group)
    expect_equal(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$group, sampleAnnot$group)

    ## Provide sample annotations in different order
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    set.seed(123)
    sampleAnnot <- sampleAnnot[sample(seq_len(nrow(sampleAnnot)),
                                      nrow(sampleAnnot)), ]
    qft1 <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                            sampleAnnot = sampleAnnot, mergeGroups = list())
    cdt1 <- SummarizedExperiment::colData(qft1)
    expect_s4_class(qft1, "QFeatures")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group_orig", "group"))
    expect_equal(cdt1$group_orig, cdt1$group)
    ordr <- match(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$sample, sampleAnnot$sample[ordr])
    expect_equal(cdt1$group, sampleAnnot$group[ordr])

    ## Extra samples in annotation
    sampleAnnot <- data.frame(sample = c(samples, paste0(samples, "_rep2")),
                              group = rep(gsub("_IP.*", "", samples), 2))
    qft1 <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                            sampleAnnot = sampleAnnot, mergeGroups = list())
    cdt1 <- SummarizedExperiment::colData(qft1)
    expect_s4_class(qft1, "QFeatures")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group_orig", "group"))
    expect_equal(cdt1$group_orig, cdt1$group)
    expect_equal(cdt1$sample, samples)
    expect_equal(cdt1$group, sampleAnnot$group[seq_len(nrow(cdt1))])

    ## Additional columns in sampleAnnot - only batch will be used
    set.seed(123)
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples),
                              batch = rep(c("b1", "b2"), c(4, 5)),
                              age = 100 * runif(length(samples)))
    qft1 <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                            sampleAnnot = sampleAnnot, mergeGroups = list())
    cdt1 <- SummarizedExperiment::colData(qft1)
    expect_s4_class(qft1, "QFeatures")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group_orig", "batch", "group"))
    expect_equal(cdt1$group_orig, cdt1$group)
    expect_equal(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$group, sampleAnnot$group)

    ## Missing samples in annotation - should fail
    sampleAnnot <- data.frame(sample = samples[1:7],
                              group = gsub("_IP.*", "", samples[1:7]))
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot, mergeGroups = list()),
                 "Some samples are missing")

    ## Include iColPattern in sample names - should fail
    sampleAnnot <- data.frame(sample = paste0("iBAQ.", samples),
                              group = gsub("_IP.*", "", samples))
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                                 sampleAnnot = sampleAnnot, mergeGroups = list()),
                 "Some samples are missing")

    ## Wrong iColPattern - should fail
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    expect_error(addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQQ\\.",
                                 sampleAnnot = sampleAnnot, mergeGroups = list()),
                 "Some samples are missing")

    ## Merge groups
    sampleAnnot <- data.frame(sample = samples,
                              group = gsub("_IP.*", "", samples))
    qft1 <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                            sampleAnnot = sampleAnnot,
                            mergeGroups = list(G1 = c("Adnp", "RBC_ctrl"),
                                               G2 = "Chd4BF"))
    cdt1 <- SummarizedExperiment::colData(qft1)
    expect_s4_class(qft1, "QFeatures")
    expect_s4_class(cdt1, "DFrame")
    expect_named(cdt1, c("sample", "group_orig", "group"))
    expect_equal(cdt1$sample, sampleAnnot$sample)
    expect_equal(cdt1$group_orig, sampleAnnot$group)
    expect_true(all(cdt1$group[cdt1$group_orig %in% c("Adnp", "RBC_ctrl")] == "G1"))
    expect_true(all(cdt1$group[cdt1$group_orig %in% c("Chd4BF")] == "G2"))
})
