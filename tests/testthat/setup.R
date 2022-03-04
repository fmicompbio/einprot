## Preparation of MQ QFeatures object
## ------------------------------------------------------------------------- ##
mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                      package = "einprot")
samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
             "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
             "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
ecol <- paste0("iBAQ.", samples)
qft_mq_initial <- QFeatures::readQFeatures(mqFile, ecol = ecol, name = "iBAQ",
                                           sep = "\t", nrows = 70)
sampleAnnot <- data.frame(sample = samples,
                          group = gsub("_IP.*", "", samples))
qft <- addSampleAnnots(qft_mq_initial, iColPattern = "^iBAQ\\.",
                       sampleAnnot = sampleAnnot, mergeGroups = list())
qft <- fixFeatureIds(qft)
qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ")
qft <- QFeatures::logTransform(qft, base = 2, i = "iBAQ", name = "log2_iBAQ_withNA")
tmp <- qft[["log2_iBAQ"]]
SummarizedExperiment::assay(tmp) <- !is.finite(SummarizedExperiment::assay(tmp))
qft <- QFeatures::addAssay(qft, tmp, name = "imputed_iBAQ")
qft <- QFeatures::zeroIsNA(qft, "iBAQ")
qft <- QFeatures::infIsNA(qft, "log2_iBAQ")
qft <- QFeatures::infIsNA(qft, "log2_iBAQ_withNA")
qft_mq_preimputation <- qft
nbr_na_mq <- QFeatures::nNA(qft_mq_preimputation, i = seq_along(qft_mq_preimputation))
set.seed(123)
qft_mq_final <- QFeatures::impute(qft_mq_preimputation, method = "MinProb", i = "log2_iBAQ")
fcoll_mq_final <- prepareFeatureCollections(
    qft = qft_mq_final, idCol = "Gene.names",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 2)
