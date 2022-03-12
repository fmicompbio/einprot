## Preparation of MQ SummarizedExperiment object
## ------------------------------------------------------------------------- ##
mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                      package = "einprot")
samples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
             "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
             "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
out <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                        includeOnlySamples = samples, nrows = 70)
sce_mq_initial <- out$sce
aName <- out$aName
sampleAnnot <- data.frame(sample = samples,
                          group = gsub("_IP.*", "", samples))
sce <- addSampleAnnots(sce_mq_initial, sampleAnnot = sampleAnnot, mergeGroups = list())
sce <- fixFeatureIds(sce, geneIdCol = "Gene.names",
                     proteinIdCol = "Majority.protein.IDs")
SummarizedExperiment::assay(sce, paste0("log2_", aName)) <-
    log2(SummarizedExperiment::assay(sce, aName))
SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA")) <-
    log2(SummarizedExperiment::assay(sce, aName))
tmp <- SummarizedExperiment::assay(sce, paste0("log2_", aName))
tmp <- !is.finite(tmp)
SummarizedExperiment::assay(sce, paste0("imputed_", aName)) <- tmp
SummarizedExperiment::assay(sce, aName)[SummarizedExperiment::assay(sce, aName) == 0] <- NA
SummarizedExperiment::assay(sce, paste0("log2_", aName))[!is.finite(SummarizedExperiment::assay(sce, paste0("log2_", aName)))] <- NA
SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA"))[
    !is.finite(SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA")))] <- NA
sce_mq_preimputation <- sce
nbr_na_mq <- QFeatures::nNA(sce)
nbr_na_mq <- lapply(nbr_na_mq, function(a) {
    a$assay <- SummarizedExperiment::assayNames(sce)[1]
    a
})
set.seed(123)
SummarizedExperiment::assay(sce, paste0("log2_", aName)) <- MsCoreUtils::impute_matrix(
    SummarizedExperiment::assay(sce, paste0("log2_", aName)), method = "MinProb"
)
sce_mq_final <- sce
fcoll_mq_final <- prepareFeatureCollections(
    sce = sce_mq_final, idCol = "Gene.names",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 2)

## Preparation of PD TMT SummarizedExperiment object
## ------------------------------------------------------------------------- ##
pdOutputFolder <- system.file("extdata", "pdtmt_example", package = "einprot")
pdResultName <- "Fig2_m23139_RTS_QC_varMods"
pdAnalysisFile <- system.file("extdata", "pdtmt_example",
                              "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
                              package = "einprot")
iColPattern = "^Abundance\\.F.+\\.Sample\\."
sampleAnnot = data.frame(
    sample = c("HIS4KO_S05", "HIS4KO_S06", "HIS4KO_S07", "HIS4KO_S08",
               "MET6KO_S01", "MET6KO_S02", "MET6KO_S03", "MET6KO_S04",
               "URA2KO_S09", "URA2KO_S10", "URA2KO_S11", "URA2KO_S12",
               "WT_S13", "WT_S14", "WT_S15", "WT_S16"),
    group = c(rep("HIS4KO", 4), rep("MET6KO", 4), rep("URA2KO", 4),
              rep("WT", 4)))
out <- importExperiment(
    inFile = file.path(pdOutputFolder, paste0(pdResultName, "_Proteins.txt")),
    iColPattern = iColPattern, nrows = 70)
sce_pd_initial <- out$sce
aName <- out$aName
sce <- addSampleAnnots(sce_pd_initial, sampleAnnot = sampleAnnot, mergeGroups = list())
sce <- fixFeatureIds(sce, geneIdCol = "Gene.Symbol",
                     proteinIdCol = "Accession")
SummarizedExperiment::assay(sce, paste0("log2_", aName)) <-
    log2(SummarizedExperiment::assay(sce, aName))
SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA")) <-
    log2(SummarizedExperiment::assay(sce, aName))
tmp <- SummarizedExperiment::assay(sce, paste0("log2_", aName))
tmp <- !is.finite(tmp)
SummarizedExperiment::assay(sce, paste0("imputed_", aName)) <- tmp
SummarizedExperiment::assay(sce, aName)[SummarizedExperiment::assay(sce, aName) == 0] <- NA
SummarizedExperiment::assay(sce, paste0("log2_", aName))[!is.finite(SummarizedExperiment::assay(sce, paste0("log2_", aName)))] <- NA
SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA"))[
    !is.finite(SummarizedExperiment::assay(sce, paste0("log2_", aName, "_withNA")))] <- NA
sce_pd_preimputation <- sce
nbr_na_pd <- QFeatures::nNA(sce)
nbr_na_pd <- lapply(nbr_na_pd, function(a) {
    a$assay <- SummarizedExperiment::assayNames(sce)[1]
    a
})
set.seed(123)
SummarizedExperiment::assay(sce, paste0("log2_", aName)) <- MsCoreUtils::impute_matrix(
    SummarizedExperiment::assay(sce, paste0("log2_", aName)), method = "MinProb"
)
sce_pd_final <- sce
fcoll_pd_final <- prepareFeatureCollections(
    sce = sce_pd_final, idCol = "Gene.Symbol",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220211_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("baker's yeast"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 2)
