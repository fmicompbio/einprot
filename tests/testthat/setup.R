## Preparation of MQ SCE object
## -----------------------------------------------------------------------------
mqFile <- system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                      package = "einprot")
mqSamples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
               "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
               "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
mqOut <- importExperiment(inFile = mqFile, iColPattern = "^iBAQ\\.",
                          includeOnlySamples = mqSamples, nrows = 150)
mqsce <- mqOut$sce
sce_mq_initial <- mqsce   ## initial object after import
mqaName <- mqOut$aName
mqSampleAnnot <- data.frame(sample = mqSamples,
                            group = gsub("_IP.*", "", mqSamples))
mqsce <- addSampleAnnots(mqsce, sampleAnnot = mqSampleAnnot)
mqsce <- fixFeatureIds(
    mqsce,
    colDefs = list(einprotId = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs")),
                   einprotLabel = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs")),
                   einprotGene = function(df) getFirstId(df, "Gene.names"),
                   einprotProtein = "Majority.protein.IDs",
                   IDsForSTRING = function(df) combineIds(df, combineCols = c("Gene.names", "Majority.protein.IDs"),
                                                          combineWhen = "missing", makeUnique = FALSE))
)
rownames(mqsce) <- rowData(mqsce)$einprotId
SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName)) <-
    log2(SummarizedExperiment::assay(mqsce, mqaName))
SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName, "_withNA")) <-
    log2(SummarizedExperiment::assay(mqsce, mqaName))
tmp <- SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName))
tmp <- !is.finite(tmp)
SummarizedExperiment::assay(mqsce, paste0("imputed_", mqaName)) <- tmp
SummarizedExperiment::assay(mqsce, mqaName)[SummarizedExperiment::assay(mqsce, mqaName) == 0] <- NA
SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName))[!is.finite(SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName)))] <- NA
SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName, "_withNA"))[
    !is.finite(SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName, "_withNA")))] <- NA
sce_mq_preimputation <- mqsce   ## fixed features, NAs for missing values
nbr_na_mq <- QFeatures::nNA(mqsce)
nbr_na_mq <- lapply(nbr_na_mq, function(a) {
    a$assay <- SummarizedExperiment::assayNames(mqsce)[1]
    a
})
set.seed(123)
SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName)) <- MsCoreUtils::impute_matrix(
    SummarizedExperiment::assay(mqsce, paste0("log2_", mqaName)), method = "MinProb"
)
sce_mq_final <- mqsce   ## final object
fcoll_mq_final <- prepareFeatureCollections(
    sce = mqsce, idCol = "Gene.names",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220323_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 1)

## Preparation of PD TMT SCE object
## -----------------------------------------------------------------------------
pdOutputFolder <- system.file("extdata", "pdtmt_example", package = "einprot")
pdResultName <- "Fig2_m23139_RTS_QC_varMods"
pdAnalysisFile <- system.file("extdata", "pdtmt_example",
                              "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
                              package = "einprot")
pdSampleAnnot = data.frame(
    sample = c("HIS4KO_S05", "HIS4KO_S06", "HIS4KO_S07", "HIS4KO_S08",
               "MET6KO_S01", "MET6KO_S02", "MET6KO_S03", "MET6KO_S04",
               "URA2KO_S09", "URA2KO_S10", "URA2KO_S11", "URA2KO_S12",
               "WT_S13", "WT_S14", "WT_S15", "WT_S16"),
    group = c(rep("HIS4KO", 4), rep("MET6KO", 4), rep("URA2KO", 4),
              rep("WT", 4)))
pdOut <- importExperiment(
    inFile = file.path(pdOutputFolder, paste0(pdResultName, "_Proteins.txt")),
    iColPattern = "^Abundance\\.F.+\\.Sample\\.", nrows = 70)
pdsce <- pdOut$sce
sce_pd_initial <- pdsce
pdaName <- pdOut$aName
pdsce <- addSampleAnnots(pdsce, sampleAnnot = pdSampleAnnot)
pdsce <- fixFeatureIds(
    pdsce,
    colDefs = list(einprotId = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession")),
                   einprotLabel = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession")),
                   einprotGene = function(df) getFirstId(df, "Gene.Symbol"),
                   einprotProtein = "Accession",
                   IDsForSTRING = function(df) combineIds(df, combineCols = c("Gene.Symbol", "Accession"),
                                                          combineWhen = "missing", makeUnique = FALSE))
)
rownames(pdsce) <- rowData(pdsce)$einprotId
SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName)) <-
    log2(SummarizedExperiment::assay(pdsce, pdaName))
SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName, "_withNA")) <-
    log2(SummarizedExperiment::assay(pdsce, pdaName))
tmp <- SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName))
tmp <- !is.finite(tmp)
SummarizedExperiment::assay(pdsce, paste0("imputed_", pdaName)) <- tmp
SummarizedExperiment::assay(pdsce, pdaName)[SummarizedExperiment::assay(pdsce, pdaName) == 0] <- NA
SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName))[!is.finite(SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName)))] <- NA
SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName, "_withNA"))[
    !is.finite(SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName, "_withNA")))] <- NA
sce_pd_preimputation <- pdsce
nbr_na_pd <- QFeatures::nNA(pdsce)
nbr_na_pd <- lapply(nbr_na_pd, function(a) {
    a$assay <- SummarizedExperiment::assayNames(pdsce)[1]
    a
})
set.seed(123)
SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName)) <- MsCoreUtils::impute_matrix(
    SummarizedExperiment::assay(pdsce, paste0("log2_", pdaName)), method = "MinProb"
)
sce_pd_final <- pdsce
fcoll_pd_final <- prepareFeatureCollections(
    sce = pdsce, idCol = "Gene.Symbol",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220323_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("baker's yeast"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 2)

## -----------------------------------------------------------------------------
## Read also the peptidegroups file
pdOut <- importExperiment(
    inFile = file.path(pdOutputFolder, paste0(pdResultName, "_PeptideGroups.txt")),
    iColPattern = "^Abundance\\.F.+\\.Sample\\.", nrows = 70)
pdsce <- pdOut$sce
sce_pd_peptide_initial <- pdsce

## Preparation of FP SCE object
## -----------------------------------------------------------------------------
fragpipeDir <- system.file("extdata", "fp_example",
                           package = "einprot")
fpSamples <- c("Adnp_IP04", "Adnp_IP05", "Adnp_IP06",
               "Chd4BF_IP07", "Chd4BF_IP08", "Chd4BF_IP09",
               "RBC_ctrl_IP01", "RBC_ctrl_IP02", "RBC_ctrl_IP03")
fpOut <- importExperiment(inFile = file.path(fragpipeDir, "combined_protein.tsv"),
                          iColPattern = "\\.MaxLFQ\\.Intensity$",
                          includeOnlySamples = fpSamples, nrows = 150)
fpsce <- fpOut$sce
sce_fp_initial <- fpsce   ## initial object after import
fpaName <- fpOut$aName
fpSampleAnnot <- data.frame(sample = fpSamples,
                            group = gsub("_IP.*", "", fpSamples))
fpsce <- addSampleAnnots(fpsce, sampleAnnot = fpSampleAnnot)
fpsce <- fixFeatureIds(
    fpsce,
    colDefs = list(einprotId = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID")),
                   einprotLabel = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID")),
                   einprotGene = function(df) getFirstId(df, "Gene"),
                   einprotProtein = "Protein.ID",
                   IDsForSTRING = function(df) combineIds(df, combineCols = c("Gene", "Protein.ID"),
                                                          combineWhen = "missing", makeUnique = FALSE))
)
rownames(fpsce) <- rowData(fpsce)$einprotId
SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName)) <-
    log2(SummarizedExperiment::assay(fpsce, fpaName))
SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName, "_withNA")) <-
    log2(SummarizedExperiment::assay(fpsce, fpaName))
tmp <- SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName))
tmp <- !is.finite(tmp)
SummarizedExperiment::assay(fpsce, paste0("imputed_", fpaName)) <- tmp
SummarizedExperiment::assay(fpsce, fpaName)[SummarizedExperiment::assay(fpsce, fpaName) == 0] <- NA
SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName))[!is.finite(SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName)))] <- NA
SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName, "_withNA"))[
    !is.finite(SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName, "_withNA")))] <- NA
sce_fp_preimputation <- fpsce   ## fixed features, NAs for missing values
nbr_na_fp <- QFeatures::nNA(fpsce)
nbr_na_fp <- lapply(nbr_na_fp, function(a) {
    a$assay <- SummarizedExperiment::assayNames(fpsce)[1]
    a
})
set.seed(123)
SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName)) <- MsCoreUtils::impute_matrix(
    SummarizedExperiment::assay(fpsce, paste0("log2_", fpaName)), method = "MinProb"
)
sce_fp_final <- fpsce   ## final object
fcoll_fp_final <- prepareFeatureCollections(
    sce = fpsce, idCol = "Gene",
    includeFeatureCollections = "complexes",
    complexDbPath = system.file("extdata", "complexes",
                                "complexdb_einprot0.5.0_20220323_orthologs.rds",
                                package = "einprot"),
    speciesInfo = getSpeciesInfo("mouse"), complexSpecies = "current",
    customComplexes = list(), minSizeToKeep = 1)


