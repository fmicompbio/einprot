library(einprot)
sampleAnnot <- read.delim(system.file("extdata", "mq_example",
                                      "1356_sampleAnnot.txt",
                                      package = "einprot"))
out <- runMaxQuantAnalysis(
    outputDir = "~/Desktop/1356_mq",
    outputBaseName = "1356",
    forceOverwrite = TRUE,
    species = "mouse",
    mqFile = system.file("extdata", "mq_example", "1356_proteinGroups.txt",
                         package = "einprot"),
    mqParameterFile = system.file("extdata", "mq_example", "1356_mqpar.xml",
                                  package = "einprot"),
    iColPattern = "^LFQ.intensity.",
    sampleAnnot = sampleAnnot,
    ctrlGroup = "RBC_ctrl",
    stringIdCol = NULL,
    includeFeatureCollections = "complexes"
)
