test_that("readProteomeDiscovererInfo works", {
    pdOutputFolder <- system.file("extdata", "pdtmt_example",
                                  package = "einprot")
    pdResultName <- "Fig2_m23139_RTS_QC_varMods"
    pdAnalysisFile <- system.file("extdata", "pdtmt_example",
                                  "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
                                  package = "einprot")

    ## Fails with wrong inputs
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = 1, pdResultName = pdResultName,
        pdAnalysisFile = pdAnalysisFile),
        "'pdOutputFolder' must be of class 'character'")
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = c(pdOutputFolder, pdOutputFolder),
        pdResultName = pdResultName,
        pdAnalysisFile = pdAnalysisFile),
        "'pdOutputFolder' must have length 1")
    expect_equal(readProteomeDiscovererInfo(
        pdOutputFolder = "missing", pdResultName = pdResultName,
        pdAnalysisFile = pdAnalysisFile),
        list())
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = 1,
        pdAnalysisFile = pdAnalysisFile),
        "'pdResultName' must be of class 'character'")
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder,
        pdResultName = c(pdResultName, pdResultName),
        pdAnalysisFile = pdAnalysisFile),
        "'pdResultName' must have length 1")
    expect_equal(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = "missing",
        pdAnalysisFile = pdAnalysisFile),
        list())
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
        pdAnalysisFile = 1),
        "'pdAnalysisFile' must be of class 'character'")
    expect_error(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
        pdAnalysisFile = c(pdAnalysisFile, pdAnalysisFile)),
        "'pdAnalysisFile' must have length 1")
    expect_equal(readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
        pdAnalysisFile = "missing"),
        list())

    ## Works with correct input
    pdi <- readProteomeDiscovererInfo(
        pdOutputFolder = pdOutputFolder, pdResultName = pdResultName,
        pdAnalysisFile = pdAnalysisFile)
    expect_type(pdi, "list")
    expect_equal(pdi$`PD result name`, pdResultName)
    expect_equal(pdi$`PD analysis file`, pdAnalysisFile)
    expect_equal(pdi$`PD version`, "2.5.0.400")
    expect_equal(pdi$`PD Processing WF`,
                 "PWF_Tribrid_TMTpro_Quan_SPS_MS3_SequestHT_Percolator")
    expect_equal(pdi$`PD Consensus WF`,
                 "CWF_Comprehensive_Enhanced Annotation_Reporter_Quan")
    expect_equal(pdi$`Search engine`, "Sequest HT")
    expect_equal(pdi$Instruments, "Orbitrap Fusion")
    expect_equal(pdi$`Sample names`, "HIS4KO_S05, HIS4KO_S06, HIS4KO_S07, HIS4KO_S08, MET6KO_S01, MET6KO_S02, MET6KO_S03, MET6KO_S04, URA2KO_S09, URA2KO_S10, URA2KO_S11, URA2KO_S12, WT_S13, WT_S14, WT_S15, WT_S16")
    expect_equal(pdi$Databases,
                 "CON_iRT_contaminants_cRAPMaxQFMI_150507.fasta; YEAST__210503.fasta")
    expect_equal(pdi$Contaminants, "CON_iRT_contaminants_cRAPMaxQFMI_150507.fasta")
    expect_equal(pdi$Enzymes, "Trypsin (Full)")
    expect_equal(pdi$`Variable modifications`, "Oxidation / +15.995 Da (M), Carbamidomethyl / +57.021 Da (C), TMTpro / +304.207 Da (K, S, T), TMTpro / +304.207 Da (N-Terminus)")
    expect_equal(pdi$`Fixed modifications`, "")
    expect_equal(pdi$`Validation method`, "PercolatorConfidenceAssignment")
    expect_equal(pdi$`Validation based on`, "Target/Decoy, q-Value")
    expect_equal(pdi$`Confidence thresholds`, "strict: 0.01, relaxed: 0.05")
})

test_that("querypdAnalysis functions work", {
    pdOutputFolder <- system.file("extdata", "pdtmt_example",
                                  package = "einprot")
    pdResultName <- "Fig2_m23139_RTS_QC_varMods"
    pdAnalysisFile <- system.file("extdata", "pdtmt_example",
                                  "Fig2_m23139_RTS_QC_varMods.pdAnalysis",
                                  package = "einprot")

    ## getContaminantsDatabaseFrompdAnalysis
    out <- getContaminantsDatabaseFrompdAnalysis(pdAnalysisFile)
    expect_type(out, "list")
    expect_equal(out$contaminantDb,
                 "CON_iRT_contaminants_cRAPMaxQFMI_150507.fasta")
    expect_equal(out$protMarkers, c("Contaminants", "YEAST",
                                    "YEAST__210503.fasta"))

    ## getSearchParametersFrompdAnalysis
    out <- getSearchParametersFrompdAnalysis(pdAnalysisFile)
    expect_type(out, "list")
    expect_equal(out$search_engine, "Sequest HT")
    expect_equal(out$fasta_db, "CON_iRT_contaminants_cRAPMaxQFMI_150507.fasta; YEAST__210503.fasta")
    expect_equal(out$enzymes, "Trypsin (Full)")
    expect_equal(out$dynamicModifications,
                 c("Oxidation / +15.995 Da (M)",
                   "Carbamidomethyl / +57.021 Da (C)",
                   "TMTpro / +304.207 Da (K, S, T)",
                   "TMTpro / +304.207 Da (N-Terminus)"))
    expect_equal(out$staticModifications, character(0))

    ## getQuantOrderFrompdAnalysis
    out <- getQuantOrderFrompdAnalysis(pdAnalysisFile)
    expect_equal(out, "MS3")

    ## getTemplateNamesFrompdAnalysis
    out <- getTemplateNamesFrompdAnalysis(pdAnalysisFile)
    expect_type(out, "list")
    expect_equal(out$Consensus,
                 "CWF_Comprehensive_Enhanced Annotation_Reporter_Quan")
    expect_equal(out$Processing,
                 "PWF_Tribrid_TMTpro_Quan_SPS_MS3_SequestHT_Percolator")

    ## getValidationInfoFrompdAnalysis
    out <- getValidationInfoFrompdAnalysis(pdAnalysisFile)
    expect_type(out, "list")
    expect_equal(out$targetFDRstrictPSM, "0.01")
    expect_equal(out$targetFDRrelaxedPSM, "0.05")
    expect_equal(out$targetFDRstrictPeptide, "0.01")
    expect_equal(out$targetFDRrelaxedPeptide, "0.05")
    expect_equal(out$validationMethod, "PercolatorConfidenceAssignment")
    expect_equal(out$validationBasedOn, "Target/Decoy, q-Value")

    ## getPSMValidationInfoFrompdAnalysis
    out <- getPSMValidationInfoFrompdAnalysis(pdAnalysisFile)
    expect_equal(out, "Percolator")

    ## getCalibrationFrompdAnalysis
    out <- getCalibrationFrompdAnalysis(pdAnalysisFile)
    expect_false(out)

    ## getQuantInfoFrompdAnalysis
    out <- getQuantInfoFrompdAnalysis(pdAnalysisFile)
    expect_equal(out$quant_mode, "Reporter Ions Quantifier")
    expect_equal(out$peptides_to_use, "Unique + Razor")
    expect_equal(out$abundance_based_on, "S/N")
    expect_equal(out$quan_value_corr, "False")
    expect_equal(out$co_isolation_thr, "50")
    expect_equal(out$ave_reporter_sn_thr, "10")
    expect_equal(out$sps_mm_pct_thr, "65")
    expect_equal(out$norm_mode, "Total Peptide Amount")
    expect_equal(out$imputation_mode, "None")
})

