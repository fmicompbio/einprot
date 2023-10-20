#' Extract information from DIA-NN logs
#'
#' Read DIA-NN log file and extract information about the run.
#'
#' @param diannLog Character scalar, the path to a log file from a DIA-NN run.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @examples
#' readDIANNInfo(system.file("extdata", "diann_example",
#'                           "diann-output.log.txt",
#'                           package = "einprot"))
#'
#' @returns A list with extracted information about the \code{DIA-NN} run.
#'
readDIANNInfo <- function(diannLog) {
    ## Based on https://github.com/vdemichev/DiaNN#command-line-reference
    .assertScalar(x = diannLog, type = "character", allowNULL = TRUE)

    if (is.null(diannLog)) {
        return(list())
    } else {
        if (!file.exists(diannLog)) {
            stop(diannLog, " doesn't exist.")
        }

        rl <- readLines(diannLog)

        ## -------------------------------------------------------------------------
        ## Version
        ## -------------------------------------------------------------------------
        diannVersion <- strsplit(rl[1], "\\(")[[1]][1]

        ## -------------------------------------------------------------------------
        ## Split command
        ## -------------------------------------------------------------------------
        cmd <- grep("--", rl, value = TRUE)[1]

        cmdspl <- strsplit(cmd, "--")[[1]]

        ## -------------------------------------------------------------------------
        ## Configuration file
        ## -------------------------------------------------------------------------
        cfg <- sub(" $", "", sub("^cfg ", "", grep("^cfg ", cmdspl,
                                                   value = TRUE)))

        ## -------------------------------------------------------------------------
        ## Mass accuracy
        ## -------------------------------------------------------------------------
        massAcc <- sub("mass-acc ", "", grep("mass-acc ", cmdspl, value = TRUE))
        massAccMS1 <- sub("mass-acc-ms1 ", "", grep("mass-acc-ms1 ", cmdspl,
                                                    value = TRUE))

        ## -------------------------------------------------------------------------
        ## Input files
        ## -------------------------------------------------------------------------
        inputFiles <- sub(" $", "", sub("^f ", "", grep("^f ", cmdspl,
                                                        value = TRUE)))

        ## -------------------------------------------------------------------------
        ## Spectral library
        ## -------------------------------------------------------------------------
        lib <- sub(" $", "", sub("^lib ", "", grep("^lib ", cmdspl, value = TRUE)))
        if (length(lib) == 1 && lib == "") lib <- c()

        genspeclib <- any(cmdspl == "gen-spec-lib ")
        if (!genspeclib) genspeclib <- c()

        ## -------------------------------------------------------------------------
        ## Sequence databases
        ## -------------------------------------------------------------------------
        fastaFiles <- sub(" $", "", sub("^fasta ", "", grep("^fasta ", cmdspl,
                                                            value = TRUE)))
        if (length(fastaFiles) != 0) {
            fastaFiles <- paste(fastaFiles, collapse = ", ")
        }

        seqdb <- grep("Library annotated with sequence", rl, value = TRUE)[1]
        seqdb <- strsplit(seqdb, "sequence database(s): ", fixed = TRUE)[[1]][2]
        if (is.na(seqdb)) seqdb <- c()

        ## -------------------------------------------------------------------------
        ## q-value
        ## -------------------------------------------------------------------------
        qval <- sub(" $", "", sub("^qvalue ", "", grep("^qvalue ", cmdspl,
                                                       value = TRUE)))

        ## -------------------------------------------------------------------------
        ## Deep learning predictor
        ## -------------------------------------------------------------------------
        dopred <- any(cmdspl == "predictor ")
        if (!dopred) {
            dopred <- noionmob <- nortpred <- c()
        } else {
            noionmob <- any(cmdspl == "dl-no-im ")
            if (!noionmob) noionmob <- c()
            nortpred <- any(cmdspl == "dl-no-rt ")
            if (!nortpred) nortpred <- c()
        }

        ## -------------------------------------------------------------------------
        ## Other specifications
        ## -------------------------------------------------------------------------
        fastasearch <- any(cmdspl == "fasta-search ")
        if (!fastasearch) fastasearch <- c()

        relaxedprotinf <- any(cmdspl == "relaxed-prot-inf ")
        if (!relaxedprotinf) relaxedprotinf <- c()

        smartprof <- any(cmdspl == "smart-profiling ")
        if (!smartprof) smartprof <- c()

        globalmasscal <- any(cmdspl == "global-mass-cal")
        if (!globalmasscal) globalmasscal <- c()

        globalnorm <- any(cmdspl == "global-norm")
        if (!globalnorm) globalnorm <- c()

        ## -------------------------------------------------------------------------
        ## Quantification mode
        ## -------------------------------------------------------------------------
        if (any(cmdspl == "peak-center ")) quantmode <- "Robust LC (peak center)"
        else if (any(cmdspl == "peak-height ")) quantmode <- "Peak height"
        else quantmode <- c()

        noifsrem <- any(cmdspl == "no-ifs-removal ")
        if (!noifsrem) noifsrem <- c()

        ## -------------------------------------------------------------------------
        ## Cleavage specification
        ## -------------------------------------------------------------------------
        cleavagespec <- sub(" $", "", sub("^cut ", "", grep("^cut ", cmdspl,
                                                            value = TRUE)))
        maxmisscleavages <- sub(" $", "", sub("^missed-cleavages ", "",
                                              grep("^missed-cleavages ", cmdspl,
                                                   value = TRUE)))

        ## -------------------------------------------------------------------------
        ## M/Z, charge, length ranges
        ## -------------------------------------------------------------------------
        fr_mz_range <- paste0(
            "[",
            sub(" $", "", sub("^min-fr-mz ", "", grep("^min-fr-mz ", cmdspl,
                                                      value = TRUE))),
            ", ",
            sub(" $", "", sub("^max-fr-mz ", "", grep("^max-fr-mz ", cmdspl,
                                                      value = TRUE))),
            "]")
        if (fr_mz_range == "[, ]") fr_mz_range <- c()

        pep_len_range <- paste0(
            "[",
            sub(" $", "", sub("^min-pep-len ", "", grep("^min-pep-len ", cmdspl,
                                                        value = TRUE))),
            ", ",
            sub(" $", "", sub("^max-pep-len ", "", grep("^max-pep-len ", cmdspl,
                                                        value = TRUE))),
            "]")
        if (pep_len_range == "[, ]") pep_len_range <- c()

        pr_mz_range <- paste0(
            "[",
            sub(" $", "", sub("^min-pr-mz ", "", grep("^min-pr-mz ", cmdspl,
                                                      value = TRUE))),
            ", ",
            sub(" $", "", sub("^max-pr-mz ", "", grep("^max-pr-mz ", cmdspl,
                                                      value = TRUE))),
            "]")
        if (pr_mz_range == "[, ]") pr_mz_range <- c()

        pr_charge_range <- paste0(
            "[",
            sub(" $", "", sub("^min-pr-charge ", "", grep("^min-pr-charge ", cmdspl,
                                                          value = TRUE))),
            ", ",
            sub(" $", "", sub("^max-pr-charge ", "", grep("^max-pr-charge ", cmdspl,
                                                          value = TRUE))),
            "]")
        if (pr_charge_range == "[, ]") pr_charge_range <- c()

        ## -------------------------------------------------------------------------
        ## Modifications
        ## -------------------------------------------------------------------------
        varmod <- sub(" $", "", sub("^var-mod ", "", grep("^var-mod ", cmdspl,
                                                          value = TRUE)))
        maxvarmods <- sub(" $", "", sub("^var-mods ", "", grep("^var-mods ", cmdspl,
                                                               value = TRUE)))
        metexcision <- any(cmdspl == "met-excision ")
        if (!metexcision) metexcision <- c()

        unimods <- sub(" $", "", sub("^unimod", "", grep("^unimod", cmdspl,
                                                         value = TRUE)))
        if (length(unimods) > 0) {
            unimods <- c("4" = "Carbamidomethylation",
                         "1" = "Acetyl (N-term)",
                         "21" = "Phospho",
                         "31" = "Oxidation (M)",
                         "121" = "K-GG")[unimods]
            modifications <- paste(unimods, collapse = ", ")
        } else {
            modifications <- c()
        }

        fixedmod <- sub(" $", "", sub("^fixed-mod ", "", grep("^fixed-mod ", cmdspl,
                                                              value = TRUE)))

        ## -------------------------------------------------------------------------
        ## Match-between-runs
        ## -------------------------------------------------------------------------
        mbr <- any(cmdspl == "reanalyse ")
        if (!mbr) mbr <- c()

        ## -------------------------------------------------------------------------
        ## Return list of values
        ## -------------------------------------------------------------------------
        l <- list("DIA-NN version" = diannVersion,
                  "DIA-NN log file" = diannLog,
                  "Configuration file" = cfg,
                  "MS2 Mass accuracy (ppm)" = massAcc,
                  "MS1 Mass accuracy (ppm)" = massAccMS1,
                  "Q-value" = qval,
                  "Perform deep learning-based prediction" = dopred,
                  "Don't perform prediction of ion mobilities" = noionmob,
                  "Don't perform prediction of retention times" = nortpred,
                  "Input files" = paste(inputFiles, collapse = ", "),
                  "Generate spectral library" = genspeclib,
                  "Spectral library" = lib,
                  "Perform in silico digest of sequence database" = fastasearch,
                  "FASTA files" = fastaFiles,
                  "Sequence databases" = seqdb,
                  "Relaxed protein inference" = relaxedprotinf,
                  "Enable smart profiling" = smartprof,
                  "Quantification mode" = quantmode,
                  "Turn off interference subtraction" = noifsrem,
                  "Cleavage specification" = cleavagespec,
                  "Max nbr of missed cleavages" = maxmisscleavages,
                  "Fragment M/Z range" = fr_mz_range,
                  "Peptide length range" = pep_len_range,
                  "Precursor M/Z range" = pr_mz_range,
                  "Precursor charge range" = pr_charge_range,
                  "Enable N-term methionine excision" = metexcision,
                  "Fixed modifications" = fixedmod,
                  "Variable modifications" = varmod,
                  "Modifications" = modifications,
                  "Max nbr of variable modifications" = maxvarmods,
                  "Match between runs" = mbr,
                  "Global mass calibration" = globalmasscal,
                  "Global normalization" = globalnorm,
                  "DIA-NN command" = cmd)
        return(l[vapply(l, length, 0) > 0])
    }
}
