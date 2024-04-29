#' Extract information from Spectronaut setup file
#'
#' Read \code{Spectronaut} setup.txt file and extract information about the
#' run.
#'
#' @param spectronautSetupFile Character scalar, the path to a
#'     \code{Spectronaut} setup file. Can be \code{NULL} (in this case, an
#'     empty list is returned).
#'
#' @author Charlotte Soneson
#' @export
#'
#' @returns A list with extracted information about the Spectronaut run.
#'
readSpectronautSetup <- function(spectronautSetupFile) {
    .assertScalar(x = spectronautSetupFile, type = "character",
                  allowNULL = TRUE)

    if (is.null(spectronautSetupFile)) {
        return(list())
    } else {
        if (!file.exists(spectronautSetupFile)) {
            stop(spectronautSetupFile, " doesn't exist.")
        }

        ## Read file and extract settings and setup part
        sf <- readLines(spectronautSetupFile)
        settings <- sf[seq(which(sf == "[BEGIN-SETTINGS]") + 1,
                           which(sf == "[END-SETTINGS]") - 1)]
        setup <- sf[seq(which(sf == "[BEGIN-SETUP]") + 1,
                        which(sf == "[END-SETUP]") - 1)]

        ## Version
        sp_version <- grep("^Spectronaut ", sf, value = TRUE)[1]

        ## Raw files
        sp_raw_files <- paste(sub("Run: ", "", grep("^Run: ", setup,
                                                    value = TRUE)),
                              collapse = ", ")

        ## Protein databases
        sp_fasta_files <- unique(sub(".*Original File: ", "",
                                     grep("Original File: ", setup,
                                          value = TRUE)))
        sp_fasta_files <- paste(sp_fasta_files, collapse = ", ")

        L <- list(
            "Spectronaut version" = sp_version,
            "Setup file" = spectronautSetupFile,
            "Raw files" = sp_raw_files,
            "Databases" = sp_fasta_files
        )
        return(L[!vapply(L, is.null, TRUE)])
    }
}
