#' Check whether pandoc and pandoc-citeproc is available
#'
#' @author Charlotte Soneson
#'
#' @param ignorePandoc logical. If TRUE, just give a warning if pandoc or
#'   pandoc-citeproc is not available. If FALSE, stop.
#'
#' @keywords internal
#' @noRd
#'
#' @return A logical(1), indicating whether pandoc can be run or not.
#'   In addition, raises either a warning or an error (depending on the
#'   value of \code{ignorePandoc}) if pandoc or pandoc-citeproc is not
#'   available.
#'
#' @importFrom rmarkdown pandoc_available pandoc_exec
#'
.checkPandoc <- function(ignorePandoc) {
    ## Initialize output to TRUE
    doRender <- TRUE

    ## First check whether pandoc is available
    if (!rmarkdown::pandoc_available()) {
        doRender <- FALSE
        ## If pandoc is not available, either give a warning or an error,
        ## depending on the value of ignorePandoc
        if (ignorePandoc) {
            ## If ignorePandoc is TRUE, just give a warning
            warning("pandoc is not available! ",
                    "HTML reports can not be generated.",
                    immediate. = TRUE)
        } else {
            ## If ignorePandoc is FALSE, stop
            stop("pandoc is not available!")
        }
    } else {
        ## If pandoc is available, check for pandoc-citeproc
        ## Only do this if the pandoc version is <2.11, since
        ## pandoc-citeproc is not included (or needed) in v2.11 and later.
        if (!rmarkdown::pandoc_available(version = "2.11")) {
            ## TRUE if the available pandoc version is not 2.11 or newer
            ## pandoc-citeproc should be found in the path, or in the
            ## same folder as the pandoc executable
            if (Sys.which("pandoc-citeproc") == "" &&
                !file.exists(file.path(dirname(rmarkdown::pandoc_exec()),
                                       "pandoc-citeproc"))) {
                doRender <- FALSE
                ## pandoc-citeproc is required, but not found
                if (ignorePandoc) {
                    ## If ignorePandoc is TRUE, just give a warning
                    warning("pandoc-citeproc is not available! ",
                            "HTML reports can not be generated.",
                            immediate. = TRUE)
                } else {
                    ## If ignorePandoc is FALSE, stop
                    stop("pandoc-citeproc is not available!")
                }
            }
        }
    }
    return(doRender)
}
