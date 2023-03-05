#' Generate config chunk to paste in Rmd template
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
.generateConfigChunk <- function(configlist, fcn = "dump") {
    .assertVector(x = configlist, type = "list")
    .assertVector(x = names(configlist), type = "character")
    .assertScalar(x = fcn, type = "character", validValues = c("deparse", "dump"))

    outstr <- paste0("```{r config, eval = TRUE}\n",
                     "## The following variables were specified as input ",
                     "arguments when calling the rendering function.\n",
                     "## They will be used in the workflow below.\n\n")

    if (fcn == "deparse") {
        for (nm in names(configlist)) {
            outstr <- paste0(outstr, nm, " <- ",
                             paste(deparse(configlist[[nm]], control = "all"),
                                   collapse = ""), "\n")
        }
        outstr <- paste0(outstr, "```")
    } else if (fcn == "dump") {
        tmpf <- tempfile()
        dump(names(configlist), file = tmpf, envir = as.environment(configlist))
        tmp <- readLines(tmpf)
        outstr <- paste(c(outstr, paste(tmp, collapse = "\n"), "```"),
                        collapse = "\n")
    }

    outstr
}
