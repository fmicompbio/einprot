#' Generate config chunk to paste in Rmd template
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
.generateConfigChunk <- function(configlist) {
    .assertVector(x = configlist, type = "list")
    .assertVector(x = names(configlist), type = "character")

    outstr <- paste0("```{r config, eval = TRUE}\n",
                     "## The following variables were specified as input ",
                     "arguments when calling the rendering function.\n",
                     "## They will be used in the workflow below.\n\n")
    for (nm in names(configlist)) {
        outstr <- paste0(outstr, nm, " <- ",
                         paste(deparse(configlist[[nm]]), collapse = ""), "\n")
    }
    outstr <- paste0(outstr, "```")

    outstr
}
