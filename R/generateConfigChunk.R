#' Generate config chunk to paste in Rmd template
#'
#' @keywords internal
#' @noRd
#' @author Charlotte Soneson
.generateConfigChunk <- function(configlist) {
    outstr <- paste0("```{r config, eval = TRUE}\n",
                     "## The following parameters were obtained from ",
                     "the input arguments to the rendering function.\n\n")
    for (nm in names(configlist)) {
        outstr <- paste0(outstr, nm, " <- ",
                         paste(deparse(configlist[[nm]]), collapse = ""), "\n")
    }
    outstr <- paste0(outstr, "```")

    outstr
}
