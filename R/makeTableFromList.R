#' Make a nice-looking table from a named list
#'
#' @param l Named list. Each element of the list must be a scalar value.
#'
#' @author Charlotte Soneson
#' @export
#'
#' @return A \code{kableExtra} table.
#'
#' @examples
#' makeTableFromList(list(first = "one", second = "two", third = 3))
#'
#' @importFrom tibble rownames_to_column
#' @importFrom dplyr %>%
#' @importFrom stats setNames
#' @importFrom kableExtra kbl kable_paper column_spec
#'
makeTableFromList <- function(l) {
    ## --------------------------------------------------------------------- ##
    ## Check arguments
    ## --------------------------------------------------------------------- ##
    .assertVector(x = l, type = "list")
    .assertVector(x = names(l), type = "character", allowNULL = FALSE)
    .assertVector(x = vapply(l, length, 0), type = "numeric",
                  validValues = c(0, 1))

    ## --------------------------------------------------------------------- ##
    ## Turn into a data.frame
    ## --------------------------------------------------------------------- ##
    df <- data.frame(
        c2 = unlist(l)
    ) %>%
        tibble::rownames_to_column("c1") %>%
        stats::setNames(NULL)

    ## --------------------------------------------------------------------- ##
    ## Format
    ## --------------------------------------------------------------------- ##
    kableExtra::kbl(df) %>%
        kableExtra::kable_paper(
            full_width = TRUE, lightable_options = c("striped"),
            html_font = "\"Trebuchet MS\", verdana, sans-serif") %>%
        kableExtra::column_spec(column = 1, width = "25em", bold = TRUE,
                                border_right = TRUE) %>%
        kableExtra::column_spec(column = 2, width = "100em")
}
