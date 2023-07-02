#' List available complex DBs
#'
#' @param dbDir Character scalar pointing to the database directory to search
#'     for complex DBs.
#'
#' @author Charlotte Soneson
#'
#' @export
#'
#' @returns A \code{data.frame} with the path to available complex databases,
#'     sorted with the most recently generated first.
#'
#' @examples
#' listComplexDBs()
#'
#' @importFrom dplyr mutate desc %>% arrange
#' @importFrom rlang .data
#'
listComplexDBs <- function(dbDir = system.file("extdata/complexes",
                                               package = "einprot")) {
    .assertScalar(x = dbDir, type = "character")

    complexdbs <- list.files(dbDir, pattern = "complexdb_einprot.*_orthologs",
                             recursive = TRUE, full.names = TRUE)
    data.frame(complexDbPath = complexdbs) %>%
        dplyr::mutate(genDate = gsub("complexdb_einprot.*_", "",
                                     gsub("(_orthologs)*.rds", "",
                                          basename(.data$complexDbPath)))) %>%
        dplyr::mutate(genDate = as.Date(.data$genDate, format = "%Y%m%d")) %>%
        dplyr::arrange(dplyr::desc(.data$genDate))
}
