## Define default filenames for files provided via einprot
## Adapt when a new version is generated

#' Constants
#'
#' Constant values used in einprot workflows, representing paths to default
#' complex database and ID conversion tables.
#'
#' @author Charlotte Soneson
#'
#' @returns Paths to built-in complex database and conversion tables.
#'
#' @name constants
NULL

#' @export
#' @rdname constants
EINPROT_COMPLEXES_FILE <- "extdata/complexes/complexdb_einprot0.9.3_20240328_orthologs.rds"
#' @export
#' @rdname constants
EINPROT_WORMBASE_CONVTABLE <- "extdata/conversion_tables/WormBaseConv_einprot0.5.0_20220211.rds"
#' @export
#' @rdname constants
EINPROT_POMBASE_CONVTABLE <- "extdata/conversion_tables/PomBaseConv_einprot0.5.0_20220211.rds"
