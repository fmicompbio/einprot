## Generate complex DB
## File name contains einprot version as well as geneation date
library(einprot)

dbdir <- tempdir()
paths <- makeComplexDB(dbDir = dbdir, customComplexTxt = NULL)
paths
file.copy(from = paths$orthPath,
          to = file.path(here::here(), "inst", "extdata", "complexes", basename(paths$orthPath)))
