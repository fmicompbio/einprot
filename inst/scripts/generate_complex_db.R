## Generate complex DB
## File name contains einprot version as well as generation date
library(einprot)

dbdir <- tempdir()
paths <- makeComplexDB(dbDir = dbdir, customComplexTxt = NULL)
paths
file.copy(from = paths$orthPath,
          to = file.path(here::here(), "inst", "extdata", "complexes", basename(paths$orthPath)))

## CYC2008 and CORUM give download issues -> use pre-downloaded files
dbdir <- tempdir()
paths <- makeComplexDB(dbDir = dbdir, customComplexTxt = NULL,
                       Cyc2008Db = read.delim("complex_db/20220522/S_cerevisiae_CYC2008_complex.tab"),
                       CorumDb = read.delim("complex_db/20240327/CORUM_allComplexes_20240327.txt"))
paths
file.copy(from = paths$orthPath,
          to = file.path(here::here(), "inst", "extdata", "complexes", basename(paths$orthPath)))
