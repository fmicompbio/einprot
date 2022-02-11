## Generate ID conversion tables
library(einprot)

wormbase <- getConvTable(type = "WormBase")
saveRDS(wormbase, file = file.path(here::here(), "inst", "extdata",
                                   "conversion_tables",
                                   paste0("WormBaseConv_einprot",
                                          utils::packageVersion("einprot"),
                                          "_", gsub("-", "", Sys.Date()), ".rds")))

pombase <- getConvTable(type = "PomBase")
saveRDS(pombase, file = file.path(here::here(), "inst", "extdata",
                                  "conversion_tables",
                                  paste0("PomBaseConv_einprot",
                                         utils::packageVersion("einprot"),
                                         "_", gsub("-", "", Sys.Date()), ".rds")))
