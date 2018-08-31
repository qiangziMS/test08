require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
require(MRMlib)



# SQlite DB ---------------------------------------------------------------
require(RSQLite)
require(dbplyr)

metaDB <- dbConnect(RSQLite::SQLite(),"~/R/DB/metaDB.sqlite")

# load from SQLite --------------------------------------------------------

hmdbInfo_list <- tbl(metaDB, "hmdbInfo") %>% as_tibble()
matchIdx_list <- tbl(metaDB, "matchIdx") %>% as_tibble()
initInfo_list <- tbl(metaDB, "initInfo") %>% as_tibble()


# saveRDS(msp_list, file = "~/R/DB/msp_list_n1_m0.1.rds")
msp_list <- readRDS(file = "~/R/DB/msp_list.rds")

cl <- makeCluster(8)
metaMsn_i <- new('metaMSn', MSn= msp_list)
traTbl <- filterMSn(object = metaMsn_i, topX = 5, type = "local", cluster = cl)

dbWriteTable(metaDB, "tarTbl-n1-m0.1-t5", traTbl, overwrite = T)

# traTbl <- tbl(metaDB, "tarTbl") %>% as_tibble()



traTbl_join <- traTbl %>% as_tibble() %>%
  left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
  left_join(., hmdbInfo_list, by = c("initialInChIKey" = "InChIKey")) %>%
  dplyr::filter(ExactMass-PrecursorMz < 2) %>%
  dplyr::filter(PrecursorMz > 50) %>%
  dplyr::filter(Polarity == "+")


dbWriteTable(metaDB, "traTbl_join-n1-m0.1-t5", traTbl_join, overwrite = T)

# ---------- LOAD DB  -----------------------------------------------------

traTbl_join <- tbl(metaDB, "traTbl_join") %>% as_tibble()


# cross ref. --------------------------------------------------------------



xDB <- xMStoDB(ms = peak0, db = traTbl_join, tol = 20)
skylineTra <- toSkyline(infoTibble = xDB, deltaMz = 12) %>% unique()


# write to skyline transition -----------------------------------------------------------------
skylineTra_2 <-
  skylineTra %>%
  mutate(`Precursor Name` =paste(`Precursor Name`, `Explicit Retention Time`, sep="--"))


write_csv(
  skylineTra_2,
    append = F,
    path = sprintf("./data/MRMtransition-3x-%s.csv", Sys.Date())
)

stopCluster(cl)




