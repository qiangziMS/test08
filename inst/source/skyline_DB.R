
# batch MRMlib ------------------------------------------------------------

<<<<<<< HEAD
mlib <- MRMlib(DIR = "./data/")
=======
mlib <- MRMlib()
>>>>>>> d44a169e6f05f3a738e62b343a555bfe41510416


# seperate ----------------------------------------------------------------



require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
require(MRMlib)
require(RSQLite)
require(dbplyr)
metaDB <- dbConnect(RSQLite::SQLite(),"~/R/DB/metaDB.sqlite")
env <- new.env()


# load raw database ---------------------------------------------------------------------------

hmdbInfo <- load("~/R/DB/DBraw/hmdbInfo.QT",env)
matchIdx <- load("~/R/DB/DBraw/initalMatch2HmdbIndex.QT", env)
initInfo <- load("~/R/DB/DBraw/initalInfo.QT", env)

hmdbInfo_list <- lapply(hmdbInfo, function(x) {
    as.name(x) %>% eval(envir = env)
    })
names(hmdbInfo_list) <- hmdbInfo

matchIdx_list <- lapply(matchIdx, function(x) {
    as.name(x) %>% eval(envir = env)
    })
names(matchIdx_list) <- matchIdx

initInfo_list <- lapply(initInfo, function(x) {
    as.name(x) %>% eval(envir = env)
    })
names(initInfo_list) <- initInfo

initInfo_list <- as_tibble(initInfo_list) %>% bind_cols(matchIdx_list,.)
hmdbInfo_list <- as_tibble(hmdbInfo_list) %>% mutate(IDX = seq(nrow(.)))

# SQlite DB ---------------------------------------------------------------

dbWriteTable(metaDB, name = 'hmdbInfo', hmdbInfo_list, overwrite = T)
dbWriteTable(metaDB, name = 'initInfo', initInfo_list, overwrite = T)

# load from SQLite --------------------------------------------------------

hmdbInfo_list <- tbl(metaDB, "hmdbInfo") %>% as_tibble()
initInfo_list <- tbl(metaDB, "initInfo") %>% as_tibble()

# load msp file -------------------------------------------------------------------------------

# file_raw <- system.file("/DB/Plant.msp", package = "MRMlib")
file_raw <- "~/R/DB/Plant.msp"
lib <- read_lines(file = file_raw)
begin_idx <- grep(pattern = "BEGIN",ignore.case = F, fixed = T,x = lib)
end_idx <- grep(pattern = "END",ignore.case = F, fixed = T,x = lib)

range_idx <-
    tibble(start = begin_idx+1,
           end = end_idx-1,
           idx = seq_along(end_idx)) %>%
    split.data.frame(f=.$idx)

# Parse .MSP database file to list ----------------------------------------

cl <- makeCluster(8)
msp_list <-
    parLapply(
        cl,
        range_idx,
        mspParser,
        lib_vct = lib,
        nbTol = 0.8,
        mz_tol = 0.01,
        mz_only = T,
        ignoreInt = F
    )


# to MSP file -------------------------------------------------------------


# check NA ----------------------------------------------------------------
msp_list <- parLapply(
  cl,
  msp_list,
  fun = function(x) {
    if (x$precursorType == ""|is.null(x$precursorType)) {
      if (x$polarity == ""|is.null(x$polarity)) {
        x$polarity <- "+"
        cat(x$formula)
      }
      x$precursorType <- sprintf("[M%sH]%s", x$polarity, x$polarity)
    }
    return(x)
  }
)

# save msp_list to file ---------------------------------------------------

saveRDS(msp_list, file = "~/R/DB/msp_list_n0.8_m0.01ig.rds")
msp_list <- readRDS(file = "~/R/DB/msp_list_n0.8_m0.01ig.rds")



# run filterMsn -----------------------------------------------------------


cl <- makeCluster(8)
metaMsn_i <- new('metaMSn', MSn= msp_list)
traTbl <- filterMSn(object = metaMsn_i, topX = 5, type = "local", a = 0.5,b=0.5,cluster = cl)



dbWriteTable(metaDB, "tarTbl-n0.8-m0.01ig-t5", traTbl, overwrite = T)
traTbl <- tbl(metaDB, "tarTbl-n0.8-m0.01ig-t5") %>% as_tibble()



traTbl_join <- traTbl %>% as_tibble() %>%
  left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
  left_join(., hmdbInfo_list, by = c("match2HmdbIndex" = "IDX")) %>%
  dplyr::filter(ExactMass-PrecursorMz < 2) %>%
  dplyr::filter(!stri_detect(Formula, regex = "[D]")) %>%
  dplyr::filter(Polarity == "+")


dbWriteTable(metaDB, "traTbl_join-n0.8-m0.01ig-t5", traTbl_join, overwrite = T)



# to skyline DB.MSP -------------------------------------------------------

tra_list <- traTbl_join %>% toSkyline(deltaMz = 12)

toSkylineDB(X = tra_list, overWrite = T)





# +++++++ START from DB +++++++ -------------------------------------------
# ---------- LOAD DB  -----------------------------------------------------

traTbl_join <- tbl(metaDB, "traTbl_join-n0.8-m0.01ig-t5") %>% as_tibble()
peak0 <- tbl(metaDB, "peak0") %>% as_tibble()

# cross ref. --------------------------------------------------------------

xDB <- xMStoDB(MS = peak0, DB = traTbl_join, tol = 5, N = 8000)
skylineTra <- toSkyline(infoTibble = xDB, deltaMz = 12) %>% unique()

# write to skyline transition -----------------------------------------------------------------
skylineTra_2 <-
  skylineTra %>%
  split.data.frame(f=.$InChiKey) %>%
  lapply(FUN = mergeOverlap, RTw = 2) %>%
  bind_rows() %>%
  mutate(`Precursor Name` = paste(`Precursor Name`, `Explicit Retention Time`, sep="--"))


write_csv(
  skylineTra_2 %>% select(-InChiKey),
    append = F,
    path = sprintf("./data/MRMtransition-%s.csv", Sys.Date())
)

# skylineTra_2$InChiKey %>% unique() %>% length()


stopCluster(cl)




