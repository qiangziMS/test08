require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
require(MRMlib)

env <- new.env()
hmdbInfo <- load("~/R/DB/DBraw/hmdbInfo.QT",env)
matchIdx <- load("~/R/DB/DBraw/initalMatch2HmdbIndex.QT", env)
initInfo <- load("~/R/DB/DBraw/initalInfo.QT", env)
# Get compounds i-~/R/---------------------------------------------

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


hmdbInfo_list <- as_tibble(hmdbInfo_list)
matchIdx_list <- as_tibble(matchIdx_list)
initInfo_list <- as_tibble(initInfo_list)


# SQlite DB ---------------------------------------------------------------
require(RSQLite)
require(dbplyr)

metaDB <- dbConnect(RSQLite::SQLite(),"~/R/DB/metaDB.sqlite")
dbWriteTable(metaDB, name = 'hmdbInfo', hmdbInfo_list)
dbWriteTable(metaDB, name = 'matchIdx', matchIdx_list)
dbWriteTable(metaDB, name = 'initInfo', initInfo_list)

# load from SQLite --------------------------------------------------------

hmdbInfo_list <- tbl(metaDB, "hmdbInfo")
matchIdx_list <- tbl(metaDB, "matchIdx")
initInfo_list <- tbl(metaDB, "initInfo")

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
        mz_tol = 0.015
    )



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

saveRDS(msp_list, file = "~/R/DB/msp_list.rds")
msp_list <- readRDS(file = "~/R/DB/msp_list.rds")


metaMsn_i <- new('metaMSn', MSn= msp_list)
traTbl <- filterMSn(object = metaMsn_i, topX = 10, type = "local", cluster = cl)

dbWriteTable(metaDB, "tarTbl", traTbl, overwrite = T)




traTbl_join <- traTbl %>% as_tibble() %>%
  left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
  left_join(., hmdbInfo_list, by = c("initialInChIKey" = "InChIKey")) %>%
  dplyr::filter(ExactMass-PrecursorMz < 2) %>%
  dplyr::filter(Polarity == "+")

# dbWriteTable(metaDB, "traTbl_join", traTbl_join)

# cross ref. --------------------------------------------------------------
xDB <- xMStoDB(ms = peak0, db = traTbl_join, tol = 5)
skylineTra <- toSkyline(infoTibble = xDB, deltaMz = 12) %>% unique()


# write to skyline transition -----------------------------------------------------------------
skylineTra_2 <-
  skylineTra %>% mutate(`Precursor Name` =
                          paste(`Precursor Name`,
                                round(`Explicit Retention Time` / 60, 1), sep="---")
                        )


write_csv(
  skylineTra_2,
    append = F,
    path = sprintf("./data/MRMtransition-%s.csv", Sys.Date())
)






