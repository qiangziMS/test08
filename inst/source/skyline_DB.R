require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
env <- new.env()
hmdbInfo <- load(system.file("DB/hmdbInfo.QT", package = "MRMlib"), env)
matchIdx <- load(system.file("DB/initalMatch2HmdbIndex.QT", package = "MRMlib"), env)
initInfo <- load(system.file("DB/initalInfo.QT", package = "MRMlib"), env)
# Get compounds idx -------------------------------------------------------

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


file_raw <- system.file("DB/Plant.msp", package = "MRMlib")
lib <- read_lines(file = file_raw)
# lib_val <- lib[]
# Name_idx <- grep(pattern = "BEGIN",ignore.case = F, fixed = T,x = lib)
begin_idx <- grep(pattern = "BEGIN",ignore.case = F, fixed = T,x = lib)
end_idx <- grep(pattern = "END",ignore.case = F, fixed = T,x = lib)

range_idx <-
    tibble(start = begin_idx+1,
           end = end_idx-1,
           idx = seq_along(end_idx)) %>%
    split.data.frame(f=.$idx)



# Parse .MSP database file to list ----------------------------------------
require(parallel)

cl <- makeCluster(8)
preFilter_idx <- parLapply(cl, range_idx, mspPreFilter, lib_vct=lib) %>% unlist()

preFilterTbl <- tibble(idx = seq_along(preFilter_idx),  splash = preFilter_idx) %>%
  left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
  left_join(., hmdbInfo_list, by = c("initialInChIKey" = "InChIKey")) %>%
  dplyr::filter(Polarity == "+") %>%
  split.data.frame(f=.$Class)


# preFilter_idx[1:3]

msp_list <-
    parLapply(
        cl,
        range_idx,
        mspParser,
        lib_vct = lib,
        nbTol = 0.8,
        mz_tol = 0.01
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

# git config --global user.email "qiangziMS@gmail.com"
# git config --global user.name "qiangziMS"


metaMsn_i <- new('metaMSn', MSn= msp_list)

traTbl <- filterMSn(metaMsn_i, topX = 5, type = "local", deltaMz = 12, newRule = NULL, cluster = cl)
# traTbl_s
traTbl_s <- traTbl %>% dplyr::filter(stri_detect(Adduct, regex = "[NA]NA"))

idx <- traTbl_s$splash %>% match(., preFilter_idx)

msp_list_NA <- msp_list[idx]

adduct <- table(traTbl$Adduct0)
# tra_dup <- traTbl_list[traTbl_list %>% duplicated() %>% which(),]
# cross ref. --------------------------------------------------------------
traTbl_list <- traTbl %>%
    left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
    left_join(., hmdbInfo_list, by = c("initialInChIKey" = "InChIKey")) %>%
    dplyr::filter(ExactMass-PrecursorMz < 2) %>%
    dplyr::filter(Polarity == "+")

xDB <- xMStoDB(ms = as_tibble(peak0),db = traTbl_list,tol = 20)


xTra <- inner_join(x = xDB, y = db0)

skylineTra <-
    toSkyline(infoTibble = xTra, deltaMz = 12) %>%
    unique() %>%
    dplyr::filter(`Precursor Adduct` == "[M+H]+")



write_csv(
    skylineTra,
    append = F,
    path = sprintf("./MRM_transition_list_%s.csv", Sys.Date())
)

stopCluster(cl)








