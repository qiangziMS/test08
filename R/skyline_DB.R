require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
env <- new.env()
hmdbInfo <- load("F://R&D/Rscript/hmdbInfo.QT", env)
matchIdx <- load("F://R&D/Rscript/initalMatch2HmdbIndex.QT", env)
initInfo <- load("F://R&D/Rscript/initalInfo.QT", env)
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


file_raw <- "E://R&D/Plant.msp"
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

# preFilter_idx[1:3]

msp_list <-
    parLapply(
        cl,
        range_idx[1:4000],
        mspParser,
        lib_vct = lib,
        nbTol = 0.8,
        mz_tol = 0.01
    )
metaMsn_i <- new('metaMSn', MSn= msp_list)

traTbl <- filterMSn(metaMsn_i, topX = 5, deltaMz = 12, cluster = 8)
# tra_dup <- traTbl_list[traTbl_list %>% duplicated() %>% which(),]
# cross ref. --------------------------------------------------------------
traTbl_list <- traTbl %>%
    left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
    left_join(., hmdbInfo_list, by = c("initialInChIKey" = "InChIKey")) %>%
    dplyr::filter(ExactMass-PrecursorMz < 2) %>%
    dplyr::filter(Polarity == "+")
    
skylineTra <- 
    toSkyline(infoTibble = traTbl_list, deltaMz = 12) %>% 
    unique() %>% 
    dplyr::filter(`Precursor Adduct` == "[M+H]+")



write_csv(
    skylineTra,
    append = F,
    path = sprintf("./MRM_transition_list_%s.csv", Sys.Date())
)

stopCluster(cl)








