require(tidyverse)
require(RSQLite)
require(dbplyr)



cpd <- readxl::read_excel(path = "~/RD/microbiome/microbiomeList.xlsx") %>%
  mutate(KEGG = stri_replace(cpd, replacement = "", regex = "cpd:"))

env <- new.env()
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
merge_info <- left_join(initInfo_list, hmdbInfo_list, by = c("match2HmdbIndex"="IDX"))
names(merge_info)

bestKEGG <- function(a,b){
  kegg_a <- stri_extract(str = a, regex = "[:alnum:]+")
  kegg_b <- stri_extract(str = b, regex = "[:alnum:]+")
  apply(X = tibble(kegg_a,kegg_b), MARGIN = 1, FUN = function(x) na.omit(x)[1])
}

merge_info <- merge_info %>% mutate(KEGGid = bestKEGG(a=KEGG,b=initialKEGG))

microCPD <-
  merge_info %>%
  filter(KEGGid %in% na.omit(cpd$KEGG)) %>%
  select(KEGGid, splash10, Names, SuperClass, Class, SubClass) %>%
  unique()

microCPD <-
  microCPD %>% left_join(x = cpd,
                         y = .,
                         by = c("KEGG" = "KEGGid"))











