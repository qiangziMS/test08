require(xcms)
require(CAMERA)
# files <- dir("E://raw/PRM/", full.names = T)
files <- dir("~//RD/mzML/20180827-3//",full.names = T)
if(Sys.info()["sysname"]=="Windows"){
  parN <- SnowParam(8)
}else{
  parN <- MulticoreParam(8)
}

xSet <-
  xcmsSet(
    files = files[1:3],
    sclass = rep("a1",3),
    method = "centWave",
    ppm = 10,
    noise = 8000,
    snthresh = 100,
    prefilter = c(5, 5000),
    BPPARAM = parN,
    fitgauss = T,
    integrate = 1,
    peakwidth = c(10, 30)
  )

xSet <-
  retcor(
    xSet,
    method = "obiwarp",
    plottype = c("none", "deviation"),
    profStep = 0.01,
    # center = NULL,
    col = NULL,
    ty = NULL,
    response = 1,
    distFunc = "cor_opt",
    gapInit = NULL,
    gapExtend = NULL,
    factorDiag = 2,
    factorGap = 2,
    localAlignment = 0,
    initPenalty = 0
  )

peakRaw <- group(xSet, method = "density", minfrac=0, bw = 10, mzwid = 0.01, minsamp=1)
peak0 <-
  xsAnnotate(peakRaw) %>%
  groupFWHM() %>%
  groupCorr() %>%
  findIsotopes(minfrac=0.1, ppm = 5) %>%
  getPeaklist() %>%
  filter(!isotopes %>% stri_detect(regex = "[M][+][:digit:]"))

dbWriteTable(metaDB, "peak0", peak0, overwrite = T)


