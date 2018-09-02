require(xcms)
require(CAMERA)
# files <- dir("E://raw/PRM/", full.names = T)
files <- dir("~/R/RD/mzML/20180827-3//",full.names = T)
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
    noise = 5000,
<<<<<<< HEAD
    snthresh = 20,
    # prefilter = c(6, 5000),
=======
    snthresh = 100,
    prefilter = c(3, 5000),
>>>>>>> aad68aad45b2fc363eb2ac85c342af92fff87292
    BPPARAM = parN,
    fitgauss = T,
    integrate = 1,
    peakwidth = c(10, 30)
  )

xSet <-
  retcor(
    xSet,
    method = "obiwarp",
<<<<<<< HEAD
    # plottype = c("none", "deviation"),
=======
    plottype = c("none", "deviation"),
>>>>>>> aad68aad45b2fc363eb2ac85c342af92fff87292
    profStep = 0.01,
    # center = NULL,
    col = NULL,
    ty = NULL,
    response = 1,
    distFunc = "cor_opt",
    gapInit = NULL,
    gapExtend = NULL,
    factorDiag = 2,
<<<<<<< HEAD
    factorGap = 1,
    localAlignment = 0,
    initPenalty = 0
  )

peakRaw <- group(xSet, method = "density", minfrac=0, max=100,bw = 30, mzwid = 0.1, minsamp=0)
peakss <- peakRaw@peaks %>% as_tibble()
xx <- peakRaw@groups %>% as_tibble() %>% filter(mzmed > 611&mzmed <612 )
xxx <- groupval(peakRaw)
pr <- peakRaw@peaks %>% as_tibble() %>% filter(mz > 610&mz <613 & rt >300&rt <400)
gidx <- peakRaw@groupidx
xxxx <- peakss[gidx[17634][[1]], ]
xx[which(xx$mzmed >= 611.158&xx$mzmed <= 611.159),]

=======
    factorGap = 2,
    localAlignment = 0,
    initPenalty = 0
  )
>>>>>>> aad68aad45b2fc363eb2ac85c342af92fff87292

peakRaw <- group(xSet, method = "density", minfrac=0, bw = 15, mzwid = 0.01, minsamp=1)
peak0 <-
  xsAnnotate(peakRaw) %>%
  groupFWHM() %>%
  groupCorr() %>%
  findIsotopes(minfrac=0.1, ppm = 5) %>%
  getPeaklist() %>%
  filter(!isotopes %>% stri_detect(regex = "[M][+][:digit:]"))

dbWriteTable(metaDB, "peak0", peak0, overwrite = T)


