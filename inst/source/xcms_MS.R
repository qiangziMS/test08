require(xcms)
# require(parallel)
# files <- dir("E://raw/PRM/", full.names = T)
files <- dir("~/RD/mzML/",full.names = T)
if(Sys.info()["sysname"]=="Windows"){
  parN <- SnowParam(8)
}else{
  parN <- MulticoreParam(8)
}

xSet <-
  xcmsSet(
    files = files[1],
    method = "centWave",
    ppm = 10,
    noise = 1000,
    snthresh = 20,
    prefilter = c(6,1000),
    BPPARAM = parN,
    fitgauss = T,
    integrate = 1,
    peakwidth = c(10, 30)
  )
peak0 <- xSet@peaks %>% as_tibble()
