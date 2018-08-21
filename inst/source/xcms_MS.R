require(xcms)
require(CAMERA)
# files <- dir("E://raw/PRM/", full.names = T)
files <- dir("~/RD/mzML/20180827-3//",full.names = T)
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
    snthresh = 20,
    prefilter = c(6, 5000),
    BPPARAM = parN,
    fitgauss = T,
    integrate = 1,
    peakwidth = c(10, 30)
  )

peakRaw <- group(xSet, method = "density", minfrac=0,bw = 10, mzwid = 0.1, minsamp=0)
peakRaws <- peakRaw@peaks %>% as_tibble()
peakRaw@groups %>% as_tibble() %>% filter(mzmed > 611&mzmed <612 & rtmed >345&rtmed <352)
peakRaw@peaks %>% as_tibble() %>% filter(mz > 611&mz <612 & rt >345&rt <352)


#Create an xsAnnotate object
xsa <- xsAnnotate(group(xSet, method = "density", minfrac=0))
#Group after RT value of the xcms grouped peak
xsaF <- groupFWHM(xsa, perfwhm=5)
# peak <- xsaF %>% getPeaklist()
#Verify grouping
xsaC <- groupCorr(xsaF)
#Annotate isotopes, could be done before groupCorr
xsaFI <- findIsotopes(xsaC, minfrac=0.1, ppm = 5)
peak <- xsaFI %>% getPeaklist()
peak0 <- peak %>% filter(!isotopes %>% stri_detect(regex = "[M][+][:digit:]"))

#Annotate adducts
rules <-
  system.file("rules/primary_adducts_pos.csv", package = "CAMERA") %>%
  read.csv()
xsaFA <- findAdducts(xsaFI, polarity="positive", rules = rules)

peak <- xsaFA %>% getPeaklist()
peak0 <- peak %>% filter(!isotopes %>% stri_detect(regex = "[M][+][:digit:]"))

xDB$ppm %>% hist(breaks=50)
