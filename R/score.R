require(xcms)
require(MassSpecWavelet)


# read mzML -----------------------------------------------------------------------------------

mzml <- dir("F://R&D/")
mrm <- function(sd) {
  a0 <- dnorm(1:10,5,sd)
  r0 <- function() rnorm(100,0,0.001)
  bi <- function() dnorm(1:20,sample(5:15,1),sample(seq(1,5,by = 0.2),1))
  ab <- function() c(r0(),bi(),a0*sample(1:5,1),bi(),r0())

  list(
    a1 = matrix(c(seq(ab(
    )), ab()), ncol = 2),
    a2 = matrix(c(seq(ab(
    )), ab()), ncol = 2),
    a3 = matrix(c(seq(ab(
    )), ab()), ncol = 2),
    a4 = matrix(c(seq(ab(
    )), ab()), ncol = 2),
    a5 = matrix(c(seq(ab(
    )), ab()), ncol = 2)
  )
}

MRMpeaks <- function(x = mrm()[[2]], scales, sn = 0) {
    x <- x[, 2]
    wCoefs <- cwt(ms = x, scales = scales)
    # wCoefs <- cbind(as.vector(x), wCoefs)
    # colnames(wCoefs) <- c(0, scales)
    localMax <- getLocalMaximumCWT(wCoefs = wCoefs)
    ridgeList <- getRidge(localMax)
    majorPeakInfo <-
      identifyMajorPeaks(
        ms = x,
        scales = scales,
        ridgeList = ridgeList,
        wCoefs = wCoefs,
        SNR.Th = sn,
        ridgeLength = 2,
        nearbyWinSize = 5,
        peakScaleRange = 5,
        minNoiseLevel = 0,
        winSize.noise = 0,
        nearbyPeak = T
      )
    if (!dev.interactive()) {
      # dev.new()
      # plot(rep(0, length(x)), type = "l", ylim = c(0, 2))
    }
    points(x,type="l",col=sample(1:10,1))
    betterPeakInfo <-
      tuneInPeakInfo(ms = x, majorPeakInfo = majorPeakInfo,maxScale = 128)
    betterPeakInfo[1:5] %>% as_tibble()
  }
# plot(rep(0,250),type="l",ylim=c(0,2))
# peakList <- lapply(mrm(3), MRMpeaks, scales = 1:12, sn = 0) %>% bind_rows()




