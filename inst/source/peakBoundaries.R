load("~/RD/Jurkat cell QC SWATH/mtpSpecs.RData")
load("~/RD/Jurkat cell QC SWATH/diaFeature.RData")
load("~/RD/Jurkat cell QC SWATH/dData.RData")

ppc.cut <- cutoff <- 0.8

peak.i <- dia.feature@peakgroup$corrected[[1]]
eic.i.spl <- mtp.spec@spectrumMatch[[1]]$cell01[[1]]
eic.i.spl.ppc <- mtp.spec@featureScore[[1]]$cell01[[1]]

EICs <- list(
  ms1 = mtp.spec@featureMatchSmooth[[1]]$cell01,
  ms2 = eic.p.s[eic.p.s.ppc >= ppc.cut])

dia.res <-
  getDIAResult(
    new('DIAResult'),
    dia.feature,
    mtp.spec,
    cutoff ,
    minfrac.vote = 0.5
  )

diaResults1 <- dia.res@result


