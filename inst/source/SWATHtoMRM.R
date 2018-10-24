library(SWATHtoMRM)
SWATHtoMRM(d.in='~/RD/Jurkat cell QC SWATH/', d.out={},
           polarity ='positive',
           adduct.list ='HILIC',
           peakwidth = c(8, 30),
           sn = 20, minfrac = 1, minfrac.vote = 0.3,
           nSlaves = 8, ppm.pd = 20,
           cutoff = 0.8,
           ppm.ms2.mtp = 30,
           int.filter.ms1.field = 'maxo',
           int.filter.ms1 = 500,
           int.filter.ms2 = 100,
           isFWHM=FALSE,
           is.plot.eic.feature = FALSE,
           is.plot.eic.spec = FALSE,
           is.ms1.only = FALSE,
           rerun = FALSE)

require(tidyverse)
require(stringi)
tra <- read_csv('~/RD/Jurkat cell QC SWATH/Results/result.csv')


transitions <- tra %>%
  transmute(
  `Precursor Name` = sprintf("%.1f-%.1f", mz.precursor, rt),
  `Precursor Charge` = 1,
  `Precursor m/z` = mz.precursor,
  `Product m/z` = mz.product,
  # `Product intensity` = int.product,
  `Product Charge` = 1,
  `Explicit Retention Time` = round(rt/60,2)
)

write_csv(transitions, "~/RD/Jurkat cell QC SWATH/Results/skyline_tra.csv")


toSearch <- function(X = tra, File = "./SWATHtoMRM-spectra.msp", overWrite = T) {
  tra_list <- X %>% unique() %>% split.data.frame(f = .$ft.idx)
  if(overWrite & file.exists(File)){
    file.remove(File)
  }
  for (i in tra_list) {
    paste(
      sprintf("BEGIN\n"),
      sprintf("id: biotree_%06d\n", i$ft.idx[1]),
      sprintf("precursorMz: %s\n", i$mz.precursor[1]),
      sprintf("retentionTime: %s\n", i$rt[1]),
      sprintf("peaks_count: %s\n", nrow(i)),
      sprintf("%.4f\t%.0f\n", i$mz.product, i$int.product) %>% paste(collapse = ""),
      sprintf("END\n"),
      "\n",
      sep = "", collaps = ""
    ) %>% write_file(File, append = T)
  }
}

toSearch(X = tra, File = "~/RD/Jurkat cell QC SWATH/Results/SWATHtoMRM-spectra.msp" )





