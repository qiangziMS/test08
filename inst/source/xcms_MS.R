xcmsPeak0 <-
  function(Files = NULL,
           DIR = "~/RD/data-20180905/",
           Class,
           SQLiteName = NULL) {
    require(xcms)
    require(CAMERA)
    # files <- dir("E://raw/PRM/", full.names = T)
    if (is.null(Files)) {
      files <- dir(DIR, full.names = T, pattern = "mzML|mzXML")
    } else{
      files <- Files
    }
    
    if (Sys.info()["sysname"] == "Windows") {
      parN <- SnowParam(8)
    } else{
      parN <- MulticoreParam(8)
    }
    
    xSet <-
      xcmsSet(
        files = files,
        sclass = rep(Class, length(files)),
        method = "centWave",
        ppm = 10,
        noise = 5000,
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
    
    peakRaw <-
      group(
        xSet,
        method = "density",
        minfrac = 0,
        bw = 10,
        mzwid = 0.01,
        minsamp = 1
      )
    peak0 <-
      xsAnnotate(peakRaw) %>%
      groupFWHM() %>%
      groupCorr() %>%
      findIsotopes(minfrac = 0.1, ppm = 5) %>%
      getPeaklist() %>%
      filter(!isotopes %>% stri_detect(regex = "[M][+][:digit:]"))
    
    if (is.null(SQLiteName)) {
      metaDB <- dbConnect(RSQLite::SQLite(), "xcmsPeak0")
    } else{
      metaDB <- dbConnect(RSQLite::SQLite(), SQLiteName)
    }
    
    dbWriteTable(metaDB, "peak0", peak0, overwrite = T)
    return(peak0)
  }

