require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
require(parallel)
# require(MRMlib)
# Get isotope proportion list ---------------------------------------------


#' Get isotope proportion range up/down
#'
#' @param iso1 isotope up
#' @param iso2 isotope down
#' @param cpdLib the database; see compoundLibraries() for a list of supported
#' databases
#' @param mzWin the mass window size for grouping compounds;
#' see massWindowSizes(compoundLibrary = "kegg") for a list of supported
#' @param mz m/z vector to be calculated
#' @param quantile_2 quantile range
#' databases for e.g. the database kegg
#' @return A list of 2-column tibbles with length equal to the multiple iso arg
#' @export NULL
#'
#' @examples NULL
get_iso_percent <-
    function(mz = mz_x,
             iso1 = 0,
             iso2 = 1,
             quantile_2 = c(0.05, 0.95),
             cpdLib = "kegg",
             mzWin = 50) {
        cQtl <-
            CAMERA::compoundQuantiles(compoundLibrary = cpdLib,
                                      massWindowSize = mzWin)
        mz[mz >= 1000] <- 999
        if (length(iso1) > 1 & length(iso2) > 1) {
            stop("only 1 iso can larger than 1 element")
        } else if (length(iso1) > 1) {
            names(iso1) <- paste(iso1, iso2, sep = "/")
            FUN <- function(iso) {
                tibble(
                    up = CAMERA::getIsotopeProportion(
                        object = cQtl,
                        isotope1 = iso,
                        isotope2 = iso2,
                        mass = mz,
                        quantile = quantile_2[1]
                    ),
                    dn = CAMERA::getIsotopeProportion(
                        object = cQtl,
                        isotope1 = iso,
                        isotope2 = iso2,
                        mass = mz,
                        quantile = quantile_2[2]
                    )
                )
            }
            lapply(iso1, FUN = FUN)
        } else if (length(iso2) >= 1) {
            names(iso2) <- paste(iso1, iso2, sep = "/")
            FUN <- function(iso) {
                tibble(
                    up = CAMERA::getIsotopeProportion(
                        object = cQtl,
                        isotope1 = iso1,
                        isotope2 = iso,
                        mass = mz,
                        quantile = quantile_2[1]
                    ),
                    dn = CAMERA::getIsotopeProportion(
                        object = cQtl,
                        isotope1 = iso1,
                        isotope2 = iso,
                        mass = mz,
                        quantile = quantile_2[2]
                    )
                )
            }
            lapply(iso2, FUN = FUN)
        }

    }



# Find_iso ----------------------------------------------------------------


#' Find isotope proportion list
#'
#' @param i The i-th product ion, i <= N
#' @param mz_x A MSn-matrix with 2 columns and N rows
#' @param mz_tol mz tolerance
#' @param predictedIso A list of j predicted isopote proportion tibble,
#'        j is the j-th isotope, length(.) == isotope number
#'
#' @return A vector of length j, indicate whether the i-th product ion
#'         is a j-th isotope
#' @export NULL
#'
#' @examples NULL

# mz_x <- spectra
# mz_tol <- 0.01


find_isotopes <- function(i, mz_x, mz_tol=0.1, binMz=0.5, mz_only = T, ignoreInt=T, predictedIso) {
        .isotopeN <- NULL
        .mz_tol <- abs(mz_tol)

        for (j in seq(predictedIso)) {
          if (.mz_tol > 1) {

            is_iso_mz <- (abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) /
                                    mz_x[seq(i), 1]) <= (.mz_tol / 1000000)
          } else if (.mz_tol <= 1) {

            is_iso_mz <- abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) <= .mz_tol

          }

          if(ignoreInt){
            is_iso_all <- is_iso_mz
            .isotopeN[j] <- any(is_iso_all)
          }else if(mz_only){
                is_iso_all <- is_iso_mz & (mz_x[seq(i), 2] > mz_x[i, 2])
                .isotopeN[j] <- any(is_iso_all)
            }else{
            is_iso_all <-
                is_iso_mz &
                mz_x[seq(i), 2] / mz_x[i, 2] >= predictedIso[[j]][[1]][seq(i)] &
                mz_x[seq(i), 2] / mz_x[i, 2] <= predictedIso[[j]][[2]][seq(i)]
            .isotopeN[j] <- any(is_iso_all)
            }
        }
        return(.isotopeN)
    }

# parse .MSP file to list -------------------------------------------------

#' Parse .MSP file to list
#'
#' @param x index of compound info. (list/vector/tibble)
#' @param lib_vct .MSP file lines in (vector), without blank lines
#' @param quantile_2 isotope proportion quantile
#' @param cpdLib library for calculate isotope proportion
#' @param nbTol merge mz tolerance
#' @param mz_tol mz tolerance of isotope
#' @param mz_only ignore isotope proportion, use mz only
#' @param ... not used
#'
#' @return list with compound info. vector and mz/intensity tibble
#' @export NULL
#'
#' @examples NULL
mspParser <- function(x = range_idx["5"],
                      quantile_2 = c(0.05,0.95),
                      cpdLib = "kegg",
                      iso1 = 0,
                      iso2 = 1:4,
                      lib_vct = lib,
                      nbTol = 0.8,
                      mz_tol = 0.01,
                      mz_only= T,
                      ignoreInt = T,...) {
  require(dplyr)
  require(tidyverse)
  require(stringi)
  require(MRMlib)

  x <- unlist(x)
    lib_val_x <- lib_vct[x[1]:x[2]]

    info_x <-
        lib_val_x %>%
        grep(
            fixed = F,
            pattern = ".+[=|:][ ]?",
            value = T
        ) %>%
        stri_replace_first(
            replacement = "",
            regex = ".+[=|:][ ]?") %>%
        as.list()

    names(info_x) <-
        lib_val_x %>%
        grep(
            fixed = F,
            pattern = ".+[=|:][ ]?",
            value = T) %>%
        lapply(., stri_extract, regex = "^[:print:]+(?==|:)")

    # Get MSn spectrum matrix
    mz_x <-
      lib_val_x %>%
      grep(fixed = F,
           pattern = "^[0-9]+",
           value = T) %>%
      stri_split(regex = "[:space:]+") %>%
      lapply(FUN = as.numeric) %>%
      do.call(rbind, .) %>%
      .[.[, 2] > 0,] %>%
      matrix(ncol = 2) %>%
      cleanMSn.iter(mz = ., tol = nbTol, echo = F)
# cat(x)

    # Get predict isotope proportion list
    if (mz_only) {
        pred_iso <- seq_along(iso1)*seq_along(iso2)
    } else{
        pred_iso <-
            get_iso_percent(
                mz = mz_x[, 1],
                cpdLib = cpdLib,
                quantile_2 = quantile_2,
                iso1 = iso1,
                iso2 = iso2
            )
    }

    # Find isotopes
    MSn_iso <-
        lapply(
            seq_along(mz_x[, 1]),
            FUN = find_isotopes,
            mz_x = mz_x,
            mz_tol = mz_tol,
            mz_only = mz_only,
            ignoreInt = ignoreInt,
            predictedIso = pred_iso
        ) %>%
        do.call(rbind, .) %>%
        as_tibble()

    names(MSn_iso) <- paste("iso", seq_along(MSn_iso), sep = "")

    mz_tbl <- as_tibble(mz_x)
    names(mz_tbl) <- c("mz", "intensity")

    info_x[['MSn']] <- bind_cols(mz_tbl, MSn_iso)
    return(info_x)
}

# Get classyFire sdf file ------------------------------------------------------

setGeneric("classyFireAPI",
           function(object, type, toDir = sprintf("./classyFire--%s.%s", Sys.Date(), type), ...) standardGeneric("classyFireAPI") )
setMethod("classyFireAPI", c(object = "list"),
          function(object, type = "sdf",
                   toDir = sprintf("./classyFire--%s.%s", Sys.Date(), type),
                   mode = "a", ...){
              lapply(object, FUN = function(msp_i){
                  url_i <- sprintf("http://classyfire.wishartlab.com/entities/%s.%s",
                                   msp_i$InChIKey, type)
                  if(RCurl::url.exists(url_i)){
                      curl::curl_download(url = url_i,
                                          destfile = toDir,
                                          mode = mode)
                  }
              })
})
setMethod("classyFireAPI", c(object = "character"),
          function(object, type = "sdf",
                   toDir = sprintf("./classyFire--%s.%s", Sys.Date(), type),
                   mode = "a", ...){
              lapply(object, FUN = function(msp_i){
                  url_i <- sprintf("http://classyfire.wishartlab.com/entities/%s.%s",
                                   msp_i, type)
                  if(RCurl::url.exists(url_i)){
                      curl::curl_download(url = url_i,
                                          destfile = toDir,
                                          mode = mode)
                  }
              })
          })


# New Generic function to fix SDF/SDFset objects --------------------------

setGeneric("fixSDFset", function(object, ...) standardGeneric("fixSDFset") )

# New method for SDF/SDFset class ---------------------------------------------
setMethod("fixSDFset", c(object = "SDF"), function(object) {
    datablock <- object@datablock[-1]
    idx <- datablock %>% grep(pattern = "> <")


    inchikey_idx <-
        datablock[idx+1] %>%
        grep(pattern = "inchikey", ignore.case = T)

    datablock[idx+1][inchikey_idx] <-
        datablock[idx+1][inchikey_idx] %>%
        stri_extract(regex = "(?<=[:alpha:]{1,8}[=]).+$")
    object@datablock <- datablock[idx + 1]
    names(object@datablock) <-
        stri_extract(str = datablock, regex = "(?<=^> <).+(?=>$)")
    object
})
setMethod("fixSDFset", c(object = "SDFset"), function(object) {
    object@SDF <-  lapply(object@SDF, FUN = fixSDFset)
    object@ID <- lapply(
        object@SDF,
        FUN = function(sdf) {
            sdf@datablock["InChIKey"]
        }
    ) %>% unlist()
    object
})



# metaMSn class definition ------------------------------------------------

setClass('metaMSn', slots = c(SDFinfo="SDFset", MSn="list"))

setGeneric("metaLink", function(object, tag, ...) standardGeneric("metaLink") )
setMethod('metaLink', c(object="metaMSn"), function(object, tag, ...){
    inchikey_idx <- datablocktag(object@SDFinfo, tag = "InChIKey")
    tagMatch_idx <- datablocktag(object@SDFinfo, tag = tag)
    names(tagMatch_idx) <- inchikey_idx
    return(tagMatch_idx)
})



# Filter product ions to MRM ------------------------------------------------------------------

#' Filter product ions to MRM
#'
#' @param object metaMSn object
#' @param topX select top X product ions
#' @param ... NOT used
#' @param type local/...
#' @param cluster
#'
#' @return Tibble of products and info (infoTibble)
#' @export NULL
#'
#' @examples NULL
setGeneric("filterMSn",function(object,topX,type,cluster,a,b,...) standardGeneric("filterMSn"))

setMethod("filterMSn", c(object = "metaMSn"),
          function(object,topX=5,type="local",cluster = cl, a=0.5, b=0.5,...) {
            if (type == "local") {
              .fun <- function(objList) {
                require(MRMlib)
                require(tidyverse)
                require(stringi)

                dplyr::filter(objList$MSn,!iso1 & !iso2 & !iso3) %>%
                  dplyr::select(mz, intensity) %>%
                  dplyr::mutate(
                    splash = objList$id,
                    iRank = rank(-intensity, ties.method = "first"),
                    Adduct0 = ifelse(is.null(objList$precursorType),
                                     NA,objList$precursorType),
                    Polarity = ifelse(is.null(objList$polarity),
                                      NA,objList$polarity),
                    ExactMass = ifelse(is.null(objList$exact_mass),
                                       NA,objList$exact_mass) %>% as.numeric(),
                    Formula = ifelse(is.null(objList$formula),
                                     NA,objList$formula),
                    CE = ifelse(is.null(objList$collisionEnergy),
                                NA,objList$collisionEnergy),
                    InstrumentType = ifelse(is.null(objList$instrument_type),
                                            NA,objList$instrument_type)
                  ) %>% as_tibble()
              }

              if (is.null(cluster)) {
                infoTibble <-
                  lapply(object@MSn, .fun) %>%
                  do.call(bind_rows, args = .) %>%
                  dplyr::filter(iRank <= 20)
              } else{
                cat("step 1 start ...\n")
                infoTibble <-
                  parallel::parLapply(cluster, object@MSn, .fun) %>%
                  do.call(bind_rows, args = .) %>%
                  dplyr::filter(iRank <= 20)
                cat("step 1 done ...\n")
              }
            }
            cat("step 2 start ...\n")

# add fragments specificity filter ------------------------------------------------------------

            tblRaw <-
              infoTibble %>%
              mutate(fMZ = round(mz, 0)) %>%
              group_by(Formula, fMZ)

            tblRank <-
              tblRaw %>%
              summarise(dupNum = n()) %>%
              left_join(tblRaw, .) %>%
              group_by(splash) %>%
              mutate(
                sRank = rank(dupNum, ties.method = "first"),
                Ranks = rank(a*iRank + b*sRank, ties.method = "first"))

            xAdduct(
              cls = cluster,
              formulas = tblRank$Formula,
              adducts = tblRank$Adduct0,
              polarities = tblRank$Polarity) %>%
              bind_cols(tblRank, .) %>%
              mutate(DeltaMz = PrecursorMz - mz, RTs = 0) %>%
              dplyr::filter(Ranks <=topX)
          })



# to skyline traTable -------------------------------------------------------------------------
#' Write to skyline traTable
#'
#' @param infoTibble DB infoTable from xMStoDB/filterMSn
#' @param deltaMz delta mz with precursor mz
#'
#' @return transition list table of skyline
#' @export NULL
#'
#' @examples NULL
setGeneric("toSkyline", function(infoTibble, deltaMz) standardGeneric("toSkyline"))
setMethod("toSkyline", c(infoTibble = "tbl"),
          function (infoTibble, deltaMz=12) {
              chargeRule <-
                  c(
                      `-` = -1,
                      `+` = 1,
                      `-2` = -2,
                      `+2` = 2,
                      `2-` = -2,
                      `2+` = 2,
                      `-3` = -3,
                      `+3` = 3,
                      `3-` = -3,
                      `3+` = 3
                  )

              bestName <- function(...){

                Args <- list(...)
                names(Args) <- paste("V", seq_along(Args),sep="")
                  as_tibble(Args) %>%
                  transmute(maxx = apply(., 1, function(x){
                    # cat(x)
                    xn <- nchar(x, type = "bytes")
                    xn[xn<1] <- NA
                    x[which.min(xn)]
                  })
                  ) %>% .$maxx
              }


              multiCETbl <-
                  infoTibble %>%
                  dplyr::filter(DeltaMz >= deltaMz) %>%
                  mutate(CE =ifelse(is.na(CE)|CE=="", "30.1", CE)) %>%
                  transmute(
                      `Note` = splash,
                      `Molecule List Name` = Class,
                      `Precursor Name` = bestName(Name, initialName),
                      `Precursor Formula` = Formula,
                      `Precursor Adduct` = Adduct1,
                      `Precursor Charge` = chargeRule[as.character(pCharge)],
                      `Precursor m/z` = PrecursorMz/`Precursor Charge`,
                      `Product m/z` = mz,
                      `Product intensity` = intensity,
                      `Product Charge` = chargeRule[as.character(Polarity)],
                      `Explicit Retention Time` = round(RTs/60,2),
                      # `Product Name` =  round(PrecursorMz,2),
                      `InChiKey` = InChIKey,
                      `Explicit Collision Energy` =
                          stri_extract_all(CE, regex = "[:digit:]+[.]?[:digit:]+") %>%
                          .[[1]] %>%
                          as.numeric() %>%
                          median()
                  ) %>%
                  dplyr::filter(stri_detect_regex(
                      InChiKey, pattern = "[:upper:]+[-][:upper:]+[-][:upper:]")
                      )

              grp.rm <- grep("Precursor Name",names(multiCETbl),fixed = T)
              multiCETbl_u <-
                multiCETbl  %>%
                group_by(!!!lapply(
                  names(multiCETbl)[- grp.rm], as.symbol)) %>%
                summarise(`Precursor Name` = first(`Precursor Name`))

              multiCETbl_u %>%
                group_by(`Note`, `InChiKey`) %>%
                  summarise(CEBest = n()) %>%
                  split.data.frame(f = .$InChiKey) %>%
                  lapply(FUN = function(tbl) {
                          tbl[which.max(tbl$CEBest), ]
                      }) %>%
                  do.call(bind_rows, args = .) %>%
                  left_join(
                      x = .,
                      y = multiCETbl_u,
                      by = c("InChiKey", "Note")) %>%
                  select(-CEBest)
          })



# set rules ---------------------------------------------------------------

xRule <- function(n) {
  neg <- c('Cl', 'CL', 'CHO2', 'HCOO','CH3COO', 'CH3CO2')
  negF <- c('Cl', 'CL', 'CHO2', 'HCOO','CH3COO', 'CH3CO2')
  negR <- c(F,F,F,F,F,F)

  pos <- c('H','Na', 'NA','K','NH4')
  posF <- c('H','Na', 'NA', 'K','NH4')
  posR <- c(F,F,F,F,F)

  neu <- c("H2O", "FA" ,'DMSO', "ACN")
  neuF <- c("H2O", "CH2O2" ,'C2H6OS', "C2NH3")
  neuR <- c(F,T,T,T)

  ruleTbl <- tibble(
    key = c(neg, pos, neu),
    keyF = c(negF, posF, neuF),
    keyR = c(negR, posR, neuR),
    sign = c(rep(-1, length(neg)), rep(1, length(pos)), rep(0, length(neu)))
  )
  if (n == 1) {
    return(ruleTbl)
  } else{
    bind_rows(
      lapply(2:n, function(x) {
        ruleTbl %>%
          mutate(
            key = paste(x, key, sep = ""),
            keyF = paste(x, keyF, sep = ""),
            keyR = keyR,
            sign = sign * x
          )
      }) %>% bind_rows(ruleTbl, .) %>%
        mutate(
          key = paste("+", key, sep = ""),
          keyF = paste("+", keyF, sep = ""),
          keyR = keyR
        ),
      lapply(2:n, function(x) {
        ruleTbl %>%
          mutate(
            key = paste(x, key, sep = ""),
            keyF = paste(x, keyF, sep = ""),
            keyR = keyR,
            sign = sign * -x
          )
      }) %>% bind_rows(ruleTbl %>% mutate(sign = sign * -1), .) %>%
        mutate(
          key = paste("-", key, sep = ""),
          keyF = paste("-", keyF, sep = ""),
          keyR = keyR
        )
    )
  }
}


# xAdduct g function --------------------------------------------------------------------------
#' calculate adduct formula and mz
#'
#' @param cls parallel cluster object
#' @param adducts vector of adduct
#' @param formulas vector of formula
#' @param polarities vector of polarity
#'
#' @return tibble with 5 columns
#' @export NULL
#'
#' @examples NULL
xAdduct <- function(cls, adducts, formulas, polarities) {
  # load elements info ------------------------------------------------------
  elsTbl <-
    read_table(
      file = system.file("DB/ICIS_elements.els", package = "MRMlib"),
      col_names = c("element", "mass", "proportion", "other")
    )
  monoMassTbl <-
    elsTbl %>% arrange(element, desc(proportion))
  monoMass <-
    monoMassTbl[!(duplicated(monoMassTbl$element)),]
  elements <- monoMass$mass
  names(elements) <- monoMass$element


  adducts <-
    adducts %>% stri_replace_all(., "", regex = "[^[:alnum:]|[-|+|\\]|\\[]]")
  adducts[adducts %in% c("", "NA", NA)] <-
    paste("M", polarities[adducts %in% c("", "NA", NA)], sep = "")

  fml <-
    stri_extract(adducts, regex = "[:alnum:]+.*[:alnum:]|[:digit:]*[M]")
  chr <-
    stri_extract_last_regex(adducts, "[:digit:]*[-|+](?=$|[\\]])")
  pol <- polarities


  rule3 <- xRule(n = 3)
  ruleR <- rule3  %>% dplyr::filter(keyR)

  formula0 <-
    formulas %>% stri_replace_all(., "", regex = "[-|+]")
  parLapply(cls, seq_along(fml), function(x) {
    require(MRMlib)
    require(tidyverse)
    require(stringi)
    M <- stri_extract(fml[x], regex = '[:digit:]*[M]')
    A <-
      stri_extract_all(fml[x], regex = '[+|-][:digit:]*[:alnum:]+')[[1]]
    if (is.na(A)[1]) {
      if (chr[x] %in% c("-", "+")) {
        chrg <- chr[x]
      } else {
        chrg <- pol[x]
      }
    } else{
      charge <-
        rule3 %>%
        filter(key %in% A | keyF %in% A) %>%
        .$sign %>%
        sum(., na.rm = T)
      if (charge == 1) {
        chrg <- "+"
      } else if (charge == -1) {
        chrg <- "-"
      } else if (charge > 1) {
        chrg <- sprintf("%s+", charge)
      } else if (charge < -1) {
        chrg <- sprintf("%s-", abs(charge))
      } else if (charge == 0) {
        chrg <- ""
      }
    }

    MAr <- MA <-
      sprintf("[%s]%s", paste(na.omit(c(M, A)), collapse = "", sep = ""), chrg)
    Ar <- A
    for (a in seq_along(ruleR$key)) {
      MAr <- stri_replace_all(MAr, ruleR$keyF[a], fixed = ruleR$key[a])
      Ar <-
        stri_replace_all(Ar, ruleR$keyF[a], fixed = ruleR$key[a])
    }
    Mr <- stri_replace(M, formula0[x], fixed = "M")

    FormulaL <-
      paste(c("+", Mr, Ar), sep="", collapse = "") %>%
      stri_replace_all(., "1", regex = "(?<=[-|+])(?=[:upper:])") %>%
      stri_split(., regex = "(?<=[:alnum:])(?=[-|+][:digit:]?)") %>% .[[1]] %>%
      stri_split_regex("(?=[:upper:])")

    countTbl <-
      lapply(FormulaL, function(y) {
        mt <-
          stri_split(y[-1], regex = "(?<=[:alpha:])(?=[:digit:]|$)") %>%
          do.call(rbind, .)
        mt[mt == ""] <- 1
        tibble(element = mt[, 1],
               number = as.numeric(mt[, 2]) * as.numeric(y[1]))
      }) %>%
      bind_rows() %>%
      group_by(element) %>%
      summarise(N = sum(number)) %>%
      mutate(mz = N * elements[element]) %>%
      dplyr::filter(N > 0)

    tibble(
      Adduct0 = fml[x],
      Adduct1 = MA,
      AdductR = MAr,
      pFormula = paste(
        countTbl$element,
        countTbl$N,
        collapse = "",
        sep = ""
      ),
      PrecursorMz = sum(countTbl$mz),
      pCharge = chrg
    )
  }) %>% bind_rows()
}


# clean MSn -----------------------------------------------------------------------------------

setGeneric("cleanMSn", function(mz, tol) standardGeneric("cleanMSn"))
setMethod("cleanMSn", c(mz = "matrix", tol = "numeric"), function(mz, tol){
    X <- lapply(seq_along(mz[,1]), function(i){
        j <- abs(mz[i,1] - mz[,1]) <= tol
        xcms::binYonX(
            mz[, 1],
            mz[, 1],
            nBins = 1,
            method = "mean",
            fromIdx = min(seq_along(j)[j]),
            toIdx = max(seq_along(j)[j])
        ) %>% unlist()
    }) %>% do.call(rbind,.) %>% as_tibble()

    Y <- lapply(seq_along(mz[,1]), function(i){
        j <- abs(mz[i,1] - mz[,1]) <= tol
        xcms::binYonX(
            mz[, 1],
            mz[, 2],
            nBins = 1,
            method = "sum",
            fromIdx = min(seq_along(j)[j]),
            toIdx = max(seq_along(j)[j])
        ) %>% unlist()
    }) %>% do.call(rbind,.) %>% as_tibble()

    left_join(X, Y, by = "x") %>%
        transmute(mz = y.x, intensity = y.y) %>%
        unique() %>% as.matrix()
})



setGeneric("cleanMSn.iter", function(mz,tol,echo = F) standardGeneric("cleanMSn.iter"))
setMethod("cleanMSn.iter", c(mz = "matrix", tol = "numeric"), function(mz, tol, echo = F) {
    i <- NULL
    while (!identical(i, mz)) {
        i <- mz
        mz <- cleanMSn(mz = mz, tol = tol)
        if(echo) cat(nrow(mz),'\t')
    }
    if(echo) cat("\n")

    return(mz)
})



# pre filter --------------------------------------------------------------

mspPreFilter <- function(x = range_idx["1"], lib_vct = lib) {
  require(MRMlib)
  require(tidyverse)
  require(stringi)

  x <- unlist(x)
  lib_val_x <- lib_vct[x[1]:x[2]]

  info_x <-
    lib_val_x %>%
    grep(fixed = F,
         pattern = ".+[=|:][ ]?",
         value = T) %>%
    stri_replace_first(replacement = "",
                       regex = ".+[=|:][ ]?") %>%
    as.list()

  names(info_x) <-
    lib_val_x %>%
    grep(fixed = F,
         pattern = ".+[=|:][ ]?",
         value = T) %>%
    lapply(., stri_extract, regex = "^[:print:]+(?==|:)")
  return(info_x$id)

}


# intersection of full-MS and Library ---------------------------------------------------------

#' intersection of full-MS and Library
#'
#' @param MS cross ms to db by ppm
#' @param DB db tibble
#' @param tol mz tolerance ppm
#' @param N segment size
#'
#' @return filtered db tibble
#' @export NULL
#'
#' @examples NULL
setGeneric("xMStoDB", function(MS,DB,tol,N) standardGeneric("xMStoDB"))
setMethod("xMStoDB", c( MS = "data.frame",DB = "data.frame",tol = "numeric",N = "numeric"),
          function(MS, DB, tol, N) {
            MS <- MS %>% as_tibble()
            db <- db0 <- DB
            segment <-
              c(seq(1, nrow(MS), by = ceiling(nrow(MS) / ceiling(nrow(MS) / N))),
                nrow(MS))
            db <- db0 %>% select(splash, PrecursorMz, pFormula) %>% unique()
            dbmz <- db$PrecursorMz
            seg <- list()
            for (i in seq_along(segment)[-1]) {
              message(sprintf("Start segment %s: %s to %s ...",i-1, segment[i-1], segment[i]))

              ms <- MS[segment[i - 1]:segment[i], ]
              msmz <- ms$mz
              mzMatrix <-
                lapply(dbmz, function(d) {
                  abs(d - msmz) / d * 1000000
                }) %>% do.call(rbind, .)
              rtMatrix <-
                matrix(
                  data = rep(ms$rt, each = length(dbmz)),
                  nrow = length(dbmz),
                  ncol = length(msmz)
                )

              mzMatrixL <- mzMatrix * (mzMatrix <= tol)
              rtMatrixL <- rtMatrix * (mzMatrix <= tol)
              Hit <-  apply((mzMatrix <= tol), 1, any, na.rm = T)
              rtList <- apply(rtMatrixL, 1, function(x) {x[x > 0]})[Hit]
              mzList <- apply(mzMatrixL, 1, function(x) {x[x > 0]})[Hit]

              dbList <- db$splash[Hit]

              seg[[i-1]] <-
                lapply(seq_along(dbList),
                       function(i) {
                         rts <- rtList[[i]]
                         mzs <- mzList[[i]]
                         tibble(splash = dbList[i],
                                # PrecursorMz = db[Hit, ][i, 2][[1]],
                                # pFormula = db[Hit, ][i, 3][[1]],
                                RTs = as.numeric(rts),
                                ppm = as.numeric(mzs)) %>% dplyr::filter(RTs > 0)
                       }) %>% bind_rows() %>% left_join(., db0 %>% select(-RTs), by = "splash")
              gc()
            }
            gc()
            return(seg %>% bind_rows() %>% unique())
          })

# summary mspList ---------------------------------------------------------


# de-duplicate transitions --------------------------------------------------------------------

#' Merge overlapping transitions
#'
#' @param traTbl transition tibble (from toSkyline() and subset by InChIKey)
#' @param RTw RT window
#'
#' @return merged traTbl
#' @export NULL
#'
#' @examples NULL
mergeOverlap <- function(traTbl, RTw = 2) {

  traTbl <- mutate(traTbl,
                   mergeRT = `Explicit Retention Time`,
                   mergeRW = RTw)
  RTs <- traTbl$`Explicit Retention Time` %>% unique() %>% sort()

  for (i in seq_along(RTs)) {

    IDX0 <- traTbl$mergeRT == RTs[i]
    IDX1 <- traTbl$mergeRT == RTs[i + 1]
    RTw0 <- traTbl$mergeRW[IDX0] %>% unique()
    RTw1 <- traTbl$mergeRW[IDX1] %>% unique()

    if (i == length(RTs)) {

      traTbl_dedup <-
        traTbl %>%
        mutate(`Explicit Retention Time` = mergeRT,
               `Explicit Retention Time Window` = mergeRW) %>%
        select(-mergeRT, -mergeRW) %>%
        unique()
      return(traTbl_dedup)

    } else{

      RT0 <- traTbl$mergeRT[IDX0] %>% unique()
      RT1 <- traTbl$mergeRT[IDX1] %>% unique()
      RT0r <- c(RT0 - RTw0, RT0 + RTw0)
      RT1r <- c(RT1 - RTw1, RT1 + RTw1)
      if (RT0r[2] >= RT1r[1]) {
        RTir <- range(c(RT0r, RT1r))
        traTbl$mergeRT[IDX0 | IDX1] <- mean(RTir)
        traTbl$mergeRW[IDX0 | IDX1] <- diff(RTir) / 2
        RTs[i:(i+1)] <- mean(RTir)
      }
    }
  }
}



# convert skylineTransition to DB.msp -------------------------------------

#' convert skylineTransition to DB.msp
#'
#' @param X retruned by traTbl_join %>% toSkyline()
#' @param File file dir to save
#' @param overWrite if overwrite exist .msp file
#'
#' @return NULL
#' @export NULL
#'
#' @examples NULL
toSkylineDB <- function(X, File = "./skyline-db.msp", overWrite = T) {
  tra_list <- X %>% unique() %>% split.data.frame(f = .$InChiKey)
  if(overWrite & file.exists(File)){
    file.remove(File)
  }
  for (i in tra_list) {
    paste(
      sprintf("Name: %s\n", i$`Precursor Name`[1]),
      sprintf("InChIKey: %s\n", i$InChiKey[1]),
      sprintf("Precursor_type: %s\n", i$`Precursor Adduct`[1]),
      sprintf("Spectrum_type: %s\n", 2),
      sprintf("PrecursorMZ: %s\n", i$`Precursor m/z`[1]),
      sprintf("Ion_mode: %s\n", ""),
      sprintf("Collision_energy: %s\n", i$`Explicit Collision Energy`[1]),
      sprintf("Formula: %s\n", i$`Precursor Formula`[1]),
      sprintf("Num Peaks: %s\n", nrow(i)),
      sprintf("%.4f\t%.4f\n", i$`Product m/z`, i$`Product intensity`) %>% paste(collapse = ""),
      "\n",
      sep = "",
      collapse = ""
    ) %>% write_file(path = File, append = T)
  }
}


# batch function MRMlib ---------------------------------------------------
setClass(
  "MRMlibrary",
  slots = c(
    RSQLite = "SQLiteConnection",
    param = "list",
    traTbl = "tbl",
    traTbl_join = "tbl",
    tra_list = "list",
    peak = "data.frame",
    xDB = "tbl",
    skylineTra = "tbl"
  )
)

MRMlib <- function(
  peak = peak0,
  paraN = 8,
  raw_msp = "~/R/DB/Plant.msp",
  infoDir = "~/R/DB/DBraw/",
  DIR = getwd(),
  SQLiteName = "20180904",
  mspParser.nbTol = 0.8,
  mspParser.mz_tol = 0.01,
  mspParser.mz_only = T,
  mspParser.ignoreInt = F,
  filterMSn.topX = 5,
  filterMSn.a = 0.5,
  filterMSn.b = 0.5,
  filterMSn.type = "local",
  xMStoDB.tol = 5,
  xMStoDB.N = 8000,
  toSkyline.deltaMz = 12,
  mergeOverlap.RTw = 2) {

  require(tidyverse)
  require(stringi)
  require(ChemmineR)
  require(MSnbase)
  require(MRMlib)
  require(RSQLite)
  require(dbplyr)

  message("New MRMlibrary class ...")
  obj <- new(Class = "MRMlibrary")
  cl <- makeCluster(paraN)

  # parameters --------------------------------------------------------------

  if (!dir.exists(DIR)) {dir.create(DIR)}

  SQLiteDir <- sprintf("%s/%s.sqlite", DIR, SQLiteName)
  tra.msp <- sprintf("%s/skyline.msp", DIR)
  tra.csv <- sprintf("%s/MRMtransition", DIR)
  msp_list.rsd <- sprintf("%s/msp_list.rds", DIR)

  obj@param <- list(
    paraN = paraN,
    raw_msp = raw_msp,
    DIR = DIR,
    SQLiteName = SQLiteName,
    mspParser.nbTol = mspParser.nbTol,
    mspParser.mz_tol = mspParser.mz_tol,
    mspParser.mz_only = mspParser.mz_only,
    mspParser.ignoreInt = mspParser.ignoreInt,
    filterMSn.topX = filterMSn.topX,
    filterMSn.a = filterMSn.a,
    filterMSn.b = filterMSn.b,
    filterMSn.type = filterMSn.type,
    xMStoDB.tol = xMStoDB.tol,
    xMStoDB.N = xMStoDB.N,
    toSkyline.deltaMz = toSkyline.deltaMz,
    mergeOverlap.RTw = mergeOverlap.RTw
  )

  # main --------------------------------------------------------------------

  obj@RSQLite <-  dbConnect(RSQLite::SQLite(), SQLiteDir)

  # load from SQLite --------------------------------------------------------
  if (db_has_table(obj@RSQLite, "hmdbInfo")) {
    hmdbInfo_list <- tbl(obj@RSQLite, "hmdbInfo") %>% as_tibble()
    initInfo_list <- tbl(obj@RSQLite, "initInfo") %>% as_tibble()
    message("load existing info table ...")
  } else{
    # load raw DB infos -------------------------------------------------------
    env <- new.env()

    hmdbInfo <-
      load(dir(
        path = infoDir,
        full.names = T,
        pattern = "hmdbInfo.QT"
      ), env)

    matchIdx <-
      load(dir(
        path = infoDir,
        full.names = T,
        pattern = "initalMatch2HmdbIndex.QT"
      ), env)

    initInfo <-
      load(dir(
        path = infoDir,
        full.names = T,
        pattern = "initalInfo.QT"
      ), env)

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

    initInfo_list <-
      as_tibble(initInfo_list) %>% bind_cols(matchIdx_list, .)
    hmdbInfo_list <-
      as_tibble(hmdbInfo_list) %>% mutate(IDX = seq(nrow(.)))

    # SQlite DB ---------------------------------------------------------------
    message("writing info lists ...")
    dbWriteTable(obj@RSQLite,
                 name = 'hmdbInfo',
                 hmdbInfo_list,
                 overwrite = T)
    dbWriteTable(obj@RSQLite,
                 name = 'initInfo',
                 initInfo_list,
                 overwrite = T)
  }

  # save msp_list to DB -----------------------------------------------------
  if (file.exists(msp_list.rsd)) {
    message("loading existing msp_list ...")
    msp_list <- readRDS(file = msp_list.rsd)
  } else{
    message("new msp_list ...")
    lib <- read_lines(file = raw_msp)
    begin_idx <-
      grep(
        pattern = "BEGIN",
        ignore.case = F,
        fixed = T,
        x = lib
      )
    end_idx <-
      grep(
        pattern = "END",
        ignore.case = F,
        fixed = T,
        x = lib
      )

    range_idx <-
      tibble(start = begin_idx + 1,
             end = end_idx - 1,
             idx = seq_along(end_idx)) %>%
      split.data.frame(f = .$idx)


    # mspParser.nbTol = 0.8
    # mspParser.mz_tol = 0.01
    # mspParser.mz_only = T
    # mspParser.ignoreInt = F

    msp_list <-
      parLapply(
        cl,
        range_idx,
        mspParser,
        lib_vct = lib,
        nbTol = mspParser.nbTol,
        mz_tol = mspParser.mz_tol,
        mz_only = mspParser.mz_only,
        ignoreInt = mspParser.ignoreInt
      )
    msp_list <- parLapply(
      cl,
      msp_list,
      fun = function(x) {
        if (x$precursorType == "" | is.null(x$precursorType)) {
          if (x$polarity == "" | is.null(x$polarity)) {
            x$polarity <- "+"
            cat(x$formula)
          }
          x$precursorType <-
            sprintf("[M%sH]%s", x$polarity, x$polarity)
        }
        return(x)
      }
    )

    message("OK, saving msp_list ...")
    saveRDS(msp_list, file = msp_list.rsd)
  }

  # traTbl ------------------------------------------------------------------
  if (db_has_table(obj@RSQLite, "traTbl")) {
    message("loading existing tarTbl ...")
    traTbl <- tbl(obj@RSQLite, "tarTbl") %>% as_tibble()
  } else {
    message("new tarTbl ...")
    metaMsn_i <- new('metaMSn', MSn = msp_list)

    # filterMSn.topX = 5
    # filterMSn.a = 0.5
    # filterMSn.b = 0.5
    # filterMSn.type = "local"

    traTbl <-
      filterMSn(
        object = metaMsn_i,
        topX = filterMSn.topX,
        type = filterMSn.type,
        a = filterMSn.a,
        b = filterMSn.b,
        cluster = cl
      )
    message("writing traTbl ...")
    dbWriteTable(obj@RSQLite, "tarTbl", traTbl, overwrite = T)
  }

  if (db_has_table(obj@RSQLite, "traTbl_join")) {
    message("loading existing traTbl_join ...")
    traTbl_join <- tbl(obj@RSQLite, "traTbl_join") %>% as_tibble()
  } else {
    message("new traTbl_join ...")
    traTbl_join <- traTbl %>% as_tibble() %>%
      left_join(., initInfo_list, by = c("splash" = "splash10")) %>%
      left_join(., hmdbInfo_list, by = c("match2HmdbIndex" = "IDX")) %>%
      dplyr::filter(ExactMass - PrecursorMz < 2) %>%
      dplyr::filter(!stri_detect(Formula, regex = "[D]")) %>%
      dplyr::filter(Polarity == "+")
    message("writing traTbl_join ...")
    dbWriteTable(obj@RSQLite, "traTbl_join", traTbl_join, overwrite = T)

  }

  if (db_has_table(obj@RSQLite, "peak") & is.null(peak)) {
    message("loading existing peak ...")
    peak <- tbl(obj@RSQLite, "peak") %>% as_tibble()
  } else{
    message("use external peak ...")
    dbWriteTable(obj@RSQLite, "peak", peak, overwrite = T)
  }

  # cross DB and DATA -------------------------------------------------------
  # xMStoDB.tol = 5
  # xMStoDB.N = 8000
  # toSkyline.deltaMz = 12
  if (db_has_table(obj@RSQLite, "xDB")) {
    message("loading existing xDB ...")
    xDB <- tbl(obj@RSQLite, "xDB") %>% as_tibble()
  } else {
    message("new xDB ...")
    xDB <-
      xMStoDB(
        MS = peak,
        DB = traTbl_join,
        tol = xMStoDB.tol,
        N = xMStoDB.N
      )

    message("writing xDB ...")
    dbWriteTable(obj@RSQLite, "xDB", xDB, overwrite = T)
  }

  skylineTra <-
    toSkyline(infoTibble = xDB, deltaMz = toSkyline.deltaMz) %>% unique()

  # write to skyline DB.msp -------------------------------------------------
  if (db_has_table(obj@RSQLite, "tra_list")) {
    message("loading existing tra_list ...")
    tra_list <- tbl(obj@RSQLite, "tra_list") %>% as_tibble()
    message("writing tra.msp file ...")
    toSkylineDB(X = tra_list,
                overWrite = T,
                File = tra.msp)
  } else {
    message("new tra_list ...")
    tra_list <- traTbl_join %>% toSkyline(deltaMz = 12)
    message("writing tra.msp file ...")
    toSkylineDB(X = tra_list,
                overWrite = T,
                File = tra.msp)
    dbWriteTable(obj@RSQLite, "tra_list", tra_list, overwrite = T)
  }

  # write to transition table -----------------------------------------------
  # mergeOverlap.RTw = 2

  skylineTra_2 <-
    skylineTra %>%
    split.data.frame(f = .$InChiKey) %>%
    lapply(FUN = mergeOverlap, RTw = mergeOverlap.RTw) %>%
    bind_rows() %>%
    mutate(`Precursor Name` = paste(`Precursor Name`,
                                    `Explicit Retention Time`,
                                    sep ="--"))

  write_csv(
    skylineTra_2 %>% select(-InChiKey),
    append = F,
    path = sprintf("%s-%s.csv", tra.csv, Sys.Date())
  )
  message("All done ...")
  stopCluster(cl)

  # return MRMlib object ----------------------------------------------------

  obj@traTbl <- traTbl
  obj@traTbl_join <- traTbl_join
  obj@tra_list <- tra_list
  obj@xDB <- xDB
  obj@peak <- peak0
  obj@skylineTra <- skylineTra
  return(obj)
}

