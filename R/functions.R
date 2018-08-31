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

find_isotopes <- function(i, mz_x, mz_tol, binMz=0.5, mz_only = T, predictedIso) {
        .isotopeN <- NULL
        .mz_tol <- abs(mz_tol)

        # if (.mz_tol > 1) {
        #     delta_type = "ppm"
        # } else if (.mz_tol <= 1) {
        #     delta_type = "abs"
        # }

        for (j in seq(predictedIso)) {
            # is_iso_mz <- switch(
            #     delta_type,
            #     ppm = (abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) /
            #         mz_x[seq(i), 1]) <= (.mz_tol / 1000000),
            #     abs = abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) <= .mz_tol
            # )
          if (.mz_tol > 1) {
            # delta_type = "ppm"
            is_iso_mz <- (abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) /
                                    mz_x[seq(i), 1]) <= (.mz_tol / 1000000)
          } else if (.mz_tol <= 1) {
            is_iso_mz <- abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) <= .mz_tol
            # delta_type = "abs"
          }


            # is_close <- abs(mz_x[i, 1] - mz_x[seq(i), 1]) <= binMz
            if(mz_only){
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
mspParser <- function(x = range_idx["1"],
                      quantile_2 = c(0.05,0.95),
                      cpdLib = "kegg",
                      iso1 = 0,
                      iso2 = 1:3,
                      lib_vct = lib,
                      nbTol = 0.8,
                      mz_tol = 0.01,
                      mz_only= T, ...) {
  require(dplyr)
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
        cleanMSn.iter(tol = nbTol,echo = T )


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


# load elements info ------------------------------------------------------



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
setGeneric("filterMSn",
           function(object,topX,type,cluster,...) standardGeneric("filterMSn"))

setMethod("filterMSn", c(object = "metaMSn"),
          function(object,topX=5,type="local",cluster = cl,...) {
            if (type == "local") {
              .fun <- function(objList) {
                require(MRMlib)
                require(tidyverse)
                require(stringi)

                dplyr::filter(objList$MSn,!iso1 & !iso2 & !iso3) %>%
                  dplyr::select(mz, intensity) %>%
                  dplyr::mutate(
                    splash = objList$id,
                    Rank = rank(-intensity, ties.method = "first"),
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
                  dplyr::filter(Rank <= topX)
              } else{
                cat("step 1 start ...\n")
                infoTibble <-
                  parallel::parLapply(cluster, object@MSn, .fun) %>%
                  do.call(bind_rows, args = .) %>%
                  dplyr::filter(Rank <= topX)
                cat("step 1 done ...\n")
              }
            }
            cat("step 2 start ...\n")

            # infoTibble <- infoTibble %>% filter(Adduct0 == "[M+2H]+")

           xAdduct(
                cls = cluster,
                formulas = infoTibble$Formula,
                adducts = infoTibble$Adduct0,
                polarities = infoTibble$Polarity
              ) %>%
              bind_cols(infoTibble, .) %>%
              mutate(DeltaMz = PrecursorMz - mz, RTs = 0)
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
                  mutate(CE =ifelse(is.na(CE)|CE=="", 30.1, CE)) %>%
                  transmute(
                      `Note` = splash,
                      `Molecule List Name` = Class,
                      `Precursor Name` = bestName(Name, initialName),
                      `Precursor Formula` = Formula,
                      `Precursor Adduct` = Adduct1,
                      `Precursor Charge` = chargeRule[as.character(pCharge)],
                      `Precursor m/z` = PrecursorMz/`Precursor Charge`,
                      `Product m/z` = mz,
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

              multiCETbl_u <-
                multiCETbl  %>%
                group_by(!!!lapply(names(multiCETbl)[-3], as.symbol)) %>%
                summarise(`Precursor Name` = first(`Precursor Name`))

              # xxx <-
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
<<<<<<< HEAD
      paste(c("+", Mr, Ar),sep="", collapse = "") %>%
=======
      paste(c("+", Mr, Ar), sep="", collapse = "") %>%
>>>>>>> 03ad0f776855de7e90fd88837a6bb4fa9389f36c
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

#' Title
#'
#' @param ms cross ms to db by ppm
#' @param db db tibble
#' @param tol mz tolerance ppm
#'
#' @return filtered db tibble
#' @export NULL
#'
#' @examples NULL
setGeneric("xMStoDB", function(ms,db,tol) standardGeneric("xMStoDB"))
setMethod("xMStoDB", c(ms = "data.frame", db = "data.frame", tol = "numeric"), function(ms,db,tol){
    ms <- ms %>% as_tibble()
    db0 <- db
    db <- db0 %>% select(splash, PrecursorMz, pFormula) %>% unique()
    dbmz <- db$PrecursorMz
    msmz <- ms$mz

    mzMatrix <-
        lapply(dbmz, function(d) {
            abs(d - msmz)/d * 1000000
        }) %>% do.call(rbind, .)
    rtMatrix <-
        matrix(data = rep(ms$rt, each=length(dbmz)),
               nrow = length(dbmz),
               ncol = length(msmz))

    mzMatrixL <- mzMatrix * (mzMatrix <= tol)

    rtMatrixL <- rtMatrix * (mzMatrix <= tol)

    Hit <-  apply((mzMatrix <= tol), 1, any, na.rm=T)

    rtList <- apply(rtMatrixL, 1, function(x){
      x[x>0]
    })[Hit]

    mzList <- apply(mzMatrixL, 1, function(x){
      x[x>0]
    })[Hit]

    dbList <- db$splash[Hit]
# xxx <-
    lapply(
      seq_along(dbList),
      function(i) {
        rts <- rtList[[i]]
        mzs <- mzList[[i]]
        tibble(
          splash = dbList[i],
          # PrecursorMz = db[Hit, ][i, 2][[1]],
          # pFormula = db[Hit, ][i, 3][[1]],
          RTs = as.numeric(rts),
          ppm = as.numeric(mzs)
        ) %>% dplyr::filter(RTs > 0)
      }
    ) %>% bind_rows() %>% left_join(., db0 %>% select(-RTs), by = "splash")
})

# summary mspList ---------------------------------------------------------
setGeneric("summary", function(x = 'list') standardGeneric("summary"))
























