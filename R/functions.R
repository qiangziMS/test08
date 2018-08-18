require(tidyverse)
require(stringi)
require(ChemmineR)
require(MSnbase)
require(parallel)
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

        if (.mz_tol > 1) {
            delta_type = "ppm"
        } else if (.mz_tol <= 1) {
            delta_type = "abs"
        }

        for (j in seq(predictedIso)) {
            is_iso_mz <- switch(
                delta_type,
                ppm = abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) /
                    mz_x[seq(i), 1] <= j*.mz_tol / 1000000,
                abs = abs(mz_x[i, 1] - mz_x[seq(i), 1] - 1.00336 * j) <= j*.mz_tol
            )

            is_close <- abs(mz_x[i, 1] - mz_x[seq(i), 1]) <= binMz
            if(mz_only){
                is_iso_all <- is_iso_mz & mz_x[seq(i), 2] > mz_x[i, 2]
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
    # source("F://R&D/Rscript/functions.R", echo = F)
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
        stri_extract(str = datablock[idx], regex = "(?<=^> <).+(?=>$)")
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

setGeneric("xFormula", function(formulaRaw, return = NULL) standardGeneric("xFormula"))
setMethod("xFormula", c(formulaRaw = "character"), function(formulaRaw, return) {
    # load elements info ------------------------------------------------------
    elsTbl <-
        read_table(
        file = system.file("DB/ICIS_elements.els", package = "MRMlib"),
        col_names = c("element", "mass", "proportion", "other"))
    monoMassTbl <-
        elsTbl %>% arrange(element, desc(proportion))
    monoMass <-
        monoMassTbl[!(duplicated(monoMassTbl$element)), ]
    elements <- monoMass$mass
    names(elements) <- monoMass$element

    # function ----------------------------------------------------------------
    patternX <-
        "(?<=[:upper:])(?=[:upper:])|(?<=[:digit:]|[+|-])(?=[:upper:])"
    cpdInfo <-
        lapply(formulaRaw,
               function(formulaRaw_i) {
                   formula_tbl <-
                       formulaRaw_i %>%
                       stri_extract(regex = "(?<=\\[|^)([:alnum:]|[-|+])+(?=\\]|$)") %>%
                       stri_split(., regex = "(?=[-|+][:digit:]?)") %>% .[[1]] %>%
                       stri_split_regex(patternX) %>%
                       lapply(
                           FUN = function(x) {
                               if (!stri_detect(x[1], regex = "^[-|+]"))
                                   x <- c("+1", x)
                               if (x[1] == "+")
                                   x[1] <- "+1"
                               if (x[1] == "-")
                                   x[1] <- "-1"
                               mt <- stri_split(x[-1],
                                                regex = "(?<=[:alpha:])(?=[:digit:]|$)") %>%
                                   do.call(rbind, .)
                               mt[mt == ""] <- 1
                               tibble(element = mt[, 1],
                                      number = as.numeric(mt[, 2]) * as.numeric(x[1]))
                           }
                       ) %>%
                       do.call(bind_rows, .) %>%
                       group_by(element) %>%
                       summarise(N = sum(number)) %>%
                       mutate(mz = N * elements[element]) %>%
                       dplyr::filter(N > 0)
                   tibble(
                       formula = paste(
                           formula_tbl$element,
                           formula_tbl$N,
                           collapse = "",
                           sep = ""
                       ),
                       mz = sum(formula_tbl$mz)
                   )
               }) %>% bind_rows()
    if (return == "formula") {
        return(cpdInfo$formula)
    } else if (return == "mz") {
        return(cpdInfo$mz)
    } else {
        return(cpdInfo)
    }
})



#' Filter product ions to MRM
#'
#' @param object metaMSn object
#' @param topX select top X product ions
#' @param ... NOT used
#' @param type local/...
#' @param newRule adducts rule
#'
#' @return Tibble of products and info
#' @export NULL
#'
#' @examples NULL
setGeneric("filterMSn",
           function(object,topX,type,newRule,cluster,...) standardGeneric("filterMSn"))

setMethod("filterMSn", c(object = "metaMSn"),
          function(object,topX=5,type="local",newRule = NULL,cluster = NULL,...) {

            if (type == "local") {
              .fun <- function(objList) {
                # source("F://R&D/Rscript/functions.R", echo = F)
                require(MRMlib)
                require(tidyverse)
                require(stringi)
                # require(parallel)

                dplyr::filter(objList$MSn,!iso1 &
                                !iso2 & !iso3) %>%
                  dplyr::select(mz, intensity) %>%
                  dplyr::mutate(
                    splash = objList$id,
                    Rank = rank(-intensity, ties.method = "first"),
                    Adduct = ifelse(
                      is.null(objList$precursorType),
                      NA,
                      xAdduct(objList$precursorType)
                    ),
                    Adduct0 = objList$precursorType,
                    Polarity = ifelse(is.null(objList$polarity),
                                      NA, objList$polarity),
                    ExactMass = ifelse(is.null(objList$exact_mass),
                                       NA, objList$exact_mass) %>% as.numeric(),
                    Formula = ifelse(is.null(objList$formula),
                                     NA, objList$formula),
                    CE = ifelse(
                      is.null(objList$collisionEnergy),
                      NA,
                      objList$collisionEnergy
                    ),
                    InstrumentType = ifelse(
                      is.null(objList$instrument_type),
                      NA,
                      objList$instrument_type
                    )
                  ) %>% as_tibble()
                # cat("----",'\n')
              }
              # object <- metaMsn_i

              if (is.null(cluster)) {
                infoTibble <-
                  lapply(object@MSn, .fun) %>%
                  do.call(bind_rows, args = .) %>%
                  dplyr::filter(Rank <= topX)
              } else{
                # clt <- parallel::makeCluster(cluster)
                cat("cluster starting ...\n")
                infoTibble <-
                  parallel::parLapply(cluster, object@MSn, .fun) %>%
                  do.call(bind_rows, args = .) %>%
                  dplyr::filter(Rank <= topX)
                # parallel::stopCluster(clt)
                cat("cluster stopping ...\n")
              }

            } else{
              tagMatchClass <-
                metaLink(object = object, tag = "Class")
              tagMatchSubclass <-
                metaLink(object = object, tag = "Subclass")
              tagMatchSuperclass <-
                metaLink(object = object, tag = "Superclass")

              NAImputation <- c(Adduct = "", CE = 30)
              infoTibble <-
                lapply(
                  object@MSn,
                  FUN = function(objList) {
                    dplyr::filter(objList$MSn,!iso1 & !iso2 & !iso3) %>%
                      dplyr::select(mz, intensity) %>%
                      dplyr::mutate(
                        Rank = rank(-intensity, ties.method = "first"),
                        InChIKey = ifelse(is.null(objList$InChIKey),
                                          NA, objList$InChIKey),
                        Name = ifelse(is.null(objList$Name),
                                      NA, objList$Name),
                        Adduct = ifelse(
                          is.null(objList$Precursor_type),
                          NAImputation[["Adduct"]],
                          objList$Precursor_type
                        ),
                        PrecursorMz = ifelse(is.null(objList$PrecursorMZ),
                                             NA,
                                             objList$PrecursorMZ) %>% as.numeric(),
                        DeltaMz = PrecursorMz - mz,
                        Polarity = ifelse(is.null(objList$Ion_mode),
                                          NA, objList$Ion_mode),
                        msLevel = ifelse(
                          is.null(objList$Spectrum_type),
                          NA,
                          objList$Spectrum_type
                        ),
                        ExactMass = ifelse(is.null(objList$ExactMass),
                                           NA,
                                           objList$ExactMass) %>% as.numeric(),
                        Formula = ifelse(is.null(objList$Formula),
                                         NA, objList$Formula),
                        CE = ifelse(
                          is.null(objList$Collision_energy),
                          NAImputation[["CE"]],
                          objList$Collision_energy
                        ) ,
                        InstrumentType = ifelse(
                          is.null(objList$Instrument_type),
                          NA,
                          objList$Instrument_type
                        ),
                        Class = ifelse(is.na(tagMatchClass[InChIKey]),
                                       "Unknown", tagMatchClass[InChIKey]),
                        Subclass = ifelse(is.na(tagMatchSubclass[InChIKey]),
                                          "Unknown",
                                          tagMatchSubclass[InChIKey]),
                        Superclass = ifelse(
                          is.na(tagMatchSuperclass[InChIKey]),
                          "Unknown",
                          tagMatchSuperclass[InChIKey]
                        )
                      ) %>% as_tibble()
                  }
                ) %>%
                do.call(bind_rows, args = .) %>%
                dplyr::filter(Rank <= topX &
                                msLevel == "MS2" &
                                DeltaMz > deltaMz)
            }

            rule <- c(FA = 'CH2O2', DMSO = "C2H6OS")
            rule <- c(rule, newRule[!(newRule %in% rule)])

            # for (a in names(rule)) {
            #   infoTibble <-
            #     infoTibble %>%
            #     mutate(Adduct0 = stri_replace_all(
            #       Adduct,
            #       rule[a],
            #       regex = sprintf("(?<=[-|+][:digit:]?)%s(?=[-|+]?)", a)
            #     ))
            # }
            # cat("replace done ...\n")
            # clt <- parallel::makeCluster(cluster)
            # system.time(
            # infoList <-
            #   infoTibble %>%
            #     split.data.frame(f = seq(nrow(.))) %>%
            #     parallel::parLapply(
            #       cl = cluster,
            #       X = .,
            #       fun = function(x) {
            #         require(tidyverse)
            #         require(stringi)
            #         require(MRMlib)
            #         x %>%
            #           mutate(Formula0 =
            #                    stri_replace_all(Formula, "", regex = "[-|+]+$")) %>%
            #           mutate(pFormula =
            #                    stri_replace_all(Adduct0,
            #                                     Formula0,
            #                                     regex = "M(?=[-|+]?)")) #%>%
            #                    # xFormula(return = "formula")) #%>%
            #           # mutate(PrecursorMz =
            #           #          xFormula(pFormula, return = "mz")) %>%
            #           # mutate(DeltaMz = PrecursorMz - mz)
            #       }
            #     )
            #   # )
            # # parallel::stopCluster(clt)
            # infoList %>% bind_rows()
            infoTibble
          })



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
                      `2+` = 2
                  )

              multiCETbl <-
                  infoTibble %>%
                  dplyr::filter(DeltaMz >= deltaMz) %>%
                  mutate(CE =ifelse(is.na(CE)|CE=="",30.1, CE)) %>%
                  transmute(
                      `Note` = splash,
                      `Molecule List Name` = Class,
                      `Precursor Name` = initialName,
                      `Precursor Formula` = Formula0,
                      `Precursor Adduct` = Adduct,
                      `Precursor Charge` = chargeRule[Polarity],
                      `Precursor m/z` = PrecursorMz,
                      `Product m/z` = mz,
                      `Product Charge` = chargeRule[Polarity],
                      `Explicit Retention Time` = RTs,
                      # `Product Name` =  round(PrecursorMz,2),
                      `InChiKey` = initialInChIKey,
                      `Explicit Collision Energy` =
                          stri_extract_all(
                              CE, regex = "[:digit:]+[.]?[:digit:]+") %>%
                          .[[1]] %>%
                          as.numeric() %>% median()
                  ) %>%
                  dplyr::filter(stri_detect_regex(
                      InChiKey, pattern = "[:upper:]+[-][:upper:]+[-][:upper:]"))

              multiCETbl  %>%
                  group_by(
                      InChiKey,
                      Note,
                      `Explicit Collision Energy`) %>%
                  summarise(CEBest = n()) %>%
                  split.data.frame(f = .$InChiKey) %>%
                  lapply(FUN = function(tbl) {
                          tbl[which.max(tbl$CEBest), ] %>% mutate(CEBest = T)
                      }) %>%
                  do.call(bind_rows, args = .) %>%
                  left_join(
                      x = multiCETbl,
                      y = .,
                      by = c("InChiKey", "Note","Explicit Collision Energy")) %>%
                  dplyr::filter(CEBest) %>% select(-CEBest)
          })




# set rules ---------------------------------------------------------------

xRule <- function(n) {
  neg <- c('Cl', 'CHO2', 'HCOO','CH3COO', 'CH3CO2')
  negF <- c('Cl', 'CHO2', 'HCOO','CH3COO', 'CH3CO2')

  pos <- c('H','Na','K','NH4')
  posF <- c('H','Na','K','NH4')

  neu <- c("H2O", "FA" ,'DMSO', "ACN")
  neuF <- c("H2O", "CH2O2" ,'C2H6OS', "CNH")

  ruleTbl <- tibble(
    key = c(neg, pos, neu),
    keyF = c(negF, posF, neuF),
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
            sign = sign * x
          )
      }) %>% bind_rows(ruleTbl, .) %>%
        mutate(
          key = paste("+", key, sep = ""),
          keyF = paste("+", keyF, sep = "")
        ),
      lapply(2:n, function(x) {
        ruleTbl %>%
          mutate(
            key = paste(x, key, sep = ""),
            keyF = paste(x, keyF, sep = ""),
            sign = sign * -x
          )
      }) %>% bind_rows(ruleTbl, .) %>%
        mutate(
          key = paste("-", key, sep = ""),
          keyF = paste("-", keyF, sep = "")
        )
    )
  }
}

# adduct <- c("[M-H]",'M-H-',"[2M+NH4-]2-","2M-2")
setGeneric("xAdduct", function(adduct) standardGeneric('xAdduct'))
setMethod('xAdduct', c("character"), function(adduct){
  adduct <- table(traTbl$Adduct0) %>% names()

 fml <- stri_extract(adduct, regex="[:alnum:]+.*[:alnum:]|[:digit:]*[M]")
 chr <- stri_extract_last_regex(adduct, ".")


 rule3 <- xRule(n = 3)

 lapply(fml, function(x){
   M <- stri_extract(x, regex = '[:digit:]*[M]')
   A <- stri_extract_all(x, regex = '[+|-][:digit:]*[:alnum:]+')[[1]]
   if(is.na(A)){
     charge <- stri_extract_last(x, regex = "[:digit:]?[-|+][:digit:]?")

   }else{
     charge <- rule3 %>% filter(key %in% A|keyF %in% A) %>% .$sign %>% sum(., na.rm = T)
     if(charge == 1) {
       chrg <- "+"
     }else if(charge == -1) {
       chrg <- "-"
     }else if(charge > 1){
       chrg <- sprinf("$s+", charge)
     }else if(charge < -1){
       chrg <- sprinf("$s-", abs(charge))
     }else if(chrg == 0)

   }


 })

 chr <- stri_extract_last(
   adduct, regex = "(?<![:alpha:]&^[M])([:digit:]?[-|+][:digit:]?)")


  stri_replace(adduct,"",regex="[-|+](?=$|[\\]])")
    a <- stri_extract_first(
        adduct, regex = "[:alnum:]+([:alnum:]|[-|+])?[:alpha:]?[:alnum:]?")
    b <- stri_extract_last(
        adduct, regex = "(?<![:alpha:]&^[M])([:digit:]?[-|+][:digit:]?)")
    sprintf('[%s]%s', a, b)
})

xAdduct(adduct = "[2M]+")
# xFormula("C13H22O4", 'mz')

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

mspPreFilter <- function(x = range_idx["1"],
                         lib_vct = lib) {
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

setGeneric("xMStoDB", function(ms,db,tol) standardGeneric("xMStoDB"))
setMethod("xMStoDB", c(ms = "data.frame", db = "data.frame", tol = "numeric"), function(ms,db,tol){
    ms <- ms %>% as_tibble()
    db0 <- db
    db <- db %>% select(splash, PrecursorMz, pFormula) %>% unique()
    dbmz <- db$PrecursorMz
    msmz <- ms$mz

    mzMatrix <-
        lapply(dbmz, function(d) {
            abs(d - msmz)/d * 1000000
        }) %>% do.call(rbind, .)


    rtMatrix <-
        matrix(data = rep(ms$rt, length(dbmz)),
               nrow = length(dbmz),
               ncol = length(msmz))

    mzMatrixL <- mzMatrix <= tol
    rtMatrixL <- rtMatrix * mzMatrixL
    Hit <-  apply(rtMatrixL, 1, sum, na.rm=T) > 0

    lapply(
      seq(sum(Hit)),
      FUN = function(i) {
        rts <- rtMatrixL[Hit, ][i, ]
        mzs <- mzMatrix[Hit, ][i, ]
        tibble(
          splash = db[Hit, ][i, 1][[1]],
          PrecursorMz = db[Hit, ][i, 2][[1]],
          pFormula = db[Hit, ][i, 3][[1]],
          RTs = as.numeric(rts),
          ppm = as.numeric(mzs)
        ) %>% dplyr::filter(RTs > 0)
      }
    ) %>%
      bind_rows()
})


# a <- 1:9
# b <- rep(5.1,9)


# summary mspList ---------------------------------------------------------

setGeneric("summary", function(x = 'list') standardGeneric("summary"))

# setMethod("summary", c(x = ))























