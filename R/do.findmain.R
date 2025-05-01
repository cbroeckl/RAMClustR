#' do.findmain
#'
#' Cluster annotation function: inference of 'M' - molecular weight of the compound giving rise to each spectrum - using the InterpretMSSpectrum::findMain function
#'
#' @param ramclustObj ramclustR object to annotate.
#' @param cmpd integer: vector defining compound numbers to annotated.  if NULL (default), all compounds
#' @param mode character: "positive" or "negative"
#' @param mzabs.error numeric: absolute mass deviation allowd, default = 0.01
#' @param ppm.error numeric: ppm mass error _added_ to mzabs.error, default = 10
#' @param ads character: vector of allowed adducts, i.e. c("[M+H]+"). if NULL, default positive mode values of H+, Na+, K+, and NH4+, as monomer, dimer, and trimer, are assigned. Negative mode include "[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-" as monomer, dimer, and trimer.
#' @param nls  character: vector of allowed neutral losses, i.e. c("[M+H-H2O]+").  if NULL, an extensive list derived from CAMERA's will be used.
#' @param scoring character: one of 'imss' , 'ramclustr', or 'auto'. default = 'auto'. see details.
#' @param plot.findmain logical: should pdf polts be generated for evaluation? detfault = TRUE. PDF saved to working.directory/spectra
#' @param writeMat logical: should individual .mat files (for MSFinder) be generated in a 'mat' subdirectory in the 'spectra' folder? default = TRUE.
#' @param writeMS logical: should individual .ms files (for Sirius) be generated in a 'ms' subdirectory in the 'spectra' folder? default = TRUE.  Note that no import functions are yet written for Sirius output.
#' @param use.z logical: if you have previously run the 'assign.z' function from ramclustR, there will be a slot reflecting the feature mass after accounting for charge (fm) - if TRUE this is used instead of feature m/z (fmz) in interpreting MS data and exporting spectra for annotation.
#' @details a partially annotated ramclustR object.  base structure is that of a standard R heirarchical clustering output, with additional slots described in ramclustR documentation (?ramclustR).  New slots added after using the interpretMSSpectrum functionality include those described below.
#' @return    $M:  The inferred molecular weight of the compound giving rise to the each spectrum
#' @return    $M.ppm:  The ppm error of all the MS signals annotated, high error values should be considered 'red flags'.
#' @return    $M.ann:  The annotated spectrum supporting the interpretation of M
#' @return    $use.findmain:  Logical vector indicating whether findmain scoring (TRUE) or ramclustR scoring (FALSE) was used to support inference of M.  By default, findmain scoring is used.  When ramclustR scoring differs from findmain scoring, the scoring metric which predicts higher M is selected.
#' @return    $M.ramclustr:  M selected using ramclustR scoring
#' @return    $M.ppm.ramclustr:  ppm error of M selected using ramclustR scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @return    $M.ann.ramclustr:  annotated spectrum supporting M using ramclustR scoring
#' @return    $M.nann.ramclustr:  number of masses annotated using ramclustR scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @return    $M.space.ramclustr:  the 'space' of scores between the best and second best ramclustR scores. Calculated as a ratio. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @return    $M.findmain:  M selected using findmain scoring
#' @return    $M.ppm.findmain:  ppm error of M selected using findmain scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @return    $M.ann.findmain:  annotated spectrum supporting M using findmain scoring
#' @return    $M.nann.findmain:  number of masses annotated using findmain scoring. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @return    $M.space.findmain:  the 'space' of scores between the best and second best findmain scores. Calculated as a ratio. Used to resolve concflicts between ramclustR and findmain M assignment when scoring = auto.
#' @references Jaeger C, ... Lisec J. Compound annotation in liquid chromatography/high-resolution mass spectrometry based metabolomics: robust adduct ion determination as a prerequisite to structure prediction in electrospray ionization mass spectra. Rapid Commun Mass Spectrom. 2017 Aug 15;31(15):1261-1266. doi: 10.1002/rcm.7905. PubMed PMID: 28499062.
#' @references Broeckling CD, Afsar FA, Neumann S, Ben-Hur A, Prenni JE. RAMClust: a novel feature clustering method enables spectral-matching-based annotation for metabolomics data. Anal Chem. 2014 Jul 15;86(14):6812-7. doi: 10.1021/ac501530d.  Epub 2014 Jun 26. PubMed PMID: 24927477.
#' @references Broeckling CD, Ganna A, Layer M, Brown K, Sutton B, Ingelsson E, Peers G, Prenni JE. Enabling Efficient and Confident Annotation of LC-MS Metabolomics Data through MS1 Spectrum and Time Prediction. Anal Chem. 2016 Sep 20;88(18):9226-34. doi: 10.1021/acs.analchem.6b02479. Epub 2016 Sep 8. PubMed PMID: 7560453.
#' @concept ramclustR
#' @concept RAMClustR
#' @concept metabolomics
#' @concept mass spectrometry
#' @concept clustering
#' @concept feature
#' @concept findMain
#' @concept interpretMSSpectrum
#' @concept xcms
#' @author Corey Broeckling
#' @export


do.findmain <- function(
    ramclustObj = NULL,
    cmpd = NULL,
    mode = "positive",
    mzabs.error = 0.005,
    ppm.error = 10,
    ads = NULL,
    nls = NULL,
    scoring = "auto",
    plot.findmain = TRUE,
    writeMat = TRUE,
    writeMS = TRUE,
    use.z = TRUE) {
  if (!requireNamespace("InterpretMSSpectrum", quietly = TRUE)) {
    stop("The use of this function requires package 'InterpretMSSpectrum'.")
  }

  if (is.null(ramclustObj)) {
    stop("must supply ramclustObj as input.  i.e. ramclustObj = RC", "\n")
  }

  score.options <- c("auto", "imss", "ramclustr")
  if (!any(scoring == score.options)) {
    stop("scoring must be set to one of 'auto', 'imss', or 'ramclustr' ")
  }
  if (use.z & any(names(ramclustObj) == "fm")) {
    use.mass <- "fm"
  } else {
    use.mass <- "fmz"
  }
  if (is.null(ads)) {
    if (grepl("p", mode)) {
      ads <- c(
        "[M+H]+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
        "[2M+H]+", "[2M+Na]+", "[2M+K]+", "[2M+NH4]+",
        "[3M+H]+", "[3M+Na]+", "[3M+K]+", "[3M+NH4]+"
      )
    }
    if (grepl("n", mode)) {
      ads <- c(
        "[M-H]-", "[M+Na-2H]-", "[M+K-2H]-", "[M+CH2O2-H]-",
        "[2M-H]-", "[2M+Na-2H]-", "[2M+K-2H]-", "[2M+CH2O2- H]-",
        "[3M-H]-", "[3M+Na-2H]-", "[3M+K-2H]-", "[3M+CH2O2- H]-"
      )
    }
    if (is.null(ads)) {
      stop("please define adducts using 'ads' or set mode to either'positive' or 'negative'")
    }
  }
  if (is.null(nls)) {
    if (grepl("p", mode)) {
      nls <- c(
        "[M+H-COCH2]+", "[M+H-C2H3NO]+", "[M+H-H2O]+",
        "[M+H-NH3]+", "[M+H-CO2]+", "[M+H-NH3-CO2]+",
        "[M+H-NH3-HCOOH]+", "[M+H-NH3-H2O]+", "[M+H-NH3-COCH2]+",
        "[M+H-S]+", "[M+H-S-NH3-HCOOH]+", "[M+H-H4O2]+",
        "[M+H-CH2]+", "[M+H-CH4O3]+", "[M+H-O]+", "[M+H-C2H2]+",
        "[M+H-C2H4]+", "[M+H-CO]+", "[M+H-C3H6]+", "[M+H-C2H4O]+",
        "[M+H-C4H6]+", "[M+H-C3H4O]+", "[M+H-C4H8]+",
        "[M+H-C5H8]+", "[M+H-C4H6O]+", "[M+H-C5H10]+",
        "[M+H-C6H12]+", "[M+H-C4H8O2]+", "[M+H-H2O-HCOOH]+",
        "[M+H-CH4]+", "[M+H-CH2O]+", "[M+H-C2H6]+",
        "[M+H-CH3OH]+", "[M+H-C3H4]+", "[M+H-C3H6O]+",
        "[M+H-CO2-C3H6]+", "[M+H-SO3]+", "[M+H-SO3-H2O]+",
        "[M+H-SO3-H2O-NH3]+", "[M+H-NH3-C3H4]+", "[M+H-H2O-CO2]+",
        "[M+H-H2O-H2O-C2H4O]+", "[M+H-NH3-CO-CO]+",
        "[M+H-NH3-CO-C2OH2]+", "[M+H-C8H6O]+", "[M+H-C8H6O-NH3]+",
        "[M+H-C8H6O-H2O]+", "[M+H-C2H2O2]+", "[M+H-C2H4O2]+",
        "[M+H-C5H8O]+", "[M+H-NH3-CO2-CH2O]+", "[M+H-N2H6-CO2-H2O]+",
        "[M+H-NH3-CO2-C3H4O]+", "[M+H-NH3-CO2-C5H8]+",
        "[M+H-C2O4H4]+", "[M+H-C2H4-CO2]+", "[M+H-C2H4-HCOOH]+",
        "[M+H-NH3-H4O2]+", "[M+H-H2O-C2H2O2]+", "[M+H-COCH2-C4H8]+",
        "[M+H-NH3-NH3-C3H4]+", "[M+H-C2H4O2-CH3OH]+",
        "[M+H-C3H6O-CH3OH]+", "[M+H-NH3-CO-COCH2-C4H6O]+",
        "[M+H-C4H6-H2O]+", "[M+H-C4H6-C2H4]+", "[M+H-C4H6-NH3-H2O]+",
        "[M+H-C4H6-COCH2]+", "[M+H-C8H12O]+", "[M+H-C3H4O-C4H6]+",
        "[M+H-C3H4O-C4H8O2]+", "[M+H-C4H8-C4H6]+", "[M+H-NH3-HCOOH-CH3OH]+",
        "[M+H-NH3-C2H6]+", "[M+H-NH3-C8H6O-CH2]+", "[M+H-NH3-C3H4-COCH2]+",
        "[M+H-C3H9N]+", "[M+H-C3H9N-C2H4O2]+", "[M+H-C6H10O7]+",
        "[M+H-C6H10O7+H2O]+", "[M+H-C6H10O7-H4O2]+",
        "[M+H-C6H12O6]+", "[M+H-C6H12O6+H2O]+", "[M+H-C6H12O6-H2O]+",
        "[M+H-C5H10O5]+", "[M+H- C5H10O5+H2O]+", "[M+H- C5H10O5-H2O]+",
        "[M+H-C2H8NO4P]+", "[M+H-H6O3-CO]+", "[M+H-C6H13NO2]+",
        "[M+H-C5H11NO2]+", "[M+H-CH3S]+", "[M+H-C8H8O2]+",
        "[M+H-C12H22O11]+", "[M+H-C12H24O12]+", "[M+H-C12H20O10]+"
      )
    }
    if (grepl("n", mode)) {
      nls <- c(
        "[M-H-NH3]-", "[M-H-H2O]-", "[M-H-COCH2]-",
        "[M-H-CO2]-", "[M-H-NH3-CO2]-", "[M-H-NH3-HCOOH]-",
        "[M-H-NH3-H2O]-", "[M-H-NH3-C2OH2]-", "[M-H-S]-",
        "[M-H-S-NH3-HCOOH]-", "[M-H-H4O2]-", "[M-H-CH2]-",
        "[M-H-O]-", "[M-H-C2H2]-", "[M-H-C2H4]-", "[M-H-CO]-",
        "[M-H-C3H6]-", "[M-H-C2H4O]-", "[M-H-C4H6]-",
        "[M-H-C3H4O]-", "[M-H-C4H8]-", "[M-H-C5H8]-",
        "[M-H-C4H6O]-", "[M-H-C5H10]-", "[M-H-C6H12]-",
        "[M-H-C4H8O2]-", "[M-H-H2O-HCOOH]-", "[M-H-CH4]-",
        "[M-H-CH2O]-", "[M-H-C2H6]-", "[M-H-CH3OH]-",
        "[M-H-C3H4]-", "[M-H-C3H6O]-", "[M-H-CO2-C3H6]-",
        "[M-H-SO3]-", "[M-H-SO3-H2O]-", "[M-H-SO3-H2O-NH3]-",
        "[M-H-NH3-C3H4]-", "[M-H-H2O-CO2]-", "[M-H-H2O-H2O-C2H4O]-",
        "[M-H-NH3-CO-CO]-", "[M-H-NH3-CO-COCH2]-", "[M-H-C8H6O]-",
        "[M-H-C8H6O-NH3]-", "[M-H-C8H6O-H2O]-", "[M-H-C2H2O2]-",
        "[M-H-C2H4O2]-", "[M-H-C5H8O]-", "[M-H-NH3-CO2-CH2O]-",
        "[M-H-NH3-CO2-NH3-H2O]-", "[M-H-NH3-CO2-C3H4O]-",
        "[M-H-NH3-CO2-C5H8]-", "[M-H-HCOOH-HCOOH]-",
        "[M-H-C2H4-CO2]-", "[M-H-C2H4-HCOOH]-", "[M-H-NH3-H2O-H2O]-",
        "[M-H-H2O-C2H2O2]-", "[M-H-C2OH2-C4H8]-", "[M-H-NH3-NH3-C3H4]-",
        "[M-H-C2H4O2-CH3OH]-", "[M-H-C3H6O-CH3OH]-",
        "[M-H-NH3-CO-COCH2-C4H6O]-", "[M-H-C4H6-H2O]-",
        "[M-H-C4H6-C2H4]-", "[M-H-C4H6-NH3-H2O]-", "[M-H-C4H6-COCH2]-",
        "[M-H-C4H6-C4H6O]-", "[M-H-C3H4O-C4H6]-", "[M-H-C3H4O-C4H8O2]-",
        "[M-H-C4H8-C4H6]-", "[M-H-NH3-HCOOH-CH3OH]-",
        "[M-H-NH3-C2H6]-", "[M-H-NH3-C8H6O-CH2]-", "[M-H-NH3-C3H4-COCH2]-",
        "[M-H-C3H9N]-", "[M-H-C3H9N-C2H4O2]-", "[M-H-C6H10O7]-",
        "[M-H-C6H10O7+H2O]-", "[M-H-C6H10O7-H4O2]-",
        "[M-H-C6H12O6]-", "[M-H-C6H12O6+H2O]-", "[M-H-C6H12O6-H2O]-",
        "[M-H-C5H10O5]-", "[M-H-C5H10O5+H2O]-", "[M-H-C5H10O5-H2O]-",
        "[M-H-C2H8NO4P]-", "[M-H-CH6O4]-", "[M-H-C6H13NO2]-",
        "[M-H-C5H11NO2]-", "[M-H-CH3S]-", "[M-H-C8H8O2]-",
        "[M-H-C12H22O11]-", "[M-H-C12H24O12]-", "[M-H-C12H20O10]-"
      )
    }
    if (is.null(nls)) {
      stop("please define neutral losses using 'nls' or set mode to either'positive' or 'negative'")
    }
  }
  M.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ppm.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ann.findmain <- as.list(rep(NA, max(ramclustObj$featclus)))
  M.nann.findmain <- rep(NA, max(ramclustObj$featclus))
  M.space.findmain <- rep(NA, max(ramclustObj$featclus))
  M.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ppm.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.rank.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.ann.ramclustr <- as.list(rep(NA, max(ramclustObj$featclus)))
  M.nann.ramclustr <- rep(NA, max(ramclustObj$featclus))
  M.space.ramclustr <- rep(NA, max(ramclustObj$featclus))
  if (is.null(cmpd)) {
    cmpd <- (1:max(ramclustObj$featclus))
  }
  for (cl in cmpd) {
    s <- data.frame(mz = ramclustObj[[use.mass]][which(ramclustObj$featclus ==
      cl)], int = ramclustObj$msint[which(ramclustObj$featclus ==
      cl)])
    s <- s[order(s$mz), ]
    out <- InterpretMSSpectrum::findMAIN(
      s,
      rules = c(ads),
      adducthyp = ads[grep("[M", ads, fixed = TRUE)],
      ionmode = mode,
      mzabs = mzabs.error,
      ppm = ppm.error
    )
    summarytable <- summary(out)

    # summarytable[which(is.na(summarytable[, "medppm"])),
    #              "medppm"] <- 2 * ppm.error
    # summarytable[which(summarytable[, "medppm"] == 0), "medppm"] <- ppm.error
    M.findmain[cl] <- summarytable[1, "neutral_mass"]
    M.ppm.findmain[cl] <- summarytable[1, "medppm"]
    M.ann.findmain[[cl]] <- out[[1]]
    M.nann.findmain[cl] <- summarytable[1, "adducts_explained"]
    M.space.findmain[cl] <- summarytable[1, "total_score"] / summarytable[2, "total_score"]
    out <- InterpretMSSpectrum::findMAIN(
      s,
      adductmz = NULL,
      ionmode = mode,
      rules = c(ads, nls),
      adducthyp = ads[grep("[M", ads, fixed = TRUE)],
      ms2spec = NULL,
      mzabs = mzabs.error,
      ppm = ppm.error,
      mainpkthr = 0.1,
      collapseResults = FALSE
    )
    summarytable <- summary(out)
    # summarytable[which(is.na(summarytable[, "medppm"])),
    #              "medppm"] <- 2 * ppm.error
    # summarytable[which(summarytable[, "medppm"] == 0), "medppm"] <- ppm.error
    summaryscores <- sapply(1:length(out), FUN = function(x) {
      is.adduct <- which(out[[x]][, "adduct"] %in% ads)
      det.adducts <- out[[x]][is.adduct, "adduct"]
      is.nl <- which(out[[x]][, "adduct"] %in% nls)
      keep <- is.adduct
      int <- sum(out[[x]][is.adduct, "int"])
      if (length(is.nl) > 0) {
        int <- sum(int, sum(0.5 * (out[[x]][is.nl, "int"])))
        keep <- c(keep, is.nl)
      }
      int <- int / sum(out[[x]][, "int"]) * length(is.adduct)
      if (any(grepl("[M+H-NH3]+", det.adducts)) & length(det.adducts ==
        1)) {
        int <- int / 2
      }
      mzerr <- (1 / summarytable[x, "medppm"])^0.5
      return(int * mzerr)
    })
    best <- which.max(summaryscores)
    if (length(best) == 0) {
      best <- 1
    }
    M.space.ramclustr[cl] <- {
      sort(summaryscores, decreasing = TRUE)[1] / sort(summaryscores,
        decreasing = TRUE
      )[2]
    }
    if (is.na(M.space.ramclustr[cl])) {
      M.space.ramclustr[cl] <- 1
    }
    M.ramclustr[cl] <- summarytable[best, "neutral_mass"]
    M.ppm.ramclustr[cl] <- summarytable[best, "medppm"]
    M.rank.ramclustr[cl] <- best
    M.ann.ramclustr[[cl]] <- out[[best]]
    M.nann.ramclustr[cl] <- summarytable[best, "adducts_explained"]
    if (10 * round(cl / 10, digits = 0) == cl) {
      cat(cl, "of", max(ramclustObj$featclus), "\n")
    }
  }
  ramclustObj$M.ramclustr <- M.ramclustr
  ramclustObj$M.ppm.ramclustr <- M.ppm.ramclustr
  ramclustObj$M.rank.ramclustr <- M.rank.ramclustr
  ramclustObj$M.ann.ramclustr <- M.ann.ramclustr
  ramclustObj$M.nann.ramclustr <- M.nann.ramclustr
  ramclustObj$M.space.ramclustr <- M.space.ramclustr
  ramclustObj$M.findmain <- M.findmain
  ramclustObj$M.ppm.findmain <- M.ppm.findmain
  ramclustObj$M.ppm.findmain[which(is.na(ramclustObj$M.ppm.findmain))] <- {
    ppm.error * 2
  }

  ramclustObj$M.ann.findmain <- M.ann.findmain
  ramclustObj$M.nann.findmain <- M.nann.findmain
  ramclustObj$use.findmain <- rep(FALSE, length(M.ppm.findmain))
  ramclustObj$M.space.findmain <- M.space.findmain
  if (scoring == "ramclustr") {
    ramclustObj$use.findmain <- rep(FALSE, length(ramclustObj$use.findmain))
  }
  if (scoring == "imss") {
    ramclustObj$use.findmain <- rep(TRUE, length(ramclustObj$use.findmain))
  }
  if (scoring == "auto") {
    resolve <- which(abs(ramclustObj$M.ramclustr - ramclustObj$M.findmain) >
      (2 * mzabs.error))
    for (i in resolve) {
      ppm.rat <- (ramclustObj$M.ppm.ramclustr[i]^0.5 / ramclustObj$M.ppm.findmain[i]^0.5)
      if (is.na(ppm.rat)) {
        ppm.rat <- 1
      }
      ppm.rat <- 1 / ppm.rat
      sel.rat <- (ramclustObj$M.space.ramclustr[i] / ramclustObj$M.space.findmain[i])
      nann.rat <- (ramclustObj$M.nann.ramclustr[i] / ramclustObj$M.nann.findmain[i])
      master.rat <- ppm.rat * sel.rat * nann.rat
      if (master.rat < 1) {
        ramclustObj$use.findmain[i] <- TRUE
      }
    }
  }

  ramclustObj$M.ann <- ramclustObj$M.ann.ramclustr
  ramclustObj$precursor.mz <- rep(NA, length(ramclustObj$M.ann))
  ramclustObj$precursor.type <- rep(NA, length(ramclustObj$M.ann))
  ramclustObj$M <- ramclustObj$M.ramclustr

  change <- which(ramclustObj$use.findmain)
  ramclustObj$M[change] <- ramclustObj$M.findmain[change]
  if (length(change) > 0) {
    for (i in change) {
      ramclustObj$M.ann[[i]] <- ramclustObj$M.ann.findmain[[i]]
    }
  }


  for (i in 1:length(ramclustObj$precursor.mz)) {
    if (is.vector(ramclustObj$M.ann[[i]])) next
    for (j in 1:length(ads)) {
      m <- which(ramclustObj$M.ann[[i]]$label == ads[[j]])
      if (length(m) == 0) next
      ramclustObj$precursor.mz[i] <- ramclustObj$M.ann[[i]][m, "mz"]
      ramclustObj$precursor.type[i] <- ads[j]
    }
  }

  if (plot.findmain) {
    cat("plotting findmain annotation results", "\n")
    if (!dir.exists("spectra")) {
      dir.create("spectra")
    }
    pdf("spectra/findmainPlots.pdf", width = 15, height = 7)
    par(mfrow = c(1, 2))
    par(xpd = TRUE)
    for (cl in cmpd) {
      InterpretMSSpectrum::PlotSpec(
        x = ramclustObj$M.ann.ramclustr[[cl]],
        txt = ramclustObj$M.ann.ramclustr[[cl]][, c("mz", "adduct")],
        cutoff = 0, masslab = 0,
        ylim = c(0, 1.1 * max(ramclustObj$M.ann.ramclustr[[cl]][, 2]))
      )
      title(main = list(
        paste(
          cl, ":", "M.ramclustr =",
          round(ramclustObj$M.ramclustr[cl], digits = 4),
          "( +/-", round(ramclustObj$M.ppm.ramclustr[cl],
            digits = 1
          ), "ppm )"
        ),
        font = if (ramclustObj$use.findmain[cl]) {
          1
        } else {
          2
        },
        col = if (ramclustObj$use.findmain[cl]) {
          1
        } else {
          2
        }
      ))
      InterpretMSSpectrum::PlotSpec(
        x = ramclustObj$M.ann.findmain[[cl]], txt = ramclustObj$M.ann.findmain[[cl]][
          ,
          c("mz", "adduct")
        ], cutoff = 0, masslab = 0,
        ylim = c(0, 1.1 * max(ramclustObj$M.ann.ramclustr[[cl]][
          ,
          2
        ]))
      )
      title(main = list(paste(
        cl, ":", "M.findmain =",
        round(ramclustObj$M.findmain[cl], digits = 4),
        "( +/-", round(ramclustObj$M.ppm.findmain[cl],
          digits = 1
        ), "ppm )"
      ), font = if (ramclustObj$use.findmain[cl]) {
        2
      } else {
        1
      }, col = if (ramclustObj$use.findmain[cl]) {
        2
      } else {
        1
      }))
    }
    dev.off()
  }

  if (writeMat) {
    if (!dir.exists("spectra")) {
      dir.create("spectra")
    }
    dir.create("spectra/mat")
    for (cl in cmpd) {
      ms <- ramclustObj$M.ann[[cl]]
      prcr <- which(ms[, "adduct"] %in% ads)
      prcr <- prcr[which.max(ms[prcr, "int"])]
      prcmz <- ms[prcr, "mz"]
      prctype <- ms[prcr, "adduct"]
      out <- paste("NAME: ", ramclustObj$cmpd[cl], "\n",
        "RETENTIONTIME: ", round(ramclustObj$clrt[cl],
          digits = 2
        ), "\n", "PRECURSORMZ: ", prcmz,
        "\n", "PRECURSORTYPE: ", prctype, "\n", "IONTYPE: ",
        mode, "\n", "SPECTRUMTYPE: Centroid", "\n",
        if ((!is.null(ramclustObj$msmsint))) {
          paste("COLLISIONENERGY: ", as.character(ramclustObj$ExpDes[[2]][which(row.names(ramclustObj$ExpDes[[2]]) ==
            "CE2"), 1]), "\n", sep = "")
        }, "MSTYPE: ", "MS1", "\n", "Num Peaks: ", nrow(ms),
        "\n",
        sep = ""
      )
      for (i in 1:nrow(ms)) {
        out <- paste(out, ms[i, 1], " ", ms[i, 2], "\n",
          sep = ""
        )
      }
      if (!is.null(ramclustObj$msmsint)) {
        do <- which(ramclustObj$featclus == cl)
        if (length(do) > 0) {
          msms <- cbind(
            mz = ramclustObj[[use.mass]][do],
            int = ramclustObj$msmsint[do]
          )
          msms <- msms[which(msms[, "mz"] <= (prcmz +
            3)), , drop = FALSE]
          msms <- msms[order(msms[, "int"], decreasing = TRUE), ,
            drop = FALSE
          ]
          if (nrow(msms) > 0) {
            out <- paste(out, "MSTYPE:", "MS2", "\n",
              "Num Peaks: ", nrow(msms), "\n",
              sep = ""
            )
            for (i in 1:nrow(msms)) {
              out <- paste(out, msms[i, 1], " ", msms[
                i,
                2
              ], "\n", sep = "")
            }
          }
        }
      } else {
        do <- which(ramclustObj$featclus == cl)
        if (length(do) > 0) {
          msms <- cbind(
            mz = ramclustObj[[use.mass]][do],
            int = ramclustObj$msint[do]
          )
          msms <- msms[which(msms[, "mz"] <= (prcmz +
            3)), , drop = FALSE]
          msms <- msms[order(msms[, "int"], decreasing = TRUE), ,
            drop = FALSE
          ]
          if (nrow(msms) > 0) {
            out <- paste(out, "MSTYPE:", "MS2", "\n",
              "Num Peaks: ", nrow(msms), "\n",
              sep = ""
            )
            for (i in 1:nrow(msms)) {
              out <- paste(out, msms[i, 1], " ", msms[
                i,
                2
              ], "\n", sep = "")
            }
          }
        }
      }
      write(out, file = paste0(
        "spectra/mat/", ramclustObj$cmpd[cl],
        ".mat"
      ))
    }
  }
  if (writeMS) {
    if (!dir.exists("spectra")) {
      dir.create("spectra")
    }
    dir.create("spectra/ms")
    for (cl in cmpd) {
      ms <- ramclustObj$M.ann[[cl]]
      m1ads <- ads[grep("[M", ads, fixed = TRUE)]
      if (length(m1ads) == 0) {
        next
      } else {
        prcr <- which(ms[, "adduct"] %in% m1ads)
      }
      prcr <- prcr[which.max(ms[prcr, "int"])]
      prcmz <- ms[prcr, "mz"]
      prctype <- ms[prcr, "adduct"]
      out <- paste(">compound ", ramclustObj$cmpd[cl],
        "\n", ">parentmass ", prcmz, "\n", ">ionization ",
        prctype, "\n", "\n",
        sep = ""
      )
      ms <- ms[which((abs(ms[, "mz"] - prcmz) < 5.5) |
        (abs(prcmz - ms[, "mz"]) < 0.2)), ]
      out <- paste(out, ">ms1peaks", "\n", sep = "")
      for (i in 1:nrow(ms)) {
        out <- paste(out, ms[i, 1], " ", ms[i, 2], "\n",
          sep = ""
        )
      }
      if (!is.null(ramclustObj$msmsint)) {
        do <- which(ramclustObj$featclus == cl)
        if (length(do) > 0) {
          msms <- cbind(
            mz = ramclustObj[[use.mass]][do],
            int = ramclustObj$msmsint[do]
          )
          msms <- msms[which(msms[, "mz"] <= (prcmz +
            3)), , drop = FALSE]
          msms <- msms[order(msms[, "int"], decreasing = TRUE), ,
            drop = FALSE
          ]
          if (nrow(msms) > 0) {
            out <- paste(out, "\n", ">collision ", ramclustObj$ExpDes$instrument[which(dimnames(ramclustObj$ExpDes$instrument)[[1]] ==
              "CE2"), 1], "\n", sep = "")
            for (i in 1:nrow(msms)) {
              out <- paste(out, msms[i, 1], " ", msms[
                i,
                2
              ], "\n", sep = "")
            }
          }
        }
      } else {
        do <- which(ramclustObj$featclus == cl)
        if (length(do) > 0) {
          msms <- cbind(
            mz = ramclustObj[[use.mass]][do],
            int = ramclustObj$msint[do]
          )
          msms <- msms[which(msms[, "mz"] <= (prcmz +
            3)), , drop = FALSE]
          msms <- msms[order(msms[, "int"], decreasing = TRUE), ,
            drop = FALSE
          ]
          if (nrow(msms) > 0) {
            out <- paste(out, "\n", ">collision ", ramclustObj$ExpDes$instrument[which(dimnames(ramclustObj$ExpDes$instrument)[[1]] ==
              "CE1"), 1], "\n", sep = "")
            for (i in 1:nrow(msms)) {
              out <- paste(out, msms[i, 1], " ", msms[
                i,
                2
              ], "\n", sep = "")
            }
          }
        }
      }
      write(out, file = paste0(
        "spectra/ms/", ramclustObj$cmpd[cl],
        ".ms"
      ))
    }
  }
  ramclustObj$history$do.findmain <- paste(
    " Molecular weight was inferred from in-source spectra (Broeckling 2016) using the do.findmain function, which calls the ",
    "interpretMSSpectrum package (Jaeger 2016). ",
    "Parameters for do.findmain were set to: ",
    "mode = ", mode, ", mzabs.error = ", mzabs.error, ", ppm.error = ",
    ppm.error, ", ads = ", paste(ads, collapse = " "), ", nls = ",
    paste(nls, collapse = " "), ", scoring = ", scoring,
    ", and use.z = ", use.z, ".",
    sep = ""
  )
  cat("finished", "\n")
  return(ramclustObj)
}
